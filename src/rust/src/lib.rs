use extendr_api::prelude::*;
use ferx_core::cancel::CancelFlag;
use ferx_core::types::*;
use nalgebra::DMatrix;
use std::path::Path;

// ---------------------------------------------------------------------------
//  R interrupt polling
// ---------------------------------------------------------------------------
//
// When R is running long native code it cannot respond to Ctrl-C — R's SIGINT
// handler only sets a flag that is checked on re-entry into R. To make our
// fit cancellable we:
//   1. Spawn the fit on a worker thread with a shared CancelFlag.
//   2. Poll on the main (R) thread via `pending_interrupt()`, which wraps
//      `R_CheckUserInterrupt` in `R_ToplevelExec` so the longjmp is contained
//      (avoids UB from skipping Rust frames).
//   3. When an interrupt is pending, flip the cancel flag and let the worker
//      exit cooperatively; ferx_core::fit returns Err("cancelled by user").
//
// Declared here rather than pulling libR-sys to keep the dependency surface
// small. These symbols are stable parts of R's public API (R.h / Rinterface.h).

extern "C" {
    fn R_CheckUserInterrupt();
    fn R_ToplevelExec(
        fun: extern "C" fn(*mut std::ffi::c_void),
        data: *mut std::ffi::c_void,
    ) -> std::ffi::c_int;
}

extern "C" fn check_interrupt_cb(_: *mut std::ffi::c_void) {
    unsafe { R_CheckUserInterrupt() };
}

/// Returns true if an R interrupt (Ctrl-C) is pending. Safe to call repeatedly
/// on the R main thread — does not longjmp because the check is wrapped in
/// `R_ToplevelExec`.
fn pending_interrupt() -> bool {
    unsafe { R_ToplevelExec(check_interrupt_cb, std::ptr::null_mut()) == 0 }
}

/// Poll interval for the interrupt-check loop. 100ms keeps Ctrl-C responsive
/// while adding negligible overhead to the worker.
const POLL_MS: u64 = 100;

/// Fit a NLME model to data.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param method Character vector of estimation methods. A single element runs
///   one stage (e.g. "focei"); multiple elements chain stages, each seeded with
///   the previous stage's converged parameters (e.g. c("saem", "focei")).
/// @param covariance Run covariance step (TRUE/FALSE)
/// @param verbose Print progress (TRUE/FALSE)
/// @param bloq_method BLOQ handling: "drop", "m3", or "" to use the model default
/// @param threads Number of rayon worker threads for the per-subject parallel
///   loops. Pass `0` (or any value `<= 0`) to leave rayon's global pool alone
///   (one worker per logical CPU). Positive values run this fit inside a
///   scoped local pool of that size.
/// @param mu_referencing Use mu-referencing for ETA initialisation (TRUE/FALSE)
/// @param sir Run SIR uncertainty estimation as a post-fit step (TRUE/FALSE).
///   Tuning knobs (`sir_samples`, `sir_resamples`, `sir_seed`) still flow
///   through `settings`; only the on/off toggle is top-level.
/// @param settings_keys Parallel vector of setting names (pre-stringified).
///   Used together with `settings_values` to pass generic estimation-method
///   options (e.g. `n_exploration`, `sir_samples`, `optimizer`,
///   `inner_maxiter`, `inner_tol`, `steihaug_max_iters`) without needing a
///   new R-Rust argument per option. Keys that duplicate a dedicated
///   argument are rejected so there is a single source of truth. Unknown
///   keys and malformed values also raise an error.
/// @param settings_values Parallel vector of setting values as strings;
///   the Rust side parses each value according to the key's expected type.
/// @return Named list with fit results
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_fit(
    model_path: &str,
    data_path: &str,
    method: Vec<String>,
    covariance: bool,
    verbose: bool,
    bloq_method: &str,
    threads: i32,
    mu_referencing: bool,
    sir: bool,
    gradient: &str,
    settings_keys: Vec<String>,
    settings_values: Vec<String>,
) -> List {
    let mut parsed =
        match ferx_core::parser::model_parser::parse_full_model_file(Path::new(model_path)) {
            Ok(p) => p,
            Err(e) => throw_r_error(format!("Error parsing model: {e}")),
        };

    // Build fit options
    let mut opts = parsed.fit_options.clone();

    // Apply generic `settings` list before the dedicated args below — the
    // dedicated args are the single source of truth for their keys and always
    // win. Reserved keys (those with a dedicated R argument) are rejected so
    // the precedence rule is explicit rather than silently clobbered.
    //
    // Data-selection keys (`ignore`, `accept`, `ignore_subjects`) are NOT
    // reserved — they flow through `apply_fit_option` which populates
    // opts.ignore_exprs / opts.accept_exprs / opts.ignore_subjects. Data
    // reading happens AFTER this loop so the SelectionFilter sees the merged
    // model-file + R-call conditions.
    if settings_keys.len() != settings_values.len() {
        throw_r_error(format!(
            "Error: settings keys/values length mismatch ({} vs {})",
            settings_keys.len(),
            settings_values.len()
        ));
    }
    // Reserved keys: these have dedicated ferx_fit() arguments. Keeping them
    // out of `settings` means there is one source of truth per value.
    // NOTE: `optimizer`, `inner_maxiter`, `inner_tol`, and
    // `steihaug_max_iters` intentionally flow through `settings` — they
    // were previously candidates for dedicated args but are fit-method
    // tuning knobs that belong alongside `n_exploration`, `sir_samples`,
    // etc.
    const RESERVED: &[&str] = &[
        "method",
        "covariance",
        "verbose",
        "bloq_method",
        "bloq",
        "threads",
        "sir",
        "gradient",
        "gradient_method",
    ];
    for (k, v) in settings_keys.iter().zip(settings_values.iter()) {
        let key = k.trim();
        if RESERVED
            .iter()
            .any(|r| r.eq_ignore_ascii_case(key))
        {
            throw_r_error(format!(
                "Error: setting `{key}` conflicts with a dedicated ferx_fit() argument — pass it via that argument instead"
            ));
        }
        match ferx_core::parser::model_parser::apply_fit_option(&mut opts, key, v) {
            Ok(true) => {}
            Ok(false) => throw_r_error(format!("Error: unknown fit setting `{key}`")),
            Err(e) => throw_r_error(format!("Error: {e}")),
        }
    }

    // Re-apply the (now settings-merged) ODE solver tolerances onto the model's
    // OdeSpec so call-time `ode_reltol` / `ode_abstol` / `ode_max_steps`
    // overrides take effect. The parser already baked the [fit_options] values;
    // this lets a `ferx_fit(settings = ...)` override win. No-op for analytical
    // models.
    parsed.model.sync_ode_solver_opts(&opts);

    // Read data via read_population_for, which handles [covariates] validation,
    // [data_selection] filters, and TTE endpoint routing in one call.
    // This replaces the previous 4-way dispatch that could not pass tte_cmts to
    // the reader (causing TTE rows to land in the Gaussian vectors instead of
    // subject.obs_records). It also fixes the pre-existing gap where the combined
    // covariates + filter case fell through to the filter-only path, losing the
    // covariate table.
    let (population, covariate_table) = {
        use ferx_core::api::read_population_for;
        use ferx_core::io::datareader::SelectionFilter;
        let filter = match SelectionFilter::from_opts(
            &opts.ignore_exprs,
            &opts.accept_exprs,
            &opts.ignore_subjects,
        ) {
            Ok(f) => f,
            Err(e) => throw_r_error(format!("Error in [data_selection]: {e}")),
        };
        let filter_opt = if filter.is_empty() { None } else { Some(&filter) };
        match read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            opts.iov_column.as_deref(),
            filter_opt,
        ) {
            Ok(result) => result,
            Err(e) => throw_r_error(format!("Error reading data: {e}")),
        }
    };

    if method.is_empty() {
        throw_r_error("Error: `method` must contain at least one estimation method");
    }
    let chain: Vec<EstimationMethod> = match method.iter().map(|m| parse_method(m)).collect() {
        Ok(v) => v,
        Err(e) => throw_r_error(format!("{e}")),
    };
    let final_method = *chain.last().unwrap();
    // The reported `interaction` flag must reflect the last *estimating* stage
    // — IMP is a diagnostic terminal stage that does not update parameters, so
    // a chain like `c("focei", "imp")` should still report `interaction = TRUE`.
    let last_estimator = chain
        .iter()
        .rev()
        .find(|m| **m != EstimationMethod::Imp)
        .copied()
        .unwrap_or(final_method);
    opts.method = final_method;
    opts.interaction = last_estimator == EstimationMethod::FoceI;
    opts.methods = if chain.len() > 1 { chain } else { Vec::new() };
    opts.run_covariance_step = covariance;
    opts.verbose = verbose;
    opts.mu_referencing = mu_referencing;
    opts.sir = sir;
    opts.threads = if threads > 0 {
        Some(threads as usize)
    } else {
        None
    };

    // Optional R-side override for BLOQ handling. Empty string → keep whatever
    // the model file specified.
    match bloq_method.trim().to_lowercase().as_str() {
        "" => {}
        "m3" => opts.bloq_method = BloqMethod::M3,
        "drop" | "none" | "ignore" => opts.bloq_method = BloqMethod::Drop,
        other => {
            rprintln!(
                "Unknown bloq_method '{}' — expected 'm3' or 'drop' (falling back to model default)",
                other
            );
        }
    }
    // Mirror onto the compiled model so likelihood functions pick it up.
    parsed.model.bloq_method = opts.bloq_method;

    // Gradient method override. Empty string → keep model/option default.
    match gradient.trim().to_lowercase().as_str() {
        "" | "auto" => opts.gradient_method = ferx_core::GradientMethod::Auto,
        "ad" | "autodiff" => opts.gradient_method = ferx_core::GradientMethod::Ad,
        "fd" | "finite" | "finite_difference" => {
            opts.gradient_method = ferx_core::GradientMethod::Fd
        }
        other => {
            rprintln!(
                "Unknown gradient method '{}' — expected 'auto', 'ad', or 'fd' (falling back to auto)",
                other
            );
            opts.gradient_method = ferx_core::GradientMethod::Auto;
        }
    }
    parsed.model.gradient_method = opts.gradient_method;

    // Install a cancellation token so Ctrl-C on the R console aborts the fit.
    let cancel = CancelFlag::new();
    opts.cancel = Some(cancel.clone());

    // Build initial parameters (initial values now live in [parameters] block)
    let init_params = parsed.model.default_params.clone();

    // Run the fit on a worker thread and poll for R interrupts on the main
    // thread. Scoped threads let us borrow `parsed.model` and `population`
    // without cloning them. The worker exits cooperatively when `cancel` is
    // set — there is a small (poll-interval bounded) tail where the worker
    // drains its last iteration before returning Err("cancelled by user").
    //
    // The worker signals completion on a channel so the main thread wakes the
    // instant the fit finishes. `POLL_MS` therefore bounds only interrupt
    // latency, not the wall time of a fast fit — a previous `is_finished()` +
    // `sleep(POLL_MS)` loop imposed a ~POLL_MS floor on every call, which
    // dominated short fits (e.g. fixed-parameter MAP, where the fit is sub-ms).
    let result = std::thread::scope(|s| {
        let (done_tx, done_rx) = std::sync::mpsc::channel::<()>();
        // Bind borrows explicitly so the `move` closure captures these Copy
        // references (and `done_tx` by value) rather than moving `parsed` /
        // `population`, which are still needed after the fit returns.
        let model_ref = &parsed.model;
        let pop_ref = &population;
        let init_ref = &init_params;
        let opts_ref = &opts;
        let handle = s.spawn(move || {
            let r = ferx_core::fit(model_ref, pop_ref, init_ref, opts_ref);
            // Ignore send errors: the receiver is only dropped once we stop
            // waiting, which happens after the worker has already finished.
            let _ = done_tx.send(());
            r
        });

        loop {
            match done_rx.recv_timeout(std::time::Duration::from_millis(POLL_MS)) {
                // Worker finished (or its sender was dropped on panic) — join
                // below collects the real result / propagates the panic.
                Ok(()) | Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => break,
                Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                    // Still running: service R interrupts, then keep waiting
                    // for the worker to drain and report (cancelled or not).
                    if pending_interrupt() {
                        cancel.cancel();
                    }
                }
            }
        }

        handle.join()
    });

    let mut result = match result {
        Ok(Ok(r)) => r,
        Ok(Err(_)) if cancel.is_cancelled() => {
            // Raise a proper R error so the user sees a clean condition
            // instead of a silent empty-list return. A small one-time leak
            // of the worker's locals is acceptable here — one R session's
            // worth of cancelled fits won't add up to anything meaningful.
            throw_r_error("ferx_fit: cancelled by user");
        }
        Ok(Err(e)) => throw_r_error(format!("Fit error: {e}")),
        Err(_) => throw_r_error("Fit error: worker thread panicked"),
    };

    // Record source-file provenance on the result so `ferx_sir(fit)` can
    // refuse to run against a tampered model or dataset later. We compute
    // the hashes here (not inside ferx_core::fit) because the R binding
    // owns the path strings the user supplied. Hashing failures are
    // non-fatal — the fit already succeeded, and missing hashes just mean
    // `ferx_sir` will fall back to a no-verification call.
    result.model_path = Some(model_path.to_string());
    result.data_path = Some(data_path.to_string());
    result.model_hash = ferx_core::io::hash::sha256_file(Path::new(model_path)).ok();
    result.data_hash = ferx_core::io::hash::sha256_file(Path::new(data_path)).ok();
    result.model_text = std::fs::read_to_string(model_path).ok();
    // Attach the covariate table built during the strict read above (None when
    // the model has no `[covariates]` block). `fit()` itself leaves this None.
    result.covariate_table = covariate_table;

    // Convert to R list
    fit_result_to_list(&result, &population, &parsed.model)
}

/// Simulate from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV (for population structure)
/// @param n_sim Number of simulations
/// @param seed Random seed
/// @param match_method Propensity-score matching method: "none" (off),
///   "optimal", "nearest", or "rank" (requires observed DV)
/// @return Data frame with ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
fn ferx_rust_simulate(
    model_path: &str,
    data_path: &str,
    n_sim: i32,
    seed: i32,
    match_method: &str,
) -> Robj {
    let match_method = match parse_match_method(match_method) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };
    // parse_full_model_file (vs parse_model_file) so iov_column is available
    // for the reader; without it, models with kappa declarations panic in
    // pk_param_fn (Eta index >= n_bsv_eta).
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) =
        match ferx_core::api::read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            iov_col.as_deref(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return ().into();
            }
        };

    let opts = ferx_core::SimulateOptions {
        seed: Some(seed as u64),
        match_method,
        // No TTE administrative horizon from the R wrapper yet (ferx-core #522);
        // `..Default` fills `horizon: None` and any future fields.
        ..Default::default()
    };
    let results = match ferx_core::simulate_with_options(
        &parsed.model,
        &population,
        &parsed.model.default_params,
        n_sim as usize,
        &opts,
    ) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Error simulating: {}", e);
            return ().into();
        }
    };

    sim_results_to_df(&results)
}

/// Simulate using fitted parameters.
///
/// @param model_path Path to .ferx model file (used for structure + default metadata)
/// @param data_path Path to NONMEM-format CSV (for population structure)
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @param n_sim Number of simulations
/// @param seed Random seed
/// @param match_method Propensity-score matching method: "none" (off),
///   "optimal", "nearest", or "rank" (requires observed DV)
/// @return Data frame with SIM, ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_simulate_from_fit(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
    n_sim: i32,
    seed: i32,
    match_method: &str,
) -> Robj {
    let match_method = match parse_match_method(match_method) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) =
        match ferx_core::api::read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            iov_col.as_deref(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return ().into();
            }
        };

    let params = match params_from_fit(&parsed.model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let opts = ferx_core::SimulateOptions {
        seed: Some(seed as u64),
        match_method,
        // No TTE administrative horizon from the R wrapper yet (ferx-core #522);
        // `..Default` fills `horizon: None` and any future fields.
        ..Default::default()
    };
    let results = match ferx_core::simulate_with_options(
        &parsed.model,
        &population,
        &params,
        n_sim as usize,
        &opts,
    ) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Error simulating: {}", e);
            return ().into();
        }
    };
    sim_results_to_df(&results)
}

/// Simulate with parameter uncertainty propagation.
///
/// For each parameter set drawn from the uncertainty distribution
/// (`method = "asymptotic"` uses the FD covariance matrix; `method = "sir"`
/// reuses SIR resamples) the per-subject eta/eps sampler runs
/// `n_sim_per_draw` times. Returns a long data frame with a leading `DRAW`
/// column so downstream code can compute prediction bands that include
/// parameter uncertainty.
///
/// Both uncertainty sources need data attached to the fit:
/// - "asymptotic" requires `cov_matrix_flat` populated (run `ferx_fit` with
///   `covariance = TRUE`).
/// - "sir" requires `sir_resamples_flat` populated (run with
///   `sir = TRUE` and `sir_keep_samples = TRUE` in `settings`).
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV (provides population structure)
/// @param theta Fitted theta vector — the ML estimate that anchors the proposal
/// @param omega_flat Row-major flattened fitted omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @param method Uncertainty method: `"asymptotic"` or `"sir"`
/// @param cov_matrix_flat Row-major flattened packed-space covariance matrix
///   (empty when not using asymptotic mode)
/// @param cov_matrix_dim Side length of the covariance matrix (0 when unused)
/// @param sir_resamples_flat Flattened SIR resample pool, row-major
///   (`n_resamples * n_packed` values; empty when not using SIR mode)
/// @param sir_resamples_n Number of SIR resamples (rows)
/// @param sir_resamples_dim Packed-parameter dimension (columns)
/// @param n_uncertainty_draws Number of parameter sets to draw from the
///   uncertainty distribution
/// @param n_sim_per_draw Number of eta/eps replicates per parameter draw
/// @param seed Random seed for reproducibility
/// @return Data frame with DRAW, SIM, ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_simulate_with_uncertainty(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
    method: &str,
    cov_matrix_flat: Vec<f64>,
    cov_matrix_dim: i32,
    sir_resamples_flat: Vec<f64>,
    sir_resamples_n: i32,
    sir_resamples_dim: i32,
    n_uncertainty_draws: i32,
    n_sim_per_draw: i32,
    seed: i32,
) -> Robj {
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();
    let (population, _) =
        match ferx_core::api::read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            iov_col.as_deref(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return ().into();
            }
        };

    // Decode the method string to the engine enum.
    let uncertainty_method = match method.trim().to_lowercase().as_str() {
        "asymptotic" | "cov" | "covariance" => {
            ferx_core::UncertaintyMethod::Asymptotic
        }
        "sir" => ferx_core::UncertaintyMethod::Sir,
        other => {
            rprintln!(
                "Unknown uncertainty method '{}' — expected 'asymptotic' or 'sir'",
                other
            );
            return ().into();
        }
    };

    let fit_result =
        match build_fit_result_for_uncertainty(
            &parsed.model,
            &theta,
            &omega_flat,
            omega_dim,
            &sigma,
            uncertainty_method,
            &cov_matrix_flat,
            cov_matrix_dim,
            &sir_resamples_flat,
            sir_resamples_n,
            sir_resamples_dim,
        ) {
            Ok(f) => f,
            Err(e) => {
                rprintln!("{}", e);
                return ().into();
            }
        };

    let opts = ferx_core::SimulateUncertaintyOptions {
        n_uncertainty_draws: n_uncertainty_draws.max(0) as usize,
        n_sim_per_draw: n_sim_per_draw.max(0) as usize,
        method: uncertainty_method,
        seed: Some(seed as u64),
    };

    match ferx_core::simulate_with_uncertainty(&parsed.model, &population, &fit_result, &opts) {
        Ok(results) => sim_results_to_df(&results),
        Err(e) => {
            rprintln!("simulate_with_uncertainty error: {}", e);
            ().into()
        }
    }
}

/// Population predictions from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @return Data frame with ID, TIME, PRED columns
/// @export
#[extendr]
fn ferx_rust_predict(
    model_path: &str,
    data_path: &str,
) -> Robj {
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) =
        match ferx_core::api::read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            iov_col.as_deref(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return ().into();
            }
        };

    let results = ferx_core::predict(&parsed.model, &population, &parsed.model.default_params);

    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let pred: Vec<f64> = results.iter().map(|r| r.pred).collect();

    data_frame!(ID = id, TIME = time, PRED = pred).into()
}

/// Population predictions using fitted parameters.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @return Data frame with ID, TIME, PRED columns
/// @export
#[extendr]
fn ferx_rust_predict_from_fit(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
) -> Robj {
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) =
        match ferx_core::api::read_population_for(
            &parsed.model,
            &parsed.covariate_decls,
            data_path,
            None,
            iov_col.as_deref(),
            None,
        ) {
            Ok(r) => r,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return ().into();
            }
        };

    let params = match params_from_fit(&parsed.model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results = ferx_core::predict(&parsed.model, &population, &params);

    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let pred: Vec<f64> = results.iter().map(|r| r.pred).collect();

    data_frame!(ID = id, TIME = time, PRED = pred).into()
}

/// Flatten survival-function predictions into an R data frame (one row per
/// subject × TTE CMT × time-grid point).
fn survival_results_to_df(results: &[ferx_core::SurvivalPredictionResult]) -> Robj {
    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let cmt: Vec<i32> = results.iter().map(|r| r.cmt as i32).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let survival: Vec<f64> = results.iter().map(|r| r.survival).collect();
    let cum_hazard: Vec<f64> = results.iter().map(|r| r.cum_hazard).collect();
    let hazard: Vec<f64> = results.iter().map(|r| r.hazard).collect();
    // Competing-risks quantities (ferx-core #501): cause-specific cumulative
    // incidence and all-cause survival. For a single TTE endpoint cif == 1 -
    // survival and survival_all == survival; across causes Σ cif + survival_all = 1.
    let cif: Vec<f64> = results.iter().map(|r| r.cif).collect();
    let survival_all: Vec<f64> = results.iter().map(|r| r.survival_all).collect();
    let median_survival: Vec<f64> = results.iter().map(|r| r.median_survival).collect();
    let mean_survival: Vec<f64> = results.iter().map(|r| r.mean_survival).collect();
    data_frame!(
        ID = id,
        CMT = cmt,
        TIME = time,
        survival = survival,
        cum_hazard = cum_hazard,
        hazard = hazard,
        cif = cif,
        survival_all = survival_all,
        median_survival = median_survival,
        mean_survival = mean_survival
    )
    .into()
}

/// Survival-function predictions for TTE endpoints from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param times Time grid at which S(t), H(t), h(t) are evaluated
/// @return Data frame with ID, CMT, TIME, survival, cum_hazard, hazard, cif,
///   survival_all, median_survival, mean_survival
/// @export
#[extendr]
fn ferx_rust_predict_survival(model_path: &str, data_path: &str, times: Vec<f64>) -> Robj {
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) = match ferx_core::api::read_population_for(
        &parsed.model,
        &parsed.covariate_decls,
        data_path,
        None,
        iov_col.as_deref(),
        None,
    ) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let results =
        ferx_core::predict_survival(&parsed.model, &population, &parsed.model.default_params, &times);
    survival_results_to_df(&results)
}

/// Survival-function predictions for TTE endpoints using fitted parameters.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param times Time grid at which S(t), H(t), h(t) are evaluated
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @return Data frame with ID, CMT, TIME, survival, cum_hazard, hazard, cif,
///   survival_all, median_survival, mean_survival
/// @export
#[extendr]
fn ferx_rust_predict_survival_from_fit(
    model_path: &str,
    data_path: &str,
    times: Vec<f64>,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
) -> Robj {
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    let (population, _) = match ferx_core::api::read_population_for(
        &parsed.model,
        &parsed.covariate_decls,
        data_path,
        None,
        iov_col.as_deref(),
        None,
    ) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&parsed.model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results = ferx_core::predict_survival(&parsed.model, &population, &params, &times);
    survival_results_to_df(&results)
}

/// Simulation-based NPDE / NPD diagnostics from fitted parameters.
///
/// Recomputes the Normalized Prediction Distribution Errors (NPDE, decorrelated
/// within subject) and Normalized Prediction Discrepancies (NPD) post-hoc, so a
/// fit produced without `npde_nsim` can still get them without re-fitting.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param theta Fitted theta vector
/// @param omega_flat Row-major flattened omega matrix
/// @param omega_dim Side length of the omega matrix
/// @param sigma Fitted sigma vector
/// @param nsim Number of Monte-Carlo replicates per subject
/// @param seed RNG seed; pass -1 for the engine default
/// @return Data frame with ID, TIME, NPDE, NPD columns
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_npde_from_fit(
    model_path: &str,
    data_path: &str,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
    nsim: i32,
    seed: i32,
) -> Robj {
    if nsim <= 0 {
        rprintln!("npde error: nsim must be a positive integer");
        return ().into();
    }

    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };
    let iov_col = parsed.fit_options.iov_column.clone();

    // Re-apply the model file's `[data_selection]` ignore/accept/ignore_subjects
    // so the population matches the one the fit was computed on. Without this the
    // NPDE decorrelation would run over a different per-subject observation set
    // than the fit, silently disagreeing with a fit-time `npde_nsim` run. (Only
    // the model-file selection is visible here; selection applied via R-side
    // `ferx_fit(settings=)` is not carried on the fit object.)
    let filter = match ferx_core::io::datareader::SelectionFilter::from_opts(
        &parsed.fit_options.ignore_exprs,
        &parsed.fit_options.accept_exprs,
        &parsed.fit_options.ignore_subjects,
    ) {
        Ok(f) => f,
        Err(e) => {
            rprintln!("Error in [data_selection]: {}", e);
            return ().into();
        }
    };
    let filter_opt = if filter.is_empty() { None } else { Some(&filter) };

    let (population, _) = match ferx_core::api::read_population_for(
        &parsed.model,
        &parsed.covariate_decls,
        data_path,
        None,
        iov_col.as_deref(),
        filter_opt,
    ) {
        Ok(r) => r,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&parsed.model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let seed_opt = if seed < 0 { None } else { Some(seed as u64) };
    let per_subject = ferx_core::stats::npde::compute_npde_npd(
        &parsed.model,
        &population,
        &params,
        nsim as usize,
        seed_opt,
    );

    // Flatten per-subject NPDE/NPD back to one row per observation. ID and TIME
    // are emitted exactly as `io::output::sdtab` builds them — numeric ID
    // (`id.parse::<f64>()`, falling back to the 1-based subject index) and the
    // raw data TIME — so the R side can align this table to `fit$sdtab`
    // positionally and assert per-row ID/TIME agreement rather than re-joining.
    let mut id: Vec<f64> = Vec::new();
    let mut time: Vec<f64> = Vec::new();
    let mut npde: Vec<f64> = Vec::new();
    let mut npd: Vec<f64> = Vec::new();
    for (si, (subj, sn)) in population.subjects.iter().zip(per_subject.iter()).enumerate() {
        let id_num = subj.id.parse::<f64>().unwrap_or(si as f64 + 1.0);
        for j in 0..subj.observations.len() {
            id.push(id_num);
            time.push(
                subj.obs_raw_times
                    .get(j)
                    .copied()
                    .unwrap_or(subj.obs_times[j]),
            );
            npde.push(sn.npde.get(j).copied().unwrap_or(f64::NAN));
            npd.push(sn.npd.get(j).copied().unwrap_or(f64::NAN));
        }
    }

    data_frame!(ID = id, TIME = time, NPDE = npde, NPD = npd).into()
}

// -- Helper: parse a single R-side propensity-matching token into an
//    Option<MatchMethod> ("none" → None disables matching). --

fn parse_match_method(
    token: &str,
) -> std::result::Result<Option<ferx_core::MatchMethod>, String> {
    match token.trim().to_lowercase().as_str() {
        "none" | "false" | "" => Ok(None),
        "optimal" | "true" => Ok(Some(ferx_core::MatchMethod::Optimal)),
        "nearest" => Ok(Some(ferx_core::MatchMethod::Nearest)),
        "rank" => Ok(Some(ferx_core::MatchMethod::Rank)),
        other => Err(format!(
            "Unknown match method '{}' — expected one of: none, optimal, nearest, rank",
            other
        )),
    }
}

// -- Helper: parse a single R-side method token into EstimationMethod --

fn parse_method(token: &str) -> std::result::Result<EstimationMethod, String> {
    let m = token.trim().to_lowercase();
    if m == "saem" {
        Ok(EstimationMethod::Saem)
    } else if m.contains("hybrid") {
        Ok(EstimationMethod::FoceGnHybrid)
    } else if m == "gn" || m.contains("gauss") {
        Ok(EstimationMethod::FoceGn)
    } else if m == "focei" || m == "foce-i" || m == "foce_i" || m.contains("interaction") {
        Ok(EstimationMethod::FoceI)
    } else if m == "foce" {
        Ok(EstimationMethod::Foce)
    } else if m == "impmap" || m == "importance_sampling_map" || m == "importance-sampling-map" {
        Ok(EstimationMethod::Impmap)
    } else if m == "imp" || m == "importance_sampling" || m == "importance-sampling" {
        Ok(EstimationMethod::Imp)
    } else if m == "impmap" || m == "importance_sampling_map" {
        Ok(EstimationMethod::Impmap)
    } else if m == "bayes" || m == "bayesian" || m == "mcmc" {
        Ok(EstimationMethod::Bayes)
    } else {
        Err(format!(
            "Unknown estimation method '{}' — expected one of: foce, focei, saem, gn, gn_hybrid, imp, impmap, bayes",
            token.trim()
        ))
    }
}

// -- Helper: materialize ModelParameters from R-side theta/omega/sigma --

fn params_from_fit(
    model: &CompiledModel,
    theta: &[f64],
    omega_flat: &[f64],
    omega_dim: i32,
    sigma: &[f64],
) -> std::result::Result<ModelParameters, String> {
    let template = &model.default_params;

    if theta.len() != template.theta.len() {
        return Err(format!(
            "Fit error: theta length {} does not match model ({} expected)",
            theta.len(),
            template.theta.len()
        ));
    }
    if sigma.len() != template.sigma.values.len() {
        return Err(format!(
            "Fit error: sigma length {} does not match model ({} expected)",
            sigma.len(),
            template.sigma.values.len()
        ));
    }
    let d = omega_dim as usize;
    if d != template.omega.dim() {
        return Err(format!(
            "Fit error: omega dim {} does not match model ({} expected)",
            d,
            template.omega.dim()
        ));
    }
    if omega_flat.len() != d * d {
        return Err(format!(
            "Fit error: omega_flat length {} does not match dim²={}",
            omega_flat.len(),
            d * d
        ));
    }

    // R/Rust both use row-major-from-column-major interchangeably here since
    // the matrix is symmetric in the intended use (covariance); reconstruct
    // via DMatrix::from_row_slice to be explicit.
    let omega_mat = DMatrix::from_row_slice(d, d, omega_flat);
    let omega = OmegaMatrix::from_matrix(
        omega_mat,
        template.omega.eta_names.clone(),
        template.omega.diagonal,
    );

    Ok(ModelParameters {
        theta: theta.to_vec(),
        theta_names: template.theta_names.clone(),
        theta_lower: template.theta_lower.clone(),
        theta_upper: template.theta_upper.clone(),
        theta_fixed: template.theta_fixed.clone(),
        omega,
        omega_fixed: template.omega_fixed.clone(),
        sigma: SigmaVector {
            values: sigma.to_vec(),
            names: template.sigma.names.clone(),
        },
        sigma_fixed: template.sigma_fixed.clone(),
        omega_iov: None,
        kappa_fixed: Vec::new(),
    })
}

// -- Helper: build a minimal FitResult from R-supplied arrays so the
//    uncertainty path can be exercised without re-running fit() from R.
//    Only the fields that simulate_with_uncertainty / fitted_params_from_result
//    actually read are populated with real values; the rest get sensible
//    defaults that won't be inspected. --
#[allow(clippy::too_many_arguments)]
fn build_fit_result_for_uncertainty(
    model: &CompiledModel,
    theta: &[f64],
    omega_flat: &[f64],
    omega_dim: i32,
    sigma: &[f64],
    method: ferx_core::UncertaintyMethod,
    cov_matrix_flat: &[f64],
    cov_matrix_dim: i32,
    sir_resamples_flat: &[f64],
    sir_resamples_n: i32,
    sir_resamples_dim: i32,
) -> std::result::Result<FitResult, String> {
    let template = &model.default_params;
    if theta.len() != template.theta.len() {
        return Err(format!(
            "Uncertainty error: theta length {} does not match model ({} expected)",
            theta.len(),
            template.theta.len()
        ));
    }
    if sigma.len() != template.sigma.values.len() {
        return Err(format!(
            "Uncertainty error: sigma length {} does not match model ({} expected)",
            sigma.len(),
            template.sigma.values.len()
        ));
    }
    let d = omega_dim as usize;
    if d != template.omega.dim() {
        return Err(format!(
            "Uncertainty error: omega dim {} does not match model ({} expected)",
            d,
            template.omega.dim()
        ));
    }
    if omega_flat.len() != d * d {
        return Err(format!(
            "Uncertainty error: omega_flat length {} does not match dim²={}",
            omega_flat.len(),
            d * d
        ));
    }
    let omega_mat = DMatrix::from_row_slice(d, d, omega_flat);

    let (covariance_matrix, sir_resamples_packed) = match method {
        ferx_core::UncertaintyMethod::Asymptotic => {
            let cd = cov_matrix_dim as usize;
            if cd == 0 || cov_matrix_flat.is_empty() {
                return Err("Asymptotic uncertainty requires a covariance matrix on the \
                    fit object — re-fit with `covariance = TRUE`."
                    .to_string());
            }
            if cov_matrix_flat.len() != cd * cd {
                return Err(format!(
                    "Uncertainty error: cov_matrix_flat length {} does not match dim²={}",
                    cov_matrix_flat.len(),
                    cd * cd
                ));
            }
            let cov = DMatrix::from_row_slice(cd, cd, cov_matrix_flat);
            (Some(cov), None)
        }
        ferx_core::UncertaintyMethod::Sir => {
            let n = sir_resamples_n as usize;
            let pd = sir_resamples_dim as usize;
            if n == 0 || pd == 0 || sir_resamples_flat.is_empty() {
                return Err("SIR uncertainty requires resamples on the fit object — \
                    re-fit with `sir = TRUE` and `sir_keep_samples = TRUE` in `settings`."
                    .to_string());
            }
            if sir_resamples_flat.len() != n * pd {
                return Err(format!(
                    "Uncertainty error: sir_resamples_flat length {} does not match \
                     n_resamples * packed_dim = {} * {} = {}",
                    sir_resamples_flat.len(),
                    n,
                    pd,
                    n * pd
                ));
            }
            let mut pool: Vec<Vec<f64>> = Vec::with_capacity(n);
            for i in 0..n {
                pool.push(sir_resamples_flat[i * pd..(i + 1) * pd].to_vec());
            }
            (None, Some(pool))
        }
    };

    Ok(default_fit_result(
        model,
        theta.to_vec(),
        omega_mat,
        sigma.to_vec(),
        covariance_matrix,
        sir_resamples_packed,
    ))
}

// Build a defaulted FitResult populated with the fields that
// fitted_params_from_result / simulate_with_uncertainty actually read. Other
// fields get neutral defaults — none of them are inspected on this path.
fn default_fit_result(
    model: &CompiledModel,
    theta: Vec<f64>,
    omega: DMatrix<f64>,
    sigma: Vec<f64>,
    covariance_matrix: Option<DMatrix<f64>>,
    sir_resamples_packed: Option<Vec<Vec<f64>>>,
) -> FitResult {
    let template = &model.default_params;
    FitResult {
        method: EstimationMethod::FoceI,
        method_chain: vec![EstimationMethod::FoceI],
        bayes: None,
        cond_dist: None,
        converged: true,
        ofv: 0.0,
        aic: 0.0,
        bic: 0.0,
        theta,
        theta_names: model.theta_names.clone(),
        eta_names: template.omega.eta_names.clone(),
        omega,
        sigma,
        sigma_names: template.sigma.names.clone(),
        error_model: model.error_model,
        covariance_matrix,
        se_theta: None,
        se_omega: None,
        se_sigma: None,
        theta_fixed: template.theta_fixed.clone(),
        omega_fixed: template.omega_fixed.clone(),
        sigma_fixed: template.sigma_fixed.clone(),
        subjects: Vec::new(),
        n_obs: 0,
        n_subjects: 0,
        n_parameters: 0,
        n_iterations: 0,
        interaction: true,
        warnings: Vec::new(),
        sir_ci_theta: None,
        sir_ci_omega: None,
        sir_ci_sigma: None,
        sir_ess: None,
        sir_resamples_packed,
        importance_sampling: None,
        impmap_trace: None,
        omega_iov: None,
        kappa_names: model.kappa_names.clone(),
        kappa_fixed: template.kappa_fixed.clone(),
        se_kappa: None,
        shrinkage_kappa: Vec::new(),
        shrinkage_kappa_by_occ: Vec::new(),
        ebe_kappas: Vec::new(),
        saem_mu_ref_m_step_evals_saved: None,
        saem_n_subjects_hmc: None,
        gradient_method_inner: String::new(),
        gradient_method_outer: String::new(),
        uses_ode_solver: model.is_ode_based(),
        n_threads_used: 1,
        nlopt_missing_algorithms: Vec::new(),
        covariance_n_evals_estimated: None,
        trace_path: None,
        ebe_convergence_warnings: 0,
        max_unconverged_subjects: 0,
        total_ebe_fallbacks: 0,
        covariance_status: CovarianceStatus::Computed,
        shrinkage_eta: Vec::new(),
        shrinkage_eps: f64::NAN,
        wall_time_secs: 0.0,
        model_name: model.name.clone(),
        ferx_version: String::new(),
        eta_param_info: Vec::new(),
        theta_transform: Vec::new(),
        sigma_types: Vec::new(),
        cov_eigenvalues: None,
        cov_condition_number: None,
        eta_log_transformed: Vec::new(),
        omega_param_corr: None,
        omega_iov_param_corr: None,
        model_path: None,
        data_path: None,
        model_hash: None,
        data_hash: None,
        dw_statistic: f64::NAN,
        iwres_lag1_r: f64::NAN,
        uses_sde: false,
        omega_init_as_sd: Vec::new(),
        sigma_init_as_sd: Vec::new(),
        kappa_init_as_sd: Vec::new(),
        warnings_structured: Vec::new(),
        model_text: None,
        theta_init: Vec::new(),
        omega_init: nalgebra::DMatrix::zeros(0, 0),
        sigma_init: Vec::new(),
        obs_time_range: None,
        final_gradient: None,
        optimizer: "auto".to_string(),
        n_starts: 1,
        multi_start_seed: None,
        saem_seed: None,
        sir_seed: None,
        imp_seed: None,
        npde_seed: None,
        bloq_method: "drop".to_string(),
        outer_maxiter: 0,
        outer_gtol: 0.0,
        inits_from_nca: None,
        covariate_names: Vec::new(),
        input_columns: Vec::new(),
        covariate_table: None,
        exclusions: None,
        #[cfg(feature = "nn")]
        neural_networks: Vec::new(),
    }
}

// -- Helper: SimulationResult slice → R data frame --

fn sim_results_to_df(results: &[ferx_core::api::SimulationResult]) -> Robj {
    let draw: Vec<i32> = results.iter().map(|r| r.draw as i32).collect();
    let sim: Vec<i32> = results.iter().map(|r| r.sim as i32).collect();
    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let ipred: Vec<f64> = results.iter().map(|r| r.ipred).collect();
    // `SimulationResult` replaced the flat `dv_sim: f64` with `outcome: SimOutcome`
    // (Gaussian vs TTE). `continuous_value()` returns the Gaussian value (NAN for
    // non-Gaussian outcomes), preserving the DV_SIM column for the existing paths.
    let dv_sim: Vec<f64> = results.iter().map(|r| r.outcome.continuous_value()).collect();

    // DRAW is leading because for legacy `simulate*` paths every row is
    // DRAW = 1 and downstream code can ignore the column; for
    // `simulate_with_uncertainty` it's the primary grouping variable for
    // computing uncertainty bands.
    data_frame!(
        DRAW = draw,
        SIM = sim,
        ID = id,
        TIME = time,
        IPRED = ipred,
        DV_SIM = dv_sim
    )
    .into()
}

// -- Helper: short label for the structural model. ODE models always report
//    "ODE" regardless of the placeholder `pk_model` value the parser stores
//    for them, so `is_ode_based()` is checked first. The PK match is
//    exhaustive over `PkModel`, so this never needs to fall back to NULL on
//    the R side; the corresponding R doc still allows `NULL` because the
//    pre-fit regex parser in `.ferx_model_type()` can fail to recognise the
//    structural form. --
fn pk_model_type_label(model: &CompiledModel) -> &'static str {
    if model.is_ode_based() {
        return "ODE";
    }
    // ferx-core #176 collapsed the split `*_iv_bolus` / `*_infusion`
    // variants into a single `*_iv` per compartment count — the bolus
    // vs infusion route is now read per dose from the RATE column.
    // The model declaration alone no longer pins one or the other, so
    // the label drops the bolus/infusion suffix.
    match model.pk_model {
        PkModel::OneCptIv => "1-cpt IV",
        PkModel::OneCptOral => "1-cpt oral",
        PkModel::TwoCptIv => "2-cpt IV",
        PkModel::TwoCptOral => "2-cpt oral",
        PkModel::ThreeCptIv => "3-cpt IV",
        PkModel::ThreeCptOral => "3-cpt oral",
    }
}

fn error_model_to_string(em: ErrorModel) -> &'static str {
    match em {
        ErrorModel::Additive => "additive",
        ErrorModel::Proportional => "proportional",
        ErrorModel::Combined => "combined",
    }
}

// Residual-error label for `model_structure`. Single-endpoint models report
// the bare error type ("proportional"); multi-endpoint (per-CMT) models report
// each compartment's error type, e.g. "per-CMT (CMT2=proportional,
// CMT3=additive)". Log-transform-both-sides (LTBS) models report
// "additive (log-transformed)", matching the R-side `.ferx_parse_structure`
// label so the fit object's `model_structure$residual` survives a `.fitrx`
// round-trip with the LTBS distinction intact.
fn residual_label(model: &CompiledModel) -> String {
    if model.log_transform {
        // LTBS (`log(DV) ~ additive` / `DV ~ log_additive`): additive error on
        // the log scale. The engine rejects LTBS with per-CMT, so this is always
        // the single-endpoint additive case.
        return "additive (log-transformed)".to_string();
    }
    match &model.error_spec {
        ErrorSpec::Single(em) => error_model_to_string(*em).to_string(),
        ErrorSpec::PerCmt(map) => {
            let mut cmts: Vec<usize> = map.keys().copied().collect();
            cmts.sort_unstable();
            let parts: Vec<String> = cmts
                .iter()
                .map(|c| format!("CMT{}={}", c, error_model_to_string(map[c].error_model)))
                .collect();
            format!("per-CMT ({})", parts.join(", "))
        }
    }
}

// -- Helper: build the model_structure sub-list from the parsed/compiled
//    model. Mirrors the shape that the R layer historically derived by
//    re-parsing the .ferx file (theta_names, model_type, iiv, iov, residual)
//    so it can be a drop-in canonical source for `ferx_model_inspect()`. --
fn model_structure_list(model: &CompiledModel) -> Robj {
    list!(
        theta_names = model.theta_names.clone(),
        model_type = pk_model_type_label(model).to_string(),
        iiv = model.eta_names.clone(),
        iov = model.kappa_names.clone(),
        residual = residual_label(model)
    )
    .into()
}

// -- Helper: convert FitResult + Population to R named list --

fn fit_result_to_list(
    result: &FitResult,
    population: &Population,
    model: &CompiledModel,
) -> List {
    // Theta
    let theta_names: Vec<String> = result.theta_names.clone();
    let theta_values: Vec<f64> = result.theta.clone();

    // Omega as flat matrix (row-major)
    let n_eta = result.omega.nrows();
    let mut omega_flat: Vec<f64> = Vec::with_capacity(n_eta * n_eta);
    for i in 0..n_eta {
        for j in 0..n_eta {
            omega_flat.push(result.omega[(i, j)]);
        }
    }

    // SE vectors
    let se_theta: Vec<f64> = result.se_theta.clone().unwrap_or_default();
    let se_omega: Vec<f64> = result.se_omega.clone().unwrap_or_default();
    let se_sigma: Vec<f64> = result.se_sigma.clone().unwrap_or_default();

    // SDTAB as data frame
    let sdtab_cols = ferx_core::io::output::sdtab(result, population);
    let sdtab = sdtab_to_dataframe(&sdtab_cols);

    // Covariate table (NULL unless the model declared a `[covariates]` block).
    let covtab = covtab_to_dataframe(result.covariate_table.as_ref());
    // Per-covariate type tag ("continuous"/"categorical") as a named character
    // vector; NULL when no `[covariates]` block. Lets R treat categoricals as
    // factors etc. (the value carried by the typed declaration).
    let covariate_types = covariate_types_robj(result.covariate_table.as_ref());

    // Warnings
    let warnings: Vec<String> = result.warnings.clone();
    // Structured warnings as four parallel vectors (severity / category /
    // message / source_method), following the same parallel-array convention
    // as sigma / sigma_names / sigma_types. Severity is lowercased so the R
    // side can match on "critical" / "warning" / "info" directly.
    // source_method is empty string when None (no chain prefix was stripped).
    let warnings_severity: Vec<String> = result
        .warnings_structured
        .iter()
        .map(|w| format!("{:?}", w.severity).to_lowercase())
        .collect();
    let warnings_category: Vec<String> = result
        .warnings_structured
        .iter()
        .map(|w| w.category.clone())
        .collect();
    let warnings_message: Vec<String> = result
        .warnings_structured
        .iter()
        .map(|w| w.message.clone())
        .collect();
    let warnings_source_method: Vec<String> = result
        .warnings_structured
        .iter()
        .map(|w| w.source_method.clone().unwrap_or_default())
        .collect();

    let method_label = result.method.label();
    let method_chain: Vec<String> = result
        .method_chain
        .iter()
        .map(|m| m.label().to_string())
        .collect();

    // SIR CIs flattened as [lo1, hi1, lo2, hi2, ...]; empty vector => not computed.
    let flatten_ci = |ci: &Option<Vec<(f64, f64)>>| -> Vec<f64> {
        ci.as_ref()
            .map(|v| v.iter().flat_map(|(lo, hi)| [*lo, *hi]).collect())
            .unwrap_or_default()
    };
    let sir_ess: f64 = result.sir_ess.unwrap_or(f64::NAN);
    let sir_ci_theta = flatten_ci(&result.sir_ci_theta);
    let sir_ci_omega = flatten_ci(&result.sir_ci_omega);
    let sir_ci_sigma = flatten_ci(&result.sir_ci_sigma);

    // SIR resamples flattened row-major as a length-(n_resamples * n_packed)
    // vector; (sir_resamples_n, sir_resamples_dim) lets R re-shape on read.
    // Empty when sir_keep_samples was off (or SIR didn't run).
    let (sir_resamples_flat, sir_resamples_n, sir_resamples_dim): (Vec<f64>, i32, i32) =
        match &result.sir_resamples_packed {
            Some(pool) if !pool.is_empty() => {
                let n = pool.len();
                let d = pool[0].len();
                let mut v = Vec::with_capacity(n * d);
                for row in pool {
                    v.extend_from_slice(row);
                }
                (v, n as i32, d as i32)
            }
            _ => (Vec::new(), 0i32, 0i32),
        };

    // Importance Sampling (IMP) diagnostics. Surfaced as a nested named list
    // when the chain ended with an `imp` stage; NULL otherwise. The flat
    // `(ids, ess_frac)` pair lets the R side rebuild a data frame without
    // hitting extendr's data-frame-in-list limitations.
    let importance_sampling: Robj = match &result.importance_sampling {
        Some(is) => {
            let kappa_treatment_str: &'static str = match is.kappa_treatment {
                KappaTreatment::NotApplicable => "not_applicable",
                KappaTreatment::FixedAtMode => "fixed_at_mode",
                KappaTreatment::Marginalized => "marginalized",
            };
            let low_ess_ids: Vec<String> =
                is.low_ess_subjects.iter().map(|(id, _)| id.clone()).collect();
            let low_ess_frac: Vec<f64> =
                is.low_ess_subjects.iter().map(|(_, f)| *f).collect();
            list!(
                minus2_log_likelihood = is.minus2_log_likelihood,
                mc_standard_error = is.mc_standard_error,
                // n_samples is usize in ferx-core; cast to f64 (rather than
                // i32) so a user-set `imp_samples > i32::MAX` can't silently
                // wrap to a negative number when bridged to R.
                n_samples = is.n_samples as f64,
                proposal_df = is.proposal_df,
                ess_min = is.ess_min,
                ess_median = is.ess_median,
                kappa_treatment = kappa_treatment_str,
                low_ess_subject_ids = low_ess_ids,
                low_ess_subject_frac = low_ess_frac,
            )
            .into()
        }
        None => ().into(),
    };

    // IMPMAP per-iteration parameter trace — flat vectors that the R side
    // reassembles into a data.frame.  NULL when `impmap_trace = false` or a
    // non-IMPMAP method was used.
    let impmap_trace: Robj = match &result.impmap_trace {
        Some(trace) => {
            let n = trace.rows.len();
            let iterations: Vec<f64> = trace.rows.iter().map(|r| r.iteration as f64).collect();
            let ofvs: Vec<f64> = trace.rows.iter().map(|r| r.ofv).collect();

            // Per-theta columns.
            let n_theta = if n > 0 { trace.rows[0].theta.len() } else { 0 };
            let mut theta_cols: Vec<Vec<f64>> = vec![Vec::with_capacity(n); n_theta];
            for row in &trace.rows {
                for (j, v) in row.theta.iter().enumerate() {
                    if j < n_theta {
                        theta_cols[j].push(*v);
                    }
                }
            }

            // Per-omega-LT columns.
            let n_omega = if n > 0 { trace.rows[0].omega_lower_tri.len() } else { 0 };
            let mut omega_cols: Vec<Vec<f64>> = vec![Vec::with_capacity(n); n_omega];
            for row in &trace.rows {
                for (j, v) in row.omega_lower_tri.iter().enumerate() {
                    if j < n_omega {
                        omega_cols[j].push(*v);
                    }
                }
            }

            // Per-sigma columns.
            let n_sigma = if n > 0 { trace.rows[0].sigma.len() } else { 0 };
            let mut sigma_cols: Vec<Vec<f64>> = vec![Vec::with_capacity(n); n_sigma];
            for row in &trace.rows {
                for (j, v) in row.sigma.iter().enumerate() {
                    if j < n_sigma {
                        sigma_cols[j].push(*v);
                    }
                }
            }

            // Flatten column data + names for the R side to unpack.
            let mut col_data: Vec<Vec<f64>> = Vec::new();
            let mut col_names: Vec<String> = Vec::new();
            col_names.push("ITERATION".to_string());
            col_data.push(iterations);
            for (j, col) in theta_cols.into_iter().enumerate() {
                col_names.push(trace.theta_names.get(j).cloned().unwrap_or_else(|| format!("THETA{}", j + 1)));
                col_data.push(col);
            }
            for (j, col) in omega_cols.into_iter().enumerate() {
                col_names.push(trace.omega_names.get(j).cloned().unwrap_or_else(|| format!("OMEGA{}", j + 1)));
                col_data.push(col);
            }
            for (j, col) in sigma_cols.into_iter().enumerate() {
                col_names.push(trace.sigma_names.get(j).cloned().unwrap_or_else(|| format!("SIGMA{}", j + 1)));
                col_data.push(col);
            }
            col_names.push("OBJ".to_string());
            col_data.push(ofvs);

            // Flatten all columns into one Vec<f64> + dimension metadata so the
            // R side can matrix(flat, nrow=n, ncol=n_cols, byrow=FALSE).
            let n_cols = col_data.len() as i32;
            let n_rows = n as i32;
            let flat: Vec<f64> = col_data.into_iter().flatten().collect();

            list!(
                flat = flat,
                n_rows = n_rows,
                n_cols = n_cols,
                col_names = col_names,
            )
            .into()
        }
        None => ().into(),
    };

    // Bayes posterior summary (`method = bayes`). Per-parameter summaries are
    // packed as parallel vectors so the R side can build a data frame without
    // hitting extendr's list-of-records limitations.
    let bayes: Robj = match &result.bayes {
        Some(b) => {
            let names: Vec<String> = b.summaries.iter().map(|s| s.name.clone()).collect();
            let mean: Vec<f64> = b.summaries.iter().map(|s| s.mean).collect();
            let sd: Vec<f64> = b.summaries.iter().map(|s| s.sd).collect();
            let q025: Vec<f64> = b.summaries.iter().map(|s| s.q025).collect();
            let median: Vec<f64> = b.summaries.iter().map(|s| s.median).collect();
            let q975: Vec<f64> = b.summaries.iter().map(|s| s.q975).collect();
            let rhat: Vec<f64> = b.summaries.iter().map(|s| s.rhat).collect();
            let ess_bulk: Vec<f64> = b.summaries.iter().map(|s| s.ess_bulk).collect();
            let ess_tail: Vec<f64> = b.summaries.iter().map(|s| s.ess_tail).collect();
            let mcse: Vec<f64> = b.summaries.iter().map(|s| s.mcse).collect();
            list!(
                n_chains = b.n_chains as f64,
                n_warmup = b.n_warmup as f64,
                n_draws_per_chain = b.n_draws_per_chain as f64,
                n_divergent = b.n_divergent as f64,
                max_rhat = b.max_rhat,
                param_names = names,
                mean = mean,
                sd = sd,
                q025 = q025,
                median = median,
                q975 = q975,
                rhat = rhat,
                ess_bulk = ess_bulk,
                ess_tail = ess_tail,
                mcse = mcse,
            )
            .into()
        }
        None => ().into(),
    };

    let trace_path: Robj = match &result.trace_path {
        Some(p) => p.clone().into(),
        None => ().into(),
    };

    // Full parameter covariance matrix (row-major flat); empty when not computed
    let (cov_matrix_flat, cov_matrix_dim): (Vec<f64>, i32) =
        match &result.covariance_matrix {
            Some(m) => {
                let n = m.nrows();
                let mut v = Vec::with_capacity(n * n);
                for i in 0..n {
                    for j in 0..n {
                        v.push(m[(i, j)]);
                    }
                }
                (v, n as i32)
            }
            None => (Vec::new(), 0i32),
        };

    let covariance_status_str = match result.covariance_status {
        CovarianceStatus::Computed => "computed",
        CovarianceStatus::Failed => "failed",
        CovarianceStatus::NotRequested => "not_requested",
        // SirFallback: FD Hessian was non-PD and `covariance_fallback = sir`
        // produced SIR-based intervals (ferx-core #245). The R side maps this
        // token to the "SirFallback" label (see persist.R).
        CovarianceStatus::SirFallback => "sir_fallback",
    };

    // Parameter-level correlation for block omega (row-major flat); empty when diagonal or absent
    let omega_param_corr_flat: Vec<f64> = match &result.omega_param_corr {
        Some(m) => {
            let n = m.nrows();
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n {
                for j in 0..n {
                    v.push(m[(i, j)]);
                }
            }
            v
        }
        None => Vec::new(),
    };
    // Parameter-level correlation for block kappa IOV (row-major flat); empty when diagonal or absent
    let omega_iov_param_corr_flat: Vec<f64> = match &result.omega_iov_param_corr {
        Some(m) => {
            let n = m.nrows();
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n {
                for j in 0..n {
                    v.push(m[(i, j)]);
                }
            }
            v
        }
        None => Vec::new(),
    };
    // Whether each BSV eta is lognormal (true) or additive/unknown (false)
    let eta_log_transformed: Vec<bool> = result.eta_log_transformed.clone();

    // IOV kappa omega (row-major flat); empty when no IOV
    let (omega_iov_flat, omega_iov_dim): (Vec<f64>, i32) = match &result.omega_iov {
        Some(m) => {
            let n = m.nrows();
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n {
                for j in 0..n {
                    v.push(m[(i, j)]);
                }
            }
            (v, n as i32)
        }
        None => (Vec::new(), 0i32),
    };
    let kappa_names: Vec<String> = result.kappa_names.clone();
    let se_kappa: Vec<f64> = result.se_kappa.clone().unwrap_or_default();
    let shrinkage_kappa: Vec<f64> = result.shrinkage_kappa.clone();

    // Per-occasion kappa shrinkage (ferx-core PR #167). Empty when no IOV or
    // when the design is unbalanced and per-occasion values are unreliable.
    // Returned as a data frame: occ (1-based) + one column per kappa name.
    let shrinkage_kappa_by_occ_df: Robj = {
        let by_occ = &result.shrinkage_kappa_by_occ;
        if by_occ.is_empty() || kappa_names.is_empty() {
            ().into()
        } else {
            let n_occ = by_occ.len();
            let n_kappa = kappa_names.len();
            let occ_col: Vec<i32> = (1..=n_occ as i32).collect();
            let mut kappa_cols: Vec<Vec<f64>> = vec![Vec::with_capacity(n_occ); n_kappa];
            for occ_sh in by_occ {
                for k in 0..n_kappa {
                    kappa_cols[k].push(if k < occ_sh.len() { occ_sh[k] } else { f64::NAN });
                }
            }
            let mut pairs: Vec<(&str, Robj)> = Vec::with_capacity(1 + n_kappa);
            pairs.push(("occ", occ_col.into()));
            for (k, name) in kappa_names.iter().enumerate() {
                pairs.push((name.as_str(), kappa_cols[k].clone().into()));
            }
            let mut df = List::from_pairs(pairs);
            df.set_class(&["data.frame"]).unwrap();
            let row_names: Vec<i32> = (1..=n_occ as i32).collect();
            df.set_attrib("row.names", row_names).unwrap();
            df.into()
        }
    };

    // Parameter transform metadata (added in ferx-core PR #54; empty vecs for
    // older binaries that don't populate these fields).
    let eta_param_types: Vec<String> = result.eta_param_info.iter().map(|info| {
        match info.param_type {
            ferx_core::types::EtaParamType::LogNormal        => "log_normal",
            ferx_core::types::EtaParamType::Additive         => "additive",
            ferx_core::types::EtaParamType::Logit            => "logit",
            ferx_core::types::EtaParamType::LogitProbability => "logit_probability",
            ferx_core::types::EtaParamType::Custom           => "custom",
        }.to_string()
    }).collect();
    // Linked theta name for each ETA (empty string when not detected).
    let eta_linked_theta: Vec<String> = result.eta_param_info.iter()
        .map(|info| info.linked_theta.clone().unwrap_or_default())
        .collect();
    let theta_transforms: Vec<String> = result.theta_transform.iter().map(|t| {
        match t {
            ferx_core::types::ThetaTransform::Identity         => "identity",
            ferx_core::types::ThetaTransform::Log              => "log",
            ferx_core::types::ThetaTransform::Logit            => "logit",
            ferx_core::types::ThetaTransform::LogitProbability => "logit_probability",
        }.to_string()
    }).collect();
    // EBE kappas as a data frame: ID, OCC, KAPPA_1, ...
    // Rows: one per (subject, occasion); empty when no IOV.
    // - ID: original subject identifier (matches sdtab).
    // - OCC: labeled occasion in first-seen order, matching the kappa[k]
    //   indexing used by the inner optimizer (see split_obs_by_occasion in
    //   ferx-core). Falls back to a 1-based positional index for any subject
    //   missing an occasion column entry.
    let ebe_kappas_df: Robj = if result.ebe_kappas.is_empty() || kappa_names.is_empty() {
        ().into()
    } else {
        let n_kappa = kappa_names.len();
        let mut ids: Vec<String> = Vec::new();
        let mut occs: Vec<i32> = Vec::new();
        let mut kappa_cols: Vec<Vec<f64>> = vec![Vec::new(); n_kappa];
        for (si, subj_kappas) in result.ebe_kappas.iter().enumerate() {
            let subj = &population.subjects[si];
            // Unique occasion labels in first-seen order (matches kappa[k]).
            let mut unique_occs: Vec<u32> = Vec::new();
            for &occ in &subj.occasions {
                if !unique_occs.contains(&occ) {
                    unique_occs.push(occ);
                }
            }
            for (oi, kappa_vec) in subj_kappas.iter().enumerate() {
                ids.push(subj.id.clone());
                let occ_label = unique_occs.get(oi).copied().unwrap_or(oi as u32 + 1);
                occs.push(occ_label as i32);
                for k in 0..n_kappa {
                    kappa_cols[k].push(if k < kappa_vec.len() { kappa_vec[k] } else { f64::NAN });
                }
            }
        }
        let n_rows = ids.len();
        let mut pairs: Vec<(&str, Robj)> = Vec::new();
        pairs.push(("ID", ids.into()));
        pairs.push(("OCC", occs.into()));
        for (k, name) in kappa_names.iter().enumerate() {
            pairs.push((name.as_str(), kappa_cols[k].clone().into()));
        }
        let mut df = List::from_pairs(pairs);
        df.set_class(&["data.frame"]).unwrap();
        let row_names: Vec<i32> = (1..=n_rows as i32).collect();
        df.set_attrib("row.names", row_names).unwrap();
        df.into()
    };

    // Per-subject EBE etas: ID + one column per BSV eta (named).
    let ebe_etas_df: Robj = build_ebe_etas(result, population);

    // Per-subject individual parameter estimates: ID + one column per
    // [individual_parameters] declaration. Computed by evaluating
    // `pk_param_fn` at the subject's BSV eta + zero kappas + covariates, so
    // each value is the subject's typical (kappa-free) parameter under the
    // model's covariate effects. Per-occasion variation lives in `ebe_kappas`.
    let individual_estimates_df: Robj = build_individual_estimates(result, population, model);

    // Pre-compute Option<Vec<f64>> fields as Robj (list! macro can't infer
    // From<Option<Vec<f64>>> because Vec<f64> has multiple ToVectorValue impls).
    let obs_time_range_robj: Robj = match result.obs_time_range {
        Some((lo, hi)) => vec![lo, hi].into(),
        None => ().into(),
    };
    let final_gradient_robj: Robj = match &result.final_gradient {
        Some(g) => g.clone().into(),
        None => ().into(),
    };

    // Data-selection exclusion summary (None → NULL; Some → named list).
    let exclusions_robj: Robj = match &result.exclusions {
        Some(ex) => list!(
            n_records_total      = ex.n_records_total as i32,
            n_obs_excluded       = ex.n_obs_excluded as i32,
            n_dose_excluded      = ex.n_dose_excluded as i32,
            n_other_excluded     = ex.n_other_excluded as i32,
            excluded_subject_ids = ex.excluded_subject_ids.clone(),
            fired_ignore         = ex.fired_ignore.clone(),
            fired_accept         = ex.fired_accept.clone(),
        ).into(),
        None => ().into(),
    };

    list!(
        converged = result.converged,
        method = method_label,
        method_chain = method_chain,
        ofv = result.ofv,
        aic = result.aic,
        bic = result.bic,
        n_subjects = result.n_subjects as i32,
        n_obs = result.n_obs as i32,
        n_parameters = result.n_parameters as i32,
        n_iterations = result.n_iterations as i32,
        theta = theta_values,
        theta_names = theta_names,
        omega = omega_flat,
        omega_dim = n_eta as i32,
        sigma = result.sigma.clone(),
        sigma_names = result.sigma_names.clone(),
        // Per-sigma component type — "proportional" or "additive". Parallel
        // to `sigma` and `sigma_names`. The R layer uses this to decide
        // whether to display CV% (proportional) or just SD/variance
        // (additive) for each component without re-deriving from
        // `model_structure$residual`. See ferx-core#57 for the YAML
        // counterpart.
        sigma_types = result
            .sigma_types
            .iter()
            .map(|t| match t {
                SigmaType::Proportional => "proportional".to_string(),
                SigmaType::Additive => "additive".to_string(),
            })
            .collect::<Vec<_>>(),
        se_theta = se_theta,
        se_omega = se_omega,
        se_sigma = se_sigma,
        sdtab = sdtab,
        warnings = warnings,
        warnings_severity = warnings_severity,
        warnings_category = warnings_category,
        warnings_message = warnings_message,
        warnings_source_method = warnings_source_method,
        sir_ess = sir_ess,
        sir_ci_theta = sir_ci_theta,
        sir_ci_omega = sir_ci_omega,
        sir_ci_sigma = sir_ci_sigma,
        sir_resamples = sir_resamples_flat,
        sir_resamples_n = sir_resamples_n,
        sir_resamples_dim = sir_resamples_dim,
        importance_sampling = importance_sampling,
        impmap_trace = impmap_trace,
        bayes = bayes,
        trace_path = trace_path,
        ebe_convergence_warnings = result.ebe_convergence_warnings as i32,
        max_unconverged_subjects = result.max_unconverged_subjects as i32,
        total_ebe_fallbacks = result.total_ebe_fallbacks as i32,
        covariance_status = covariance_status_str,
        shrinkage_eta = result.shrinkage_eta.clone(),
        shrinkage_eps = result.shrinkage_eps,
        wall_time_secs = result.wall_time_secs,
        model_name = result.model_name.clone(),
        ferx_version = result.ferx_version.clone(),
        cov_matrix = cov_matrix_flat,
        cov_matrix_dim = cov_matrix_dim,
        cov_eigenvalues = result.cov_eigenvalues.clone().unwrap_or_default(),
        cov_condition_number = result.cov_condition_number.unwrap_or(f64::NAN),
        omega_iov = omega_iov_flat,
        omega_iov_dim = omega_iov_dim,
        kappa_names = kappa_names,
        se_kappa = se_kappa,
        shrinkage_kappa = shrinkage_kappa,
        shrinkage_kappa_by_occ = shrinkage_kappa_by_occ_df,
        ebe_kappas = ebe_kappas_df,
        ebe_etas = ebe_etas_df,
        individual_estimates = individual_estimates_df,
        model_structure = model_structure_list(model),
        omega_param_corr = omega_param_corr_flat,
        omega_iov_param_corr = omega_iov_param_corr_flat,
        eta_log_transformed = eta_log_transformed,
        eta_param_types = eta_param_types,
        eta_linked_theta = eta_linked_theta,
        eta_names = result.eta_names.clone(),
        theta_transforms = theta_transforms,
        // Resolved gradient methods reported by the engine. `gradient_method_inner`
        // tracks the per-subject EBE gradient (what the user's `gradient` argument
        // controls — Auto resolves to AD or FD at fit time based on build and model
        // shape); `gradient_method_outer` tracks the population-level optimizer's
        // gradient (always FD or N/A today). Strings come from
        // `GradientMethodKind::as_str`: "Enzyme AD", "finite differences", "N/A".
        gradient_method_inner = result.gradient_method_inner.clone(),
        gradient_method_outer = result.gradient_method_outer.clone(),
        // Source-file provenance. `ferx_sir(fit)` reads these to re-parse
        // the model + data and verify SHA-256 hashes; empty strings mean
        // the fit didn't carry them (e.g. older binaries, or a path the
        // hasher couldn't read).
        model_path = result.model_path.clone().unwrap_or_default(),
        data_path = result.data_path.clone().unwrap_or_default(),
        model_hash = result.model_hash.clone().unwrap_or_default(),
        data_hash = result.data_hash.clone().unwrap_or_default(),
        uses_sde = result.uses_sde,
        saem_n_subjects_hmc = result.saem_n_subjects_hmc.map(|n| n as i32),
        dw_statistic = result.dw_statistic,
        iwres_lag1_r = result.iwres_lag1_r,
        // Per-parameter flag: true when the user wrote `(sd)` in the model
        // file for that omega/sigma/kappa line. R printers can use these to
        // annotate estimates with "[initial specified as SD]".
        omega_init_as_sd = result.omega_init_as_sd.clone(),
        sigma_init_as_sd = result.sigma_init_as_sd.clone(),
        kappa_init_as_sd = result.kappa_init_as_sd.clone(),
        // `[covariate_nn]` blocks from the model file (Phase A M1 of ferx-core's
        // DCM plan). One R sub-list per NN; empty list when the `nn` feature is
        // off or no NN blocks are declared. Inspectable in R as
        // `fit$neural_networks[[1]]$shape`, etc.
        neural_networks = neural_networks_list(result),
        // Engine diagnostics (Step 3 — ferx-core feature-parity).
        // saem_mu_ref_m_step_evals_saved: M-step OFV evaluations avoided by
        //   mu-referencing; NULL (NA) for non-SAEM fits.
        // covariance_n_evals_estimated: estimated OFV calls needed for the
        //   covariance step on large models; NULL when not computed.
        // n_threads_used: actual thread count the engine used.
        // nlopt_missing_algorithms: NLopt algorithm names not available on
        //   this platform (empty on most builds).
        saem_mu_ref_m_step_evals_saved = result.saem_mu_ref_m_step_evals_saved.map(|x| x as f64),
        covariance_n_evals_estimated   = result.covariance_n_evals_estimated.map(|x| x as f64),
        n_threads_used                 = result.n_threads_used as i32,
        nlopt_missing_algorithms       = result.nlopt_missing_algorithms.clone(),
        // ── runlog fields (ferx-core#172) ──────────────────────────────────
        // model_text: verbatim .ferx source; NULL for in-memory fits.
        model_text = result.model_text.clone(),
        // theta_init / omega_init / sigma_init: initial estimates before optimisation.
        theta_init = result.theta_init.clone(),
        omega_init = {
            let m = &result.omega_init;
            let n = m.nrows();
            let mut v = Vec::with_capacity(n * n);
            for i in 0..n { for j in 0..n { v.push(m[(i, j)]); } }
            v
        },
        omega_init_dim = result.omega_init.nrows() as i32,
        sigma_init = result.sigma_init.clone(),
        // obs_time_range: c(min_time, max_time) or NULL.
        obs_time_range = obs_time_range_robj,
        // final_gradient: gradient at best-OFV point; NULL for BOBYQA/BFGS/GN/SAEM.
        final_gradient = final_gradient_robj,
        // ── run-settings fields (ferx-core#172 Step 7) ────────────────────
        // optimizer: human-readable label string.
        optimizer_label   = result.optimizer.clone(),
        n_starts          = result.n_starts as i32,
        multi_start_seed  = result.multi_start_seed.map(|s| s as f64),
        saem_seed         = result.saem_seed.map(|s| s as f64),
        sir_seed_used     = result.sir_seed.map(|s| s as f64),
        imp_seed          = result.imp_seed.map(|s| s as f64),
        bloq_method_label = result.bloq_method.clone(),
        outer_maxiter     = result.outer_maxiter as i32,
        outer_gtol        = result.outer_gtol,
        inits_from_nca    = result.inits_from_nca.clone(),
        // covariate_names: character vector of non-standard data columns.
        covariate_names   = result.covariate_names.clone(),
        // input_columns: all CSV header names in file order (analogous to NONMEM $INPUT).
        input_columns     = result.input_columns.clone(),
        // exclusions: NULL when no [data_selection] filter was active; otherwise a named
        // list with counts and fired-condition strings (mirrors ExclusionSummary).
        exclusions        = exclusions_robj,
        // covtab: per-input-row covariate table (ID, TIME, EVID + declared
        // covariate columns); NULL unless the model declares a [covariates] block.
        covtab            = covtab,
        // covariate_types: named char vector covariate -> "continuous"/"categorical".
        covariate_types   = covariate_types
    )
}

/// Build the R-side `neural_networks` list from `result.neural_networks`.
/// Empty list when the `nn` feature is off (each NN's info doesn't exist
/// on FitResult in that build). Each sub-list mirrors `NeuralNetworkInfo`
/// plus a `weights` vector of the trained values for downstream rendering
/// (e.g. a `print()` method or a heatmap).
#[cfg(feature = "nn")]
fn neural_networks_list(result: &FitResult) -> Robj {
    let mut out: Vec<Robj> = Vec::with_capacity(result.neural_networks.len());
    let mut names: Vec<String> = Vec::with_capacity(result.neural_networks.len());
    for nn in &result.neural_networks {
        let weights: Vec<f64> = result.theta
            [nn.weights_offset..nn.weights_offset + nn.n_weights]
            .to_vec();
        let entry = list!(
            name = nn.name.clone(),
            shape = nn.shape.iter().map(|s| *s as i32).collect::<Vec<i32>>(),
            hidden_activation = nn.hidden_activation.clone(),
            output_activation = nn.output_activation.clone(),
            n_weights = nn.n_weights as i32,
            weights_offset = nn.weights_offset as i32,
            input_names = nn.input_names.clone(),
            output_names = nn.output_names.clone(),
            weights = weights,
        );
        names.push(nn.name.clone());
        out.push(entry.into());
    }
    let mut list = List::from_values(out);
    let _ = list.set_names(&names);
    list.into()
}

/// Empty list when the `nn` feature isn't built. Keeps the R-side schema
/// stable (callers always see `fit$neural_networks` as a list) so they can
/// `length(fit$neural_networks)` without branching on build features.
#[cfg(not(feature = "nn"))]
fn neural_networks_list(_result: &FitResult) -> Robj {
    let v: Vec<Robj> = Vec::new();
    List::from_values(v).into()
}

// -- Helper: build per-subject EBE eta data frame --
//
// Returns a data.frame with one row per subject and columns ID + one column
// per BSV eta (named after the model's eta declarations). Returns NULL when
// there are no subjects or the model declares no etas.
fn build_ebe_etas(result: &FitResult, population: &Population) -> Robj {
    let n_subj = result.subjects.len();
    if n_subj == 0 {
        return ().into();
    }
    let n_eta = if result.subjects.is_empty() {
        0
    } else {
        result.subjects[0].eta.len()
    };
    if n_eta == 0 {
        return ().into();
    }

    let eta_names: Vec<String> = if result.eta_names.len() == n_eta {
        result.eta_names.clone()
    } else {
        (0..n_eta).map(|i| format!("ETA{}", i + 1)).collect()
    };

    let mut ids: Vec<String> = Vec::with_capacity(n_subj);
    let mut eta_cols: Vec<Vec<f64>> = (0..n_eta).map(|_| Vec::with_capacity(n_subj)).collect();

    for (si, sr) in result.subjects.iter().enumerate() {
        let subj = &population.subjects[si];
        ids.push(subj.id.clone());
        for k in 0..n_eta {
            eta_cols[k].push(sr.eta.get(k).copied().unwrap_or(f64::NAN));
        }
    }

    let mut pairs: Vec<(&str, Robj)> = Vec::new();
    pairs.push(("ID", ids.into()));
    for k in 0..n_eta {
        pairs.push((eta_names[k].as_str(), eta_cols[k].clone().into()));
    }
    let mut df = List::from_pairs(pairs);
    df.set_class(&["data.frame"]).unwrap();
    let row_names: Vec<i32> = (1..=n_subj as i32).collect();
    df.set_attrib("row.names", row_names).unwrap();
    df.into()
}

// -- Helper: build per-subject individual parameter data frame --
//
// Returns a data.frame with one row per subject and columns ID + one column
// per `[individual_parameters]` declaration (in declaration order). Each cell
// is the value produced by evaluating `pk_param_fn` at the subject's EBE
// eta + zero kappas + covariates. Returns NULL when there are no subjects or
// the model declares no individual parameters.
fn build_individual_estimates(
    result: &FitResult,
    population: &Population,
    model: &CompiledModel,
) -> Robj {
    let n_subj = result.subjects.len();
    if n_subj == 0 {
        return ().into();
    }
    let n_eta = model.n_eta;
    let n_kappa = model.n_kappa;
    let indiv_names = &model.indiv_param_names;
    let n_indiv = indiv_names.len();
    if n_indiv == 0 {
        return ().into();
    }

    let mut ids: Vec<String> = Vec::with_capacity(n_subj);
    let mut param_cols: Vec<Vec<f64>> =
        (0..n_indiv).map(|_| Vec::with_capacity(n_subj)).collect();

    let is_ode = model.is_ode_based();
    let mut eta_buf: Vec<f64> = vec![0.0; n_eta + n_kappa];

    for (si, sr) in result.subjects.iter().enumerate() {
        let subj = &population.subjects[si];
        ids.push(subj.id.clone());

        for k in 0..n_eta {
            eta_buf[k] = sr.eta.get(k).copied().unwrap_or(0.0);
        }
        for k in 0..n_kappa {
            eta_buf[n_eta + k] = 0.0;
        }
        let pk = (model.pk_param_fn)(&result.theta, &eta_buf, &subj.covariates);
        for i in 0..n_indiv {
            // Analytical models route via pk_indices; ODE models write
            // sequentially into slots 0..n_indiv.
            let slot = if is_ode {
                i
            } else {
                model.pk_indices.get(i).copied().unwrap_or(i)
            };
            let v = pk.values.get(slot).copied().unwrap_or(f64::NAN);
            param_cols[i].push(v);
        }
    }

    let mut pairs: Vec<(&str, Robj)> = Vec::new();
    pairs.push(("ID", ids.into()));
    for i in 0..n_indiv {
        pairs.push((indiv_names[i].as_str(), param_cols[i].clone().into()));
    }
    let mut df = List::from_pairs(pairs);
    df.set_class(&["data.frame"]).unwrap();
    let row_names: Vec<i32> = (1..=n_subj as i32).collect();
    df.set_attrib("row.names", row_names).unwrap();
    df.into()
}

fn sdtab_to_dataframe(cols: &[(String, Vec<f64>)]) -> Robj {
    if cols.is_empty() {
        return ().into();
    }

    let mut pairs: Vec<(&str, Robj)> = Vec::new();
    for (name, values) in cols {
        pairs.push((name.as_str(), values.clone().into()));
    }

    // Build data frame via list + class attribute
    let mut df = List::from_pairs(pairs);
    df.set_class(&["data.frame"]).unwrap();

    let n_rows = cols[0].1.len();
    let row_names: Vec<i32> = (1..=n_rows as i32).collect();
    df.set_attrib("row.names", row_names).unwrap();

    df.into()
}

/// Convert a ferx-core `CovariateTable` into an R data.frame, or `NULL` when the
/// model declared no `[covariates]` block. Columns: `ID` (character), `TIME`
/// (numeric), `EVID` (integer), then one numeric column per declared covariate.
/// One row per input dataset record (including dose / EVID rows). Missing
/// covariate values are returned as `NA`.
fn covtab_to_dataframe(table: Option<&ferx_core::CovariateTable>) -> Robj {
    let Some(table) = table else {
        return ().into();
    };

    let n_rows = table.rows.len();
    let ids: Vec<String> = table.rows.iter().map(|r| r.id.clone()).collect();
    let times: Vec<f64> = table.rows.iter().map(|r| r.time).collect();
    let evids: Vec<i32> = table.rows.iter().map(|r| r.evid as i32).collect();

    let mut pairs: Vec<(&str, Robj)> = vec![
        ("ID", ids.into()),
        ("TIME", times.into()),
        ("EVID", evids.into()),
    ];
    // One column per declared covariate, in declaration order (parallel to
    // each row's `values`). Checked indexing (`get`) so a hypothetical short
    // row can't panic across the FFI boundary, and missing values (`NaN`) map
    // to R's `NA_real_` rather than surfacing as `NaN`.
    for (ci, name) in table.names.iter().enumerate() {
        let col: Doubles = table
            .rows
            .iter()
            .map(|r| match r.values.get(ci).copied() {
                Some(v) if !v.is_nan() => Rfloat::from(v),
                _ => Rfloat::na(),
            })
            .collect();
        pairs.push((name.as_str(), col.into()));
    }

    let mut df = List::from_pairs(pairs);
    df.set_class(&["data.frame"]).unwrap();
    let row_names: Vec<i32> = (1..=n_rows as i32).collect();
    df.set_attrib("row.names", row_names).unwrap();

    df.into()
}

/// Per-covariate type tag as a named R character vector (covariate name ->
/// "continuous"/"categorical"), or `NULL` when the model declared no
/// `[covariates]` block. Carries the typed declaration through to R (e.g. so a
/// categorical covariate can be treated as a factor downstream) — the value the
/// numeric `covtab` columns alone don't convey.
fn covariate_types_robj(table: Option<&ferx_core::CovariateTable>) -> Robj {
    let Some(table) = table else {
        return ().into();
    };
    let labels: Vec<String> = table.kinds.iter().map(|k| k.label().to_string()).collect();
    let names: Vec<&str> = table.names.iter().map(|s| s.as_str()).collect();
    let mut robj: Robj = labels.into();
    robj.set_attrib("names", names).unwrap();
    robj
}

/// Returns TRUE if the Rust library was compiled with Enzyme autodiff support.
///
/// The Enzyme autodiff path was retired in ferx-core
/// (FeRx-NLME/ferx-core#381) in favour of hand-rolled analytic sensitivities,
/// so this is always FALSE. Retained for backwards compatibility with callers
/// that probe it.
#[extendr]
fn ferx_rust_autodiff_enabled() -> bool {
    false
}

/// Validate a .ferx model file (and optionally its dataset) without fitting.
///
/// Runs the parser plus every data-independent check, and — when `data_path`
/// is non-empty — the data-dependent checks too (covariates, per-CMT
/// scaling/error-model, steady-state, lagtime). Returns the full
/// `CheckReport` flattened to R-friendly parallel vectors.
///
/// @param model_path Path to .ferx model file
/// @param data_path  Path to NONMEM CSV, or `""` to skip data-dependent checks
/// @return Named list with
///   * `ok` (logical) — true when no error-severity diagnostics
///   * `model` (chr) — model file stem
///   * `data` (chr) — data path used, or `""`
///   * `severity` / `code` / `message` / `block` / `line` / `suggestion`
///     — parallel vectors, one entry per diagnostic
/// @export
#[extendr]
fn ferx_rust_validate_model(model_path: &str, data_path: &str) -> List {
    let data_opt: Option<&str> = if data_path.is_empty() { None } else { Some(data_path) };
    let report = ferx_core::validate_model_file(model_path, data_opt);

    let n = report.diagnostics.len();
    let mut severity: Vec<String> = Vec::with_capacity(n);
    let mut code:     Vec<String> = Vec::with_capacity(n);
    let mut message:  Vec<String> = Vec::with_capacity(n);
    let mut block:    Vec<String> = Vec::with_capacity(n);
    let mut line:     Vec<i32>    = Vec::with_capacity(n);
    let mut suggestion: Vec<String> = Vec::with_capacity(n);
    for d in &report.diagnostics {
        severity.push(match d.severity {
            ferx_core::Severity::Error => "error".to_string(),
            ferx_core::Severity::Warning => "warning".to_string(),
        });
        code.push(d.code.clone());
        message.push(d.message.clone());
        block.push(d.block.clone().unwrap_or_default());
        // 0 sentinel for "no line" — R side maps to NA_integer_.
        line.push(d.line.map(|l| l as i32).unwrap_or(0));
        suggestion.push(d.suggestion.clone().unwrap_or_default());
    }

    list!(
        ok = report.valid,
        model = report.model,
        data = report.data.unwrap_or_default(),
        severity = severity,
        code = code,
        message = message,
        block = block,
        line = line,
        suggestion = suggestion,
    )
}

/// Derive NCA-based starting values from the data without running a fit.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param method One of "nca", "nca_sweep", "nca_ebe"
/// @return Named list with `theta_names`, `theta`, `theta_fixed`, `eta_names`,
///   `omega` (row-major flattened matrix), `omega_dim`, `method`, and
///   `warnings`. Returns an empty list on a parse/read error.
/// @export
#[extendr]
fn ferx_rust_inits_from_nca(model_path: &str, data_path: &str, method: &str) -> List {
    let parsed = match ferx_core::parser::model_parser::parse_full_model_file(Path::new(model_path))
    {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return List::new(0);
        }
    };

    let iov_col = parsed.fit_options.iov_column.clone();
    let population =
        match ferx_core::read_nonmem_csv(Path::new(data_path), None, iov_col.as_deref()) {
            Ok(p) => p,
            Err(e) => {
                rprintln!("Error reading data: {}", e);
                return List::new(0);
            }
        };

    let nca_method = match method.trim().to_lowercase().as_str() {
        "nca" => ferx_core::NcaInit::Nca,
        "" | "true" | "sweep" | "nca_sweep" => ferx_core::NcaInit::Sweep,
        "ebe" | "nca_ebe" => ferx_core::NcaInit::Ebe,
        other => {
            rprintln!(
                "Unknown inits_from_nca method '{}' — expected 'nca', 'nca_sweep', or 'nca_ebe'",
                other
            );
            return List::new(0);
        }
    };
    let method_label = match nca_method {
        ferx_core::NcaInit::Nca => "nca",
        ferx_core::NcaInit::Sweep => "nca_sweep",
        ferx_core::NcaInit::Ebe => "nca_ebe",
    };

    let suggested = ferx_core::inits_from_nca(&parsed.model, &population, nca_method);
    let params = &suggested.params;

    // Omega as a row-major flattened matrix (same convention as fit_result_to_list).
    let n_eta = params.omega.dim();
    let mut omega_flat: Vec<f64> = Vec::with_capacity(n_eta * n_eta);
    for i in 0..n_eta {
        for j in 0..n_eta {
            omega_flat.push(params.omega.matrix[(i, j)]);
        }
    }

    list!(
        theta_names = params.theta_names.clone(),
        theta = params.theta.clone(),
        theta_fixed = params.theta_fixed.clone(),
        eta_names = params.omega.eta_names.clone(),
        omega = omega_flat,
        omega_dim = n_eta as i32,
        method = method_label,
        warnings = suggested.warnings.clone()
    )
}

/// Standalone SIR — run Sampling Importance Resampling against an existing fit.
///
/// The R wrapper `ferx_sir()` flattens the fit list into the primitives this
/// binding expects: parameter estimates, per-subject EBE etas, the asymptotic
/// covariance matrix, the original fit's OFV, and the source paths + hashes
/// captured at fit time. We rebuild a minimal `FitResult` Rust-side (only the
/// fields `ferx_core::run_sir` reads), then call the standalone wrapper which
/// re-parses model + data from the paths, verifies hashes, and runs SIR.
///
/// @param model_path Path to the `.ferx` model file as recorded on the fit.
/// @param data_path Path to the NONMEM CSV as recorded on the fit.
/// @param model_hash Expected SHA-256 of the model file; empty string disables the check.
/// @param data_hash Expected SHA-256 of the data file; empty string disables the check.
/// @param ofv Original fit's OFV (= 2 * nll).
/// @param interaction TRUE for a FOCEI fit, FALSE for FOCE. Controls the inner-loop NLL.
/// @param theta Vector of theta point estimates.
/// @param omega_flat Row-major flattened omega matrix.
/// @param omega_dim Dimension of the omega matrix.
/// @param sigma Vector of sigma point estimates.
/// @param cov_matrix_flat Row-major flattened parameter covariance matrix.
/// @param cov_matrix_dim Dimension of the covariance matrix.
/// @param eta_hats_flat Row-major flattened per-subject EBE etas (n_subjects × n_eta).
/// @param n_subjects Number of subjects.
/// @param sir_samples Number of proposal samples (M).
/// @param sir_resamples Number of resamples (m); must be <= M.
/// @param sir_seed Random seed; pass -1 for the engine default.
/// @param sir_keep_samples When TRUE, retains the resampled packed parameter vectors.
/// @param verbose When TRUE, the engine prints progress to stderr.
/// @return Named list with `sir_ess`, `sir_ci_theta`, `sir_ci_omega`, `sir_ci_sigma`,
///   `sir_resamples`, `sir_resamples_n`, `sir_resamples_dim`.
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_sir(
    model_path: &str,
    data_path: &str,
    model_hash: &str,
    data_hash: &str,
    ofv: f64,
    interaction: bool,
    theta: Vec<f64>,
    omega_flat: Vec<f64>,
    omega_dim: i32,
    sigma: Vec<f64>,
    cov_matrix_flat: Vec<f64>,
    cov_matrix_dim: i32,
    eta_hats_flat: Vec<f64>,
    n_subjects: i32,
    sir_samples: i32,
    sir_resamples: i32,
    sir_seed: i32,
    sir_keep_samples: bool,
    verbose: bool,
) -> Robj {
    // All error paths in this binding throw an R condition (via
    // `throw_r_error`) rather than printing to stderr + returning NULL.
    // The latter pattern (used by `ferx_rust_fit`) loses the engine
    // message to stderr — callers using `tryCatch()` or
    // `expect_error(..., regexp = ...)` see only a generic
    // "backend returned no result" from the R wrapper. Throwing
    // propagates the actual message (e.g. "hash mismatch") into the R
    // condition, which is what test code expects and what users want.
    let parsed = match ferx_core::parse_full_model_file(Path::new(model_path)) {
        Ok(p) => p,
        Err(e) => throw_r_error(&format!(
            "ferx_sir: error parsing model at {}: {}",
            model_path, e
        )),
    };
    let model = &parsed.model;
    let template = &model.default_params;

    let n_theta = theta.len();
    let n_sigma = sigma.len();
    let n_eta = omega_dim as usize;
    let n_subj = n_subjects.max(0) as usize;
    let n_packed = cov_matrix_dim as usize;

    if n_theta != template.theta.len() {
        throw_r_error(&format!(
            "ferx_sir: theta length {} does not match model ({} expected)",
            n_theta,
            template.theta.len()
        ));
    }
    if n_sigma != template.sigma.values.len() {
        throw_r_error(&format!(
            "ferx_sir: sigma length {} does not match model ({} expected)",
            n_sigma,
            template.sigma.values.len()
        ));
    }
    if n_eta != template.omega.dim() {
        throw_r_error(&format!(
            "ferx_sir: omega dim {} does not match model ({} expected)",
            n_eta,
            template.omega.dim()
        ));
    }
    if omega_flat.len() != n_eta * n_eta {
        throw_r_error(&format!(
            "ferx_sir: omega_flat length {} does not match dim^2 = {}",
            omega_flat.len(),
            n_eta * n_eta
        ));
    }
    if n_packed == 0 || cov_matrix_flat.len() != n_packed * n_packed {
        throw_r_error(&format!(
            "ferx_sir: cov_matrix is missing or malformed (dim={}, len={}). \
             Re-fit with `covariance = TRUE`.",
            n_packed,
            cov_matrix_flat.len()
        ));
    }
    if eta_hats_flat.len() != n_subj * n_eta {
        throw_r_error(&format!(
            "ferx_sir: eta_hats_flat length {} does not match n_subjects * n_eta = {}",
            eta_hats_flat.len(),
            n_subj * n_eta
        ));
    }

    let omega_mat = DMatrix::from_row_slice(n_eta, n_eta, &omega_flat);
    let cov_mat = DMatrix::from_row_slice(n_packed, n_packed, &cov_matrix_flat);

    // Build SubjectResult vec with only `eta` populated (the only field
    // ferx_core::run_sir reads off subjects). IDs are synthesised because
    // SIR doesn't consume them; the original IDs live on the R fit list.
    let mut subjects: Vec<SubjectResult> = Vec::with_capacity(n_subj);
    for i in 0..n_subj {
        let mut eta = nalgebra::DVector::<f64>::zeros(n_eta);
        for k in 0..n_eta {
            eta[k] = eta_hats_flat[i * n_eta + k];
        }
        subjects.push(SubjectResult {
            id: format!("{}", i + 1),
            eta,
            ipred: Vec::new(),
            pred: Vec::new(),
            iwres: Vec::new(),
            cwres: Vec::new(),
            ofv_contribution: 0.0,
            cens: Vec::new(),
            n_obs: 0,
            extra_columns: Vec::new(),
            per_obs_tad: Vec::new(),
            // PR #207 (ferx-core) added this field; the SIR path never reads it.
            compartment_states: Vec::new(),
            // PR #377 (ferx-core) added NPDE/NPD; the SIR path never reads them.
            npde: Vec::new(),
            npd: Vec::new(),
        });
    }

    // Skeleton FitResult — only the fields ferx_core::run_sir actually
    // reads are populated; everything else gets a neutral default.
    let fit = FitResult {
        method: if interaction {
            EstimationMethod::FoceI
        } else {
            EstimationMethod::Foce
        },
        method_chain: vec![if interaction {
            EstimationMethod::FoceI
        } else {
            EstimationMethod::Foce
        }],
        bayes: None,
        cond_dist: None,
        converged: true,
        ofv,
        aic: 0.0,
        bic: 0.0,
        theta: theta.clone(),
        theta_names: template.theta_names.clone(),
        eta_names: template.omega.eta_names.clone(),
        omega: omega_mat,
        sigma: sigma.clone(),
        sigma_names: template.sigma.names.clone(),
        error_model: model.error_model,
        covariance_matrix: Some(cov_mat),
        se_theta: None,
        se_omega: None,
        se_sigma: None,
        theta_fixed: template.theta_fixed.clone(),
        omega_fixed: template.omega_fixed.clone(),
        sigma_fixed: template.sigma_fixed.clone(),
        subjects,
        n_obs: 0,
        n_subjects: n_subj,
        n_parameters: n_packed,
        n_iterations: 0,
        interaction,
        warnings: Vec::new(),
        sir_ci_theta: None,
        sir_ci_omega: None,
        sir_ci_sigma: None,
        sir_ess: None,
        sir_resamples_packed: None,
        importance_sampling: None,
        impmap_trace: None,
        omega_iov: None,
        kappa_names: model.kappa_names.clone(),
        kappa_fixed: template.kappa_fixed.clone(),
        se_kappa: None,
        shrinkage_kappa: Vec::new(),
        shrinkage_kappa_by_occ: Vec::new(),
        ebe_kappas: Vec::new(),
        saem_mu_ref_m_step_evals_saved: None,
        saem_n_subjects_hmc: None,
        gradient_method_inner: String::new(),
        gradient_method_outer: String::new(),
        uses_ode_solver: model.is_ode_based(),
        n_threads_used: 1,
        nlopt_missing_algorithms: Vec::new(),
        covariance_n_evals_estimated: None,
        trace_path: None,
        ebe_convergence_warnings: 0,
        max_unconverged_subjects: 0,
        total_ebe_fallbacks: 0,
        covariance_status: CovarianceStatus::Computed,
        shrinkage_eta: Vec::new(),
        shrinkage_eps: f64::NAN,
        wall_time_secs: 0.0,
        model_name: model.name.clone(),
        ferx_version: String::new(),
        eta_param_info: Vec::new(),
        theta_transform: Vec::new(),
        sigma_types: Vec::new(),
        cov_eigenvalues: None,
        cov_condition_number: None,
        eta_log_transformed: Vec::new(),
        omega_param_corr: None,
        omega_iov_param_corr: None,
        model_path: Some(model_path.to_string()),
        data_path: Some(data_path.to_string()),
        model_hash: if model_hash.is_empty() {
            None
        } else {
            Some(model_hash.to_string())
        },
        data_hash: if data_hash.is_empty() {
            None
        } else {
            Some(data_hash.to_string())
        },
        dw_statistic: f64::NAN,
        iwres_lag1_r: f64::NAN,
        uses_sde: false,
        omega_init_as_sd: Vec::new(),
        sigma_init_as_sd: Vec::new(),
        kappa_init_as_sd: Vec::new(),
        warnings_structured: Vec::new(),
        model_text: None,
        theta_init: Vec::new(),
        omega_init: nalgebra::DMatrix::zeros(0, 0),
        sigma_init: Vec::new(),
        obs_time_range: None,
        final_gradient: None,
        optimizer: "auto".to_string(),
        n_starts: 1,
        multi_start_seed: None,
        saem_seed: None,
        sir_seed: None,
        imp_seed: None,
        npde_seed: None,
        bloq_method: "drop".to_string(),
        outer_maxiter: 0,
        outer_gtol: 0.0,
        inits_from_nca: None,
        covariate_names: Vec::new(),
        input_columns: Vec::new(),
        covariate_table: None,
        exclusions: None,
        #[cfg(feature = "nn")]
        neural_networks: Vec::new(),
    };

    let mut opts = FitOptions::default();
    opts.sir_samples = sir_samples.max(0) as usize;
    opts.sir_resamples = sir_resamples.max(0) as usize;
    opts.sir_seed = if sir_seed < 0 {
        None
    } else {
        Some(sir_seed as u64)
    };
    opts.sir_keep_samples = sir_keep_samples;
    opts.interaction = interaction;
    opts.verbose = verbose;

    // ferx_core::run_sir verifies hashes (when set), re-parses model + data,
    // reconstructs the inner ModelParameters, and runs SIR. We pass None for
    // model/population so that path goes through and the integrity check
    // fires — that is the whole point of having hashes here.
    let new_fit = match ferx_core::run_sir(&fit, None, None, &opts) {
        Ok(f) => f,
        Err(e) => throw_r_error(&format!("ferx_sir: {}", e)),
    };

    let flatten_ci = |ci: &Option<Vec<(f64, f64)>>| -> Vec<f64> {
        ci.as_ref()
            .map(|v| v.iter().flat_map(|(lo, hi)| [*lo, *hi]).collect())
            .unwrap_or_default()
    };
    let (sir_resamples_flat, sir_resamples_n, sir_resamples_dim): (Vec<f64>, i32, i32) =
        match &new_fit.sir_resamples_packed {
            Some(pool) if !pool.is_empty() => {
                let n = pool.len();
                let d = pool[0].len();
                let mut v = Vec::with_capacity(n * d);
                for row in pool {
                    v.extend_from_slice(row);
                }
                (v, n as i32, d as i32)
            }
            _ => (Vec::new(), 0i32, 0i32),
        };

    list!(
        sir_ess = new_fit.sir_ess.unwrap_or(f64::NAN),
        sir_ci_theta = flatten_ci(&new_fit.sir_ci_theta),
        sir_ci_omega = flatten_ci(&new_fit.sir_ci_omega),
        sir_ci_sigma = flatten_ci(&new_fit.sir_ci_sigma),
        sir_resamples = sir_resamples_flat,
        sir_resamples_n = sir_resamples_n,
        sir_resamples_dim = sir_resamples_dim
    )
    .into()
}

/// Prepare a FREM (Full Random Effects Model) dataset and model file.
///
/// Transforms a base model + dataset into a FREM model that treats covariates
/// as additional dependent variables. The covariance structure of an extended
/// omega matrix captures covariate-parameter relationships implicitly.
///
/// @param model_path Path to base .ferx model file
/// @param data_path Path to NONMEM-format CSV
/// @param covariates Character vector of covariate column names
/// @param categorical_covariates Character vector of categorical covariate
///   names (subset of `covariates`). These are binarized (one-hot encoded)
///   before FREM transformation. Pass an empty vector if none.
/// @param output_model_path Optional path for the output .ferx model file.
///   If empty, defaults to `<model_dir>/<model_stem>_frem.ferx`.
/// @param output_data_path Optional path for the output CSV file.
///   If empty, defaults to `<model_dir>/<model_stem>_frem_data.csv`.
/// @return Named list with model_path, data_path, covariate_means,
///   covariate_variances, fremtype_map, n_total_etas
/// @export
#[extendr]
fn ferx_rust_prepare_frem(
    model_path: &str,
    data_path: &str,
    covariates: Vec<String>,
    categorical_covariates: Vec<String>,
    output_model_path: &str,
    output_data_path: &str,
) -> List {
    let cat_opt: Option<Vec<String>> = if categorical_covariates.is_empty() {
        None
    } else {
        Some(categorical_covariates)
    };
    let out_model: Option<&Path> = if output_model_path.is_empty() {
        None
    } else {
        Some(Path::new(output_model_path))
    };
    let out_data: Option<&Path> = if output_data_path.is_empty() {
        None
    } else {
        Some(Path::new(output_data_path))
    };

    let result = match ferx_core::prepare_frem(
        Path::new(model_path),
        Path::new(data_path),
        &covariates,
        cat_opt.as_deref(),
        out_model,
        out_data,
        None, // missing value indicator (default: -99)
    ) {
        Ok(r) => r,
        Err(e) => throw_r_error(format!("Error in prepare_frem: {e}")),
    };

    let mean_names: Vec<String> = result.covariate_means.iter().map(|(n, _)| n.clone()).collect();
    let mean_vals: Vec<f64> = result.covariate_means.iter().map(|(_, v)| *v).collect();
    let var_names: Vec<String> = result.covariate_variances.iter().map(|(n, _)| n.clone()).collect();
    let var_vals: Vec<f64> = result.covariate_variances.iter().map(|(_, v)| *v).collect();
    let ft_names: Vec<String> = result.fremtype_map.iter().map(|(n, _)| n.clone()).collect();
    let ft_vals: Vec<i32> = result.fremtype_map.iter().map(|(_, v)| *v as i32).collect();

    let warnings: Vec<String> = result.warnings.clone();

    list!(
        model_path = result.model_path.to_string_lossy().to_string(),
        data_path = result.data_path.to_string_lossy().to_string(),
        covariate_mean_names = mean_names,
        covariate_mean_values = mean_vals,
        covariate_var_names = var_names,
        covariate_var_values = var_vals,
        fremtype_names = ft_names,
        fremtype_values = ft_vals,
        n_total_etas = result.n_total_etas as i32,
        warnings = warnings
    )
}

extendr_module! {
    mod ferx;
    fn ferx_rust_fit;
    fn ferx_rust_simulate;
    fn ferx_rust_simulate_from_fit;
    fn ferx_rust_simulate_with_uncertainty;
    fn ferx_rust_predict;
    fn ferx_rust_predict_from_fit;
    fn ferx_rust_predict_survival;
    fn ferx_rust_predict_survival_from_fit;
    fn ferx_rust_npde_from_fit;
    fn ferx_rust_sir;
    fn ferx_rust_autodiff_enabled;
    fn ferx_rust_validate_model;
    fn ferx_rust_inits_from_nca;
    fn ferx_rust_prepare_frem;
}
