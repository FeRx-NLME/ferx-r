use extendr_api::prelude::*;
use ferx_nlme::cancel::CancelFlag;
use ferx_nlme::types::*;
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
//      exit cooperatively; ferx_nlme::fit returns Err("cancelled by user").
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
        match ferx_nlme::parser::model_parser::parse_full_model_file(Path::new(model_path)) {
            Ok(p) => p,
            Err(e) => {
                rprintln!("Error parsing model: {}", e);
                return List::new(0);
            }
        };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return List::new(0);
        }
    };

    // Build fit options
    let mut opts = parsed.fit_options.clone();

    // Apply generic `settings` list before the dedicated args below — the
    // dedicated args are the single source of truth for their keys and always
    // win. Reserved keys (those with a dedicated R argument) are rejected so
    // the precedence rule is explicit rather than silently clobbered.
    if settings_keys.len() != settings_values.len() {
        rprintln!(
            "Error: settings keys/values length mismatch ({} vs {})",
            settings_keys.len(),
            settings_values.len()
        );
        return List::new(0);
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
            rprintln!(
                "Error: setting `{}` conflicts with a dedicated ferx_fit() argument — pass it via that argument instead",
                key
            );
            return List::new(0);
        }
        match ferx_nlme::parser::model_parser::apply_fit_option(&mut opts, key, v) {
            Ok(true) => {}
            Ok(false) => {
                rprintln!("Error: unknown fit setting `{}`", key);
                return List::new(0);
            }
            Err(e) => {
                rprintln!("Error: {}", e);
                return List::new(0);
            }
        }
    }

    if method.is_empty() {
        rprintln!("Error: `method` must contain at least one estimation method");
        return List::new(0);
    }
    let chain: Vec<EstimationMethod> = match method.iter().map(|m| parse_method(m)).collect() {
        Ok(v) => v,
        Err(e) => {
            rprintln!("{}", e);
            return List::new(0);
        }
    };
    let final_method = *chain.last().unwrap();
    opts.method = final_method;
    opts.interaction = final_method == EstimationMethod::FoceI;
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
        "" | "auto" => opts.gradient_method = ferx_nlme::GradientMethod::Auto,
        "ad" | "autodiff" => opts.gradient_method = ferx_nlme::GradientMethod::Ad,
        "fd" | "finite" | "finite_difference" => {
            opts.gradient_method = ferx_nlme::GradientMethod::Fd
        }
        other => {
            rprintln!(
                "Unknown gradient method '{}' — expected 'auto', 'ad', or 'fd' (falling back to auto)",
                other
            );
            opts.gradient_method = ferx_nlme::GradientMethod::Auto;
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
    let result = std::thread::scope(|s| {
        let handle = s.spawn(|| ferx_nlme::fit(&parsed.model, &population, &init_params, &opts));

        while !handle.is_finished() {
            if pending_interrupt() {
                cancel.cancel();
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(POLL_MS));
        }

        handle.join()
    });

    let result = match result {
        Ok(Ok(r)) => r,
        Ok(Err(_)) if cancel.is_cancelled() => {
            // Raise a proper R error so the user sees a clean condition
            // instead of a silent empty-list return. A small one-time leak
            // of the worker's locals is acceptable here — one R session's
            // worth of cancelled fits won't add up to anything meaningful.
            throw_r_error("ferx_fit: cancelled by user");
        }
        Ok(Err(e)) => {
            rprintln!("Fit error: {}", e);
            return List::new(0);
        }
        Err(_) => {
            rprintln!("Fit error: worker thread panicked");
            return List::new(0);
        }
    };

    // Convert to R list
    fit_result_to_list(&result, &population)
}

/// Simulate from a NLME model.
///
/// @param model_path Path to .ferx model file
/// @param data_path Path to NONMEM-format CSV (for population structure)
/// @param n_sim Number of simulations
/// @param seed Random seed
/// @return Data frame with ID, TIME, IPRED, DV_SIM columns
/// @export
#[extendr]
fn ferx_rust_simulate(
    model_path: &str,
    data_path: &str,
    n_sim: i32,
    seed: i32,
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::simulate_with_seed(
        &model,
        &population,
        &model.default_params,
        n_sim as usize,
        seed as u64,
    );

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
) -> Robj {
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results =
        ferx_nlme::simulate_with_seed(&model, &population, &params, n_sim as usize, seed as u64);
    sim_results_to_df(&results)
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
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::predict(&model, &population, &model.default_params);

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
    let model = match ferx_nlme::parse_model_file(Path::new(model_path)) {
        Ok(m) => m,
        Err(e) => {
            rprintln!("Error parsing model: {}", e);
            return ().into();
        }
    };

    let population = match ferx_nlme::read_nonmem_csv(Path::new(data_path), None) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("Error reading data: {}", e);
            return ().into();
        }
    };

    let params = match params_from_fit(&model, &theta, &omega_flat, omega_dim, &sigma) {
        Ok(p) => p,
        Err(e) => {
            rprintln!("{}", e);
            return ().into();
        }
    };

    let results = ferx_nlme::predict(&model, &population, &params);

    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let pred: Vec<f64> = results.iter().map(|r| r.pred).collect();

    data_frame!(ID = id, TIME = time, PRED = pred).into()
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
    } else {
        Err(format!(
            "Unknown estimation method '{}' — expected one of: foce, focei, saem, gn, gn_hybrid",
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
    })
}

// -- Helper: SimulationResult slice → R data frame --

fn sim_results_to_df(results: &[ferx_nlme::api::SimulationResult]) -> Robj {
    let sim: Vec<i32> = results.iter().map(|r| r.sim as i32).collect();
    let id: Vec<String> = results.iter().map(|r| r.id.clone()).collect();
    let time: Vec<f64> = results.iter().map(|r| r.time).collect();
    let ipred: Vec<f64> = results.iter().map(|r| r.ipred).collect();
    let dv_sim: Vec<f64> = results.iter().map(|r| r.dv_sim).collect();

    data_frame!(
        SIM = sim,
        ID = id,
        TIME = time,
        IPRED = ipred,
        DV_SIM = dv_sim
    )
    .into()
}

// -- Helper: convert FitResult + Population to R named list --

fn fit_result_to_list(result: &FitResult, population: &Population) -> List {
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
    let sdtab_cols = ferx_nlme::io::output::sdtab(result, population);
    let sdtab = sdtab_to_dataframe(&sdtab_cols);

    // Warnings
    let warnings: Vec<String> = result.warnings.clone();

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

    let trace_path: Robj = match &result.trace_path {
        Some(p) => p.clone().into(),
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
        se_theta = se_theta,
        se_omega = se_omega,
        se_sigma = se_sigma,
        sdtab = sdtab,
        warnings = warnings,
        sir_ess = sir_ess,
        sir_ci_theta = sir_ci_theta,
        sir_ci_omega = sir_ci_omega,
        sir_ci_sigma = sir_ci_sigma,
        trace_path = trace_path
    )
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

/// Returns TRUE if the Rust library was compiled with the `autodiff` feature
/// (Enzyme toolchain), FALSE otherwise.
#[extendr]
fn ferx_rust_autodiff_enabled() -> bool {
    cfg!(feature = "autodiff")
}

extendr_module! {
    mod ferx;
    fn ferx_rust_fit;
    fn ferx_rust_simulate;
    fn ferx_rust_simulate_from_fit;
    fn ferx_rust_predict;
    fn ferx_rust_predict_from_fit;
    fn ferx_rust_autodiff_enabled;
}
