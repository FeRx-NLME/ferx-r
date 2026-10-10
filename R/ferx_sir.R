#' Run SIR against an existing fit
#'
#' Run a Sampling Importance Resampling (SIR) uncertainty step against a
#' fit that was produced earlier - useful when the original fit was
#' expensive and you want to add SIR without re-estimating, or when working
#' with a fit loaded from a `.fitrx` bundle.
#'
#' `ferx_sir()` re-uses the fit's asymptotic covariance matrix as the SIR
#' proposal distribution and the per-subject empirical Bayes ETAs as
#' warm-starts for the inner loop. The returned fit is the input with
#' `sir_ess`, `sir_ci_theta`, `sir_ci_omega`, `sir_ci_sigma`, `sir_ci_kappa`
#' (IOV models only) and, when `sir_keep_samples = TRUE`, `sir_resamples` /
#' `sir_resamples_n` / `sir_resamples_dim` populated.
#'
#' The fit is rebuilt at its fitted estimates, including the IOV (kappa)
#' covariance `fit$omega_iov`, under the inner loop of its last estimation
#' method (FOCE or FOCEI).
#'
#' ## Settings
#'
#' A fit made with `ferx_fit(..., sir = TRUE)`, or returned by an earlier
#' `ferx_sir()`, records the settings its SIR ran under as `fit$sir_settings`:
#' the sample sizes, the seed, the proposal degrees of freedom (`df`), the
#' `scale`, `keep_samples`, and under `scoring` the inner-loop and ODE settings
#' each draw was scored with. Every argument left at `NULL` takes the recorded
#' value, so `ferx_sir(fit)` repeats the fit's SIR: its intervals and effective
#' sample size are `identical()` to the in-fit ones, also after a
#' [ferx_save_fit()] / [ferx_load_fit()] round trip. The exception is a
#' normal proposal (`sir_df = Inf`): a `.fitrx` bundle cannot record it yet
#' (ferx-core#1819), so `ferx_save_fit()` warns and saves no SIR settings, and
#' the loaded fit runs as one with no record. An argument you pass
#' replaces that one setting and keeps the rest. The returned fit's
#' `sir_settings` records the run it made.
#'
#' A fit with no SIR record (SIR never ran on it, or it was saved before the
#' record existed) runs with the engine's defaults for every argument left at
#' `NULL`: 1000 samples, 250 resamples, seed 12345, `sir_df = 5`,
#' `sir_scale = "packed"`, and default inner-loop and ODE settings
#' (tolerances written in the model file's `[fit_options]` are kept). On such
#' a fit, settings passed to `ferx_fit()` through `settings =` (an
#' `inner_maxiter`, an `ode_reltol`) are not used: the engine reads them only
#' from a SIR record (ferx-core #1806).
#'
#' ## Data selection
#'
#' The data are re-read with the selection the fit was made with, recorded as
#' `fit$reader_settings`: a `ferx_fit(ignore =, accept =, ignore_ids =)`
#' selection applies as well as the model file's `[data_selection]`, and the
#' engine checks the rows it reads against `fit$population_fingerprint`. A fit
#' made before these were recorded is refused when its `fit$exclusions` show a
#' clause, or its `fit$call_settings` an `iov_column`, that the model file does
#' not state, rather than scored on other data.
#'
#' ## Integrity check
#'
#' `ferx_fit()` records the model and data file paths plus SHA-256 hashes
#' on the returned fit (`fit$model_path`, `fit$data_path`,
#' `fit$model_hash`, `fit$data_hash`). `ferx_sir()` re-reads those files
#' and verifies the hashes; if either file changed since the fit, the call
#' is a **hard error**. The point of running SIR against the original fit
#' is to refine its uncertainty estimate, which is meaningless against a
#' modified model or dataset.
#'
#' Edge case: if hashing failed at fit time (e.g. permission flip between
#' parse and hash) the corresponding `fit$*_hash` is `NA` and the
#' integrity check silently passes on that side. This is rare in
#' practice - hashing would have to fail while the parse just
#' succeeded - but if you require integrity verification, check
#' `!is.na(fit$model_hash) && !is.na(fit$data_hash)` before relying on
#' the protection.
#'
#' @param fit A `ferx_fit` object produced by [ferx_fit()] or
#'   [ferx_load_fit()].
#' @param sir_samples Number of proposal samples drawn from the asymptotic
#'   distribution. Higher values give tighter weights at proportional cost.
#'   `NULL` (the default) takes the fit's recorded value, else 1000.
#' @param sir_resamples Number of resampled vectors. Must be `<= sir_samples`.
#'   `NULL` takes the fit's recorded value, else 250.
#' @param sir_seed Integer RNG seed for reproducibility. `NULL` takes the
#'   seed the fit's SIR recorded, else the engine's built-in seed (12345).
#' @param sir_keep_samples When `TRUE`, retain the resampled packed
#'   parameter vectors on the returned fit. Required for
#'   [ferx_simulate_with_uncertainty()] with `method = "sir"`. `NULL` takes
#'   the fit's recorded value, else `FALSE`.
#' @param verbose When `TRUE`, the engine prints progress to stderr.
#'   Default `FALSE`.
#' @param sir_scale The parameter scale SIR's importance-sampling target is
#'   flat on (ferx-core #1723). `"packed"` (the default, and the only scale
#'   before ferx-core #1723) is flat on the optimizer's packed scale (log-sd
#'   for a variance); a variance the data cannot bound away from zero then has
#'   a likelihood shelf down to the parameter box floor, so its SIR lower limit
#'   tracks the box. `"natural"` is flat on the reported scale (Omega / kappa
#'   variances, sigma as a variance, theta as declared), the PsN SIR
#'   convention: its lower limits no longer depend on the box, but a variance
#'   informed by few groups gets a heavy upper tail. `"natural"` is refused for
#'   a model with `prior(...)`. The same option is
#'   `settings = list(sir_scale = ...)` in [ferx_fit()] and `sir_scale` in the
#'   model's `[fit_options]`. `NULL` (the default) takes the fit's recorded
#'   scale, else `"packed"`; this argument does not read the model file's
#'   value.
#' @param sir_df Degrees of freedom of the Student-t SIR proposal, at least 1
#'   (`Inf` is a multivariate-normal proposal). The same option is
#'   `settings = list(sir_df = ...)` in [ferx_fit()]. `NULL` (the default)
#'   takes the fit's recorded value, else 5.
#'
#' @return The input `fit`, augmented with `sir_ess`, `sir_ci_theta`,
#'   `sir_ci_omega`, `sir_ci_sigma`, `sir_ci_kappa` (one row per IOV kappa
#'   variance, named by `kappa_names`; NULL without IOV), `sir_seed_used`
#'   (the seed this run resampled with: `sir_seed`, or 12345 when it is
#'   `NULL`), `sir_settings` (the settings this run used; see [ferx_fit()]),
#'   and (when the run kept them) `sir_resamples` / `sir_resamples_n` /
#'   `sir_resamples_dim`.
#'   Any warnings the SIR step emitted are appended to `fit$warnings` and to
#'   `fit$warnings_structured` under the `sir` category - in particular the
#'   proposal diagnostics: a covariance that is rank-deficient beyond its
#'   `FIX`ed parameters, or a proposal direction shrunk to keep draws inside
#'   the parameter bounds. Both name the parameters involved and mean the same
#'   thing - those directions are not identified by the data, and their SIR
#'   intervals understate the uncertainty. A run whose effective sample size
#'   is below 100 adds a `SIR:` low-ESS warning naming the draw the intervals
#'   hinge on and, under `sir_scale = "packed"`, every Omega / kappa variance
#'   the data do not bound away from zero, whose SIR lower limit then reflects
#'   the parameter box (ferx-core #1723). This run's `SIR:`
#'   lines replace the ones an earlier SIR run (inline `sir = TRUE`, or a
#'   previous `ferx_sir()`) left on the fit, as do its `SIR failed:` lines, so
#'   the warnings describe the intervals on the returned fit. See
#'   [ferx_get_warnings()].
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE)
#'
#' # Run SIR post-hoc (equivalent to sir = TRUE at fit time)
#' fit <- ferx_sir(fit, sir_samples = 2000, sir_resamples = 500, sir_seed = 42)
#'
#' # Effective sample size - closer to sir_resamples means good coverage
#' fit$sir_ess
#'
#' # 95% CI for fixed-effect thetas
#' fit$sir_ci_theta
#'
#' # 95% CI for IIV omega diagonal
#' fit$sir_ci_omega
#'
#' # 95% CI for residual error sigma
#' fit$sir_ci_sigma
#'
#' # Below ESS 100 a "SIR:" warning names the variances whose lower limit
#' # tracks the parameter box; the natural scale removes that dependence
#' grep("^SIR:", fit$warnings, value = TRUE)
#' fit_nat <- ferx_sir(fit, sir_seed = 42, sir_scale = "natural")
#' fit_nat$sir_ci_omega
#'
#' # Retain resamples for downstream uncertainty simulation
#' fit2 <- ferx_sir(fit, sir_samples = 2000, sir_resamples = 500,
#'                  sir_keep_samples = TRUE)
#' sims <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit2,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 5,
#'   method = "sir"
#' )
#' head(sims)
#' }
#'
#' @seealso [ferx_fit()] for the inline SIR option (`sir = TRUE`),
#'   [ferx_simulate_with_uncertainty()] for downstream consumption of the
#'   retained resamples.
#' @family fitting
#' @export
ferx_sir <- function(fit,
                     sir_samples = NULL,
                     sir_resamples = NULL,
                     sir_seed = NULL,
                     sir_keep_samples = NULL,
                     verbose = FALSE,
                     sir_scale = NULL,
                     sir_df = NULL) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }

  # Validate the SIR knobs up-front so user errors surface as clean R
  # conditions rather than as the engine's "weighted sampler failed"
  # or the binding's i32-clamped 0 after a bad cast.
  is_positive_count <- function(x) {
    length(x) == 1L && is.finite(x) && x > 0 && x == as.integer(x)
  }
  if (!is.null(sir_samples) && !is_positive_count(sir_samples)) {
    stop("`sir_samples` must be NULL or a single positive integer (got: ",
         paste(format(sir_samples), collapse = ", "), ").")
  }
  if (!is.null(sir_resamples) && !is_positive_count(sir_resamples)) {
    stop("`sir_resamples` must be NULL or a single positive integer (got: ",
         paste(format(sir_resamples), collapse = ", "), ").")
  }
  if (!is.null(sir_seed)) {
    if (!is_positive_count(sir_seed) && !(length(sir_seed) == 1L &&
                                          is.finite(sir_seed) &&
                                          sir_seed == as.integer(sir_seed) &&
                                          sir_seed >= 0)) {
      stop("`sir_seed` must be NULL or a single non-negative integer.")
    }
  }
  if (!is.null(sir_keep_samples) &&
        !(is.logical(sir_keep_samples) && length(sir_keep_samples) == 1L &&
            !is.na(sir_keep_samples))) {
    stop("`sir_keep_samples` must be NULL, TRUE or FALSE.")
  }
  if (!is.null(sir_scale)) sir_scale <- match.arg(sir_scale, c("packed", "natural"))
  if (!is.null(sir_df) &&
        !(is.numeric(sir_df) && length(sir_df) == 1L && !is.na(sir_df) && sir_df >= 1)) {
    stop("`sir_df` must be NULL or a single number >= 1 (`Inf` for a normal proposal).")
  }

  model_path <- fit$model_path
  if (is.null(model_path) || is.na(model_path) || !nzchar(model_path)) {
    stop(
      "ferx_sir: fit has no recorded model_path; cannot locate the model file. ",
      "Re-fit via ferx_fit(model, data) so the path is recorded, or pass ",
      "the model and data explicitly via ferx_rust_sir()."
    )
  }
  if (!file.exists(model_path)) {
    stop("ferx_sir: model file not found at ", model_path)
  }

  data_path <- fit$data_path
  if (is.null(data_path) || is.na(data_path) || !nzchar(data_path)) {
    stop(
      "ferx_sir: fit has no recorded data_path; cannot locate the data file. ",
      "Re-fit via ferx_fit(model, data) so the path is recorded."
    )
  }
  if (!file.exists(data_path)) {
    stop("ferx_sir: data file not found at ", data_path)
  }

  if (is.null(fit$cov_matrix)) {
    stop(
      "ferx_sir: fit has no covariance matrix to use as the SIR proposal. ",
      "Re-fit with `covariance = TRUE` (and verify the cov step converged)."
    )
  }

  # The explicit arguments, edited into the fit's SIR record (or the engine's
  # default record when the fit has none); every argument left NULL keeps the
  # recorded value, so `ferx_sir(fit)` repeats the fit's SIR (#472).
  explicit <- list(
    samples = if (!is.null(sir_samples)) as.integer(sir_samples),
    resamples = if (!is.null(sir_resamples)) as.integer(sir_resamples),
    seed = if (!is.null(sir_seed)) as.numeric(sir_seed),
    df = if (!is.null(sir_df)) as.numeric(sir_df),
    scale = sir_scale,
    keep_samples = sir_keep_samples
  )
  explicit <- explicit[!vapply(explicit, is.null, logical(1L))]
  record_args <- .ferx_fit_record_args(fit, sir = explicit)
  resolved <- record_args$sir_settings %||% .ferx_default_settings("sir_settings")
  if (resolved$resamples > resolved$samples) {
    # Name where each value came from: an argument the caller did not pass is
    # the fit's recorded value, or the engine default on a fit with none.
    source_of <- function(arg) {
      if (!is.null(arg)) ""
      else if (!is.null(fit$sir_settings)) ", recorded on the fit"
      else ", the default"
    }
    stop("`sir_resamples` (", resolved$resamples, source_of(sir_resamples),
         ") must be <= `sir_samples` (", resolved$samples, source_of(sir_samples), ")")
  }

  # Build the flat eta_hats matrix. `ebe_etas` is a data frame: first column
  # ID, subsequent columns are one per ETA. We strip ID and flatten row-major
  # so each subject contributes `n_eta` contiguous values.
  ebes <- fit$ebe_etas
  n_eta <- nrow(fit$omega)
  if (is.null(n_eta)) n_eta <- 0L

  if (n_eta == 0L) {
    # Fixed-effects-only (naive-pooled) fit - ferx-core #989. No inner
    # empirical-Bayes problem means no EBEs to warm-start from, and
    # `fit$ebe_etas` is NULL rather than an empty data frame. SIR itself still
    # applies: it resamples the theta/sigma block, which is the whole parameter
    # vector here.
    eta_hats_flat <- numeric(0)
  } else {
    if (is.null(ebes) || nrow(ebes) == 0L) {
      stop(
        "ferx_sir: fit$ebe_etas is empty, but fit$omega is ", n_eta, "x", n_eta,
        " so this fit should carry per-subject EBEs. Cannot warm-start the ",
        "inner loop."
      )
    }
    eta_cols <- setdiff(names(ebes), c("ID", "ofv_contribution", "n_obs"))
    if (length(eta_cols) != n_eta) {
      stop(
        "ferx_sir: fit$ebe_etas has ", length(eta_cols), " ETA columns (",
        paste(eta_cols, collapse = ", "),
        ") but fit$omega is ", n_eta, "x", n_eta, ". ",
        "The EBE table and the omega matrix must agree on n_eta - was ",
        "this fit object hand-edited or assembled from incompatible parts?"
      )
    }
    eta_mat <- as.matrix(ebes[, eta_cols, drop = FALSE])
    storage.mode(eta_mat) <- "double"
    eta_hats_flat <- as.numeric(t(eta_mat))  # row-major
  }
  # The fit's own subject IDs, in fit order (#468); their count is the
  # subject count, and in the n_eta > 0 branch they come from the same
  # `ebe_etas` rows as `eta_mat`.
  subject_ids <- .ferx_fit_subject_ids(fit, "ferx_sir")

  # The Rust binding wants row-major matrices and treats empty hash strings
  # as "no integrity check needed". Pass the recorded hashes through; the
  # Rust wrapper enforces equality when non-empty.
  #
  # Hash plumbing has three states we need to handle:
  #   1. Non-empty hex string: forward it; Rust enforces equality.
  #   2. NULL (older binary that didn't populate the field): forward "";
  #      Rust skips the check. (This is the "as-given" semantics.)
  #   3. NA_character_ (Rust ran but sha256_file failed; the post-fit
  #      `.ok()` in api.rs converts the Err to None, which the R wrapper
  #      stores as NA): we MUST coerce to "" here. `%||%` only catches
  #      NULL, and NA at the FFI boundary stringifies to "NA", which
  #      compares unequal to any real digest and would trigger a
  #      spurious "hash mismatch" error.
  # (`.ferx_hash_arg()` in zzz.R implements the NULL/NA -> "" coercion.)
  model_hash_arg <- .ferx_hash_arg(fit$model_hash)
  data_hash_arg <- .ferx_hash_arg(fit$data_hash)
  if (!nzchar(model_hash_arg) || !nzchar(data_hash_arg)) {
    warning(
      "ferx_sir: one or more file hashes are missing on the fit; ",
      "the integrity check will be skipped for the affected file(s). ",
      "Re-fit (or load a newer .fitrx) to enable hash verification."
    )
  }

  omega_flat <- as.numeric(t(fit$omega))
  cov_flat <- as.numeric(t(fit$cov_matrix))

  # The fitted kappa (IOV) covariance: without it the engine would rebuild the
  # fit at the model file's initial kappa and resample around the wrong centre,
  # so the glue refuses a kappa fit that lacks it (#465).
  iov_args <- .ferx_omega_iov_args(fit)
  binding_args <- .ferx_fit_binding_args(fit)
  reader_args <- .ferx_fit_reader_args(fit)
  .ferx_refuse_unrecorded_selection(
    fit, "ferx_sir",
    "Refit with `ferx_fit(..., sir = TRUE)`, or refit and call ferx_sir() on the new fit."
  )
  # A refusal reaches the caller as `ferx_engine_error` with ferx-core's code
  # when the engine gave it one (#498); see `.ferx_engine_call()`.
  raw <- .ferx_engine_call(ferx_rust_sir(
    model_path = model_path,
    data_path = data_path,
    model_hash = model_hash_arg,
    data_hash = data_hash_arg,
    ofv = as.numeric(fit$ofv),
    ofv_prior = .ferx_ofv_prior(fit),
    # `fit$interaction` is not plumbed to R, so reading it gave FALSE on every
    # fit and a FOCEI fit was resampled under the FOCE inner loop.
    interaction = .ferx_fit_interaction(fit),
    theta = as.numeric(fit$theta),
    omega_flat = omega_flat,
    omega_dim = nrow(fit$omega),
    sigma = as.numeric(fit$sigma),
    omega_iov_flat = iov_args$omega_iov_flat,
    omega_iov_dim = iov_args$omega_iov_dim,
    residual_rho = .ferx_residual_rho_vec(fit),
    cov_matrix_flat = cov_flat,
    cov_matrix_dim = nrow(fit$cov_matrix),
    eta_hats_flat = eta_hats_flat,
    subject_ids = subject_ids,
    # The occasion rule the fit derived its occasions with; without it the
    # engine falls back to the model file's and refuses a rule passed only
    # through `settings =` (#512).
    iov_occasion = as.character(fit$iov_occasion %||% ""),
    # The reader settings the fit read its data with and the fingerprint of
    # what it read: the engine re-reads the fit's own rows and verifies them
    # (#462).
    reader_settings = reader_args$reader_settings,
    population_fingerprint = reader_args$population_fingerprint,
    # How the fit was scored, and its exact packed estimate (#472).
    scoring_settings = record_args$scoring_settings,
    sir_settings = record_args$sir_settings,
    packed_estimate = record_args$packed_estimate,
    verbose = isTRUE(verbose),
    fit_bindings = binding_args$fit_bindings
  ), model_path, data_path)

  # All error paths inside `ferx_rust_sir` raise an R condition and never
  # return `NULL`, so we don't need to test for that here.

  # Merge results onto the fit object. Mirrors the post-processing
  # `ferx_fit()` applies to its inline SIR output so downstream code sees
  # the same shapes regardless of which entry point produced them.
  #
  # Non-finite ESS shouldn't happen in the standalone path (the engine
  # would have thrown an error before returning), but if it does, warn
  # before discarding it so users know the run was degenerate.
  if (!is.finite(raw$sir_ess)) {
    warning(
      "ferx_sir: effective sample size is not finite (got ", raw$sir_ess,
      "). The proposal distribution may be a poor match for the true ",
      "uncertainty - increase `sir_samples`, or reconsider the underlying fit."
    )
    fit$sir_ess <- NULL
  } else {
    fit$sir_ess <- raw$sir_ess
  }

  ci <- .ferx_sir_ci_matrices(
    raw$sir_ci_theta, raw$sir_ci_omega, raw$sir_ci_sigma,
    names(fit$theta), fit$eta_names, fit$sigma_names
  )
  fit$sir_ci_theta <- ci$sir_ci_theta
  fit$sir_ci_omega <- ci$sir_ci_omega
  fit$sir_ci_sigma <- ci$sir_ci_sigma
  fit$sir_ci_kappa <- .ferx_sir_ci_kappa(raw$sir_ci_kappa, fit$kappa_names)
  # The seed and settings this run used (ferx-core #1767), not the input fit's.
  fit$sir_seed_used <- raw$sir_seed_used
  fit$sir_settings <- raw$sir_settings

  # A run that did not keep its draws drops the input fit's: they belong to
  # an earlier run, and ferx_save_fit() would otherwise write them beside
  # this run's intervals (#549).
  keep <- isTRUE(fit$sir_settings$keep_samples)
  fit$sir_resamples <- if (keep) as.numeric(raw$sir_resamples)
  fit$sir_resamples_n <- if (keep) as.integer(raw$sir_resamples_n)
  fit$sir_resamples_dim <- if (keep) as.integer(raw$sir_resamples_dim)

  # Append any SIR-step warnings to the structured warnings table. The engine
  # can emit warnings during the SIR run (e.g. ESS collapse, proposal issues);
  # they would otherwise be silently dropped since only ferx_fit() calls the
  # assembler. We append rather than replace so pre-SIR warnings are preserved.
  #
  # Two sources feed the table: rows derived R-side by the assembler, and the
  # engine's own SIR-step messages, which the binding returns as a flat
  # character vector with no severity/category (ferx-core#1021 added the
  # proposal-conditioning diagnostics that make this matter: a rank-deficient
  # proposal, or one shrunk to keep draws inside the parameter bounds, qualifies
  # the CIs this very call just wrote onto the fit). Fold them in as
  # sir-category rows, exactly as ferx_covariance() does for its own.
  #
  # This run's `SIR: ` lines replace an earlier SIR run's, and its success
  # replaces an in-fit `SIR failed: ` line, as ferx-core's own run_sir does
  # (#1723): otherwise a re-run with `sir_scale = "natural"` - the low-ESS
  # warning's own advice - would print the old run's "ESS 3.5" beside the new
  # intervals. `SIR requested ...` and `SIR fallback ...` lines describe a
  # different step and stay.
  stale_sir <- function(msg) startsWith(msg, "SIR: ") | startsWith(msg, "SIR failed: ")
  if (length(fit$warnings) > 0L) {
    fit$warnings <- fit$warnings[!stale_sir(as.character(fit$warnings))]
  }
  if (is.data.frame(fit$warnings_structured) && nrow(fit$warnings_structured) > 0L) {
    fit$warnings_structured <- fit$warnings_structured[
      !stale_sir(as.character(fit$warnings_structured$message)), , drop = FALSE]
  }
  if (length(raw$warnings) > 0L) {
    fit$warnings <- unique(c(fit$warnings, as.character(raw$warnings)))
  }
  assembled <- .ferx_assemble_structured_warnings(raw, fit)
  engine_rows <- if (length(raw$warnings) > 0L) {
    data.frame(
      severity      = "warning",
      category      = "sir",
      message       = as.character(raw$warnings),
      source_method = "",
      stringsAsFactors = FALSE
    )
  } else {
    assembled[0, , drop = FALSE]
  }
  sir_warnings_df <- rbind(assembled, engine_rows)
  if (nrow(sir_warnings_df) > 0L) {
    existing <- fit$warnings_structured
    fit$warnings_structured <- if (is.data.frame(existing) && nrow(existing) > 0L)
      unique(rbind(existing, sir_warnings_df))
    else
      unique(sir_warnings_df)
  }

  fit
}
