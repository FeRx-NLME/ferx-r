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
#' [ferx_save_fit()] / [ferx_load_fit()] round trip, a normal proposal
#' (`sir_df = Inf`) included (ferx-core#1819). An argument you pass
#' replaces that one setting and keeps the rest. The returned fit's
#' `sir_settings` records the run it made. The scoring half includes the
#' gradient route the fit took (`gradient_method`), so a fit made with
#' `gradient = "fd"` stays on finite differences here.
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
#' @param sir_seed Integer RNG seed for reproducibility, a whole number from 0
#'   to 2^53 (it need not fit an R integer). `NULL` takes the
#'   seed the fit's SIR recorded, else the engine's built-in seed (12345).
#' @param sir_keep_samples When `TRUE`, retain the resampled packed
#'   parameter vectors on the returned fit. Required for
#'   [ferx_simulate_with_uncertainty()] with `method = "sir"`. `NULL` takes
#'   the fit's recorded value, else `FALSE`.
#' @param threads Worker threads the step runs on. `NULL` (the default)
#'   takes the width the fit ran on, `fit$n_threads_used`; `0` asks for the
#'   engine default; a positive integer pins the width. It does not change
#'   the result, only the resources the step uses. The returned fit's
#'   `n_threads_used` reports the width this step ran on.
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
#' @return A new fit with new bytes, the input `fit` being left as it was: the
#'   step runs on `fit$bundle`, so a change made to the fit's list by hand does
#'   not reach it, and a fit with no bundle, such as one from an `.rds` written
#'   by an older ferx, is refused with an error of class
#'   `ferx_no_bundle_error`: refit it. The new fit has `sir_ess`, `sir_ci_theta`,
#'   `sir_ci_omega`, `sir_ci_sigma`, `sir_ci_kappa` (one row per IOV kappa
#'   variance, named by `kappa_names`; NULL without IOV), `sir_seed_used`
#'   (the seed this run resampled with: `sir_seed`, or 12345 when it is
#'   `NULL`), `sir_settings` (the settings this run used; see [ferx_fit()]),
#'   `n_threads_used` (the width this run ran on), and (when the run kept
#'   them) `sir_resamples` / `sir_resamples_n` / `sir_resamples_dim`.
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
                     sir_df = NULL,
                     threads = NULL) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }
  threads_arg <- .ferx_step_threads(threads, fit)

  # Validate the SIR knobs up-front so user errors surface as clean R
  # conditions rather than as the engine's "weighted sampler failed"
  # or the binding's i32-clamped 0 after a bad cast.
  # The range is tested before the whole-number test: `as.integer()` of a
  # count beyond `.Machine$integer.max` is NA, and an NA in `if` stops with
  # R's "missing value where TRUE/FALSE needed" instead of ours (#536).
  is_whole <- function(x, lower, upper) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower &&
      x <= upper && x == floor(x)
  }
  if (!is.null(sir_samples) && !is_whole(sir_samples, 1, .Machine$integer.max)) {
    stop("`sir_samples` must be NULL or a single positive integer (got: ",
         paste(format(sir_samples), collapse = ", "), ").")
  }
  if (!is.null(sir_resamples) && !is_whole(sir_resamples, 1, .Machine$integer.max)) {
    stop("`sir_resamples` must be NULL or a single positive integer (got: ",
         paste(format(sir_resamples), collapse = ", "), ").")
  }
  # The engine's seed is a u64, and the seed travels to it as a double, so any
  # whole number a double holds exactly (up to 2^53) is the seed that runs.
  if (!is.null(sir_seed) && !is_whole(sir_seed, 0, 2^53)) {
    stop("`sir_seed` must be NULL or a single non-negative integer no larger ",
         "than 2^53 (got: ", paste(format(sir_seed), collapse = ", "), ").")
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
  bundle <- .ferx_fit_bundle(fit, "ferx_sir")
  # Before the model-file checks, so a loaded fit without a layout is refused
  # for that reason even when its bundle had no model file (#531).
  .ferx_refuse_unavailable_levels(fit)

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

  data_path <- .ferx_fit_data_path(fit)
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

  # The engine checks the files against the hashes the bytes record; a fit that
  # records none is not checked (`.ferx_hash_arg()` in zzz.R reads NULL and NA
  # alike).
  if (!nzchar(.ferx_hash_arg(fit$model_hash)) || !nzchar(.ferx_hash_arg(fit$data_hash))) {
    warning(
      "ferx_sir: one or more file hashes are missing on the fit; ",
      "the integrity check will be skipped for the affected file(s). ",
      "Re-fit (or load a newer .fitrx) to enable hash verification."
    )
  }

  .ferx_refuse_unrecorded_selection(
    fit, "ferx_sir",
    "Refit with `ferx_fit(..., sir = TRUE)`, or refit and call ferx_sir() on the new fit."
  )
  # A refusal reaches the caller as `ferx_engine_error` with ferx-core's code
  # when the engine gave it one (#498); see `.ferx_engine_call()`.
  new_bundle <- .ferx_engine_call(ferx_rust_sir(
    bundle = bundle,
    model_path = model_path,
    data_path = data_path,
    sir_settings = record_args$sir_settings,
    threads = threads_arg,
    verbose = isTRUE(verbose)
  ), model_path, data_path)

  # A new fit from the new bytes: `fit` is left as it was.
  new_fit <- .ferx_fit_with_bundle(fit, new_bundle)
  ess <- new_fit$sir_ess
  if (length(ess) == 1L && !is.finite(ess)) {
    warning(
      "ferx_sir: effective sample size is not finite (got ", ess,
      "). The proposal distribution may be a poor match for the true ",
      "uncertainty - increase `sir_samples`, or reconsider the underlying fit."
    )
  }
  new_fit
}
