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
#' `sir_ess`, `sir_ci_theta`, `sir_ci_omega`, `sir_ci_sigma` (and, when
#' `sir_keep_samples = TRUE`, `sir_resamples` / `sir_resamples_n` /
#' `sir_resamples_dim`) populated.
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
#'   Default 1000.
#' @param sir_resamples Number of resampled vectors. Must be `<= sir_samples`.
#'   Default 250.
#' @param sir_seed Integer RNG seed for reproducibility. `NULL` (the default)
#'   uses the engine's built-in seed.
#' @param sir_keep_samples When `TRUE`, retain the resampled packed
#'   parameter vectors on the returned fit. Required for
#'   [ferx_simulate_with_uncertainty()] with `method = "sir"`. Default
#'   `FALSE`.
#' @param verbose When `TRUE`, the engine prints progress to stderr.
#'   Default `FALSE`.
#'
#' @return The input `fit`, augmented with `sir_ess`, `sir_ci_theta`,
#'   `sir_ci_omega`, `sir_ci_sigma`, and (when requested) `sir_resamples` /
#'   `sir_resamples_n` / `sir_resamples_dim`.
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
                     sir_samples = 1000L,
                     sir_resamples = 250L,
                     sir_seed = NULL,
                     sir_keep_samples = FALSE,
                     verbose = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }

  # Validate the SIR knobs up-front so user errors surface as clean R
  # conditions rather than as the engine's "weighted sampler failed"
  # or the binding's i32-clamped 0 after a bad cast.
  is_positive_count <- function(x) {
    length(x) == 1L && is.finite(x) && x > 0 && x == as.integer(x)
  }
  if (!is_positive_count(sir_samples)) {
    stop("`sir_samples` must be a single positive integer (got: ",
         paste(format(sir_samples), collapse = ", "), ").")
  }
  if (!is_positive_count(sir_resamples)) {
    stop("`sir_resamples` must be a single positive integer (got: ",
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

  if (sir_resamples > sir_samples) {
    stop("`sir_resamples` (", sir_resamples,
         ") must be <= `sir_samples` (", sir_samples, ")")
  }

  # Build the flat eta_hats matrix. `ebe_etas` is a data frame: first column
  # ID, subsequent columns are one per ETA. We strip ID and flatten row-major
  # so each subject contributes `n_eta` contiguous values.
  ebes <- fit$ebe_etas
  if (is.null(ebes) || nrow(ebes) == 0L) {
    stop(
      "ferx_sir: fit$ebe_etas is empty; cannot warm-start the inner loop. ",
      "This fit appears to be missing per-subject EBEs."
    )
  }
  eta_cols <- setdiff(names(ebes), c("ID", "ofv_contribution", "n_obs"))
  if (length(eta_cols) == 0L) {
    stop("ferx_sir: fit$ebe_etas has no ETA columns.")
  }
  n_eta <- nrow(fit$omega)
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

  raw <- ferx_rust_sir(
    model_path = model_path,
    data_path = data_path,
    model_hash = model_hash_arg,
    data_hash = data_hash_arg,
    ofv = as.numeric(fit$ofv),
    interaction = isTRUE(fit$interaction),
    theta = as.numeric(fit$theta),
    omega_flat = omega_flat,
    omega_dim = nrow(fit$omega),
    sigma = as.numeric(fit$sigma),
    cov_matrix_flat = cov_flat,
    cov_matrix_dim = nrow(fit$cov_matrix),
    eta_hats_flat = eta_hats_flat,
    n_subjects = nrow(eta_mat),
    sir_samples = as.integer(sir_samples),
    sir_resamples = as.integer(sir_resamples),
    sir_seed = if (is.null(sir_seed)) -1L else as.integer(sir_seed),
    sir_keep_samples = isTRUE(sir_keep_samples),
    verbose = isTRUE(verbose)
  )

  # All error paths inside `ferx_rust_sir` throw an R condition (via
  # `throw_r_error`) and never return `NULL`, so we don't need to test
  # for that here.

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

  reshape_ci <- function(v, row_names = NULL) {
    if (length(v) == 0L) return(NULL)
    matrix(
      v,
      ncol = 2L,
      byrow = TRUE,
      dimnames = list(row_names, c("lower", "upper"))
    )
  }
  fit$sir_ci_theta <- reshape_ci(raw$sir_ci_theta, names(fit$theta))
  en <- fit$eta_names
  eta_row_names <- if (!is.null(en) && length(en) == n_eta) en else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
  fit$sir_ci_omega <- reshape_ci(raw$sir_ci_omega, eta_row_names)
  sn <- fit$sigma_names
  sig_names <- if (!is.null(sn) && length(sn) == length(fit$sigma)) sn else paste0("SIGMA(", seq_along(fit$sigma), ")")
  fit$sir_ci_sigma <- reshape_ci(raw$sir_ci_sigma, sig_names)

  if (isTRUE(sir_keep_samples)) {
    fit$sir_resamples <- as.numeric(raw$sir_resamples)
    fit$sir_resamples_n <- as.integer(raw$sir_resamples_n)
    fit$sir_resamples_dim <- as.integer(raw$sir_resamples_dim)
  }

  # Append any SIR-step warnings to the structured warnings table. The engine
  # can emit warnings during the SIR run (e.g. ESS collapse, proposal issues);
  # they would otherwise be silently dropped since only ferx_fit() calls the
  # assembler. We append rather than replace so pre-SIR warnings are preserved.
  sir_warnings_df <- .ferx_assemble_structured_warnings(raw, fit)
  if (nrow(sir_warnings_df) > 0L) {
    existing <- fit$warnings_structured
    fit$warnings_structured <- if (is.data.frame(existing) && nrow(existing) > 0L)
      rbind(existing, sir_warnings_df)
    else
      sir_warnings_df
  }

  fit
}
