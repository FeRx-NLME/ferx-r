#' Run SIR against an existing fit
#'
#' Run a Sampling Importance Resampling (SIR) uncertainty step against a
#' fit that was produced earlier — useful when the original fit was
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
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE)
#' fit <- ferx_sir(fit, sir_samples = 2000, sir_resamples = 500, sir_seed = 42)
#' fit$sir_ess
#' fit$sir_ci_theta
#' }
#'
#' @seealso [ferx_fit()] for the inline SIR option (`sir = TRUE`),
#'   [ferx_simulate_with_uncertainty()] for downstream consumption of the
#'   retained resamples.
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
  eta_mat <- as.matrix(ebes[, eta_cols, drop = FALSE])
  storage.mode(eta_mat) <- "double"
  eta_hats_flat <- as.numeric(t(eta_mat))  # row-major

  # The Rust binding wants row-major matrices and treats empty hash strings
  # as "no integrity check needed". Pass the recorded hashes through; the
  # Rust wrapper enforces equality.
  omega_flat <- as.numeric(t(fit$omega))
  cov_flat <- as.numeric(t(fit$cov_matrix))

  raw <- ferx_rust_sir(
    model_path = model_path,
    data_path = data_path,
    model_hash = fit$model_hash %||% "",
    data_hash = fit$data_hash %||% "",
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

  if (is.null(raw) || length(raw) == 0L) {
    stop("ferx_sir: backend returned no result (see messages above).")
  }

  # Merge results onto the fit object. Mirrors the post-processing
  # `ferx_fit()` applies to its inline SIR output so downstream code sees
  # the same shapes regardless of which entry point produced them.
  fit$sir_ess <- if (is.finite(raw$sir_ess)) raw$sir_ess else NULL

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
  n_eta <- nrow(fit$omega)
  eta_row_names <- paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
  fit$sir_ci_omega <- reshape_ci(raw$sir_ci_omega, eta_row_names)
  sig_names <- paste0("SIGMA(", seq_along(fit$sigma), ")")
  fit$sir_ci_sigma <- reshape_ci(raw$sir_ci_sigma, sig_names)

  if (isTRUE(sir_keep_samples)) {
    fit$sir_resamples <- as.numeric(raw$sir_resamples)
    fit$sir_resamples_n <- as.integer(raw$sir_resamples_n)
    fit$sir_resamples_dim <- as.integer(raw$sir_resamples_dim)
  }

  fit
}
