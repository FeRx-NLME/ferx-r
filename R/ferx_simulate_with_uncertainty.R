#' Simulate with parameter-uncertainty propagation
#'
#' For each parameter set drawn from the uncertainty distribution, the
#' per-subject random-effect / residual-error simulator runs
#' \code{n_sim_per_draw} times. The result includes both individual
#' variability (etas, epsilons) and parameter uncertainty - useful for
#' uncertainty-aware VPCs, dose-recommendation intervals, and any analysis
#' where treating the ML estimates as fixed would understate variability.
#'
#' Two uncertainty sources are supported:
#' \describe{
#'   \item{\code{method = "asymptotic"}}{Multivariate normal around the ML
#'     estimate in the engine's packed (log-theta, Cholesky-omega, log-sigma)
#'     parameter space, using \code{fit$cov_matrix}. Requires \code{fit}
#'     to come from a \code{ferx_fit()} call with \code{covariance = TRUE}.}
#'   \item{\code{method = "sir"}}{Sample with replacement from
#'     \code{fit$sir_resamples}. Requires the fit to have been run with
#'     \code{sir = TRUE} and \code{sir_keep_samples = TRUE} (passed via
#'     \code{settings}).}
#' }
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure).
#'   The \code{DV} column may be left empty (\code{.} / \code{NA}) on the
#'   sampling rows - the DV is what the simulation produces, so an empty cell
#'   means "simulate here" (a placeholder value is not needed). Rows marked
#'   \code{MDV = 1} are excluded, as always. Kept empty-DV records are counted
#'   in the \code{simulation_warnings} attribute and re-emitted as an R warning:
#'   \code{ferx_fit()} skips those same records, so simulated rows at those
#'   times have no counterpart in a fit's \code{sdtab} (do not overlay the two,
#'   e.g. in a VPC).
#' @param fit A \code{ferx_fit} result. Must carry either \code{cov_matrix}
#'   (asymptotic) or \code{sir_resamples} (SIR) depending on \code{method}.
#' @param n_uncertainty_draws Number of parameter sets to draw from the
#'   uncertainty distribution
#' @param n_sim_per_draw Number of eta/eps replicates per parameter draw
#' @param method Either \code{"asymptotic"} (default) or \code{"sir"}
#' @param seed Random seed for reproducibility
#'
#' @return A data.frame with columns: DRAW, SIM, ID, TIME, CMT, IPRED, DV_SIM,
#'   OBSERVED. \code{CMT} is the observation compartment; \code{OBSERVED} is
#'   \code{NA} throughout (this path is Gaussian-only - a drug-driven ODE-TTE
#'   endpoint is not supported here). Row count:
#'   \code{n_uncertainty_draws * n_sim_per_draw * n_obs}.
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE)
#'
#' # Asymptotic (default): fast, MVN draws around the ML estimate
#' sims <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10
#' )
#' head(sims)         # SIM, ID, TIME, IPRED, DV_SIM
#' range(sims$DV_SIM) # sanity check on simulated range
#'
#' # SIR method: requires sir = TRUE + sir_keep_samples = TRUE at fit time
#' fit_sir <- ferx_fit(ex$model, ex$data,
#'                     covariance = TRUE, sir = TRUE,
#'                     settings = list(sir_keep_samples = TRUE))
#' sims_sir <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit_sir,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10,
#'   method = "sir"
#' )
#'
#' # Summarise: 90% prediction interval per time point
#' pi90 <- aggregate(DV_SIM ~ TIME, data = sims,
#'                   FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
#' }
#'
#' @family simulation
#' @export
ferx_simulate_with_uncertainty <- function(model, data, fit,
                                           n_uncertainty_draws = 100L,
                                           n_sim_per_draw = 1L,
                                           method = c("asymptotic", "sir"),
                                           seed = 42L) {
  stopifnot(file.exists(model), file.exists(data))
  method <- match.arg(method)

  # Check the cheap, user-facing arg constraints before the fit-object
  # validation so callers passing `n_uncertainty_draws = 0` see the obvious
  # error rather than a follow-on complaint about `cov_matrix` / `sir_resamples`.
  if (n_uncertainty_draws < 1L) {
    stop("`n_uncertainty_draws` must be >= 1.")
  }
  if (n_sim_per_draw < 1L) {
    stop("`n_sim_per_draw` must be >= 1.")
  }

  fit_pieces <- validate_fit_for_params(fit)
  unc_pieces <- validate_fit_for_uncertainty(fit, method)

  res <- ferx_rust_simulate_with_uncertainty(
    model_path           = normalizePath(model),
    data_path            = normalizePath(data),
    theta                = fit_pieces$theta,
    omega_flat           = fit_pieces$omega_flat,
    omega_dim            = fit_pieces$omega_dim,
    sigma                = fit_pieces$sigma,
    omega_iov_flat       = fit_pieces$omega_iov_flat,
    omega_iov_dim        = fit_pieces$omega_iov_dim,
    method               = method,
    cov_matrix_flat      = unc_pieces$cov_matrix_flat,
    cov_matrix_dim       = unc_pieces$cov_matrix_dim,
    sir_resamples_flat   = unc_pieces$sir_resamples_flat,
    sir_resamples_n      = unc_pieces$sir_resamples_n,
    sir_resamples_dim    = unc_pieces$sir_resamples_dim,
    n_uncertainty_draws  = as.integer(n_uncertainty_draws),
    n_sim_per_draw       = as.integer(n_sim_per_draw),
    seed                 = as.integer(seed)
  )

  # Same `simulation_warnings` channel `ferx_simulate()` uses - here it carries
  # the design-point count (a kept empty-DV record; see the `data` note above),
  # which otherwise diverges silently from what `ferx_fit()` scored.
  .ferx_surface_sim_warnings(res, "ferx_simulate_with_uncertainty")
}

# Internal: pull the uncertainty payload out of a ferx_fit result, validate it
# matches the requested method, and return flat representations for FFI. The
# empty branches give the Rust side something to ignore - it switches on
# `method` and only inspects the relevant arrays.
validate_fit_for_uncertainty <- function(fit, method) {
  if (method == "asymptotic") {
    cov <- fit$cov_matrix
    if (is.null(cov) || length(cov) == 0L) {
      stop("`fit$cov_matrix` is empty - re-fit with `covariance = TRUE` for ",
           "asymptotic uncertainty.")
    }
    # `cov` is a square R matrix after `process_fit_result()`. Flatten
    # row-major to match the engine's `DMatrix::from_row_slice` reader on
    # the Rust side - same convention used for `omega` above. The
    # transpose is mathematically a no-op for any well-formed covariance
    # (which must be symmetric), but the explicit `t()` keeps the FFI
    # contract uniform with `omega_flat` and prevents a future footgun
    # if the Rust reader ever changes.
    if (!is.matrix(cov) || nrow(cov) != ncol(cov)) {
      stop("`fit$cov_matrix` must be a square matrix.")
    }
    list(
      cov_matrix_flat    = as.numeric(t(cov)),
      cov_matrix_dim     = as.integer(nrow(cov)),
      sir_resamples_flat = numeric(0),
      sir_resamples_n    = 0L,
      sir_resamples_dim  = 0L
    )
  } else {
    resamples <- fit$sir_resamples
    n <- fit$sir_resamples_n
    d <- fit$sir_resamples_dim
    if (is.null(resamples) || length(resamples) == 0L ||
        is.null(n) || n == 0L || is.null(d) || d == 0L) {
      stop("`fit$sir_resamples` is empty - re-fit with `sir = TRUE` and ",
           "`sir_keep_samples = TRUE` in `settings` for SIR uncertainty.")
    }
    list(
      cov_matrix_flat    = numeric(0),
      cov_matrix_dim     = 0L,
      sir_resamples_flat = as.numeric(resamples),
      sir_resamples_n    = as.integer(n),
      sir_resamples_dim  = as.integer(d)
    )
  }
}
