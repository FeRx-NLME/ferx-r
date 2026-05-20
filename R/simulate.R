#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error. When a \code{fit} is supplied, the fitted theta, omega,
#' and sigma replace the model file's initial values - which is the usual flow
#' after \code{\link{ferx_fit}} (e.g. for posterior-predictive checks or VPCs).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure: doses, obs times)
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#' @param fit Optional \code{ferx_fit} result. When provided, simulation uses
#'   \code{fit$theta}, \code{fit$omega}, and \code{fit$sigma} instead of the
#'   model file's initial values.
#'
#' @return A data.frame with columns: SIM, ID, TIME, IPRED, DV_SIM
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' sim <- ferx_simulate(ex$model, ex$data, n_sim = 10L, seed = 1L, fit = fit)
#' head(sim)
#'
#' @export
ferx_simulate <- function(model, data, n_sim = 1L, seed = 42L, fit = NULL) {
  stopifnot(file.exists(model), file.exists(data))

  if (is.null(fit)) {
    return(ferx_rust_simulate(
      model_path = normalizePath(model),
      data_path = normalizePath(data),
      n_sim = as.integer(n_sim),
      seed = as.integer(seed)
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_simulate_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed)
  )
}

#' Population predictions from a NLME model
#'
#' Computes population-level predictions (eta = 0) for all subjects.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV
#' @param fit Optional \code{ferx_fit} result. When provided, predictions use
#'   \code{fit$theta} instead of the model file's initial estimate for theta.
#'
#' @return A data.frame with columns: ID, TIME, PRED
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' preds <- ferx_predict(ex$model, ex$data, fit = fit)
#' head(preds)
#'
#' @export
ferx_predict <- function(model, data, fit = NULL) {
  stopifnot(file.exists(model), file.exists(data))

  if (is.null(fit)) {
    return(ferx_rust_predict(
      model_path = normalizePath(model),
      data_path = normalizePath(data)
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_predict_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma
  )
}

# Internal: pull theta/omega/sigma out of a ferx_fit result for FFI.
# Flattens omega row-major for the Rust side.
validate_fit_for_params <- function(fit) {
  if (!is.list(fit) || is.null(fit$theta) || is.null(fit$omega) || is.null(fit$sigma)) {
    stop("`fit` must be a ferx_fit result with theta, omega, and sigma components.")
  }
  theta <- as.numeric(fit$theta)
  sigma <- as.numeric(fit$sigma)
  omega <- fit$omega
  if (!is.matrix(omega) || nrow(omega) != ncol(omega)) {
    stop("`fit$omega` must be a square matrix.")
  }
  list(
    theta = theta,
    omega_flat = as.numeric(t(omega)),  # row-major
    omega_dim = as.integer(nrow(omega)),
    sigma = sigma
  )
}

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
#' @param data Path to a NONMEM-format CSV (provides population structure)
#' @param fit A \code{ferx_fit} result. Must carry either \code{cov_matrix}
#'   (asymptotic) or \code{sir_resamples} (SIR) depending on \code{method}.
#' @param n_uncertainty_draws Number of parameter sets to draw from the
#'   uncertainty distribution
#' @param n_sim_per_draw Number of eta/eps replicates per parameter draw
#' @param method Either \code{"asymptotic"} (default) or \code{"sir"}
#' @param seed Random seed for reproducibility
#'
#' @return A data.frame with columns: DRAW, SIM, ID, TIME, IPRED, DV_SIM.
#'   Row count: \code{n_uncertainty_draws * n_sim_per_draw * n_obs}.
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

  ferx_rust_simulate_with_uncertainty(
    model_path           = normalizePath(model),
    data_path            = normalizePath(data),
    theta                = fit_pieces$theta,
    omega_flat           = fit_pieces$omega_flat,
    omega_dim            = fit_pieces$omega_dim,
    sigma                = fit_pieces$sigma,
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
