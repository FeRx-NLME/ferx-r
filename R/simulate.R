#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error. When a \code{fit} is supplied, the fitted theta, omega,
#' and sigma replace the model file's initial values — which is the usual flow
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
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Simulate at the model file's initial estimates:
#' sim0 <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42)
#'
#' # Fit, then simulate at the fitted estimates (typical VPC flow):
#' fit <- ferx_fit(ex$model, ex$data)
#' sim <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42, fit = fit)
#' }
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
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' preds <- ferx_predict(ex$model, ex$data)
#'
#' fit <- ferx_fit(ex$model, ex$data)
#' preds_fitted <- ferx_predict(ex$model, ex$data, fit = fit)
#' }
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
