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
