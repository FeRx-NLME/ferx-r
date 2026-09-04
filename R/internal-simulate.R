# Internal: pull theta/omega/sigma/omega_iov out of a ferx_fit result for FFI.
# Flattens the matrices row-major for the Rust side.
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
  # Fitted IOV (kappa) covariance. A `kappa` model needs it downstream: the engine
  # draws (simulate) or conditions on (predict / npde) one kappa vector per occasion
  # from it, and dropping it panicked ferx-core's simulate path (#1019). NULL for a
  # non-IOV fit, which maps to an empty vector + dim 0 ("no IOV") on the Rust side.
  omega_iov <- fit$omega_iov
  if (!is.null(omega_iov) &&
      (!is.matrix(omega_iov) || nrow(omega_iov) != ncol(omega_iov))) {
    stop("`fit$omega_iov` must be a square matrix (or NULL for a model without IOV).")
  }
  list(
    theta = theta,
    omega_flat = as.numeric(t(omega)),  # row-major
    omega_dim = as.integer(nrow(omega)),
    sigma = sigma,
    omega_iov_flat = if (is.null(omega_iov)) numeric(0) else as.numeric(t(omega_iov)),
    omega_iov_dim = if (is.null(omega_iov)) 0L else as.integer(nrow(omega_iov)),
    # Fitted `block_sigma` residual correlations, in model declaration order.
    # A plain (non-FIX) block estimates rho, so passing them is what keeps the
    # engine from rebuilding this fit at the model file's declared correlation;
    # empty for a model that declares none.
    residual_rho = .ferx_residual_rho_vec(fit)
  )
}

# The fitted residual correlations as a bare numeric vector for the FFI.
# Reads the tidy frame ferx_fit() builds, and stays empty for a fit that
# predates it (an older .fitrx bundle) so the engine keeps the declared values.
.ferx_residual_rho_vec <- function(fit) {
  rc <- fit$residual_correlations
  if (is.data.frame(rc) && nrow(rc) > 0L) return(as.numeric(rc$rho))
  numeric(0)
}
