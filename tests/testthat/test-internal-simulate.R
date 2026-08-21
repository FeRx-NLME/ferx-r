
# ---- header from test-simulate-validate.R ----
# Tests for the internal FFI-payload validators in simulate.R. These guard the
# shapes pulled out of a ferx_fit before handing them to the Rust simulator.

.validate_params      <- getFromNamespace("validate_fit_for_params",      "ferx")
.validate_uncertainty <- getFromNamespace("validate_fit_for_uncertainty", "ferx")










test_that("validate_fit_for_params requires theta/omega/sigma", {
  expect_error(.validate_params(list()), "theta, omega, and sigma")
})
test_that("validate_fit_for_params requires a square omega matrix", {
  expect_error(
    .validate_params(list(theta = 1, omega = 1, sigma = 1)),  # omega not a matrix
    "square matrix"
  )
})
test_that("validate_fit_for_params flattens a valid fit row-major", {
  fit <- list(theta = c(1, 2), omega = matrix(1:4, 2, 2), sigma = 0.05)
  out <- .validate_params(fit)
  expect_identical(out$omega_flat, c(1, 3, 2, 4))
})

# ferx-core #1019: the fitted IOV covariance must reach the engine, or the
# occasion-aware kappa draw has nothing to sample from (it used to panic).
test_that("validate_fit_for_params reports no IOV for a fit without omega_iov", {
  fit <- list(theta = 1, omega = matrix(0.1, 1, 1), sigma = 0.05)
  out <- .validate_params(fit)
  expect_identical(out$omega_iov_flat, numeric(0))
  expect_identical(out$omega_iov_dim, 0L)
})
test_that("validate_fit_for_params flattens omega_iov row-major", {
  fit <- list(
    theta = 1, omega = matrix(0.1, 1, 1), sigma = 0.05,
    omega_iov = matrix(1:4, 2, 2)  # column-major fill: rows are (1, 3), (2, 4)
  )
  out <- .validate_params(fit)
  expect_identical(out$omega_iov_flat, c(1, 3, 2, 4))
  expect_identical(out$omega_iov_dim, 2L)
})
test_that("validate_fit_for_params requires a square omega_iov", {
  fit <- list(
    theta = 1, omega = matrix(0.1, 1, 1), sigma = 0.05,
    omega_iov = matrix(1:6, 2, 3)
  )
  expect_error(.validate_params(fit), "square matrix")
})
