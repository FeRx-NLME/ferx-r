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
  fit <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.2), 2, 2),
              sigma = 0.05)
  out <- .validate_params(fit)
  expect_identical(out$omega_flat, c(0.1, 0, 0, 0.2))
})

test_that("validate_fit_for_uncertainty (asymptotic) needs a non-empty cov matrix", {
  expect_error(.validate_uncertainty(list(cov_matrix = NULL), "asymptotic"),
               "cov_matrix` is empty")
})

test_that("validate_fit_for_uncertainty (asymptotic) requires a square cov matrix", {
  expect_error(
    .validate_uncertainty(list(cov_matrix = matrix(1:6, 2, 3)), "asymptotic"),
    "square matrix"
  )
})

test_that("validate_fit_for_uncertainty (asymptotic) flattens a valid cov matrix", {
  cov <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  out <- .validate_uncertainty(list(cov_matrix = cov), "asymptotic")
  expect_identical(out$cov_matrix_dim, 2L)
  expect_identical(out$cov_matrix_flat, as.numeric(t(cov)))
})

test_that("validate_fit_for_uncertainty (SIR) errors when resamples are empty", {
  expect_error(.validate_uncertainty(list(), "sir"), "sir_resamples` is empty")
})

test_that("validate_fit_for_uncertainty (SIR) passes through resamples", {
  fit <- list(sir_resamples = c(1, 2, 3, 4), sir_resamples_n = 2L,
              sir_resamples_dim = 2L)
  out <- .validate_uncertainty(fit, "sir")
  expect_identical(out$sir_resamples_n, 2L)
  expect_identical(out$sir_resamples_flat, c(1, 2, 3, 4))
})
