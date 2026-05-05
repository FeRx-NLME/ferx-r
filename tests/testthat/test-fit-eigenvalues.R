# Tests for the R-layer handling of eigenvalues and condition_number fields
# that are pre-computed by the Rust backend and exposed via the FFI.
#
# The Rust side (ferx-nlme) is tested for computation correctness
# (fixed-param exclusion, Inf for non-positive eigenvalue) in api.rs
# tests_cov_diagnostics. These tests cover the R-side sentinel-to-NULL
# conversion and the high-condition-number warning injection.

make_fit <- function(cov_eigenvalues = numeric(0),
                     cov_condition_number = NaN,
                     warnings = character(0)) {
  # Minimal ferx_fit-like list with just the fields the R layer touches.
  structure(
    list(
      cov_eigenvalues    = cov_eigenvalues,
      cov_condition_number = cov_condition_number,
      warnings           = warnings,
      # other fields the R conversion code may touch
      cov_matrix         = NULL,
      cov_matrix_dim     = NULL,
      ebe_etas           = NULL,
      covariance_status  = "computed",
      omega_iov_dim      = NULL,
      omega_iov          = NULL,
      kappa_names        = character(0),
      kappa_fixed        = logical(0),
      se_kappa           = numeric(0),
      shrinkage_kappa    = numeric(0),
      ebe_kappas         = list(),
      model_structure    = NULL
    ),
    class = "ferx_fit_raw"  # not a real ferx_fit — we test the helper directly
  )
}

.apply_cov_sentinels <- function(result) {
  result$eigenvalues <- if (length(result$cov_eigenvalues) == 0L) NULL else result$cov_eigenvalues
  result$condition_number <- if (is.nan(result$cov_condition_number)) NULL else result$cov_condition_number
  result$cov_eigenvalues <- NULL
  result$cov_condition_number <- NULL
  if (!is.null(result$condition_number) && is.finite(result$condition_number) &&
        result$condition_number > 1000) {
    result$warnings <- c(
      result$warnings,
      sprintf(
        "High condition number (%.1f) — parameter space may be ill-conditioned",
        result$condition_number
      )
    )
  }
  result
}

test_that("empty cov_eigenvalues and NaN cov_condition_number become NULL", {
  r <- .apply_cov_sentinels(make_fit())
  expect_null(r$eigenvalues)
  expect_null(r$condition_number)
})

test_that("populated eigenvalues and finite condition_number are forwarded as-is", {
  r <- .apply_cov_sentinels(make_fit(
    cov_eigenvalues = c(1.5, 0.8, 0.3),
    cov_condition_number = 5.0
  ))
  expect_equal(r$eigenvalues, c(1.5, 0.8, 0.3))
  expect_equal(r$condition_number, 5.0)
})

test_that("Inf condition_number is forwarded as Inf (not NULL)", {
  r <- .apply_cov_sentinels(make_fit(
    cov_eigenvalues = c(2.5, -0.5),
    cov_condition_number = Inf
  ))
  expect_true(is.infinite(r$condition_number))
  expect_false(is.null(r$condition_number))
})

test_that("condition_number > 1000 appends a warning", {
  r <- .apply_cov_sentinels(make_fit(
    cov_eigenvalues = c(1500.0, 1.5),
    cov_condition_number = 1000.1
  ))
  expect_true(any(grepl("High condition number", r$warnings)))
})

test_that("condition_number <= 1000 does not append a warning", {
  r <- .apply_cov_sentinels(make_fit(
    cov_eigenvalues = c(1.2, 1.0),
    cov_condition_number = 1000.0
  ))
  expect_false(any(grepl("High condition number", r$warnings)))
})

test_that("Inf condition_number does not trigger the high-condition-number warning", {
  r <- .apply_cov_sentinels(make_fit(
    cov_eigenvalues = c(2.5, -0.5),
    cov_condition_number = Inf
  ))
  expect_false(any(grepl("High condition number", r$warnings)))
})
