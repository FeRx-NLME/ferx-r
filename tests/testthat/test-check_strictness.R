# A synthetic ferx_fit carrying only what check_strictness() reads. Defaults
# describe a clean fit: converged, covariance computed, nothing near a bound,
# well away from the initial estimates.
strict_fit <- function(...) {
  base <- list(
    converged = TRUE,
    covariance_status = "computed",
    cov_matrix = diag(2),
    cov_condition_number = 12,
    max_abs_correlation = 0.3,
    near_boundary = FALSE,
    boundary_estimate_message = "",
    stalled_at_init = FALSE,
    sir_ess = NA_real_
  )
  overrides <- list(...)
  base[names(overrides)] <- overrides
  structure(base, class = "ferx_fit")
}

test_that("a clean fit passes every default gate", {
  v <- check_strictness(strict_fit())
  expect_identical(names(v), c("passed", "failures", "skipped"))
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 0L)
})

test_that("a live fit is checkable and reports the same shape", {
  v <- check_strictness(warfarin_fit_cov())
  expect_type(v$passed, "logical")
  expect_length(v$passed, 1L)
  expect_type(v$failures, "character")
  expect_type(v$skipped, "character")
  expect_identical(v$passed, length(v$failures) == 0L)
})

test_that("convergence is gated, and the gate can be turned off", {
  fit <- strict_fit(converged = FALSE)
  v <- check_strictness(fit)
  expect_false(v$passed)
  expect_match(v$failures, "did not converge", all = FALSE)

  expect_true(check_strictness(fit, require_converged = FALSE)$passed)
})

test_that("require_covariance names why the step delivered nothing", {
  cases <- list(
    list(fit = strict_fit(covariance_status = "not_requested", cov_matrix = NULL,
                          cov_condition_number = NaN),
         pattern = "not requested"),
    list(fit = strict_fit(covariance_status = "failed", cov_matrix = NULL,
                          cov_condition_number = NaN),
         pattern = "covariance step failed"),
    list(fit = strict_fit(covariance_status = "computed", cov_matrix = NULL,
                          cov_condition_number = NaN),
         pattern = "stored no matrix"),
    list(fit = strict_fit(covariance_status = "sir_fallback", cov_matrix = NULL,
                          cov_condition_number = NaN, sir_ess = 1.2),
         pattern = "SIR fallback collapsed")
  )
  for (case in cases) {
    v <- check_strictness(case$fit, require_covariance = TRUE)
    expect_false(v$passed, info = case$pattern)
    expect_match(v$failures, case$pattern, all = FALSE, info = case$pattern)
    # Off by default: the same fit passes when uncertainty is not required.
    expect_true(check_strictness(case$fit)$passed, info = case$pattern)
  }
})

test_that("a SIR fallback that delivered intervals satisfies require_covariance", {
  fit <- strict_fit(covariance_status = "sir_fallback", cov_matrix = NULL,
                    cov_condition_number = NaN, max_abs_correlation = NaN,
                    sir_ess = 40)
  v <- check_strictness(fit, require_covariance = TRUE)
  expect_true(v$passed)
  # It stores no covariance matrix, so the two threshold gates have no input.
  expect_length(v$skipped, 2L)
  expect_match(v$skipped, "condition number", all = FALSE)
  expect_match(v$skipped, "parameter correlation", all = FALSE)
})

test_that("the status label is read in either casing", {
  # A fresh fit reports "computed"; one rebuilt from a bundle reports
  # "Computed". Both have to reach the same gate.
  lower <- strict_fit(covariance_status = "not_requested", cov_matrix = NULL,
                      cov_condition_number = NaN)
  upper <- strict_fit(covariance_status = "NotRequested", cov_matrix = NULL,
                      cov_condition_number = NaN)
  expect_equal(
    check_strictness(lower, require_covariance = TRUE)$failures,
    check_strictness(upper, require_covariance = TRUE)$failures
  )
})

test_that("the threshold gates fail, pass and switch off", {
  cn <- strict_fit(cov_condition_number = 5000)
  expect_false(check_strictness(cn)$passed)
  expect_match(check_strictness(cn)$failures, "condition number", all = FALSE)
  expect_true(check_strictness(cn, max_condition_number = 1e5)$passed)
  expect_true(check_strictness(cn, max_condition_number = NULL)$passed)
  expect_true(check_strictness(cn, max_condition_number = NA)$passed)

  r <- strict_fit(max_abs_correlation = 0.99)
  expect_false(check_strictness(r)$passed)
  expect_match(check_strictness(r)$failures, "correlation", all = FALSE)
  expect_true(check_strictness(r, max_correlation = 0.999)$passed)
  expect_true(check_strictness(r, max_correlation = NULL)$passed)
})

test_that("a gate with no input is skipped, not failed", {
  fit <- strict_fit(cov_matrix = NULL, cov_condition_number = NaN,
                    max_abs_correlation = NaN, stalled_at_init = NA)
  v <- check_strictness(fit)
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 3L)
  expect_match(v$skipped, "init stall", all = FALSE)

  # Making the covariance step mandatory is what turns the untestable case
  # into a failure.
  fit$covariance_status <- "failed"
  expect_false(check_strictness(fit, require_covariance = TRUE)$passed)
})

test_that("a boundary estimate fails and carries the engine's message", {
  fit <- strict_fit(
    near_boundary = TRUE,
    boundary_estimate_message = "TVCL is at its upper bound (100)"
  )
  v <- check_strictness(fit)
  expect_false(v$passed)
  expect_match(v$failures, "pinned to a declared bound", all = FALSE)
  expect_match(v$failures, "TVCL", all = FALSE)
  expect_true(check_strictness(fit, reject_on_boundary = FALSE)$passed)
})

test_that("a run that never left its initial estimates fails", {
  fit <- strict_fit(stalled_at_init = TRUE)
  v <- check_strictness(fit)
  expect_false(v$passed)
  expect_match(v$failures, "stalled at the initial estimates", all = FALSE)
  expect_true(check_strictness(fit, reject_init_stall = FALSE)$passed)
})

test_that("every gate off passes anything", {
  worst <- strict_fit(
    converged = FALSE, covariance_status = "failed", cov_matrix = NULL,
    cov_condition_number = 1e9, max_abs_correlation = 0.999,
    near_boundary = TRUE, boundary_estimate_message = "TVCL at bound",
    stalled_at_init = TRUE
  )
  v <- check_strictness(
    worst,
    require_converged = FALSE, require_covariance = FALSE,
    max_condition_number = NULL, max_correlation = NULL,
    reject_on_boundary = FALSE, reject_init_stall = FALSE
  )
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 0L)

  # ... and with the gates back on, every one that tripped is reported, not
  # just the first. The condition number is not among them: this fit stores no
  # covariance matrix, so that gate has no input and is skipped.
  v <- check_strictness(worst)
  expect_length(v$failures, 4L)
  expect_length(v$skipped, 1L)
  expect_match(v$skipped, "condition number", all = FALSE)
})

test_that("check_strictness() validates its arguments", {
  expect_error(check_strictness(list(converged = TRUE)), "must be a ferx_fit")
  expect_error(
    check_strictness(strict_fit(), max_condition_number = c(1, 2)),
    "single number"
  )
  expect_error(
    check_strictness(strict_fit(), max_correlation = "high"),
    "single number"
  )
})

test_that("the strictness ingredients survive a .fitrx round trip", {
  fit <- warfarin_fit_cov()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)
  ferx_save_fit(fit, tmp)
  back <- ferx_load_fit(tmp)

  expect_equal(back$max_abs_correlation, fit$max_abs_correlation)
  expect_equal(back$near_boundary, fit$near_boundary)
  expect_equal(back$boundary_estimate_message, fit$boundary_estimate_message)
  expect_equal(back$stalled_at_init, fit$stalled_at_init)
  expect_equal(check_strictness(back), check_strictness(fit))
})
