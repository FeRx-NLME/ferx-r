# Synthetic fit carrying only what check_strictness() reads. Defaults describe
# a clean, converged, well-conditioned fit so each test can flip one thing.
strict_fit <- function(...) {
  base <- list(
    converged = TRUE,
    covariance_status = "computed",
    cov_matrix = matrix(c(1, 0.1, 0.1, 1), 2, 2),
    condition_number = 12,
    max_abs_correlation = 0.1,
    estimate_near_boundary = FALSE,
    stalled_at_init = FALSE,
    sir_ess = NULL
  )
  structure(utils::modifyList(base, list(...)), class = "ferx_fit")
}

test_that("a clean fit passes every default gate", {
  v <- check_strictness(strict_fit())
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 0L)
})

test_that("turning every gate off passes and skips nothing", {
  broken <- strict_fit(converged = FALSE, condition_number = 1e9,
                       max_abs_correlation = 0.99,
                       estimate_near_boundary = TRUE, stalled_at_init = TRUE)
  v <- check_strictness(broken, require_converged = FALSE,
                        require_covariance = FALSE,
                        max_condition_number = NULL, max_correlation = NULL,
                        reject_on_boundary = FALSE, reject_init_stall = FALSE)
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 0L)
})

test_that("each gate reports its own failure", {
  v <- check_strictness(strict_fit(converged = FALSE))
  expect_false(v$passed)
  expect_match(v$failures, "did not converge")

  v <- check_strictness(strict_fit(condition_number = 5000))
  expect_match(v$failures, "condition number")

  v <- check_strictness(strict_fit(max_abs_correlation = 0.99))
  expect_match(v$failures, "correlation")

  v <- check_strictness(strict_fit(estimate_near_boundary = TRUE))
  expect_match(v$failures, "pinned to a declared bound")

  v <- check_strictness(strict_fit(stalled_at_init = TRUE))
  expect_match(v$failures, "stalled at the initial estimates")
})

test_that("a threshold gate with no input is skipped, not failed", {
  v <- check_strictness(strict_fit(condition_number = NULL,
                                   max_abs_correlation = NA_real_,
                                   stalled_at_init = NA))
  expect_true(v$passed)
  expect_length(v$failures, 0L)
  expect_length(v$skipped, 3L)
  expect_match(paste(v$skipped, collapse = " "), "condition number")
  expect_match(paste(v$skipped, collapse = " "), "parameter correlation")
  expect_match(paste(v$skipped, collapse = " "), "init stall")
})

test_that("the boundary failure names the offending parameters when known", {
  ws <- data.frame(
    severity = "warning", category = "boundary_estimate",
    message = "TVCL is at its lower bound", stringsAsFactors = FALSE
  )
  v <- check_strictness(strict_fit(estimate_near_boundary = TRUE,
                                   warnings_structured = ws))
  expect_match(v$failures, "TVCL is at its lower bound", fixed = TRUE)
})

test_that("require_covariance reads both status spellings", {
  # ferx_fit() carries the engine token, ferx_load_fit() a CamelCase label.
  for (token in c("not_requested", "NotRequested")) {
    v <- check_strictness(strict_fit(covariance_status = token,
                                     cov_matrix = NULL),
                          require_covariance = TRUE)
    expect_match(v$failures, "not requested")
  }
  for (token in c("failed", "Failed")) {
    v <- check_strictness(strict_fit(covariance_status = token,
                                     cov_matrix = NULL),
                          require_covariance = TRUE)
    expect_match(v$failures, "covariance step failed")
  }
})

test_that("require_covariance rejects a computed step that stored no matrix", {
  v <- check_strictness(strict_fit(cov_matrix = NULL),
                        require_covariance = TRUE)
  expect_match(v$failures, "stored no matrix")
})

test_that("a SIR fallback counts as uncertainty only when it delivered", {
  delivered <- strict_fit(covariance_status = "sir_fallback",
                          cov_matrix = NULL, sir_ess = 40)
  expect_true(check_strictness(delivered, require_covariance = TRUE)$passed)

  collapsed <- strict_fit(covariance_status = "SirFallback",
                          cov_matrix = NULL, sir_ess = 1.2)
  v <- check_strictness(collapsed, require_covariance = TRUE)
  expect_match(v$failures, "SIR fallback collapsed")

  unjudgeable <- strict_fit(covariance_status = "sir_fallback",
                            cov_matrix = NULL, sir_ess = NULL)
  v <- check_strictness(unjudgeable, require_covariance = TRUE)
  expect_match(v$failures, "no effective sample size")
})

test_that("gate thresholds are validated", {
  expect_error(check_strictness(strict_fit(), max_correlation = c(0.5, 0.9)),
               "single number")
  expect_error(check_strictness(list(converged = TRUE)), "ferx_fit")
})

test_that("a real fit carries the fields the gates read", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  expect_type(fit$bic_inputs, "list")
  expect_false(is.null(fit$estimate_near_boundary))
  expect_length(fit$stalled_at_init, 1L)
  expect_length(fit$max_abs_correlation, 1L)
  v <- check_strictness(fit)
  expect_type(v$passed, "logical")
  expect_type(v$failures, "character")
  expect_type(v$skipped, "character")
})
