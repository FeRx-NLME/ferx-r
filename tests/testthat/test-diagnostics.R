# check_diagnostics() — Tier 1

test_that("check_diagnostics returns a named list with autocorrelation and shrinkage", {
  fit <- warfarin_fit()
  d   <- check_diagnostics(fit)
  expect_type(d, "list")
  expect_named(d, c("autocorrelation", "shrinkage"))
})

test_that("check_diagnostics$autocorrelation is NULL or a data frame with required columns", {
  fit <- warfarin_fit()
  d   <- check_diagnostics(fit)
  if (!is.null(d$autocorrelation)) {
    expect_s3_class(d$autocorrelation, "data.frame")
    expect_true(all(c("dw_statistic", "lag1_r", "flag") %in% names(d$autocorrelation)))
    expect_equal(nrow(d$autocorrelation), 1L)
    expect_type(d$autocorrelation$flag, "character")
  }
})

test_that("check_diagnostics$shrinkage is NULL or a data frame with required columns", {
  fit <- warfarin_fit()
  d   <- check_diagnostics(fit)
  if (!is.null(d$shrinkage)) {
    expect_s3_class(d$shrinkage, "data.frame")
    expect_true(all(c("param", "type", "shrinkage", "shrinkage_pct") %in% names(d$shrinkage)))
    expect_true(all(d$shrinkage$type %in% c("eta", "eps")))
  }
})

test_that("dw_statistic and iwres_lag1_r on fit are numeric scalars or NA", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$dw_statistic) && length(fit$dw_statistic) == 1L)
  expect_true(is.numeric(fit$iwres_lag1_r) && length(fit$iwres_lag1_r) == 1L)
})

test_that("check_diagnostics works on a fit with no shrinkage data", {
  fit <- warfarin_fit()
  fit2 <- fit
  fit2$shrinkage_eta <- NULL
  fit2$shrinkage_eps <- NULL
  d <- check_diagnostics(fit2)
  expect_null(d$shrinkage)
})

test_that("check_diagnostics works on a fit with no autocorrelation data", {
  fit <- warfarin_fit()
  fit2 <- fit
  fit2$dw_statistic <- NA_real_
  fit2$iwres_lag1_r <- NA_real_
  d <- check_diagnostics(fit2)
  expect_null(d$autocorrelation)
})
