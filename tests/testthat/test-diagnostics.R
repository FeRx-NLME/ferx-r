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

test_that("check_diagnostics flag values match .dw_label output", {
  fit <- warfarin_fit()
  fit$dw_statistic <- 1.2
  fit$iwres_lag1_r <- 0.4
  d <- check_diagnostics(fit)
  expect_equal(d$autocorrelation$flag, "positive autocorrelation")

  fit$dw_statistic <- 2.8
  d <- check_diagnostics(fit)
  expect_equal(d$autocorrelation$flag, "negative autocorrelation")

  fit$dw_statistic <- 2.0
  d <- check_diagnostics(fit)
  expect_equal(d$autocorrelation$flag, "no autocorrelation")
})

test_that("check_diagnostics handles multi-sigma shrinkage_eps vector", {
  fit <- warfarin_fit()
  fit$shrinkage_eps  <- c(0.10, 0.15)
  fit$sigma_names    <- c("EPS_PROP", "EPS_ADD")
  fit$shrinkage_eta  <- NULL
  d <- check_diagnostics(fit)
  expect_equal(nrow(d$shrinkage), 2L)
  expect_equal(d$shrinkage$param, c("EPS_PROP", "EPS_ADD"))
  expect_equal(d$shrinkage$shrinkage_pct, c(10, 15))
})

test_that("check_diagnostics uses EPS1/EPS2 fallback when sigma_names absent for multi-sigma", {
  fit <- warfarin_fit()
  fit$shrinkage_eps <- c(0.08, 0.12)
  fit$sigma_names   <- NULL
  fit$shrinkage_eta <- NULL
  d <- check_diagnostics(fit)
  expect_equal(d$shrinkage$param, c("EPS1", "EPS2"))
})

test_that("positive DW < 1.5 emits message; negative DW > 2.5 emits message; no message otherwise", {
  fit <- list(dw_statistic = 1.2, iwres_lag1_r = 0.4, uses_sde = FALSE)
  expect_message(ferx:::.ferx_emit_dw_message(fit), "Positive IWRES autocorrelation")

  fit$dw_statistic <- 2.8
  expect_message(ferx:::.ferx_emit_dw_message(fit), "Negative IWRES autocorrelation")

  fit$dw_statistic <- 2.0
  expect_silent(ferx:::.ferx_emit_dw_message(fit))
})
