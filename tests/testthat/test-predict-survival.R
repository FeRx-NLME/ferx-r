test_that("ferx_predict_survival returns a survival data.frame for a TTE model", {
  # Minimal standalone exponential TTE model (lambda_pop = 0.1) + tiny dataset.
  model <- tempfile(fileext = ".ferx")
  data <- tempfile(fileext = ".csv")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "",
    "[event_model]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), model)
  writeLines(c(
    "ID,TIME,DV,EVID,CMT,MDV",
    "1,7.2,1,0,2,0",
    "2,24.0,0,0,2,0",
    "3,3.1,1,0,2,0"
  ), data)

  preds <- ferx_predict_survival(model, data, times = c(0, 6, 12, 24))

  expect_s3_class(preds, "data.frame")
  expect_true(all(
    c(
      "ID", "CMT", "TIME", "survival", "cum_hazard", "hazard",
      "median_survival", "mean_survival"
    ) %in% names(preds)
  ))
  expect_equal(unique(preds$CMT), 2L)

  # S(t) in [0, 1]; S(0) = 1; monotone non-increasing within a subject.
  expect_true(all(preds$survival >= 0 & preds$survival <= 1 + 1e-9))
  s0 <- preds$survival[preds$TIME == 0]
  expect_true(all(abs(s0 - 1) < 1e-9))
  s_subj1 <- preds$survival[preds$ID == "1"][order(preds$TIME[preds$ID == "1"])]
  expect_true(all(diff(s_subj1) <= 1e-9))

  # Exponential population median = log(2) / lambda = log(2) / 0.1 ~= 6.93.
  expect_equal(
    unique(preds$median_survival[preds$ID == "1"]),
    log(2) / 0.1,
    tolerance = 0.05
  )
})

test_that("ferx_predict_survival validates its arguments", {
  expect_error(ferx_predict_survival("does_not_exist.ferx", "nope.csv", times = 1))
  # Non-finite times must be rejected up front, not forwarded to the engine as
  # silent NaN survival rows.
  m <- tempfile(fileext = ".ferx")
  d <- tempfile(fileext = ".csv")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "[event_model]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), m)
  writeLines(c("ID,TIME,DV,EVID,CMT,MDV", "1,7.2,1,0,2,0"), d)
  expect_error(ferx_predict_survival(m, d, times = c(0, NA, 12)))
  expect_error(ferx_predict_survival(m, d, times = c(0, Inf)))
})

test_that("ferx_predict_survival uses fitted theta when a fit is supplied (from_fit path)", {
  # Exercises the validate_fit_for_params -> ferx_rust_predict_survival_from_fit
  # marshalling (theta / omega / sigma), which the population-path tests don't.
  model <- tempfile(fileext = ".ferx")
  data <- tempfile(fileext = ".csv")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "",
    "[event_model]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), model)
  writeLines(c(
    "ID,TIME,DV,EVID,CMT,MDV",
    "1,7.2,1,0,2,0",
    "2,24.0,0,0,2,0",
    "3,3.1,1,0,2,0"
  ), data)

  fit   <- ferx_fit(model, data, method = "focei")
  preds <- ferx_predict_survival(model, data, times = c(0, 6, 12), fit = fit)

  expect_s3_class(preds, "data.frame")
  expect_true(all(
    c("ID", "CMT", "TIME", "survival", "median_survival") %in% names(preds)
  ))
  # Valid survival function on the fitted parameters: S in [0, 1], S(0) = 1.
  expect_true(all(preds$survival >= 0 & preds$survival <= 1 + 1e-9))
  expect_true(all(abs(preds$survival[preds$TIME == 0] - 1) < 1e-9))
  expect_true(all(is.finite(preds$median_survival) & preds$median_survival > 0))
})

test_that("ferx_predict_survival returns an empty frame for a model with no TTE endpoint", {
  # A PK-only (non-[event_model]) model has no TTE endpoint, so predict_survival
  # yields zero rows rather than erroring.
  ex    <- ferx_example("one_cpt_iv")
  preds <- ferx_predict_survival(ex$model, ex$data, times = c(0, 6, 12))

  expect_s3_class(preds, "data.frame")
  expect_equal(nrow(preds), 0L)
})
