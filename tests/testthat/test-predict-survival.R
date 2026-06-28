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
      "cif", "survival_all", "median_survival", "mean_survival"
    ) %in% names(preds)
  ))
  expect_equal(unique(preds$CMT), 2L)

  # S(t) in [0, 1]; S(0) = 1; monotone non-increasing within a subject.
  expect_true(all(preds$survival >= 0 & preds$survival <= 1 + 1e-9))
  s0 <- preds$survival[preds$TIME == 0]
  expect_true(all(abs(s0 - 1) < 1e-9))
  s_subj1 <- preds$survival[preds$ID == "1"][order(preds$TIME[preds$ID == "1"])]
  expect_true(all(diff(s_subj1) <= 1e-9))

  # Single endpoint: cif reduces to 1 - survival and survival_all == survival.
  expect_true(all(abs(preds$cif - (1 - preds$survival)) < 1e-9))
  expect_true(all(abs(preds$survival_all - preds$survival) < 1e-9))

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

test_that("ferx_predict_survival reports competing-risks CIF with sum(cif) + survival_all = 1", {
  # Two cause-specific exponential hazards (CMT 2, CMT 3) sharing a frailty.
  model <- tempfile(fileext = ".ferx")
  data  <- tempfile(fileext = ".csv")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA_A(0.10, 0.001, 10.0)",
    "  theta TVLAMBDA_B(0.06, 0.001, 10.0)",
    "  omega ETA_F ~ 0.09",
    "",
    "[event_model cause_a]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA_A * exp(ETA_F)",
    "",
    "[event_model cause_b]",
    "  cmt    = 3",
    "  family = exponential",
    "  scale  = TVLAMBDA_B * exp(ETA_F)"
  ), model)
  # Cause-specific layout: one row per cause CMT per subject (DV=1 the observed
  # cause, DV=0 the other censored at the same time).
  writeLines(c(
    "ID,TIME,DV,EVID,CMT,MDV",
    "1,5.0,1,0,2,0",
    "1,5.0,0,0,3,0",
    "2,8.0,0,0,2,0",
    "2,8.0,1,0,3,0",
    "3,14.0,0,0,2,0",
    "3,14.0,0,0,3,0"
  ), data)

  times <- c(0, 4, 10, 20)
  preds <- ferx_predict_survival(model, data, times = times)

  expect_s3_class(preds, "data.frame")
  expect_true(all(c("cif", "survival_all") %in% names(preds)))
  # Both cause CMTs are predicted.
  expect_setequal(unique(preds$CMT), c(2L, 3L))
  # CIF in [0, 1] and the all-cause survival is no greater than any cause-specific
  # survival (Σ_j H_j >= H_k).
  expect_true(all(preds$cif >= -1e-9 & preds$cif <= 1 + 1e-9))
  expect_true(all(preds$survival_all <= preds$survival + 1e-9))

  # Partition invariant: at each (ID, TIME), Σ_k cif_k + survival_all = 1.
  cif_sum <- aggregate(cif ~ ID + TIME, data = preds, FUN = sum)
  s_all   <- unique(preds[, c("ID", "TIME", "survival_all")])
  merged  <- merge(cif_sum, s_all, by = c("ID", "TIME"))
  expect_true(all(abs(merged$cif + merged$survival_all - 1) < 1e-9))
})

test_that("joint PK-TTE: bundled pktte_joint fits and predicts an ODE-accumulated hazard", {
  # Drug-driven hazard h = H0 * exp(BETA * Cc) carried as a cumulative-hazard ODE
  # state and estimated jointly with the PK (ferx-core #564). Exercises both the
  # joint fit and the ODE-hazard read in ferx_predict_survival through the wrapper.
  ex  <- ferx_example("pktte_joint")
  fit <- ferx_fit(ex$model, ex$data, method = "focei")
  expect_s3_class(fit, "ferx_fit")
  expect_true(is.finite(fit$ofv))

  surv <- ferx_predict_survival(ex$model, ex$data, times = c(0, 6, 12, 24), fit = fit)
  expect_s3_class(surv, "data.frame")
  expect_true(all(
    c("ID", "CMT", "TIME", "survival", "cum_hazard", "hazard") %in% names(surv)
  ))
  # Event endpoint is on CMT 3 (PK is on CMT 2).
  expect_equal(unique(surv$CMT), 3L)
  # Valid survival function: S in [0, 1], S(0) = 1, monotone non-increasing; the
  # drug-driven hazard is strictly positive.
  expect_true(all(surv$survival >= 0 & surv$survival <= 1 + 1e-9))
  expect_true(all(abs(surv$survival[surv$TIME == 0] - 1) < 1e-9))
  expect_true(all(surv$hazard >= 0))
  id1 <- surv$ID[1]
  s1  <- surv$survival[surv$ID == id1][order(surv$TIME[surv$ID == id1])]
  expect_true(all(diff(s1) <= 1e-9))
})
