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

# ferx-core #1261: an ODE-accumulated hazard that reads a dose-time anchor. The
# hazard readout re-evaluated the ODE right-hand side with the bare PK parameter
# array, which has no TAD / TAFD slots, so `hazard` came back NaN and a fit scored
# the subject at the 2e20 rejection sentinel. The hazard, doses and event follow
# ferx-core's own regression test: H0 = 0.1, KT = 0.01, doses at t = 0 and
# t = 12, one exact event at t = 30. After the second dose TAD (time after the
# last dose) and TAFD (time after the first dose) differ, so the two anchors
# cannot stand in for each other.
# `hazard = NULL` writes the PK-only sibling: same [odes], doses and PK rows, no
# [event_model] and no event row.
write_tad_hazard_fixture <- function(hazard = "H0 * (1.0 + KT * TAD)") {
  tte <- !is.null(hazard)
  model <- tempfile(fileext = ".ferx")
  data  <- tempfile(fileext = ".csv")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(1.0, 0.01, 100.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    if (tte) c("  theta TVH0(0.1, 0.001, 10.0)", "  theta TVKT(0.01, 0.0001, 1.0)"),
    "  omega ETA_CL ~ 0.09",
    "  sigma PROP_ERR ~ 0.1 (sd)",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    if (tte) c("  H0 = TVH0", "  KT = TVKT"),
    "",
    "[structural_model]",
    "  ode(obs_cmt=central, states=[central])",
    "",
    "[odes]",
    "  d/dt(central) = -CL / V * central",
    "",
    if (tte) c("[event_model]", "  cmt    = 3", paste("  hazard =", hazard), ""),
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ), model)
  writeLines(c(
    "ID,TIME,DV,EVID,AMT,CMT,MDV",
    "1,0,.,1,100,1,1",
    "1,2,80,0,.,1,0",
    "1,12,.,1,100,1,1",
    "1,14,105,0,.,1,0",
    "1,24,40,0,.,1,0",
    if (tte) "1,30,1,0,.,3,0"
  ), data)
  list(model = model, data = data)
}

test_that("ferx_predict_survival reads TAD / TAFD in an ODE-accumulated hazard", {
  # h(t) = H0 * (1 + KT * anchor). At t = 6, 18, 30 TAD is 6, 6, 18 and TAFD is
  # 6, 18, 30; H(t) integrates h piecewise across the second dose. Measured error
  # against these closed forms is ~1e-16 for both columns.
  expected <- list(
    TAD  = list(hazard = c(0.106, 0.106, 0.118), cum_hazard = c(0.618, 1.890, 3.234)),
    TAFD = list(hazard = c(0.106, 0.118, 0.130), cum_hazard = c(0.618, 1.962, 3.450))
  )
  for (anchor in names(expected)) {
    fx   <- write_tad_hazard_fixture(paste0("H0 * (1.0 + KT * ", anchor, ")"))
    surv <- ferx_predict_survival(fx$model, fx$data, times = c(6, 18, 30))
    surv <- surv[order(surv$TIME), ]
    expect_equal(surv$hazard, expected[[anchor]]$hazard,
                 tolerance = 1e-6, label = paste(anchor, "hazard"))
    expect_equal(surv$cum_hazard, expected[[anchor]]$cum_hazard,
                 tolerance = 1e-6, label = paste(anchor, "cum_hazard"))
  }
})

test_that("ferx_fit scores a joint PK-TTE subject whose ODE-accumulated hazard reads TAD", {
  # The PK-only sibling shares the Gaussian term, so the difference is the TTE
  # term alone: the exact event at t = 30 adds H(30) - log h(30) =
  # 3.234 - log(0.118) to the NLL, twice that to the OFV. The measured gap to that
  # closed form is 2.8e-4 on the OFV (2.5e-5 with the random effect removed, so
  # most of it is the inner-loop EBE solve); the relative tolerance below is about
  # 1e-3 absolute.
  ofv_at_init <- function(fx) {
    ferx_fit(fx$model, fx$data, method = "focei", covariance = FALSE,
             verbose = FALSE, settings = list(maxiter = 0))$ofv
  }
  ofv_joint <- ofv_at_init(write_tad_hazard_fixture())
  ofv_pk    <- ofv_at_init(write_tad_hazard_fixture(hazard = NULL))
  expect_equal(ofv_joint - ofv_pk, 2 * (3.234 - log(0.118)), tolerance = 1e-4)
})

# -- Engine findings on the survival predictor (#433) --------------------------
# Both glue entry points relay the `warnings` of core's `predict_survival_diag`
# (the non-fit bundle: parse, reader-through-filter, model/data, ODE solver) on
# the `simulation_warnings` channel `ferx_predict()` uses.

write_exp_tte <- function(extra = character()) {
  model <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    extra,
    "  omega ETA_LAMBDA ~ 0.09",
    "",
    "[event_model]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), model)
  model
}

write_exp_tte_data <- function() {
  data <- tempfile(fileext = ".csv")
  writeLines(c(
    "ID,TIME,DV,EVID,CMT,MDV",
    "1,7.2,1,0,2,0",
    "2,24.0,0,0,2,0",
    "3,3.1,1,0,2,0"
  ), data)
  data
}

unused_theta_note <- "theta 'UNUSED' is declared"

test_that("ferx_predict_survival relays the engine's parse findings (#433)", {
  model <- write_exp_tte("  theta UNUSED(1.0, 0.1, 10.0)")
  expect_warning(
    preds <- ferx_predict_survival(model, write_exp_tte_data(), times = c(0, 6)),
    unused_theta_note, fixed = TRUE
  )
  w <- attr(preds, "simulation_warnings", exact = TRUE)
  expect_type(w, "character")
  expect_true(any(grepl(unused_theta_note, w, fixed = TRUE)))
})

test_that("ferx_predict_survival(fit = ) relays the engine's findings too (#433)", {
  # The separate glue entry point (`ferx_rust_predict_survival_from_fit`).
  model <- write_exp_tte("  theta UNUSED(1.0, 0.1, 10.0)")
  data <- write_exp_tte_data()
  fit <- suppressWarnings(ferx_fit(model, data, method = "foce", covariance = FALSE))
  expect_warning(
    preds <- ferx_predict_survival(model, data, times = c(0, 6), fit = fit),
    unused_theta_note, fixed = TRUE
  )
  expect_true(any(grepl(unused_theta_note,
                        attr(preds, "simulation_warnings", exact = TRUE),
                        fixed = TRUE)))
})

test_that("a clean ferx_predict_survival run carries an empty attribute (#433)", {
  # A compartment-free TTE model on dose-free data: the reader's moot
  # `W_NO_DOSES` / `W_CMT_DEFAULTED` stay behind core's filter.
  expect_no_warning(
    preds <- ferx_predict_survival(write_exp_tte(), write_exp_tte_data(),
                                   times = c(0, 6))
  )
  expect_identical(attr(preds, "simulation_warnings", exact = TRUE), character(0))
})

# -- type = "individual": each subject's curve at its EBEs (ferx-core #1885) ---
#
# On the bundled exponential model the hazard is TVLAMBDA * exp(ETA_LAMBDA), so
# a subject's individual cumulative hazard is the population one times
# exp(eta_hat): an oracle that needs no second implementation of the curve.
tte_exp_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("tte_exponential")
      fit <<- suppressWarnings(ferx_fit(ex$model, ex$data, method = "focei",
                                        covariance = FALSE, verbose = FALSE))
    }
    fit
  }
})

test_that("type = \"individual\" scales each subject's H(t) by exp(eta_hat)", {
  ex <- ferx_example("tte_exponential")
  fit <- tte_exp_fit()
  times <- c(0, 5, 20)
  pop <- ferx_predict_survival(ex$model, ex$data, times = times, fit = fit)
  ind <- ferx_predict_survival(ex$model, ex$data, times = times, fit = fit,
                               type = "individual")

  expect_identical(names(ind), names(pop))
  expect_identical(ind$ID, pop$ID)
  expect_identical(ind$TIME, pop$TIME)
  eta <- setNames(fit$ebe_etas$ETA_LAMBDA, as.character(fit$ebe_etas$ID))
  expect_gt(stats::sd(eta), 0)
  later <- ind$TIME > 0
  expect_equal(ind$cum_hazard[later] / pop$cum_hazard[later],
               unname(exp(eta[as.character(ind$ID[later])])), tolerance = 1e-10)
  expect_equal(ind$survival, exp(-ind$cum_hazard), tolerance = 1e-12)
  # Subjects now differ from one another; the population curve is shared.
  at20 <- ind$TIME == 20
  expect_gt(length(unique(round(ind$survival[at20], 12))), 1L)
  expect_length(unique(round(pop$survival[pop$TIME == 20], 12)), 1L)
})

test_that("type defaults to \"population\" and needs a fit for \"individual\"", {
  ex <- ferx_example("tte_exponential")
  fit <- tte_exp_fit()
  expect_identical(
    ferx_predict_survival(ex$model, ex$data, times = c(0, 10), fit = fit),
    ferx_predict_survival(ex$model, ex$data, times = c(0, 10), fit = fit,
                          type = "population")
  )
  expect_error(
    ferx_predict_survival(ex$model, ex$data, times = c(0, 10), type = "individual"),
    "type = \"individual\" needs `fit`", fixed = TRUE
  )
  expect_error(
    ferx_predict_survival(ex$model, ex$data, times = c(0, 10), fit = fit,
                          type = "posthoc"),
    "should be one of"
  )
})

test_that("type = \"individual\" takes a subset of the fitted subjects, refuses a new one", {
  ex <- ferx_example("tte_exponential")
  fit <- tte_exp_fit()
  rows <- utils::read.csv(ex$data)
  ids <- unique(rows$ID)[1:3]
  subset <- tempfile(fileext = ".csv")
  utils::write.csv(rows[rows$ID %in% ids, ], subset, row.names = FALSE)
  ind <- ferx_predict_survival(ex$model, subset, times = c(0, 10), fit = fit,
                               type = "individual")
  expect_setequal(unique(as.character(ind$ID)), as.character(ids))

  extra <- rows[rows$ID == ids[1], ]
  extra$ID <- 999999L
  stranger <- tempfile(fileext = ".csv")
  utils::write.csv(rbind(rows[rows$ID %in% ids, ], extra), stranger, row.names = FALSE)
  e <- expect_error(
    ferx_predict_survival(ex$model, stranger, times = c(0, 10), fit = fit,
                          type = "individual")
  )
  expect_match(conditionMessage(e), "E_SUPPLIED_ETAS|999999")
})

test_that("type = \"individual\" on a fit without random effects is the population curve", {
  model <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "",
    "[event_model]",
    "  cmt    = 2",
    "  family = exponential",
    "  scale  = TVLAMBDA"
  ), model)
  data <- write_exp_tte_data()
  fit <- suppressWarnings(ferx_fit(model, data, method = "foce",
                                   covariance = FALSE, verbose = FALSE))
  expect_null(fit$ebe_etas)
  expect_identical(
    ferx_predict_survival(model, data, times = c(0, 6), fit = fit, type = "individual"),
    ferx_predict_survival(model, data, times = c(0, 6), fit = fit)
  )
})
