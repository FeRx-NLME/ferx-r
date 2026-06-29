# Joint PK-TTE event-time simulation through ferx_simulate() (ferx-core #564,
# Slice 2.2). The bundled `pktte_joint` model has a drug-driven hazard
# (`[event_model] hazard = H0 * exp(BETA * (central / V))`, event on CMT 3), so
# simulate() samples event times by integrating the augmented hazard ODE until
# the cumulative hazard reaches -log(U). A finite `horizon` is required.

test_that("ferx_simulate samples drug-driven TTE event times given a horizon", {
  ex <- ferx_example("pktte_joint")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, horizon = 24)

  expect_s3_class(sim, "data.frame")
  expect_true(all(c("CMT", "TIME", "IPRED", "DV_SIM", "OBSERVED") %in% names(sim)))

  # Event rows: OBSERVED is non-NA (1 = event, 0 = right-censored), on the event
  # CMT (3), with TIME inside the horizon.
  ev <- sim[!is.na(sim$OBSERVED), ]
  expect_gt(nrow(ev), 0)
  expect_true(all(ev$OBSERVED %in% c(0, 1)))
  expect_true(all(ev$CMT == 3))
  expect_true(all(ev$TIME <= 24 + 1e-6))
  # An observed event strictly precedes the horizon; a censor sits at it.
  expect_true(all(ev$TIME[ev$OBSERVED == 1] < 24))
  expect_true(all(abs(ev$TIME[ev$OBSERVED == 0] - 24) < 1e-6))
  # Event rows carry no Gaussian prediction.
  expect_true(all(is.na(ev$DV_SIM)))

  # Continuous PK rows are still present, with OBSERVED = NA and finite DV_SIM.
  cont <- sim[is.na(sim$OBSERVED), ]
  expect_gt(nrow(cont), 0)
  expect_true(all(is.finite(cont$DV_SIM)))
})

test_that("an ODE-accumulated TTE model without a horizon does not simulate", {
  # ferx-core rejects drug-driven event-time sampling without a finite horizon
  # (the hazard can vanish, so there is no implicit observation window); the
  # simulate glue surfaces that as a NULL result rather than fabricating rows.
  ex <- ferx_example("pktte_joint")
  res <- suppressWarnings(
    ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  )
  expect_null(res)
})

test_that("ferx_simulate from a fitted joint PK-TTE model honours `horizon`", {
  # The `fit` path threads `horizon` the same way as the default-parameter path;
  # exercise it end-to-end with a cheap, iteration-capped fit (a converged fit is
  # not needed - any valid theta/omega/sigma drives the event-time sampler).
  ex <- ferx_example("pktte_joint")
  fit <- ferx_fit(ex$model, ex$data, method = "focei", covariance = FALSE,
                  settings = list(maxiter = 3L))
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 2L,
                       fit = fit, horizon = 24)

  expect_s3_class(sim, "data.frame")
  ev <- sim[!is.na(sim$OBSERVED), ]
  expect_gt(nrow(ev), 0)
  expect_true(all(ev$CMT == 3))
  expect_true(all(ev$OBSERVED %in% c(0, 1)))
  expect_true(all(ev$TIME <= 24 + 1e-6))
})

test_that("ferx_simulate validates the horizon argument", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate(ex$model, ex$data, horizon = -1),
    "finite positive"
  )
  expect_error(
    ferx_simulate(ex$model, ex$data, horizon = c(1, 2)),
    "single finite positive"
  )
  expect_error(
    ferx_simulate(ex$model, ex$data, horizon = "24"),
    "finite positive"
  )
})
