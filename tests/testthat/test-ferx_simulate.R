
# ---- header from test-simulate-tte.R ----
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

test_that("ferx_simulate returns a data frame", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_s3_class(sim, "data.frame")
})
test_that("ferx_simulate has required columns: SIM, ID, TIME, IPRED, DV_SIM", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_true(all(c("SIM", "ID", "TIME", "IPRED", "DV_SIM") %in% names(sim)))
})
test_that("n_sim = 3 produces exactly 3 unique SIM values", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 3L, seed = 1L)
  expect_equal(length(unique(sim$SIM)), 3L)
})
test_that("n_sim = 1 nrow matches observation rows in source data", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  dat <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(sim), n_obs)
})
test_that("same seed produces identical output", {
  ex   <- ferx_example("warfarin")
  sim1 <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 42L)
  sim2 <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 42L)
  expect_equal(sim1, sim2)
})
test_that("different seeds produce different DV_SIM", {
  ex   <- ferx_example("warfarin")
  sim1 <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  sim2 <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 99L)
  expect_false(identical(sim1$DV_SIM, sim2$DV_SIM))
})
test_that("DV_SIM is finite numeric with no NAs", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_true(is.numeric(sim$DV_SIM))
  expect_true(all(is.finite(sim$DV_SIM)))
})
test_that("simulate-from-fit returns data frame with required columns", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, fit = warfarin_fit())
  expect_s3_class(sim, "data.frame")
  expect_true(all(c("SIM", "ID", "TIME", "IPRED", "DV_SIM") %in% names(sim)))
})
test_that("simulate-from-fit IPRED values are finite", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, fit = warfarin_fit())
  expect_true(all(is.finite(sim$IPRED)))
})
test_that("simulate-from-fit produces different IPRED than population simulation", {
  ex      <- ferx_example("warfarin")
  sim_pop <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)
  sim_ind <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, fit = warfarin_fit())
  expect_false(identical(sim_pop$IPRED, sim_ind$IPRED))
})
test_that("each match method returns the same shape as the unmatched simulation", {
  ex  <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 1L)
  for (m in list(TRUE, "optimal", "nearest", "rank")) {
    simm <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 1L, match = m)
    expect_s3_class(simm, "data.frame")
    expect_true(all(c("SIM", "ID", "TIME", "IPRED", "DV_SIM") %in% names(simm)),
                info = paste("method", m))
    expect_equal(nrow(simm), nrow(sim), info = paste("method", m))
    expect_true(all(is.finite(simm$DV_SIM)), info = paste("method", m))
  }
})
test_that("each match method is reproducible and differs from the unmatched path", {
  ex        <- ferx_example("warfarin")
  unmatched <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 7L, match = FALSE)
  for (m in c("optimal", "nearest", "rank")) {
    a <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 7L, match = m)
    b <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 7L, match = m)
    expect_equal(a, b, info = paste("method", m))
    expect_false(identical(a$DV_SIM, unmatched$DV_SIM), info = paste("method", m))
  }
})
test_that("match = TRUE is equivalent to match = \"optimal\"", {
  ex <- ferx_example("warfarin")
  a  <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 3L, match = TRUE)
  b  <- ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 3L, match = "optimal")
  expect_equal(a, b)
})
test_that("match methods work with a fit", {
  ex <- ferx_example("warfarin")
  for (m in c("optimal", "nearest", "rank")) {
    simm <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L,
                          fit = warfarin_fit(), match = m)
    expect_s3_class(simm, "data.frame")
    expect_true(all(is.finite(simm$DV_SIM)), info = paste("method", m))
  }
})
test_that("ferx_simulate rejects an unknown match method", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, match = "bogus"),
    "match"
  )
})
test_that("ferx_simulate errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_simulate("no_such_model.ferx", ex$data))
})
test_that("ferx_simulate errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_simulate(ex$model, "no_such_data.csv"))
})

# ---- Simulation diagnostics channel (ferx-core #762/#763) ------------------
# ferx-core surfaces per-subject simulation diagnostics (a degenerate or
# pathological hazard that would otherwise censor a subject with no signal)
# via simulate_with_options_diag; the Rust glue attaches them as a
# `simulation_warnings` attribute and ferx_simulate emits them as an R warning.

test_that("ferx_simulate attaches a (clean-run empty) simulation_warnings attribute", {
  ex <- ferx_example("pktte_joint")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, horizon = 24)
  # The channel is always present so downstream code can rely on it; a clean run
  # yields a zero-length character vector, not NULL.
  w <- attr(sim, "simulation_warnings", exact = TRUE)
  expect_false(is.null(w))
  expect_type(w, "character")
  expect_length(w, 0)
})

test_that(".ferx_surface_sim_warnings emits attached diagnostics as one warning, returns object", {
  fake <- structure(
    data.frame(ID = 1L, TIME = 0),
    simulation_warnings = c(
      "W_TTE_DEGENERATE_HAZARD: subject 3 (CMT 2)",
      "W_RTTE_DEGENERATE: subject 7"
    )
  )
  expect_warning(out <- .ferx_surface_sim_warnings(fake), "2 simulation diagnostics")
  # Returned unchanged, attribute preserved for programmatic access.
  expect_identical(out, fake)
  expect_length(attr(out, "simulation_warnings"), 2L)
})

test_that(".ferx_surface_sim_warnings is silent with no diagnostics (empty or NULL)", {
  clean <- structure(data.frame(ID = 1L), simulation_warnings = character(0))
  expect_no_warning(.ferx_surface_sim_warnings(clean))
  expect_no_warning(.ferx_surface_sim_warnings(NULL))
})
