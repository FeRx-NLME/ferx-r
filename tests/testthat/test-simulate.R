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
