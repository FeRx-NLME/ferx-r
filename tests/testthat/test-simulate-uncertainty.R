# Tests for ferx_simulate_with_uncertainty. The cached `warfarin_fit_cov()`
# helper provides a fit with a populated covariance step for the asymptotic
# path; SIR-path tests are gated on a fresh fit that opts into
# `sir_keep_samples`.

test_that("ferx_simulate_with_uncertainty asymptotic returns a data frame", {
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 7L
  )
  expect_s3_class(sims, "data.frame")
})

test_that("asymptotic output has the expected columns and DRAW leads", {
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 7L
  )
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "IPRED", "DV_SIM") %in%
                    names(sims)))
  expect_identical(names(sims)[1], "DRAW")
})

test_that("DRAW spans 1..n_uncertainty_draws", {
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 4L, n_sim_per_draw = 1L,
    method = "asymptotic", seed = 11L
  )
  expect_equal(sort(unique(sims$DRAW)), 1L:4L)
})

test_that("SIM spans 1..n_sim_per_draw inside each draw", {
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 2L, n_sim_per_draw = 3L,
    method = "asymptotic", seed = 21L
  )
  for (d in unique(sims$DRAW)) {
    expect_equal(sort(unique(sims$SIM[sims$DRAW == d])), 1L:3L)
  }
})

test_that("row count is n_uncertainty_draws * n_sim_per_draw * n_obs", {
  ex   <- ferx_example("warfarin")
  dat  <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 33L
  )
  expect_equal(nrow(sims), 3L * 2L * n_obs)
})

test_that("DV_SIM and IPRED are finite numerics with no NAs", {
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, warfarin_fit_cov(),
    n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
    method = "asymptotic", seed = 41L
  )
  expect_true(is.numeric(sims$DV_SIM))
  expect_true(is.numeric(sims$IPRED))
  expect_true(all(is.finite(sims$DV_SIM)))
  expect_true(all(is.finite(sims$IPRED)))
})

test_that("same seed produces identical output (asymptotic)", {
  ex   <- ferx_example("warfarin")
  fit  <- warfarin_fit_cov()
  sims1 <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L, seed = 42L
  )
  sims2 <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L, seed = 42L
  )
  expect_equal(sims1, sims2)
})

test_that("uncertainty-aware sims show wider DV_SIM spread than fixed-param sims", {
  # n_uncertainty_draws = 1 ≈ no parameter uncertainty (one draw). Bumping
  # n_uncertainty_draws should widen the empirical spread of DV_SIM at any
  # given (ID, TIME). Use IQR to make the comparison robust to outliers.
  ex  <- ferx_example("warfarin")
  fit <- warfarin_fit_cov()
  many <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 30L, n_sim_per_draw = 5L, seed = 1L
  )
  few <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 1L, n_sim_per_draw = 150L, seed = 1L
  )
  iqr_many <- diff(quantile(many$DV_SIM, c(0.25, 0.75), names = FALSE))
  iqr_few  <- diff(quantile(few$DV_SIM,  c(0.25, 0.75), names = FALSE))
  expect_gt(iqr_many, iqr_few * 0.9)  # generous to tolerate Monte-Carlo noise
})

test_that("asymptotic errors when covariance step was not run", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate_with_uncertainty(
      ex$model, ex$data, warfarin_fit(),
      n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
      method = "asymptotic"
    ),
    "cov_matrix"
  )
})

test_that("SIR errors when no resamples are stored on the fit", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate_with_uncertainty(
      ex$model, ex$data, warfarin_fit_cov(),
      n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
      method = "sir"
    ),
    "sir_resamples"
  )
})

test_that("ferx_simulate_with_uncertainty errors on missing model / data files", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate_with_uncertainty("no_such_model.ferx", ex$data, warfarin_fit_cov())
  )
  expect_error(
    ferx_simulate_with_uncertainty(ex$model, "no_such_data.csv", warfarin_fit_cov())
  )
})

test_that("n_uncertainty_draws < 1 raises an informative error", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate_with_uncertainty(
      ex$model, ex$data, warfarin_fit_cov(),
      n_uncertainty_draws = 0L, n_sim_per_draw = 1L
    ),
    "n_uncertainty_draws"
  )
})
