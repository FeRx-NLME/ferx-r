# Tests for ferx_simulate_with_uncertainty. The cached `warfarin_fit_cov()`
# helper provides a fit with `covariance = TRUE`; in CI's no-autodiff build
# the covariance step occasionally fails to converge with maxiter = 30 (same
# limitation as the existing tests in test-fit.R, which already use
# `skip_if(is.null(fit$cov_matrix), ...)`). The asymptotic happy-path tests
# below follow that pattern. SIR-path tests are gated on `sir_resamples`,
# which the fixture never populates, so the negative SIR test runs
# unconditionally.

cov_skip <- "covariance step did not converge — skipping"

test_that("ferx_simulate_with_uncertainty asymptotic returns a data frame", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 7L
  )
  expect_s3_class(sims, "data.frame")
})

test_that("asymptotic output has the expected columns and DRAW leads", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 7L
  )
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "IPRED", "DV_SIM") %in%
                    names(sims)))
  expect_identical(names(sims)[1], "DRAW")
})

test_that("DRAW spans 1..n_uncertainty_draws", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 4L, n_sim_per_draw = 1L,
    method = "asymptotic", seed = 11L
  )
  expect_equal(sort(unique(sims$DRAW)), 1L:4L)
})

test_that("SIM spans 1..n_sim_per_draw inside each draw", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 2L, n_sim_per_draw = 3L,
    method = "asymptotic", seed = 21L
  )
  for (d in unique(sims$DRAW)) {
    expect_equal(sort(unique(sims$SIM[sims$DRAW == d])), 1L:3L)
  }
})

test_that("row count is n_uncertainty_draws * n_sim_per_draw * n_obs", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex    <- ferx_example("warfarin")
  dat   <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  sims  <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 3L, n_sim_per_draw = 2L,
    method = "asymptotic", seed = 33L
  )
  expect_equal(nrow(sims), 3L * 2L * n_obs)
})

test_that("DV_SIM and IPRED are finite numerics with no NAs", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
    method = "asymptotic", seed = 41L
  )
  expect_true(is.numeric(sims$DV_SIM))
  expect_true(is.numeric(sims$IPRED))
  expect_true(all(is.finite(sims$DV_SIM)))
  expect_true(all(is.finite(sims$IPRED)))
})

test_that("same seed produces identical output (asymptotic)", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex    <- ferx_example("warfarin")
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
  # given (ID, TIME). Both calls use the same total sample size (150) and
  # the same seed, so the only difference is whether the population
  # parameters are perturbed between draws — which means `many` must show a
  # strictly wider IQR than `few`. A small safety margin (`* 1.02`)
  # guards against ULP-level ties without weakening the directional claim.
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex   <- ferx_example("warfarin")
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
  expect_gt(iqr_many, iqr_few * 1.02)
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
    ferx_simulate_with_uncertainty("no_such_model.ferx", ex$data, warfarin_fit_cov()),
    "file.exists"
  )
  expect_error(
    ferx_simulate_with_uncertainty(ex$model, "no_such_data.csv", warfarin_fit_cov()),
    "file.exists"
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

test_that("n_sim_per_draw < 1 raises an informative error", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_simulate_with_uncertainty(
      ex$model, ex$data, warfarin_fit_cov(),
      n_uncertainty_draws = 2L, n_sim_per_draw = 0L
    ),
    "n_sim_per_draw"
  )
})
