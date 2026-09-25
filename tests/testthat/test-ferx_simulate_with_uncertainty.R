
# ---- header from test-simulate-uncertainty.R ----
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
test_that("uncertainty-aware sims label each parameter draw and vary across them", {
  # Structural check on the uncertainty pipeline: the API contract says one
  # DRAW index per parameter set, and the engine must actually perturb state
  # between draws (theta and/or eta). A statistical IQR-widening claim was
  # tried here previously but is underpowered - parameter SEs contribute a
  # sub-percent share of total DV_SIM variance versus BSV + residual error,
  # so the directional comparison flips with Monte-Carlo noise.
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex <- ferx_example("warfarin")
  sims <- ferx_simulate_with_uncertainty(
    ex$model, ex$data, fit,
    n_uncertainty_draws = 8L, n_sim_per_draw = 1L, seed = 1L
  )
  expect_true("DRAW" %in% names(sims))
  expect_equal(length(unique(sims$DRAW)), 8L)
  # At a single (ID, TIME), IPRED across DRAWs would be identical if the
  # uncertainty layer were a no-op (same theta, same eta seeds per draw).
  # A non-zero spread proves it perturbs at least one of them.
  first <- sims[sims$ID == sims$ID[1] & sims$TIME == sims$TIME[1], ]
  expect_gt(stats::sd(first$IPRED), 0)
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

# ---- header from test-simulate-validate.R ----
# Tests for the internal FFI-payload validators in simulate.R. These guard the
# shapes pulled out of a ferx_fit before handing them to the Rust simulator.

.validate_params      <- getFromNamespace("validate_fit_for_params",      "ferx")
.validate_uncertainty <- getFromNamespace("validate_fit_for_uncertainty", "ferx")










test_that("validate_fit_for_uncertainty (asymptotic) needs a non-empty cov matrix", {
  expect_error(.validate_uncertainty(list(cov_matrix = NULL), "asymptotic"),
               "cov_matrix` is empty")
})
test_that("validate_fit_for_uncertainty (asymptotic) requires a square cov matrix", {
  expect_error(
    .validate_uncertainty(list(cov_matrix = matrix(1:6, 2, 3)), "asymptotic"),
    "square matrix"
  )
})
test_that("validate_fit_for_uncertainty (asymptotic) flattens a valid cov matrix", {
  cov <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  out <- .validate_uncertainty(list(cov_matrix = cov), "asymptotic")
  expect_identical(out$cov_matrix_dim, 2L)
  expect_identical(out$cov_matrix_flat, as.numeric(t(cov)))
})
test_that("validate_fit_for_uncertainty (SIR) errors when resamples are empty", {
  expect_error(.validate_uncertainty(list(), "sir"), "sir_resamples` is empty")
})
test_that("validate_fit_for_uncertainty (SIR) passes through resamples", {
  fit <- list(sir_resamples = c(1, 2, 3, 4), sir_resamples_n = 2L,
              sir_resamples_dim = 2L)
  out <- .validate_uncertainty(fit, "sir")
  expect_identical(out$sir_resamples_n, 2L)
  expect_identical(out$sir_resamples_flat, c(1, 2, 3, 4))
})
