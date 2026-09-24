
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

# ---- logit_probability draws (#373) ----
# The engine draws in its packed space, where a theta with a non-negative
# lower bound is `log(theta)`, so a `logit_probability` theta (declared on
# (0, 1)) is drawn log-normally with no ceiling at 1. The wrapper cannot change
# the draw, but it warns with the expected share that reaches 1.

.logit_prob_share <- getFromNamespace(".ferx_logit_probability_draw_share", "ferx")
.warn_logit_prob  <- getFromNamespace(".ferx_warn_logit_probability_draws", "ferx")

logit_prob_fit <- function(theta_f = 0.7923592, sd_f = 0.3136337) {
  # Shape of the bundled `bioavailability` fit the issue reports: THETA_F is
  # theta 4, and its packed SD is the relative SE `se / estimate`.
  theta <- c(TVCL = 5.71, TVV = 56.7, TVKA = 1.50, THETA_F = theta_f)
  cov <- diag(c(0.307, 0.302, 0.049, sd_f, 0.2)^2)
  dimnames(cov) <- list(c(names(theta), "ETA_CL"), c(names(theta), "ETA_CL"))
  list(
    theta = theta,
    theta_transforms = c(TVCL = "identity", TVV = "identity",
                         TVKA = "identity", THETA_F = "logit_probability"),
    cov_matrix = cov
  )
}

test_that("draw share for a logit_probability theta is P(log-normal draw >= 1)", {
  share <- .logit_prob_share(logit_prob_fit())
  expect_named(share, "THETA_F")
  # log(0.7923592) / 0.3136337 = -0.742: about 23% of draws reach 1.
  expect_equal(unname(share), pnorm(log(0.7923592) / 0.3136337))
  expect_gt(share, 0.2)
  expect_lt(share, 0.25)
})

test_that("draw share ignores other transforms and degenerate entries", {
  fit <- logit_prob_fit()
  fit$theta_transforms[["THETA_F"]] <- "logit"
  expect_length(.logit_prob_share(fit), 0L)

  # FIXed theta: no variance, no draws leave it.
  fit <- logit_prob_fit(sd_f = 0)
  expect_length(.logit_prob_share(fit), 0L)

  # No covariance matrix, or no transforms (an old fit): nothing to report.
  fit <- logit_prob_fit()
  fit$cov_matrix <- NULL
  expect_length(.logit_prob_share(fit), 0L)
  fit <- logit_prob_fit()
  fit$theta_transforms <- NULL
  expect_length(.logit_prob_share(fit), 0L)
})

test_that("the warning names the theta and share, and is classed", {
  w <- expect_warning(
    .warn_logit_prob(logit_prob_fit(), "ferx_simulate_with_uncertainty"),
    class = "ferx_logit_probability_draws"
  )
  expect_match(conditionMessage(w), "logit_probability theta THETA_F 22.9%.", fixed = TRUE)
  expect_match(conditionMessage(w), "method = \"sir\"", fixed = TRUE)
  expect_match(conditionMessage(w), "inv_logit(LOGIT_F + ETA_F)", fixed = TRUE)
})

test_that("no warning when the expected share is negligible", {
  # A tight estimate far from 1: pnorm(log(0.3) / 0.05) is ~0.
  expect_no_warning(
    .warn_logit_prob(logit_prob_fit(theta_f = 0.3, sd_f = 0.05), "f")
  )
})

test_that("asymptotic simulation warns for an exposed logit_probability theta", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex <- ferx_example("warfarin")
  # Relabel TVCL (~0.13, bounds 0.001-10) as a logit_probability theta and
  # widen its packed SD to 1: pnorm(log(0.13) / 1) is ~2%, and the engine's
  # draws stay well inside TVCL's box, so the simulation itself runs normally.
  fit$theta_transforms[[1]] <- "logit_probability"
  fit$cov_matrix[1, 1] <- 1
  expect_warning(
    sims <- ferx_simulate_with_uncertainty(
      ex$model, ex$data, fit,
      n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
      method = "asymptotic", seed = 7L
    ),
    class = "ferx_logit_probability_draws"
  )
  expect_s3_class(sims, "data.frame")
})

test_that("asymptotic simulation without a logit_probability theta does not warn", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  ex <- ferx_example("warfarin")
  expect_no_warning(
    ferx_simulate_with_uncertainty(
      ex$model, ex$data, fit,
      n_uncertainty_draws = 2L, n_sim_per_draw = 1L,
      method = "asymptotic", seed = 7L
    ),
    class = "ferx_logit_probability_draws"
  )
})
