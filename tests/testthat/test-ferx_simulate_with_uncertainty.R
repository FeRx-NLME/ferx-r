
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

# ---- logit_probability draws (#373, ferx-core #1548) ----
# A `logit_probability` theta - used as inv_logit(logit(THETA_F) + ETA_F) -
# is drawn on the logit scale: centre logit(theta_hat), SD by the delta method
# (packed SD / (1 - theta_hat) for a log-packed theta). Before ferx-core #1548
# the draw was log-normal: 23% of it sat above 1, and was either rejected
# (upper bound 0.999: truncated, median pulled to ~0.72) or clamped to F = 1
# (upper bound above 1).
#
# The fixture fit pins every other parameter (packed SD 1e-6) and shrinks both
# omegas to the smallest variance the packed box allows, so within a draw every
# subject carries F almost exactly, and IPRED is proportional to it. Dividing
# each draw's IPRED by a reference run at theta_hat recovers the drawn F.

lp_est <- 0.7923592 # THETA_F on the bundled bioavailability fit
lp_sd  <- 0.3136337 # its packed (log-scale) SD, i.e. the relative SE

lp_fit <- function(sd_f = lp_sd) {
  theta <- c(TVCL = 5.7131442, TVV = 56.6724035, TVKA = 1.5000921,
             THETA_F = lp_est)
  # Packed layout: 4 thetas, 2 omega diagonals (log Cholesky), 1 sigma.
  pn  <- c(names(theta), "ETA_CL", "ETA_F", "PROP_ERR")
  cov <- diag(c(rep(1e-6, 3), sd_f, rep(1e-6, 3))^2)
  dimnames(cov) <- list(pn, pn)
  list(
    theta      = theta,
    omega      = diag(c(1e-5, 1e-5), 2, 2,
                      names = list(c("ETA_CL", "ETA_F"), c("ETA_CL", "ETA_F"))),
    sigma      = c(PROP_ERR = 0.15),
    cov_matrix = cov
  )
}

# The bundled model with THETA_F's upper bound replaced by `upper`.
lp_model <- function(upper) {
  ex  <- ferx_example("bioavailability")
  src <- readLines(ex$model, warn = FALSE)
  i   <- grep("theta THETA_F(", src, fixed = TRUE)
  src[i] <- sprintf("  theta THETA_F(0.70, 0.001, %s)", format(upper))
  path <- tempfile(fileext = ".ferx")
  writeLines(src, path)
  normalizePath(path)
}

# One drawn F per uncertainty draw.
lp_drawn_f <- function(model, n_draws = 400L) {
  ex  <- ferx_example("bioavailability")
  ref <- ferx_simulate_with_uncertainty(
    model, ex$data, lp_fit(sd_f = 1e-6),
    n_uncertainty_draws = 1L, n_sim_per_draw = 1L, seed = 3L
  )
  sims <- ferx_simulate_with_uncertainty(
    model, ex$data, lp_fit(),
    n_uncertainty_draws = n_draws, n_sim_per_draw = 1L, seed = 11L
  )
  keep <- ref$IPRED > 0
  vapply(split(sims$IPRED, sims$DRAW), function(ipred) {
    lp_est * stats::median(ipred[keep] / ref$IPRED[keep])
  }, numeric(1))
}

expect_logit_normal_draws <- function(f) {
  expect_true(all(f > 0 & f < 1))
  # No point mass at F = 1: P(F > 0.99) is ~1.5% under the logit-normal draw,
  # >= 23% under the old log-normal one when the upper bound is above 1.
  expect_lt(mean(f > 0.99), 0.06)
  # Logit-normal: the median is theta_hat (the truncated draw's was ~0.72) ...
  expect_lt(abs(stats::median(f) - lp_est), 0.04)
  # ... and the logit-scale IQR is 2 * qnorm(0.75) * the delta-method SD.
  y   <- stats::qlogis(f)
  iqr <- unname(diff(stats::quantile(y, c(0.25, 0.75))))
  expect_equal(iqr, 2 * stats::qnorm(0.75) * lp_sd / (1 - lp_est),
               tolerance = 0.15)
}

test_that("asymptotic logit_probability draws are logit-normal (upper 0.999)", {
  ex <- ferx_example("bioavailability")
  expect_logit_normal_draws(lp_drawn_f(ex$model))
})

test_that("asymptotic logit_probability draws are logit-normal (upper above 1)", {
  expect_logit_normal_draws(lp_drawn_f(lp_model(5)))
})
