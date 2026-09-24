
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
# (0, 1)) is drawn log-normally with no ceiling at 1. Past a declared upper
# bound <= 1 the draw is rejected (truncation); with an upper bound > 1 it is
# simulated and `logit()` clamps it to a probability of 1. The wrapper cannot
# change the draw (ferx-core #1548), but it reports it.

.lp_share   <- getFromNamespace(".ferx_logit_probability_draw_share", "ferx")
.lp_warn    <- getFromNamespace(".ferx_warn_logit_probability_draws", "ferx")
.lp_packing <- getFromNamespace(".ferx_theta_packing", "ferx")

# Estimate and packed SD of THETA_F on the bundled `bioavailability` fit, as
# the issue reports them (the SD is the relative SE `se / estimate`).
lp_est <- 0.7923592
lp_sd  <- 0.3136337

lp_packing <- function(upper = 0.999, lower = 0.001) {
  list(names = c("TVCL", "TVV", "TVKA", "THETA_F"),
       lower = c(0.1, 5, 0.05, lower),
       upper = c(50, 500, 20, upper),
       transform = c("identity", "identity", "identity", "logit_probability"))
}

lp_fit <- function(theta_f = lp_est, sd_f = lp_sd, other_sd = 1e-6) {
  theta <- c(TVCL = 5.7131442, TVV = 56.6724035, TVKA = 1.5000921,
             THETA_F = theta_f)
  # Packed layout: 4 thetas, 2 omega diagonals (log Cholesky), 1 sigma.
  pn <- c(names(theta), "ETA_CL", "ETA_F", "PROP_ERR")
  cov <- diag(c(rep(other_sd, 3), sd_f, rep(other_sd, 3))^2)
  dimnames(cov) <- list(pn, pn)
  list(
    theta = theta,
    # The smallest variances the packed box allows (log L_ii >= -6), so every
    # subject carries F almost exactly and IPRED is proportional to it.
    omega = diag(c(1e-5, 1e-5), 2, 2,
                 names = list(c("ETA_CL", "ETA_F"), c("ETA_CL", "ETA_F"))),
    sigma = c(PROP_ERR = 0.15),
    theta_transforms = c(TVCL = "identity", TVV = "identity",
                         TVKA = "identity", THETA_F = "logit_probability"),
    cov_matrix = cov
  )
}

test_that("asymptotic share, declared upper <= 1: draws past it are rejected", {
  ex <- .lp_share(lp_fit(), lp_packing(upper = 0.999), "asymptotic")
  expect_identical(ex$theta, "THETA_F")
  expect_identical(ex$outcome, "rejected")
  expect_equal(ex$share, pnorm((log(lp_est) - log(0.999)) / lp_sd))
  # An upper bound of exactly 1 still rejects: log(1) = 0 is the packed bound.
  ex1 <- .lp_share(lp_fit(), lp_packing(upper = 1), "asymptotic")
  expect_identical(ex1$outcome, "rejected")
  expect_equal(ex1$share, pnorm(log(lp_est) / lp_sd))
})

test_that("asymptotic share, declared upper > 1: draws in (1, upper] are clamped", {
  ex <- .lp_share(lp_fit(), lp_packing(upper = 5), "asymptotic")
  expect_identical(ex$outcome, "clamped")
  expect_equal(ex$share, pnorm(log(lp_est) / lp_sd) -
                 pnorm((log(lp_est) - log(5)) / lp_sd))
  # ~22.9% on the bioavailability fit, as the issue predicts.
  expect_gt(ex$share, 0.2)
  expect_lt(ex$share, 0.25)
  # The parser's default upper bound (1e9): the whole tail above 1 is clamped.
  exd <- .lp_share(lp_fit(), lp_packing(upper = 1e9), "asymptotic")
  expect_equal(exd$share, pnorm(log(lp_est) / lp_sd), tolerance = 1e-12)
})

test_that("asymptotic share uses the natural scale for an identity-packed theta", {
  # A negative lower bound packs the theta as itself, so cov is absolute.
  ex <- .lp_share(lp_fit(sd_f = 0.2), lp_packing(upper = 5, lower = -1),
                  "asymptotic")
  expect_equal(ex$share, pnorm((lp_est - 1) / 0.2) - pnorm((lp_est - 5) / 0.2))
})

test_that("draw share ignores other transforms and degenerate entries", {
  pk <- lp_packing(upper = 5)
  pk$transform[4] <- "logit"
  expect_equal(nrow(.lp_share(lp_fit(), pk, "asymptotic")), 0L)
  # FIXed theta: no variance, no draws leave it.
  expect_equal(nrow(.lp_share(lp_fit(sd_f = 0), lp_packing(5), "asymptotic")), 0L)
  # No covariance matrix, or a model that could not be read: nothing to report.
  fit <- lp_fit()
  fit$cov_matrix <- NULL
  expect_equal(nrow(.lp_share(fit, lp_packing(5), "asymptotic")), 0L)
  expect_equal(nrow(.lp_share(lp_fit(), NULL, "asymptotic")), 0L)
})

test_that("SIR share is counted in the pool, and only clamping can show", {
  fit <- lp_fit()
  pool <- matrix(0, nrow = 10, ncol = 7)
  pool[, 4] <- log(c(rep(0.7, 7), 1.2, 1.5, 3))   # 3 of 10 above 1
  fit$sir_resamples <- as.numeric(t(pool))
  fit$sir_resamples_n <- 10L
  fit$sir_resamples_dim <- 7L
  ex <- .lp_share(fit, lp_packing(upper = 5), "sir")
  expect_identical(ex$outcome, "clamped")
  expect_equal(ex$share, 0.3)
  expect_identical(ex$n, 10L)
  # With the upper bound at or below 1 the pool stays inside it.
  expect_equal(nrow(.lp_share(fit, lp_packing(upper = 0.999), "sir")), 0L)
  # No pool: nothing to count.
  expect_equal(nrow(.lp_share(lp_fit(), lp_packing(upper = 5), "sir")), 0L)
})

test_that("the warning is classed and worded for the bound that applies", {
  w <- expect_warning(
    .lp_warn(lp_fit(), lp_packing(0.999), "asymptotic", "fn"),
    class = "ferx_logit_probability_draws"
  )
  m <- conditionMessage(w)
  expect_match(m, "THETA_F: about 23.0% of draws exceed the declared upper bound 0.999 and are rejected", fixed = TRUE)
  expect_match(m, "ferx-core #1548", fixed = TRUE)
  expect_match(m, "inv_logit(LOGIT_F + ETA_F)", fixed = TRUE)
  expect_no_match(m, "method = \"sir\"", fixed = TRUE)

  w <- expect_warning(
    .lp_warn(lp_fit(), lp_packing(5), "asymptotic", "fn"),
    class = "ferx_logit_probability_draws"
  )
  expect_match(conditionMessage(w),
               "THETA_F: about 22.9% of draws land above 1 (declared upper bound 5)",
               fixed = TRUE)

  fit <- lp_fit()
  pool <- matrix(0, nrow = 4, ncol = 7)
  pool[, 4] <- log(c(0.7, 0.8, 0.6, 1.5))
  fit$sir_resamples <- as.numeric(t(pool))
  fit$sir_resamples_n <- 4L
  fit$sir_resamples_dim <- 7L
  w <- expect_warning(.lp_warn(fit, lp_packing(5), "sir", "fn"),
                      class = "ferx_logit_probability_draws")
  expect_match(conditionMessage(w),
               "THETA_F: 1 of 4 pooled SIR draws (25.0%) are above 1", fixed = TRUE)
})

test_that("no warning when the expected share is negligible", {
  # A tight estimate far from 1: pnorm(log(0.3) / 0.05) is ~0.
  expect_no_warning(
    .lp_warn(lp_fit(theta_f = 0.3, sd_f = 0.05), lp_packing(5),
             "asymptotic", "fn")
  )
})

test_that("theta packing reports the declared bounds and transforms", {
  ex <- ferx_example("bioavailability")
  pk <- .lp_packing(ex$model)
  expect_identical(pk$names, c("TVCL", "TVV", "TVKA", "THETA_F"))
  expect_identical(pk$transform[4], "logit_probability")
  expect_equal(pk$lower[4], 0.001)
  expect_equal(pk$upper[4], 0.999)
  expect_null(.lp_packing(tempfile(fileext = ".ferx")))
})

# ---- engine regression: what the draws of THETA_F actually are ----
# A hand-built fit on the bundled bioavailability model: every parameter but
# THETA_F carries (almost) no uncertainty and the etas are at the packed
# floor, so each draw's IPRED is its F times the point-estimate profile. That
# recovers the drawn F from the engine's output. These tests pin the engine
# behaviour the warning describes; they fail - on purpose - once ferx-core
# #1548 draws the theta on the logit scale, and the warning must go then.

lp_one_subject_data <- function() {
  ex <- ferx_example("bioavailability")
  d <- read.csv(ex$data)
  path <- tempfile(fileext = ".csv")
  write.csv(d[d$ID == d$ID[1], ], path, row.names = FALSE, quote = FALSE)
  path
}

lp_model_with_upper <- function(upper) {
  ex <- ferx_example("bioavailability")
  src <- readLines(ex$model)
  src <- sub("theta THETA_F(0.70, 0.001, 0.999)",
             sprintf("theta THETA_F(0.70, 0.001, %s)", upper), src, fixed = TRUE)
  path <- tempfile(fileext = ".ferx")
  writeLines(src, path)
  path
}

# Drawn F per DRAW: the median IPRED ratio to the point-estimate simulation.
lp_drawn_f <- function(sims, model, data, fit) {
  ref <- suppressWarnings(ferx_simulate(model, data, n_sim = 1L, fit = fit))
  ref <- ref[is.finite(ref$IPRED) & ref$IPRED > 0, c("TIME", "IPRED")]
  m <- merge(sims[, c("DRAW", "TIME", "IPRED")], ref, by = "TIME",
             suffixes = c("", "_ref"))
  r <- tapply(m$IPRED / m$IPRED_ref, m$DRAW, stats::median)
  lp_est * as.numeric(r)
}

test_that("engine: declared upper 0.999 truncates the draws below it", {
  model <- ferx_example("bioavailability")$model
  data  <- lp_one_subject_data()
  fit   <- lp_fit()
  expect_warning(
    sims <- ferx_simulate_with_uncertainty(model, data, fit,
      n_uncertainty_draws = 300L, n_sim_per_draw = 1L, seed = 5L),
    class = "ferx_logit_probability_draws"
  )
  f <- lp_drawn_f(sims, model, data, fit)
  expect_length(f, 300L)
  # In range, but truncated: the draws are log-normal cut at 0.999, not the
  # logit-normal ones the fit implies. The truncated log-normal's median sits
  # below the estimate by the rejected share.
  expect_true(all(f < 0.999 + 0.01))
  expect_lt(mean(f > lp_est), 0.45)
})

test_that("engine: declared upper > 1 simulates the tail above 1 as F = 1", {
  model <- lp_model_with_upper(5)
  data  <- lp_one_subject_data()
  fit   <- lp_fit()
  expect_warning(
    sims <- ferx_simulate_with_uncertainty(model, data, fit,
      n_uncertainty_draws = 300L, n_sim_per_draw = 1L, seed = 5L),
    class = "ferx_logit_probability_draws"
  )
  f <- lp_drawn_f(sims, model, data, fit)
  # No draw exceeds 1: the clamp caps F there ...
  expect_true(all(f < 1.01))
  # ... and the share piled up at F = 1 is the log-normal tail the warning
  # reports (22.9%, binomial SD ~2.4% over 300 draws).
  predicted <- .lp_share(fit, .lp_packing(model), "asymptotic")$share
  expect_equal(mean(f > 0.995), predicted, tolerance = 0.08 / predicted)
})

test_that("engine: a SIR pool draw above 1 is simulated as F = 1", {
  model <- lp_model_with_upper(5)
  data  <- lp_one_subject_data()
  fit   <- lp_fit()
  x_hat <- c(log(fit$theta), log(sqrt(diag(fit$omega))), log(fit$sigma))
  pool  <- rbind(x_hat, x_hat)
  pool[2, 4] <- log(1.5)                       # THETA_F = 1.5 in the pool
  fit$sir_resamples <- as.numeric(t(pool))
  fit$sir_resamples_n <- 2L
  fit$sir_resamples_dim <- ncol(pool)
  expect_warning(
    sims <- ferx_simulate_with_uncertainty(model, data, fit,
      n_uncertainty_draws = 40L, n_sim_per_draw = 1L, method = "sir",
      seed = 5L),
    class = "ferx_logit_probability_draws"
  )
  f <- lp_drawn_f(sims, model, data, fit)
  # Draws come from two pool rows: F = 0.79 or the clamped 1.5 -> 1.
  expect_true(all(abs(f - lp_est) < 0.01 | abs(f - 1) < 0.01))
  expect_true(any(abs(f - 1) < 0.01))
})

test_that("a model that never reaches 1 does not warn", {
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
