# ferx_sir() and ferx_covariance() score the population ferx_fit() scored.
#
# ferx-core #1790 (#1783): the standalone steps re-read the data but did not
# prepare it as `fit()` does. On a `log(DV) ~ ...` (LTBS) fit DV was not
# log-transformed, so the step reported `Computed` with standard errors
# hundreds of times too large; on a fit whose occasions come from an
# `iov_occasion` rule no occasions were derived, so the covariance step failed
# and SIR kept one draw. Both are measured against the in-fit step, which ran on
# the prepared population.
#
# #512: the rule is taken from the fit (`fit$iov_occasion`), so a rule passed
# only through `ferx_fit(settings = )` reaches the steps, and survives a
# `.fitrx` round trip. Without it the engine falls back to the model file,
# which has none, and refuses.
#
# Reddened by the pin before #1790 (e2f9f641): P1 and P2 report SEs far from
# the fit's and an SIR effective sample size of about 1. P3/P4 are reddened by a
# glue whose skeleton drops the rule (`iov_occasion: None` in `fit_skeleton`).

pp_scratch <- function(ext, env = parent.frame()) {
  withr::local_tempfile(fileext = ext, .local_envir = env)
}

# warfarin_iov without its OCC column, so occasions can only come from a
# rule; the rule `dose` gives the two occasions the column did.
pp_iov_data <- function(env = parent.frame()) {
  rows <- utils::read.csv(ferx_example("warfarin_iov")$data, na.strings = ".")
  rows$OCC <- NULL
  path <- pp_scratch(".csv", env)
  utils::write.csv(rows, path, row.names = FALSE, quote = FALSE, na = ".")
  path
}

# warfarin_iov with its `iov_column` line replaced by `line` (NULL: dropped).
pp_iov_model <- function(line, env = parent.frame()) {
  text <- readLines(ferx_example("warfarin_iov")$model)
  at <- grep("iov_column", text, fixed = TRUE)
  text <- if (is.null(line)) text[-at] else replace(text, at, line)
  path <- pp_scratch(".ferx", env)
  writeLines(text, path)
  path
}

pp_fit <- function(model, data, ...) {
  suppressWarnings(ferx_fit(model, data, covariance = TRUE, verbose = FALSE, ...))
}

pp_sir <- function(fit) {
  ferx_sir(fit, sir_samples = 300L, sir_resamples = 100L, sir_seed = 1L)
}

# The standalone step agrees with the in-fit one, and SIR is not degenerate.
pp_expect_scored <- function(fit) {
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  cov <- ferx_covariance(fit)
  expect_identical(cov$covariance_status, "computed")
  expect_true(all(is.finite(cov$se_theta)))
  expect_equal(cov$se_theta, fit$se_theta, tolerance = 1e-4)
  sir <- pp_sir(fit)
  expect_true(all(is.finite(sir$sir_ci_theta)))
  # Measured at 826d3bb9: 105 (LTBS) and 109 (iov_occasion) of 100 resamples
  # from 300 draws; 1 before #1790 on the IOV fit.
  expect_gt(sir$sir_ess, 30)
}

test_that("P1: ferx_covariance() / ferx_sir() on an LTBS fit score log(DV)", {
  ex <- ferx_example("warfarin_ltbs")
  pp_expect_scored(pp_fit(ex$model, ex$data))
})

test_that("P2: ferx_covariance() / ferx_sir() derive a model-file iov_occasion's occasions", {
  fit <- pp_fit(pp_iov_model("  iov_occasion = dose"), pp_iov_data())
  expect_identical(fit$iov_occasion, "dose")
  pp_expect_scored(fit)
})

test_that("P3: an iov_occasion passed only through settings = reaches the steps (#512)", {
  data <- pp_iov_data()
  fit <- pp_fit(pp_iov_model(NULL), data, sir = TRUE,
                settings = list(iov_occasion = "dose", sir_samples = 300L,
                                sir_resamples = 100L, sir_seed = 1L))
  expect_identical(fit$iov_occasion, "dose")
  pp_expect_scored(fit)
  # The same SIR as the fit's own, draw for draw.
  expect_identical(pp_sir(fit)$sir_ci_theta, fit$sir_ci_theta)
  # And the same fit as the one that read its occasions from the column.
  col <- pp_fit(ferx_example("warfarin_iov")$model, ferx_example("warfarin_iov")$data)
  expect_identical(col$iov_occasion, "column")
  expect_equal(fit$ofv, col$ofv, tolerance = 1e-8)

  # Without the rule the engine has nothing to derive occasions with.
  bare <- fit
  bare$iov_occasion <- NULL
  expect_error(ferx_covariance(bare), "records no IOV occasion rule", fixed = TRUE)
})

test_that("P4: the rule survives ferx_save_fit() / ferx_load_fit() (#512)", {
  data <- pp_iov_data()
  for (rule in c("dose", "time(120.1)")) {
    fit <- pp_fit(pp_iov_model(NULL), data, settings = list(iov_occasion = rule))
    expect_identical(fit$iov_occasion, rule)
    path <- pp_scratch(".fitrx")
    unlink(path)
    ferx_save_fit(fit, path)
    loaded <- ferx_load_fit(path)
    expect_identical(loaded$iov_occasion, rule, info = rule)
    # Only the covariance comparison needs a converged step; a skip here would
    # also drop the next rule's round trip. The bundle carries the packed
    # estimate (ferx-core #1815), so the loaded fit repeats the in-memory step.
    if (!is.null(fit$cov_matrix)) {
      expect_identical(ferx_covariance(loaded)$se_theta,
                       ferx_covariance(fit)$se_theta, info = rule)
    }
  }
})

test_that("P5: the .fitrx spelling of the rule is ferx-core's", {
  to <- ferx:::.fitrx_iov_occasion_to_wire
  from <- ferx:::.fitrx_iov_occasion_from_wire
  expect_null(to(NULL))
  expect_identical(to("column"), "column")
  expect_identical(to("dose"), "per_dose")
  expect_identical(to("time(24, 48.5)"), list(time_windows = list(24, 48.5)))
  # One edge stays an array on the wire.
  expect_identical(
    as.character(jsonlite::toJSON(to("time(24)"), auto_unbox = TRUE)),
    '{"time_windows":[24]}'
  )
  expect_identical(from("per_dose"), "dose")
  expect_identical(from(list(time_windows = list(24, 48.5))), "time(24, 48.5)")
  expect_identical(from(list(time_windows = list(0.1 + 0.2))), "time(0.30000000000000004)")
  # No exponent, as Rust's `{}` writes none: round edges are where R's
  # format() would switch to scientific notation.
  expect_identical(from(list(time_windows = list(1e5, 2e6))), "time(100000, 2000000)")
  expect_identical(from(list(time_windows = list(1e-5))), "time(0.00001)")
  expect_error(to("weekly"), "not \"column\"", fixed = TRUE)
  expect_error(from("weekly"), "is not \"column\"", fixed = TRUE)
})

# #523 / #519: the from-fit predict / simulate / npde paths read the data with
# readers that fill occasions from `iov_column` only, so under a derived rule
# every subject came back unlabelled and the engine ran each occasion at
# kappa = 0, with no message. Each path is measured against the twin that reads
# the same occasions from the OCC column, at the same estimates: before the fix
# the simulated IPRED differed by up to 8.5 and the NPDE by up to 3.4.

# The column twin, fitted with a covariance step, and a fit made under `rule`
# from the column-less data, carrying the twin's estimates so that only the
# occasions can differ between the two.
pp_twins <- function(rule, env = parent.frame()) {
  col <- pp_fit(ferx_example("warfarin_iov")$model, ferx_example("warfarin_iov")$data)
  data <- pp_iov_data(env)
  model <- pp_iov_model(NULL, env)
  fit <- pp_fit(model, data, settings = list(iov_occasion = rule))
  expect_identical(fit$iov_occasion, rule)
  for (k in c("theta", "omega", "sigma", "omega_iov", "cov_matrix")) fit[[k]] <- col[[k]]
  list(col = col, fit = fit, model = model, data = data)
}

test_that("P6: from-fit ferx_simulate() derives the fit's occasions (#523)", {
  for (rule in c("dose", "time(120.1)")) {
    tw <- pp_twins(rule)
    sim <- ferx_simulate(tw$model, tw$data, fit = tw$fit, n_sim = 2L, seed = 7L)
    ref <- ferx_simulate(ferx_example("warfarin_iov")$model,
                         ferx_example("warfarin_iov")$data,
                         fit = tw$col, n_sim = 2L, seed = 7L)
    # Under the window rule a handful of records land a few ULPs from the
    # column twin's (measured: 2e-15 relative), so not identical.
    expect_equal(sim$IPRED, ref$IPRED, tolerance = 1e-10, info = rule)
    expect_equal(sim$DV_SIM, ref$DV_SIM, tolerance = 1e-10, info = rule)
  }
})

test_that("P7: ferx_calc_npde() derives the fit's occasions (#523)", {
  tw <- pp_twins("dose")
  npde <- ferx_calc_npde(tw$fit, nsim = 200L, seed = 3L)
  ref <- ferx_calc_npde(tw$col, nsim = 200L, seed = 3L)
  expect_identical(npde$sdtab$NPDE, ref$sdtab$NPDE)
  expect_identical(npde$sdtab$NPD, ref$sdtab$NPD)
})

test_that("P8: ferx_simulate_with_uncertainty() derives the fit's occasions (#523)", {
  tw <- pp_twins("dose")
  skip_if(is.null(tw$col$cov_matrix), "covariance step did not converge - skipping")
  sim <- ferx_simulate_with_uncertainty(tw$model, tw$data, tw$fit,
                                        n_uncertainty_draws = 3L, seed = 11L)
  ref <- ferx_simulate_with_uncertainty(ferx_example("warfarin_iov")$model,
                                        ferx_example("warfarin_iov")$data,
                                        tw$col, n_uncertainty_draws = 3L, seed = 11L)
  expect_identical(sim$IPRED, ref$IPRED)
})

test_that("P9: from-fit ferx_predict() runs on the fit's occasions (#523)", {
  # PRED is at eta = kappa = 0, so the occasions cannot move it: this pins
  # that the path takes the rule without refusing, not a numeric change.
  tw <- pp_twins("dose")
  pred <- ferx_predict(tw$model, tw$data, fit = tw$fit)
  ref <- ferx_predict(ferx_example("warfarin_iov")$model,
                      ferx_example("warfarin_iov")$data, fit = tw$col)
  expect_identical(pred$PRED, ref$PRED)
})

test_that("P10: a model-file iov_occasion reaches ferx_simulate() without a fit (#523)", {
  model <- pp_iov_model("  iov_occasion = dose")
  sim <- ferx_simulate(model, pp_iov_data(), n_sim = 2L, seed = 5L)
  ref <- ferx_simulate(ferx_example("warfarin_iov")$model,
                       ferx_example("warfarin_iov")$data, n_sim = 2L, seed = 5L)
  expect_identical(sim$IPRED, ref$IPRED)
})

test_that("P11: a fit that records no rule falls back to the model file's, or is refused (#523)", {
  tw <- pp_twins("dose")
  legacy <- tw$fit
  legacy$iov_occasion <- NULL
  # No rule in the model file and no OCC column: nothing to derive occasions with.
  expect_error(ferx_simulate(tw$model, tw$data, fit = legacy, n_sim = 1L, seed = 1L),
               "records no IOV occasion rule", fixed = TRUE)
  expect_error(ferx_calc_npde(legacy, nsim = 50L, seed = 1L),
               "records no IOV occasion rule", fixed = TRUE)
  # The model file's rule is used, and the result says it was a guess.
  file_model <- pp_iov_model("  iov_occasion = dose")
  expect_warning(
    sim <- ferx_simulate(file_model, tw$data, fit = legacy, n_sim = 1L, seed = 1L),
    "derived with the model file's `iov_occasion = dose`", fixed = TRUE
  )
  ref <- suppressWarnings(
    ferx_simulate(tw$model, tw$data, fit = tw$fit, n_sim = 1L, seed = 1L)
  )
  expect_identical(sim$IPRED, ref$IPRED)
})

test_that("P12: from-fit ferx_predict() reads fit$iov_occasion (#523 review r1)", {
  # PRED cannot show the occasions (P9), so the path's use of the rule is
  # pinned by what it does with a rule it cannot use, and with none at all.
  tw <- pp_twins("dose")
  bad <- tw$fit
  bad$iov_occasion <- "weekly"
  expect_error(ferx_predict(tw$model, tw$data, fit = bad),
               "ferx_predict: `fit$iov_occasion`", fixed = TRUE)
  legacy <- tw$fit
  legacy$iov_occasion <- NULL
  expect_error(ferx_predict(tw$model, tw$data, fit = legacy),
               "records no IOV occasion rule", fixed = TRUE)
})

test_that("P13: from-fit ferx_predict_survival() reads fit$iov_occasion (#523 review r1)", {
  # No bundled TTE model declares kappa, so the rule cannot move a survival
  # curve here; the path is pinned by refusing a rule it cannot parse.
  model <- pp_scratch(".ferx")
  data <- pp_scratch(".csv")
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
  fit <- suppressWarnings(ferx_fit(model, data, method = "focei", verbose = FALSE))
  expect_identical(fit$iov_occasion, "column")
  expect_no_error(ferx_predict_survival(model, data, times = c(0, 6), fit = fit))
  fit$iov_occasion <- "weekly"
  expect_error(ferx_predict_survival(model, data, times = c(0, 6), fit = fit),
               "ferx_predict_survival: `fit$iov_occasion`", fixed = TRUE)
})

test_that("P14: a model-file rule that overrides iov_column says so without a fit (#523 review r1)", {
  text <- readLines(ferx_example("warfarin_iov")$model)
  at <- grep("iov_column", text, fixed = TRUE)
  model <- pp_scratch(".ferx")
  writeLines(append(text, "  iov_occasion = dose", after = at), model)
  ex <- ferx_example("warfarin_iov")
  msg <- "the iov_column dataset column is ignored"
  expect_warning(ferx_simulate(model, ex$data, n_sim = 1L, seed = 1L), msg, fixed = TRUE)
  expect_warning(ferx_predict(model, ex$data), msg, fixed = TRUE)
  # Without the rule the column is the source, and nothing is said.
  expect_no_warning(ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L))
})
