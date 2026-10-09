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
    # also drop the next rule's round trip. A bundle carries no packed estimate
    # (ferx-core#1815), so the in-memory side drops it too (#511).
    if (!is.null(fit$cov_matrix)) {
      in_memory <- fit
      in_memory$packed_estimate <- NULL
      expect_identical(ferx_covariance(loaded)$se_theta,
                       ferx_covariance(in_memory)$se_theta, info = rule)
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
