# The legacy reader-settings guard after #530: `.ferx_refuse_unrecorded_selection()`
# on fits whose bytes record no `reader_settings`.

# -- Gap 2: first-match attribution -------------------------------------------

test_that("the call's clause is credited with rows the file's clause also drops (#530)", {
  fit <- rs_accept_fit()
  expect_identical(fit$exclusions$fired_ignore, "ignore: TIME > 110")
  expect_length(fit$exclusions$fired_accept, 0L)
})

test_that("a fit with a fingerprint and the file's population is not refused (#530)", {
  base <- rs_accept_fit()
  skip_if(is.null(base$cov_matrix), "covariance step did not converge - skipping")
  fit <- rs_fingerprint_only(base)
  expect_null(fit$reader_settings)
  expect_false(is.null(fit$population_fingerprint))
  out <- rs_steps(fit)
  for (step in names(out)) {
    expect_false(inherits(out[[step]], "error"),
                 info = paste(step, if (inherits(out[[step]], "error")) conditionMessage(out[[step]])))
  }
  # Same rows as the fit: the standalone covariance is the in-fit one.
  expect_lt(rs_rel(out$covariance$se_theta, base$se_theta), 1e-8)
})

test_that("a fit with a fingerprint and another population is still refused, naming the clause (#530)", {
  fit <- rs_fingerprint_only(rs_fit(ignore = rs_ignore))
  out <- rs_steps(fit)
  for (step in names(out)) {
    e <- out[[step]]
    expect_s3_class(e, "error")
    expect_match(conditionMessage(e), "`ignore: EVID == 0 && DV < 1.0`", fixed = TRUE,
                 info = step)
    expect_match(conditionMessage(e), "predates the record of the data selection",
                 fixed = TRUE, info = step)
  }
})

test_that("with no fingerprint the first-match label still decides (#530)", {
  # Nothing to compare the re-read population with, so the guard refuses.
  fit <- rs_legacy(rs_accept_fit())
  e <- tryCatch(suppressWarnings(ferx_covariance(fit)), error = function(e) e)
  expect_s3_class(e, "error")
  expect_match(conditionMessage(e), "`ignore: TIME > 110`", fixed = TRUE)
})

test_that("a `data =` override is not compared with the fit's fingerprint (#530)", {
  fit <- rs_fingerprint_only(rs_accept_fit())
  data <- withr::local_tempfile(fileext = ".csv")
  file.copy(ferx_example("warfarin")$data, data)
  e <- tryCatch(suppressWarnings(ferx_calc_npde(fit, data = data, nsim = 20L, seed = 5L)),
                error = function(e) e)
  # A copy of the fit's data is still another file: npde has no fingerprint to
  # check it against, so the guard names the clause.
  expect_s3_class(e, "error")
  expect_match(conditionMessage(e), "`ignore: TIME > 110`", fixed = TRUE)
})

# -- Gap 1: engine-written bundles --------------------------------------------
#
# Format version 2 (ferx-core #2695) is one bundle for both writers, and R's fit
# bytes are core's: the engine writes `reader_settings`, the fingerprint and
# `exclusions` into fit.json for every fit that has them, and a bundle of an
# older format does not load. So no loadable bundle reads as "no selection"
# while it had one.

test_that("a fit's bytes carry reader_settings, fingerprint and exclusions in core's fit.json (#530)", {
  fit <- rs_fit(ignore = rs_ignore)
  wire <- fitrx_fit_json(fit$bundle)
  expect_true(all(c("reader_settings", "population_fingerprint", "exclusions") %in% names(wire)))
  expect_identical(unlist(wire$exclusions$fired_ignore), "ignore: EVID == 0 && DV < 1.0")
  loaded <- ferx_load_fit({
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    path
  })
  expect_false(is.null(loaded$reader_settings))
  expect_identical(loaded$exclusions$fired_ignore, fit$exclusions$fired_ignore)
})
