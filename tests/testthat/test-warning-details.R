# A warning's structured `details` payload (ferx-core #2268): `fit$warning_details`
# live, beside its row in ferx_get_warnings(as_df = TRUE), through
# ferx_save_fit() / ferx_load_fit() in ferx-core's own `.fitrx` slot, and through
# the post-hoc steps. warfarin's fit carries one payload, the Durbin-Watson
# statistic of its `dw_autocorrelation` warning; the oracle for the bundle is
# ferx-core's own loader and writer, not R reading back what R wrote.

wd_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
    }
    fit
  }
})

wd_dw <- function(fit) {
  Filter(function(e) identical(e$category, "dw_autocorrelation"), fit$warning_details)
}

test_that("a fit carries each engine warning's details, matched to its message", {
  skip_on_cran()
  fit <- wd_fit()
  dw <- wd_dw(fit)
  expect_length(dw, 1L)
  e <- dw[[1L]]
  expect_setequal(names(e$details), c("durbin_watson", "iwres_lag1_autocorr"))
  expect_true(e$message %in% fit$warnings)
  # The statistic the message quotes, to its printed two decimals.
  expect_match(e$message, sprintf("Durbin-Watson = %.2f", e$details$durbin_watson),
               fixed = TRUE)
  # Only entries with a payload are kept; every one has the WarningEntry fields.
  for (x in fit$warning_details) {
    expect_true(all(c("severity", "category", "message", "details") %in% names(x)))
  }
})

test_that("ferx_get_warnings() puts the details beside their row", {
  skip_on_cran()
  fit <- wd_fit()
  df <- ferx_get_warnings(fit, as_df = TRUE)
  expect_named(df, c("severity", "category", "message", "source_method", "details"))
  i <- which(df$category == "dw_autocorrelation")
  expect_length(i, 1L)
  expect_identical(df$details[[i]], wd_dw(fit)[[1L]]$details)
  expect_true(all(vapply(df$details[-i], is.null, logical(1L))))
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("details: durbin_watson, iwres_lag1_autocorr", out, fixed = TRUE)))
  # A fit with no warning at all still gets the column.
  bare <- fit
  bare <- unsealed_fit(bare)
  bare$warnings <- character()
  bare$warnings_structured <- NULL
  bare$warning_details <- NULL
  df0 <- ferx_get_warnings(bare, as_df = TRUE)
  expect_identical(nrow(df0), 0L)
  expect_true("details" %in% names(df0))
})

test_that("the details survive ferx_save_fit() / ferx_load_fit() in ferx-core's slot", {
  skip_on_cran()
  fit <- wd_fit()
  # ferx-core's key, one WarningEntry object per payload.
  wire <- fitrx_fit_json(fit$bundle)
  expect_length(wire$warning_details, length(fit$warning_details))
  # The bundle loads to the fit's own list, bit for bit.
  path <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, path))
  expect_identical(suppressWarnings(ferx_load_fit(path))$warning_details,
                   fit$warning_details)
})

test_that("a bundle with an unknown warning category still loads", {
  skip_on_cran()
  fit <- wd_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, path))
  # A category this ferx does not know (a newer ferx-core): core's loader keeps
  # the bundle, and the message stays.
  newer <- fitrx_edit_bundle(path, function(w) {
    w$warning_details[[1L]]$category <- "a_future_code"
    w
  })
  back <- suppressWarnings(ferx_load_fit(newer))
  expect_identical(back$warnings, fit$warnings)
})

test_that("a whole-number float in the details reloads as a float", {
  # `details` is an untyped serde_json::Value, which keeps `50` an integer and
  # `50.0` a float; a loaded `threshold_pct = 50.0` must not come back `50L`,
  # which is not identical() (#646 review).
  json <- paste0('{"severity":"Warning","category":"inflated_rse","message":"m",',
                 '"details":{"threshold_pct":50.0,"n":3,"x":2.5}}')
  # A live entry is the engine's own serialisation of it.
  live <- ferx:::.ferx_parse_warning_details(ferx:::ferx_rust_warning_details_canonical(json))
  expect_identical(live[[1L]]$details$threshold_pct, 50)
  expect_identical(live[[1L]]$details$n, 3L)

  # And one in a bundle, spelled as the engine spells it.
  skip_on_cran()
  fit <- wd_fit()
  path <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, path))
  staging <- withr::local_tempdir()
  utils::unzip(path, files = "fit.json", exdir = staging)
  json <- file.path(staging, "fit.json")
  text <- readLines(json, warn = FALSE)
  edited <- sub('("durbin_watson": *)[-+0-9.eE]+', "\\150.0", text)
  expect_false(identical(edited, text))
  writeLines(edited, json)
  back <- suppressWarnings(ferx_load_fit(fitrx_replace_entry(path, "fit.json", json)))
  expect_identical(wd_dw(back)[[1L]]$details$durbin_watson, 50)
})

test_that("the post-hoc steps keep the details in step with the warnings", {
  skip_on_cran()
  fit <- wd_fit()
  cv <- suppressWarnings(ferx_covariance(fit))
  expect_identical(wd_dw(cv), wd_dw(fit))
  sir <- suppressWarnings(ferx_sir(fit, sir_samples = 200L, sir_resamples = 50L,
                                   sir_seed = 1L))
  expect_identical(wd_dw(sir), wd_dw(fit))
})
