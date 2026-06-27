# Tests for the model-file vs. call-time [fit_options] conflict detection
# (#31 / PR #66). Covers the parser ferx:::.ferx_parse_model_fit_options() and the
# warning helper ferx:::.ferx_warn_fit_option_conflicts() in isolation, plus a small
# end-to-end check via ferx_fit() against the bundled warfarin example.

# ---------------------------------------------------------------------------
# .ferx_parse_model_fit_options
# ---------------------------------------------------------------------------

write_fit_options_model <- function(lines) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)",
               "[fit_options]", lines), path)
  path
}

test_that("ferx:::.ferx_parse_model_fit_options() returns named char vec of raw values", {
  path <- write_fit_options_model(c("  method = focei", "  maxiter = 300"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "focei", maxiter = "300"))
})

test_that("ferx:::.ferx_parse_model_fit_options() lower-cases keys, preserves value case", {
  path <- write_fit_options_model(c("  Method = FOCE", "  Covariance = TRUE"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(names(opts), c("method", "covariance"))
  expect_equal(unname(opts), c("FOCE", "TRUE"))
})

test_that("ferx:::.ferx_parse_model_fit_options() strips comments via .ferx_extract_blocks()", {
  path <- write_fit_options_model(c(
    "  method = foce  # the FOCE method",
    "  # a full-line comment",
    "  maxiter = 300 // also a comment"
  ))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "foce", maxiter = "300"))
})

test_that("ferx:::.ferx_parse_model_fit_options() preserves embedded `=` in values", {
  path <- write_fit_options_model("  custom = a=b=c")
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(unname(opts["custom"]), "a=b=c")
})

test_that("ferx:::.ferx_parse_model_fit_options() returns empty char vec when block absent", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)"), path)
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_length(opts, 0L)
  expect_type(opts, "character")
})

test_that("ferx:::.ferx_parse_model_fit_options() returns empty char vec when block is empty", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)",
               "[fit_options]"), path)
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_length(opts, 0L)
})

test_that("ferx:::.ferx_parse_model_fit_options() ignores malformed lines (no `=`)", {
  path <- write_fit_options_model(c("  method = foce", "  garbage_no_equals"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "foce"))
})

# ---------------------------------------------------------------------------
# .ferx_warn_fit_option_conflicts
# ---------------------------------------------------------------------------

empty_settings <- function() list(keys = character(0), values = character(0))

test_that("no warning when model_file_opts is empty", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = setNames(character(0), character(0)),
      dedicated_explicit = list(method = "focei"),
      settings_parts     = empty_settings()
    )
  )
})

test_that("no warning when explicit dedicated arg matches model file", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "foce"),
      dedicated_explicit = list(method = "foce"),
      settings_parts     = empty_settings()
    )
  )
})

test_that("warning when explicit dedicated arg disagrees with model file", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "foce"),
      dedicated_explicit = list(method = "focei"),
      settings_parts     = empty_settings()
    ),
    "method = foce.*overrides it with `focei`"
  )
})

test_that("no warning when dedicated arg is NULL (default accepted, not explicit)", {
  # This is the regression test for the spurious-warning bug: a user who does
  # not pass `method` should not be warned that they "overrode" the model file.
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "focei"),
      dedicated_explicit = list(method = NULL),  # user did not pass `method`
      settings_parts     = empty_settings()
    )
  )
})

test_that("comparison is case-insensitive (TRUE vs true does not warn)", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(covariance = "TRUE"),
      dedicated_explicit = list(covariance = "true"),
      settings_parts     = empty_settings()
    )
  )
})

test_that("settings list key conflict warns", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(maxiter = "300"),
      dedicated_explicit = list(),
      settings_parts     = list(keys = "maxiter", values = "500")
    ),
    "maxiter = 300.*overrides it with `500`"
  )
})

test_that("settings list key matching model file does not warn", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(maxiter = "300"),
      dedicated_explicit = list(),
      settings_parts     = list(keys = "maxiter", values = "300")
    )
  )
})

test_that("aliased model-file keys are matched (bloq → bloq_method)", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(bloq = "drop"),
      dedicated_explicit = list(bloq_method = "m3"),
      settings_parts     = empty_settings()
    ),
    "bloq = drop.*overrides it with `m3`"
  )
})

test_that("aliased model-file keys are matched (gradient_method → gradient)", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(gradient_method = "ad"),
      dedicated_explicit = list(gradient = "fd"),
      settings_parts     = empty_settings()
    ),
    "gradient_method = ad.*overrides it with `fd`"
  )
})

# ---------------------------------------------------------------------------
# End-to-end: ferx_fit() integration with match.call()-based detection
# ---------------------------------------------------------------------------

# Helper: copy warfarin example to a temp dir and rewrite [fit_options] so the
# fit runs in just a few iterations. The model file's `method = foce` and
# `covariance = true` lines are preserved so we can test conflict detection
# against them without paying for a full estimation.
make_fast_warfarin <- function() {
  ex <- ferx_example("warfarin")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  model_path <- file.path(tmpdir, "warfarin.ferx")
  data_path  <- file.path(tmpdir, "warfarin.csv")
  file.copy(ex$data, data_path)
  src <- readLines(ex$model)
  src <- sub("^(\\s*maxiter\\s*=\\s*)\\d+", "\\15", src)
  writeLines(src, model_path)
  list(model = model_path, data = data_path, dir = tmpdir)
}

test_that("ferx_fit() does not warn when defaults are accepted (regression for #66 review)", {
  # The model file specifies method = foce, covariance = true. Calling
  # ferx_fit() without overriding either should NOT warn — match.call() must
  # distinguish "user accepted the default" from "user explicitly passed the
  # default".
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  # Use expect_no_warning() rather than expect_silent() — ferx_fit() may print
  # informational messages (e.g. "Mu-referencing detected for..."), which are
  # unrelated to the conflict-detection logic under test.
  expect_no_warning(
    fit <- suppressMessages(ferx_fit(
      fx$model, fx$data,
      verbose = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ))
  )
  expect_s3_class(fit, "ferx_fit")
  expect_equal(fit$model_file_settings$method, "foce")
})

test_that("ferx_fit() warns when an explicit dedicated arg overrides the model file", {
  fx <- make_fast_warfarin()  # model file: method = foce
  on.exit(unlink(fx$dir, recursive = TRUE))

  expect_warning(
    ferx_fit(
      fx$model, fx$data,
      method = "focei",
      verbose = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ),
    "method = foce.*overrides it with `focei`"
  )
})

test_that("ferx_fit() runs the model file's method when none is passed (#558)", {
  # The model file sets `method = foce`. Before #558 the R-side default
  # method = "focei" silently overrode it; now passing no method keeps FOCE.
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressMessages(ferx_fit(
    fx$model, fx$data,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  ))
  expect_equal(fit$method, "FOCE")
})

test_that("ferx_fit() method argument overrides the model file's method (#558)", {
  fx <- make_fast_warfarin()  # model file: method = foce
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressWarnings(suppressMessages(ferx_fit(
    fx$model, fx$data,
    method = "focei",
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  )))
  expect_equal(fit$method, "FOCEI")
})

# Helper: like make_fast_warfarin(), but force `covariance = false` in the model
# file so we can test that the R-side `covariance` default no longer overrides it.
make_fast_warfarin_no_cov <- function() {
  fx <- make_fast_warfarin()
  src <- readLines(fx$model)
  src <- src[!grepl("^\\s*covariance\\s*=", src)]
  src <- sub("(\\[fit_options\\])", "\\1\n  covariance = false", src)
  writeLines(src, fx$model)
  fx
}

test_that("ferx_fit() honours the model file's `covariance = false` when none is passed (#558)", {
  # Before extending #558, the R default `covariance = TRUE` ran the covariance
  # step even though the model file disabled it. Passing no `covariance` must now
  # keep the model file's setting.
  fx <- make_fast_warfarin_no_cov()
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressMessages(ferx_fit(
    fx$model, fx$data,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  ))
  expect_equal(fit$covariance_status, "not_requested")
})

test_that("ferx_fit() covariance argument overrides the model file (#558)", {
  fx <- make_fast_warfarin_no_cov()  # model file: covariance = false
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressWarnings(suppressMessages(ferx_fit(
    fx$model, fx$data,
    covariance = TRUE,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  )))
  expect_true(fit$covariance_status %in% c("computed", "failed", "sir_fallback"))
})
