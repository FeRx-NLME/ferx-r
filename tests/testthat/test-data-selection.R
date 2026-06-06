# Tests for [data_selection] / ferx_selection() / ferx_fit(ignore=) support.

# Local alias avoids the ::: operator.
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ferx_data objects are data.frame subclasses; attributes are NOT accessible
# via `$` (which looks at columns). Use attr() for metadata.
.sel_attr <- function(x, name) attr(x, name, exact = TRUE)

# Cached fit for warfarin_data_selection (FD gradient, no Enzyme required).
warfarin_sel_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_data_selection")
      fit <<- ferx_fit(ex$model, ex$data, verbose = FALSE,
                       settings = list(maxiter = 30L))
    }
    fit
  }
})

# ---------------------------------------------------------------------------
# 1. Plain warfarin fit has NULL exclusions
# ---------------------------------------------------------------------------
test_that("plain warfarin fit has NULL exclusions", {
  skip_on_cran()
  fit <- warfarin_fit()
  expect_null(fit$exclusions)
})

# ---------------------------------------------------------------------------
# 2. [data_selection] model populates exclusions
# ---------------------------------------------------------------------------
test_that("data_selection fit has positive n_obs_excluded", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  ex  <- fit$exclusions
  expect_false(is.null(ex))
  expect_true(ex$n_obs_excluded > 0L)
})

test_that("data_selection exclusions total is consistent", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  ex  <- fit$exclusions
  expect_equal(
    ex$n_records_total,
    ex$n_obs_excluded + ex$n_dose_excluded + ex$n_other_excluded +
      fit$n_obs + fit$n_subjects,
    tolerance = 0L
  )
})

test_that("exclusions$fired_ignore is non-empty for data_selection model", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  expect_true(length(fit$exclusions$fired_ignore) > 0L)
})

# ---------------------------------------------------------------------------
# 3. ferx_fit(data, ignore = ...) produces exclusions matching model-file path
# ---------------------------------------------------------------------------
test_that("ferx_fit(ignore=) produces same exclusion count as [data_selection] model", {
  skip_on_cran()
  ex     <- ferx_example("warfarin")         # base model, no data_selection
  ex_sel <- ferx_example("warfarin_data_selection")  # has [data_selection] ignore = DV < 1.0

  # R-arg path: base model + explicit ignore arg
  fit_arg <- ferx_fit(
    ex$model, ex$data,
    ignore  = "DV < 1.0",
    verbose = FALSE,
    settings = list(maxiter = 5L)
  )
  # Model-file path: data_selection model, same dataset
  fit_mod <- ferx_fit(
    ex_sel$model, ex_sel$data,
    verbose = FALSE,
    settings = list(maxiter = 5L)
  )

  expect_equal(
    fit_arg$exclusions$n_obs_excluded,
    fit_mod$exclusions$n_obs_excluded
  )
})

# ---------------------------------------------------------------------------
# 4. ferx_selection() pure-R preview
# ---------------------------------------------------------------------------
test_that("ferx_selection() returns ferx_data S3 class", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
  expect_s3_class(sel, "ferx_data")
  expect_s3_class(sel, "data.frame")
})

test_that("ferx_selection() preview n_obs_excluded matches fit$exclusions", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
  fit <- warfarin_sel_fit()

  expect_equal(
    .sel_attr(sel, "exclusions")$n_obs_excluded,
    fit$exclusions$n_obs_excluded
  )
  expect_equal(
    .sel_attr(sel, "exclusions")$n_records_total,
    fit$exclusions$n_records_total
  )
})

test_that("ferx_selection() retained + excluded rows total original nrow", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")

  expect_equal(nrow(sel) + nrow(.sel_attr(sel, "excluded_rows")), nrow(df))
})

test_that("ferx_selection() n_obs_excluded matches manual count", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")

  expected_excl <- sum(df$DV < 1.0 & df$EVID == 0L & (is.na(df$MDV) | df$MDV == 0L),
                       na.rm = TRUE)
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, expected_excl)
})

test_that("ferx_selection() source_path is set when path supplied", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
  expect_equal(.sel_attr(sel, "source_path"), ex$data)
})

test_that("ferx_selection() fired_ignore contains the matched rule", {
  skip_on_cran()
  ex    <- ferx_example("warfarin")
  sel   <- ferx_selection(ex$data, ignore = "DV < 1.0")
  fired <- .sel_attr(sel, "exclusions")$fired_ignore
  expect_true(any(grepl("DV < 1.0", fired, fixed = TRUE)))
})

test_that("ferx_selection() returns 0 excluded for condition that never fires", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_selection(ex$data, ignore = "DV < 0.05")
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, 0L)
  expect_equal(length(.sel_attr(sel, "exclusions")$fired_ignore), 0L)
})

# ---------------------------------------------------------------------------
# 5. ferx_data passed directly to ferx_fit()
# ---------------------------------------------------------------------------
test_that("ferx_fit accepts ferx_data as data argument", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
  fit <- ferx_fit(ex$model, sel, verbose = FALSE, settings = list(maxiter = 5L))
  expect_false(is.null(fit$exclusions))
  expect_true(fit$exclusions$n_obs_excluded > 0L)
})

# ---------------------------------------------------------------------------
# 6. print.ferx_fit() shows DATA SELECTION block
# ---------------------------------------------------------------------------
test_that("print.ferx_fit shows DATA SELECTION when exclusions present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("DATA SELECTION", out, fixed = TRUE)))
  expect_true(any(grepl("Obs excl", out, fixed = TRUE)))
})

test_that("print.ferx_fit does NOT show DATA SELECTION for plain fit", {
  skip_on_cran()
  fit <- warfarin_fit()
  out <- capture.output(print(fit))
  expect_false(any(grepl("DATA SELECTION", out, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 7. ferx_runlog() shows exclusion line
# ---------------------------------------------------------------------------
test_that("ferx_runlog shows exclusion line when exclusions present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_true(any(grepl("Data selection", out, fixed = TRUE)))
})

test_that("ferx_runlog does NOT show exclusion line for plain fit", {
  skip_on_cran()
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_false(any(grepl("Data selection", out, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 8. persist round-trip preserves exclusions
# ---------------------------------------------------------------------------
test_that("ferx_save/ferx_load round-trip preserves exclusions", {
  skip_on_cran()
  fit  <- warfarin_sel_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$exclusions))
  expect_equal(loaded$exclusions$n_obs_excluded,   fit$exclusions$n_obs_excluded)
  expect_equal(loaded$exclusions$n_records_total,  fit$exclusions$n_records_total)
  expect_equal(loaded$exclusions$fired_ignore,     fit$exclusions$fired_ignore)
  expect_equal(loaded$exclusions$excluded_subject_ids, fit$exclusions$excluded_subject_ids)
})

test_that("ferx_save/ferx_load preserves NULL exclusions", {
  skip_on_cran()
  fit  <- warfarin_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_null(loaded$exclusions)
})

# ---------------------------------------------------------------------------
# 9. ferx_selection_excluded() helpers
# ---------------------------------------------------------------------------
test_that("ferx_selection_excluded.ferx_data returns a data.frame with .exclude_reason", {
  skip_on_cran()
  ex   <- ferx_example("warfarin")
  sel  <- ferx_selection(ex$data, ignore = "DV < 1.0")
  excl <- ferx_selection_excluded(sel)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  expect_true(".exclude_reason" %in% names(excl))
})
