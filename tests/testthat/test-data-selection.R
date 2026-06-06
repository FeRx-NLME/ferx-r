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
# 9. print.ferx_data
# ---------------------------------------------------------------------------
test_that("print.ferx_data shows retention summary header", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
  out <- capture.output(print(sel))
  expect_true(any(grepl("ferx_data", out, fixed = TRUE)))
  expect_true(any(grepl("retained", out, fixed = TRUE)))
  expect_true(any(grepl("excluded", out, fixed = TRUE)))
  expect_true(any(grepl("Fired ignore", out, fixed = TRUE)))
})

test_that("print.ferx_data shows no excluded message when filter does not fire", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_selection(ex$data, ignore = "DV < 0.05")
  out <- capture.output(print(sel))
  expect_true(any(grepl("ferx_data", out, fixed = TRUE)))
  expect_false(any(grepl("ferx_selection_excluded", out, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 10. ferx_selection_excluded() helpers
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

# ---------------------------------------------------------------------------
# 11. accept= and ignore_ids= coverage for ferx_selection()
# ---------------------------------------------------------------------------
test_that("ferx_selection() accept= excludes records failing the condition", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  # keep only obs with DV > 5; MDV=1 / dose rows (DV=0) also fail this accept
  sel <- ferx_selection(ex$data, accept = "DV > 5")
  excl <- ferx_selection_excluded(sel)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  # retained rows all satisfy DV > 5 (for non-dose rows)
  fired <- .sel_attr(sel, "exclusions")$fired_accept
  expect_true(any(grepl("DV > 5", fired, fixed = TRUE)))  # accept: DV > 5
  # round-trip: retained + excluded == original
  expect_equal(nrow(sel) + nrow(excl), nrow(df))
})

test_that("ferx_selection() ignore_ids= excludes entire subject", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_selection(ex$data, ignore_ids = 1L)
  excl <- ferx_selection_excluded(sel)

  # subject 1 rows are all excluded
  expect_true(all(excl$ID == 1L))
  # subject 1 is in excluded_subject_ids
  expect_true("1" %in% .sel_attr(sel, "exclusions")$excluded_subject_ids)
  # retained data has no subject 1
  expect_false(1L %in% sel$ID)
  # total count preserved
  expect_equal(nrow(sel) + nrow(excl), nrow(df))
})

test_that("ferx_selection() fired_ignore format matches fit$exclusions format", {
  skip_on_cran()
  ex    <- ferx_example("warfarin_data_selection")
  sel   <- ferx_selection(ex$data, ignore = "DV < 1.0")
  fired <- .sel_attr(sel, "exclusions")$fired_ignore
  fit   <- warfarin_sel_fit()

  # Both R preview and Rust engine use "ignore: <expr>" format consistently.
  expect_true("ignore: DV < 1.0" %in% fired)
  expect_true("ignore: DV < 1.0" %in% fit$exclusions$fired_ignore)
  expect_equal(fired, fit$exclusions$fired_ignore)
})

# ---------------------------------------------------------------------------
# 12. ferx_selection() data.frame input
# ---------------------------------------------------------------------------
test_that("ferx_selection() accepts a data.frame and sets source_path to NULL", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_selection(df, ignore = "DV < 1.0")
  expect_s3_class(sel, "ferx_data")
  expect_null(.sel_attr(sel, "source_path"))
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, 2L)
})

test_that("ferx_fit errors clearly when ferx_data has no source_path", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_selection(df, ignore = "DV < 1.0")
  expect_error(
    ferx_fit(ex$model, sel),
    regexp = "no source file path"
  )
})

# ---------------------------------------------------------------------------
# 13. ferx_selection_excluded.ferx_fit replays fired conditions
# ---------------------------------------------------------------------------
test_that("ferx_selection_excluded.ferx_fit returns excluded rows with .exclude_reason", {
  skip_on_cran()
  fit  <- warfarin_sel_fit()
  excl <- ferx_selection_excluded(fit)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  expect_true(".exclude_reason" %in% names(excl))
  # row count matches fit$exclusions count
  expect_equal(nrow(excl), fit$exclusions$n_obs_excluded +
                 fit$exclusions$n_dose_excluded + fit$exclusions$n_other_excluded)
})

# ---------------------------------------------------------------------------
# 14. ferx_selection_excluded.ferx_fit NULL exclusions path
# ---------------------------------------------------------------------------
test_that("ferx_selection_excluded.ferx_fit messages when no data-selection rules", {
  skip_on_cran()
  fit <- warfarin_fit()   # plain fit; fit$exclusions is NULL
  expect_message(ferx_selection_excluded(fit), "No data-selection rules")
})

# ---------------------------------------------------------------------------
# 15. ferx_fit() type validation stops for ignore/accept/ignore_ids
# ---------------------------------------------------------------------------
test_that("ferx_fit errors when ignore is not character", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, ignore = 123L, verbose = FALSE),
    regexp = "character vector"
  )
})

test_that("ferx_fit errors when accept is not character", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, accept = 42.0, verbose = FALSE),
    regexp = "character vector"
  )
})

test_that("ferx_fit errors when ignore_ids is not numeric or character", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, ignore_ids = list("1"), verbose = FALSE),
    regexp = "numeric or character vector"
  )
})

# ---------------------------------------------------------------------------
# 16. ferx_fit() accept= and ignore_ids= settings_parts code paths
# ---------------------------------------------------------------------------
test_that("ferx_fit with accept= runs and populates fit$exclusions", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  # DV < 100 keeps all warfarin obs; filter is still active -> exclusions struct returned
  fit <- ferx_fit(ex$model, ex$data, accept = "DV < 100", verbose = FALSE,
                  settings = list(maxiter = 3L))
  # exclusions is non-NULL whenever the filter is active
  expect_false(is.null(fit$exclusions))
})

test_that("ferx_fit with ignore_ids= excludes a subject and populates excluded_subject_ids", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  fit <- ferx_fit(ex$model, ex$data, ignore_ids = 1L, verbose = FALSE,
                  settings = list(maxiter = 3L))
  expect_false(is.null(fit$exclusions))
  expect_true("1" %in% fit$exclusions$excluded_subject_ids)
})

# ---------------------------------------------------------------------------
# 17. print.ferx_fit fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------
test_that("print.ferx_fit shows fired_accept and excluded_subject_ids when present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  # Augment the exclusions in-place to exercise the display branches.
  fit$exclusions$fired_accept          <- c("accept: DV > -999")
  fit$exclusions$excluded_subject_ids  <- c("3", "5")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Fired accept", out, fixed = TRUE)))
  expect_true(any(grepl("Subjects excluded", out, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# 18. ferx_runlog fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------
test_that("ferx_runlog shows fired_accept and excluded_subject_ids when present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  fit$exclusions$fired_accept         <- c("accept: DV > -999")
  fit$exclusions$excluded_subject_ids <- c("3")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_true(grepl("Fired accept", out, fixed = TRUE))
  expect_true(grepl("Excl. subjs", out, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# 19. .cmp() operator coverage (pure R, via ferx_selection() with inline df)
# ---------------------------------------------------------------------------

# Tiny inline dataset: 4 obs records, predictable values.
.sel_test_df <- function() {
  data.frame(
    ID   = c(1L, 1L, 2L, 2L),
    TIME = c(0, 1, 0, 1),
    DV   = c(0.5, 2.0, 3.0, 5.5),
    EVID = c(0L, 0L, 0L, 0L),
    AMT  = c(0, 0, 0, 0),
    MDV  = c(0L, 0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
}

test_that("ferx_selection() == operator excludes matching row", {
  skip_on_cran()
  sel <- ferx_selection(.sel_test_df(), ignore = "DV == 2.0")
  expect_equal(nrow(ferx_selection_excluded(sel)), 1L)
  expect_equal(ferx_selection_excluded(sel)$DV, 2.0)
})

test_that("ferx_selection() != operator excludes non-matching rows", {
  skip_on_cran()
  sel <- ferx_selection(.sel_test_df(), ignore = "DV != 2.0")
  excl <- ferx_selection_excluded(sel)
  expect_true(nrow(excl) > 0L)
  expect_true(all(excl$DV != 2.0))
})

test_that("ferx_selection() <= operator excludes rows at or below threshold", {
  skip_on_cran()
  sel  <- ferx_selection(.sel_test_df(), ignore = "DV <= 0.5")
  excl <- ferx_selection_excluded(sel)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$DV, 0.5)
})

test_that("ferx_selection() >= operator excludes rows at or above threshold", {
  skip_on_cran()
  sel  <- ferx_selection(.sel_test_df(), ignore = "DV >= 5.5")
  excl <- ferx_selection_excluded(sel)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$DV, 5.5)
})

# ---------------------------------------------------------------------------
# 20. .parse_filter_expr() coverage: character rhs, &&, !found
# ---------------------------------------------------------------------------

test_that("ferx_selection() handles character rhs (quoted string)", {
  skip_on_cran()
  df <- data.frame(
    ID = c("A", "B", "A"), DV = c(1, 2, 3),
    EVID = 0L, AMT = 0, MDV = 0L,
    stringsAsFactors = FALSE
  )
  sel  <- ferx_selection(df, ignore = "ID == 'B'")
  excl <- ferx_selection_excluded(sel)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$ID, "B")
})

test_that("ferx_selection() handles && (AND) multi-part expression", {
  skip_on_cran()
  # Only the row where DV >= 3 AND ID == 2 should be excluded.
  sel  <- ferx_selection(.sel_test_df(), ignore = "DV >= 3.0 && ID == 2")
  excl <- ferx_selection_excluded(sel)
  expect_true(nrow(excl) > 0L)
  expect_true(all(excl$ID == 2L))
  expect_true(all(excl$DV >= 3.0))
})

test_that("ferx_selection() ignores unparseable expression (returns nothing excluded)", {
  skip_on_cran()
  # No operator -> .parse_filter_expr returns NULL -> .eval_expr returns FALSE -> no exclusion
  sel <- ferx_selection(.sel_test_df(), ignore = "UNPARSEABLE_EXPR")
  expect_equal(nrow(ferx_selection_excluded(sel)), 0L)
  expect_equal(length(.sel_attr(sel, "exclusions")$fired_ignore), 0L)
})
