
# ---- header from test-data-selection.R ----
# Tests for [data_selection] / ferx_apply_selection() / ferx_fit(ignore=) support.

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

# ---------------------------------------------------------------------------
# 2. [data_selection] model populates exclusions
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 3. ferx_fit(data, ignore = ...) produces exclusions matching model-file path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4. ferx_apply_selection() pure-R preview
# ---------------------------------------------------------------------------







# ---------------------------------------------------------------------------
# 5. ferx_data passed directly to ferx_fit()
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6. print.ferx_fit() shows DATA SELECTION block
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 7. ferx_runlog() shows exclusion line
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 8. persist round-trip preserves exclusions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 9. print.ferx_data
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 10. ferx_apply_selection(excluded = TRUE) on a ferx_data
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 11. accept= and ignore_ids= coverage for ferx_apply_selection()
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 12. ferx_apply_selection() data.frame input
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 13. ferx_apply_selection(fit, excluded = TRUE) replays fired conditions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 14. ferx_apply_selection(fit, excluded = TRUE) NULL exclusions path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 15. ferx_fit() type validation stops for ignore/accept/ignore_ids
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 16. ferx_fit() accept= and ignore_ids= settings_parts code paths
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 17. print.ferx_fit fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 18. ferx_runlog fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 19. .cmp() operator coverage (pure R, via ferx_apply_selection() with inline df)
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





# ---------------------------------------------------------------------------
# 20. .parse_filter_expr() coverage: character rhs, &&, !found
# ---------------------------------------------------------------------------





test_that("plain warfarin fit has NULL exclusions", {
  skip_on_cran()
  fit <- warfarin_fit()
  expect_null(fit$exclusions)
})
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
test_that("ferx_apply_selection() returns ferx_data S3 class", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  expect_s3_class(sel, "ferx_data")
  expect_s3_class(sel, "data.frame")
})
test_that("ferx_apply_selection() preview n_obs_excluded matches fit$exclusions", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
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
test_that("ferx_apply_selection() retained + excluded rows total original nrow", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")

  expect_equal(nrow(sel) + nrow(.sel_attr(sel, "excluded_rows")), nrow(df))
})
test_that("ferx_apply_selection() n_obs_excluded matches manual count", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")

  expected_excl <- sum(df$DV < 1.0 & df$EVID == 0L & (is.na(df$MDV) | df$MDV == 0L),
                       na.rm = TRUE)
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, expected_excl)
})
test_that("ferx_apply_selection() source_path is set when path supplied", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  expect_equal(.sel_attr(sel, "source_path"), ex$data)
})
test_that("ferx_apply_selection() fired_ignore contains the matched rule", {
  skip_on_cran()
  ex    <- ferx_example("warfarin")
  sel   <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  fired <- .sel_attr(sel, "exclusions")$fired_ignore
  expect_true(any(grepl("DV < 1.0", fired, fixed = TRUE)))
})
test_that("ferx_apply_selection() returns 0 excluded for condition that never fires", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 0.05")
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, 0L)
  expect_equal(length(.sel_attr(sel, "exclusions")$fired_ignore), 0L)
})
test_that("ferx_fit accepts ferx_data as data argument", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  fit <- ferx_fit(ex$model, sel, verbose = FALSE, settings = list(maxiter = 5L))
  expect_false(is.null(fit$exclusions))
  expect_true(fit$exclusions$n_obs_excluded > 0L)
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
test_that("print.ferx_data shows retention summary header", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  out <- capture.output(print(sel))
  expect_true(any(grepl("ferx_data", out, fixed = TRUE)))
  expect_true(any(grepl("retained", out, fixed = TRUE)))
  expect_true(any(grepl("excluded", out, fixed = TRUE)))
  expect_true(any(grepl("Fired ignore", out, fixed = TRUE)))
})
test_that("print.ferx_data shows no excluded message when filter does not fire", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  sel <- ferx_apply_selection(ex$data, ignore = "DV < 0.05")
  out <- capture.output(print(sel))
  expect_true(any(grepl("ferx_data", out, fixed = TRUE)))
  expect_false(any(grepl("excluded = TRUE", out, fixed = TRUE)))
})
test_that("ferx_apply_selection(sel, excluded = TRUE) returns a data.frame with .exclude_reason", {
  skip_on_cran()
  ex   <- ferx_example("warfarin")
  sel  <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  expect_true(".exclude_reason" %in% names(excl))
})
test_that("ferx_apply_selection() accept= excludes records failing the condition", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  # keep only obs with DV > 5; MDV=1 / dose rows (DV=0) also fail this accept
  sel <- ferx_apply_selection(ex$data, accept = "DV > 5")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  # retained rows all satisfy DV > 5 (for non-dose rows)
  fired <- .sel_attr(sel, "exclusions")$fired_accept
  expect_true(any(grepl("DV > 5", fired, fixed = TRUE)))  # accept: DV > 5
  # round-trip: retained + excluded == original
  expect_equal(nrow(sel) + nrow(excl), nrow(df))
})
test_that("ferx_apply_selection() ignore_ids= excludes entire subject", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_apply_selection(ex$data, ignore_ids = 1L)
  excl <- ferx_apply_selection(sel, excluded = TRUE)

  # subject 1 rows are all excluded
  expect_true(all(excl$ID == 1L))
  # subject 1 is in excluded_subject_ids
  expect_true("1" %in% .sel_attr(sel, "exclusions")$excluded_subject_ids)
  # retained data has no subject 1
  expect_false(1L %in% sel$ID)
  # total count preserved
  expect_equal(nrow(sel) + nrow(excl), nrow(df))
})
test_that("ferx_apply_selection() fired_ignore format matches fit$exclusions format", {
  skip_on_cran()
  ex    <- ferx_example("warfarin_data_selection")
  sel   <- ferx_apply_selection(ex$data, ignore = "DV < 1.0")
  fired <- .sel_attr(sel, "exclusions")$fired_ignore
  fit   <- warfarin_sel_fit()

  # Both R preview and Rust engine use "ignore: <expr>" format consistently.
  expect_true("ignore: DV < 1.0" %in% fired)
  expect_true("ignore: DV < 1.0" %in% fit$exclusions$fired_ignore)
  expect_equal(fired, fit$exclusions$fired_ignore)
})
test_that("ferx_apply_selection() accepts a data.frame and sets source_path to NULL", {
  skip_on_cran()
  ex  <- ferx_example("warfarin")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_apply_selection(df, ignore = "DV < 1.0")
  expect_s3_class(sel, "ferx_data")
  expect_null(.sel_attr(sel, "source_path"))
  expect_equal(.sel_attr(sel, "exclusions")$n_obs_excluded, 2L)
})
test_that("ferx_fit errors clearly when ferx_data has no source_path", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_data_selection")
  df  <- utils::read.csv(ex$data, stringsAsFactors = FALSE, check.names = FALSE)
  sel <- ferx_apply_selection(df, ignore = "DV < 1.0")
  expect_error(
    ferx_fit(ex$model, sel),
    regexp = "no source file path"
  )
})
test_that("ferx_apply_selection(fit, excluded = TRUE) returns excluded rows with .exclude_reason", {
  skip_on_cran()
  fit  <- warfarin_sel_fit()
  excl <- ferx_apply_selection(fit, excluded = TRUE)
  expect_s3_class(excl, "data.frame")
  expect_true(nrow(excl) > 0L)
  expect_true(".exclude_reason" %in% names(excl))
  # row count matches fit$exclusions count
  expect_equal(nrow(excl), fit$exclusions$n_obs_excluded +
                 fit$exclusions$n_dose_excluded + fit$exclusions$n_other_excluded)
})
test_that("ferx_apply_selection(fit, excluded = TRUE) messages when no data-selection rules", {
  skip_on_cran()
  fit <- warfarin_fit()   # plain fit; fit$exclusions is NULL
  expect_message(ferx_apply_selection(fit, excluded = TRUE), "No data-selection rules")
})
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
test_that("ferx_apply_selection() == operator excludes matching row", {
  skip_on_cran()
  sel <- ferx_apply_selection(.sel_test_df(), ignore = "DV == 2.0")
  expect_equal(nrow(ferx_apply_selection(sel, excluded = TRUE)), 1L)
  expect_equal(ferx_apply_selection(sel, excluded = TRUE)$DV, 2.0)
})
test_that("ferx_apply_selection() != operator excludes non-matching rows", {
  skip_on_cran()
  sel <- ferx_apply_selection(.sel_test_df(), ignore = "DV != 2.0")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_true(nrow(excl) > 0L)
  expect_true(all(excl$DV != 2.0))
})
test_that("ferx_apply_selection() <= operator excludes rows at or below threshold", {
  skip_on_cran()
  sel  <- ferx_apply_selection(.sel_test_df(), ignore = "DV <= 0.5")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$DV, 0.5)
})
test_that("ferx_apply_selection() >= operator excludes rows at or above threshold", {
  skip_on_cran()
  sel  <- ferx_apply_selection(.sel_test_df(), ignore = "DV >= 5.5")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$DV, 5.5)
})
test_that("ferx_apply_selection() handles character rhs (quoted string)", {
  skip_on_cran()
  df <- data.frame(
    ID = c("A", "B", "A"), DV = c(1, 2, 3),
    EVID = 0L, AMT = 0, MDV = 0L,
    stringsAsFactors = FALSE
  )
  sel  <- ferx_apply_selection(df, ignore = "ID == 'B'")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_equal(nrow(excl), 1L)
  expect_equal(excl$ID, "B")
})
test_that("ferx_apply_selection() handles && (AND) multi-part expression", {
  skip_on_cran()
  # Only the row where DV >= 3 AND ID == 2 should be excluded.
  sel  <- ferx_apply_selection(.sel_test_df(), ignore = "DV >= 3.0 && ID == 2")
  excl <- ferx_apply_selection(sel, excluded = TRUE)
  expect_true(nrow(excl) > 0L)
  expect_true(all(excl$ID == 2L))
  expect_true(all(excl$DV >= 3.0))
})
test_that("ferx_apply_selection() ignores unparseable expression (returns nothing excluded)", {
  skip_on_cran()
  # No operator AND not a lone identifier (a bare token like `C` is now the valid
  # IGNORE=C shorthand) -> .parse_filter_expr returns NULL -> .eval_expr returns
  # FALSE -> no exclusion.
  sel <- ferx_apply_selection(.sel_test_df(), ignore = "UNPARSEABLE EXPR HERE")
  expect_equal(nrow(ferx_apply_selection(sel, excluded = TRUE)), 0L)
  expect_equal(length(.sel_attr(sel, "exclusions")$fired_ignore), 0L)
})

# ---- header from test-selection-helpers.R ----
# Unit tests for the pure-R data-selection helpers in selection.R:
# .parse_filter_expr(), .cmp() and .get_col_sel(). No model fit required.

# ---------------------------------------------------------------------------
# .parse_filter_expr
# ---------------------------------------------------------------------------








# ---------------------------------------------------------------------------
# .cmp
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# .get_col_sel
# ---------------------------------------------------------------------------


test_that(".parse_filter_expr parses a single numeric comparison and lower-cases the column", {
  res <- ferx:::.parse_filter_expr("AMT == 0")
  expect_length(res, 1L)
  expect_identical(res[[1L]]$col, "amt")
  expect_identical(res[[1L]]$op, "==")
  expect_identical(res[[1L]]$rhs, 0)            # numeric coercion
})
test_that(".parse_filter_expr strips quotes and keeps the value as a string", {
  res <- ferx:::.parse_filter_expr("MDV != '.'")
  expect_identical(res[[1L]]$op, "!=")
  expect_identical(res[[1L]]$rhs, ".")          # quotes removed, stays character
})
test_that(".parse_filter_expr handles && chains and prefers >= over >", {
  res <- ferx:::.parse_filter_expr("TIME >= 24 && CMT > 1")
  expect_length(res, 2L)
  expect_identical(res[[1L]]$op, ">=")          # two-char operator matched first
  expect_identical(res[[1L]]$rhs, 24)
  expect_identical(res[[2L]]$op, ">")
  expect_identical(res[[2L]]$col, "cmt")
})
test_that(".parse_filter_expr returns NULL on malformed clauses", {
  expect_null(ferx:::.parse_filter_expr("no operator here"))
  expect_null(ferx:::.parse_filter_expr("AMT == "))   # empty rhs
  expect_null(ferx:::.parse_filter_expr("== 5"))       # empty lhs
})
test_that(".parse_filter_expr expands a bare identifier to COL == COL (IGNORE=C)", {
  # NONMEM `IGNORE=C` shorthand: a lone column name drops rows whose column
  # holds the literal label. Mirrors ferx-core; without it the preview would
  # report zero exclusions while the Rust fit drops the flagged rows.
  res <- ferx:::.parse_filter_expr("C")
  expect_length(res, 1L)
  expect_identical(res[[1L]]$col, "c")     # column lowercased for lookup
  expect_identical(res[[1L]]$op, "==")
  expect_identical(res[[1L]]$rhs, "C")     # RHS case preserved to match raw cell
})
test_that(".parse_filter_expr treats non-finite numeric spellings as label strings", {
  # `NaN` / `Inf` parse as doubles but could never match numerically, so they
  # are label strings (matching ferx-core's filter_expr parser).
  res_nan <- ferx:::.parse_filter_expr("FLAG == NaN")
  expect_identical(res_nan[[1L]]$rhs, "NaN")
  res_inf <- ferx:::.parse_filter_expr("FLAG == Inf")
  expect_identical(res_inf[[1L]]$rhs, "Inf")
  # A finite literal is still numeric.
  res_num <- ferx:::.parse_filter_expr("FLAG == 70")
  expect_identical(res_num[[1L]]$rhs, 70)
})
test_that(".parse_filter_expr rejects an ordered comparison against a non-numeric RHS", {
  # ferx-core errors on `BW < abc` (ordered op needs a numeric value); the lenient
  # R preview mirrors that by returning NULL (no exclusion) instead of a silent
  # lexical string comparison that the Rust fit never performs.
  expect_null(ferx:::.parse_filter_expr("BW < abc"))
  expect_null(ferx:::.parse_filter_expr("BW >= abc"))
  expect_null(ferx:::.parse_filter_expr("BW < 'abc'"))   # quoted is no different
  # A whole && clause is dropped if any ordered sub-expression is non-numeric.
  expect_null(ferx:::.parse_filter_expr("EVID == 0 && BW < abc"))
  # Ordered comparison against a numeric RHS is still parsed.
  res <- ferx:::.parse_filter_expr("BW < 48")
  expect_identical(res[[1L]]$op, "<")
  expect_identical(res[[1L]]$rhs, 48)
})
test_that(".parse_filter_expr expands a bare identifier regardless of ignore/accept context", {
  # The parser is context-free: ferx-core routes both ignore and accept clauses
  # through the same FilterClause::parse, so a bare token expands to COL == COL in
  # either. ferx_apply_selection() applies it on both paths; this locks that intent.
  res <- ferx:::.parse_filter_expr("STUDY")
  expect_identical(res[[1L]]$col, "study")
  expect_identical(res[[1L]]$op, "==")
  expect_identical(res[[1L]]$rhs, "STUDY")
})
test_that(".cmp compares numeric columns and treats NA lhs as FALSE", {
  expect_identical(ferx:::.cmp(c("1", "2", "3"), "==", 2), c(FALSE, TRUE, FALSE))
  expect_identical(ferx:::.cmp(c("1", NA, "3"), ">", 1), c(FALSE, FALSE, TRUE))
  expect_identical(ferx:::.cmp(c("1", "2", "3"), "<=", 2), c(TRUE, TRUE, FALSE))
})
test_that(".cmp compares character columns when the rhs is a string", {
  expect_identical(ferx:::.cmp(c("a", "b"), "==", "b"), c(FALSE, TRUE))
  expect_identical(ferx:::.cmp(c("a", "b"), "!=", "b"), c(TRUE, FALSE))
})
test_that(".cmp returns all FALSE for NA / NaN rhs", {
  expect_identical(ferx:::.cmp(c("1", "2"), "==", NA), c(FALSE, FALSE))
  expect_identical(ferx:::.cmp(c("1", "2"), "==", NaN), c(FALSE, FALSE))
})
test_that(".get_col_sel matches column names case-insensitively, else NULL", {
  df <- data.frame(AMT = c(10, 0), TIME = c(0, 1))
  expect_identical(ferx:::.get_col_sel(df, "amt"), c(10, 0))
  expect_identical(ferx:::.get_col_sel(df, "time"), c(0, 1))
  expect_null(ferx:::.get_col_sel(df, "dv"))
})
