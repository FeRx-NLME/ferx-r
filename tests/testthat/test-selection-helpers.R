# Unit tests for the pure-R data-selection helpers in selection.R:
# .parse_filter_expr(), .cmp() and .get_col_sel(). No model fit required.

# ---------------------------------------------------------------------------
# .parse_filter_expr
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
  # either. ferx_selection() applies it on both paths; this locks that intent.
  res <- ferx:::.parse_filter_expr("STUDY")
  expect_identical(res[[1L]]$col, "study")
  expect_identical(res[[1L]]$op, "==")
  expect_identical(res[[1L]]$rhs, "STUDY")
})

# ---------------------------------------------------------------------------
# .cmp
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# .get_col_sel
# ---------------------------------------------------------------------------
test_that(".get_col_sel matches column names case-insensitively, else NULL", {
  df <- data.frame(AMT = c(10, 0), TIME = c(0, 1))
  expect_identical(ferx:::.get_col_sel(df, "amt"), c(10, 0))
  expect_identical(ferx:::.get_col_sel(df, "time"), c(0, 1))
  expect_null(ferx:::.get_col_sel(df, "dv"))
})
