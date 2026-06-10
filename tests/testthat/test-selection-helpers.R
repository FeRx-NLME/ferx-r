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
