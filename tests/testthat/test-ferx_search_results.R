# Reading a run's `candidates.csv` back. The file is the engine's format, so
# the fixtures below are written exactly as `search::output` writes it: an
# empty cell wherever there is no value.

candidates_csv <- function(dir, file = "candidates.csv") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, file)
  writeLines(c(
    paste(ferx:::ferx_rust_search_table_columns(), collapse = ","),
    # A winner: fitted, converged, passed the gate.
    'c1,,ab12,"COVARIATE(CL,WT,pow)",412.500000,400.250000,true,true,,,12.500000,,,,false',
    # Excluded by strictness: fitted, but the gate says no.
    'c2,c1,cd34,"COVARIATE(V1,WT,pow)",418.000000,405.000000,false,false,not converged,condition number,9.000000,,,,false',
    # Never fitted: no criterion, no OFV, an error that a resume would retry.
    'c3,c1,ef56,"COVARIATE(Q,WT,pow)",,,,false,,,0.100000,ODE solver failed,true,,false',
    # A duplicate reused from the journal.
    'c4,c1,ab12,"COVARIATE(CL,WT,pow)",412.500000,400.250000,true,true,,,0.000000,,,c1,true'
  ), path)
  path
}

test_that("a run directory reads back with typed columns", {
  dir <- tempfile()
  candidates_csv(dir)
  res <- ferx_search_results(dir)

  expect_s3_class(res, "data.frame")
  expect_equal(names(res), ferx:::ferx_rust_search_table_columns())
  expect_equal(nrow(res), 4L)

  expect_type(res$criterion, "double")
  expect_type(res$ofv, "double")
  expect_type(res$seconds, "double")
  expect_type(res$converged, "logical")
  expect_type(res$passed, "logical")
  expect_type(res$reused, "logical")
  expect_type(res$retryable, "logical")

  expect_equal(res$id, c("c1", "c2", "c3", "c4"))
  expect_equal(res$criterion[1], 412.5)
  expect_equal(res$converged, c(TRUE, FALSE, NA, TRUE))
  expect_equal(res$passed, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(res$reused, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("the engine's empty cells become NA, never NaN or an empty string", {
  dir <- tempfile()
  candidates_csv(dir)
  res <- ferx_search_results(dir)

  # A candidate that never fitted has no criterion and no OFV.
  expect_true(is.na(res$criterion[3]))
  expect_true(is.na(res$ofv[3]))
  expect_false(any(is.nan(res$criterion)))

  # `parent` is empty for the root, `retryable` empty for a row that did not
  # fail, `duplicate_of` empty for a row that is not a duplicate.
  expect_true(is.na(res$parent[1]))
  expect_equal(res$retryable, c(NA, NA, TRUE, NA))
  expect_equal(res$duplicate_of[4], "c1")
  expect_true(is.na(res$duplicate_of[1]))

  # Strictness reasons survive as text, `;`-joined by the engine.
  expect_true(is.na(res$failures[1]))
  expect_match(res$failures[2], "not converged")
})

test_that("a cancelled run's partial table is readable", {
  dir <- tempfile()
  candidates_csv(dir, "candidates.partial.csv")

  res <- ferx_search_results(dir)
  expect_equal(nrow(res), 4L)
  expect_true(attr(res, "partial"))
  expect_match(attr(res, "path"), "candidates\\.partial\\.csv$")

  # A complete table beside it wins, unless the partial one is demanded.
  candidates_csv(dir)
  expect_false(attr(ferx_search_results(dir), "partial"))
  expect_true(attr(ferx_search_results(dir, partial = TRUE), "partial"))
})

test_that("a file path can be passed directly", {
  dir  <- tempfile()
  path <- candidates_csv(dir)
  res  <- ferx_search_results(path)
  expect_equal(nrow(res), 4L)
  expect_false(attr(res, "partial"))
})

test_that("a missing or wrong-shaped table is an error saying so", {
  dir <- tempfile()
  dir.create(dir)
  expect_error(ferx_search_results(dir), "No candidate table")
  expect_error(ferx_search_results(dir, partial = TRUE), "No partial candidate")

  candidates_csv(dir, "candidates.partial.csv")
  expect_error(ferx_search_results(dir, partial = FALSE), "No candidate table")

  writeLines(c("id,ofv", "c1,400"), file.path(dir, "candidates.csv"))
  expect_error(ferx_search_results(dir, partial = FALSE), "missing columns")

  expect_error(ferx_search_results(1), "single path")
  expect_error(ferx_search_results(dir, partial = "yes"), "TRUE, FALSE, or NULL")
})

test_that("the column list comes from the engine", {
  cols <- ferx:::ferx_rust_search_table_columns()
  expect_length(cols, 15L)
  expect_equal(cols[1:4], c("id", "parent", "hash", "features"))
  expect_true(all(c("criterion", "passed", "failures", "reused") %in% cols))
})
