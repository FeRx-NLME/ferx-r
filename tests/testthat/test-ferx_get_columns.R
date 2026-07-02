test_that("path input returns correct column names invisibly", {
  ex   <- ferx_example("warfarin")
  cols <- ferx_get_columns(ex$data)
  expect_type(cols, "character")
  expect_identical(cols, c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "RATE", "MDV"))
})
test_that("ferx_example() list input works", {
  ex   <- ferx_example("warfarin")
  cols <- ferx_get_columns(ex)
  expect_identical(cols, c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "RATE", "MDV"))
})
test_that("covariate dataset shows covariates / other group", {
  ex   <- ferx_example("two_cpt_oral_cov")
  cols <- ferx_get_columns(ex)
  expect_true(all(c("WT", "CRCL") %in% cols))
  expect_identical(cols[9:10], c("WT", "CRCL"))
})
test_that("output is invisible", {
  ex <- ferx_example("warfarin")
  expect_invisible(ferx_get_columns(ex$data))
})
test_that("output is printed to console", {
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_get_columns(ex$data))
  expect_true(any(grepl("Required NONMEM", out)))
  expect_true(any(grepl("Optional NONMEM", out)))
  expect_true(any(grepl("\\[1\\] ID", out)))
})
test_that("error on missing file", {
  expect_error(ferx_get_columns("/no/such/file.csv"), regexp = "not found")
})
test_that("error on wrong type", {
  expect_error(ferx_get_columns(42L), regexp = "file path")
})
test_that("error on NA path", {
  expect_error(ferx_get_columns(NA_character_), regexp = "file path")
})
test_that("error on length-2 character vector", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_get_columns(c(ex$data, ex$data)), regexp = "file path")
})
test_that("ferx_fit input resolves correctly", {
  fit  <- warfarin_fit()
  cols <- ferx_get_columns(fit)
  expect_type(cols, "character")
  expect_true("ID" %in% cols)
  expect_true("DV" %in% cols)
})
test_that("ferx_fit with NA data_path errors informatively", {
  fit            <- warfarin_fit()
  fit$data_path  <- NA_character_
  expect_error(ferx_get_columns(fit), regexp = "valid data_path")
})
test_that("list with invalid $data errors informatively", {
  expect_error(ferx_get_columns(list(data = 123L)), regexp = "\\$data")
})
test_that("empty file errors informatively", {
  tmp <- tempfile(fileext = ".csv")
  file.create(tmp)
  on.exit(unlink(tmp))
  expect_error(ferx_get_columns(tmp), regexp = "appears to be empty")
})
test_that("whitespace-only header errors informatively", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("   ", tmp)
  on.exit(unlink(tmp))
  expect_error(ferx_get_columns(tmp), regexp = "appears to be empty")
})
test_that("long column names wrap at narrow console width", {
  tmp <- tempfile(fileext = ".csv")
  long_cols <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT",
                 "COVARIATE_ONE", "COVARIATE_TWO", "COVARIATE_THREE",
                 "COVARIATE_FOUR", "COVARIATE_FIVE")
  writeLines(paste(long_cols, collapse = ","), tmp)
  on.exit(unlink(tmp))
  old_width <- options(width = 50L)
  on.exit(options(old_width), add = TRUE)
  out <- capture.output(cols <- ferx_get_columns(tmp))
  expect_identical(cols, long_cols)
  expect_true(length(out) > 5L)
})
test_that("single-column file prints 'column' not 'columns'", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("ID", tmp)
  on.exit(unlink(tmp))
  out <- capture.output(ferx_get_columns(tmp))
  expect_true(any(grepl("^1 column$", out)))
})
