test_that("known name returns list with $model and $data", {
  ex <- ferx_example("warfarin")
  expect_true(is.list(ex))
  expect_true(all(c("model", "data") %in% names(ex)))
})
test_that("$model path exists on disk", {
  ex <- ferx_example("warfarin")
  expect_true(file.exists(ex$model))
})
test_that("$data path exists on disk", {
  ex <- ferx_example("warfarin")
  expect_true(file.exists(ex$data))
})
test_that("$model has .ferx extension", {
  ex <- ferx_example("warfarin")
  expect_true(grepl("\\.ferx$", ex$model))
})
test_that("$data has .csv extension", {
  ex <- ferx_example("warfarin")
  expect_true(grepl("\\.csv$", ex$data))
})
test_that("no-arg call returns character vector of available names", {
  res <- ferx_example()
  expect_type(res, "character")
  expect_true(length(res) > 0L)
  expect_true("warfarin" %in% res)
})
test_that("unknown name errors with available names in message", {
  expect_error(ferx_example("not_a_real_example"), regexp = "warfarin")
})
test_that("$data path exists on disk for every bundled example", {
  # Every example -- whether it ships its own dataset or resolves to one via the
  # .data_aliases map in ferx_example() -- must point at a CSV that exists.
  # Driven off ferx_example() rather than a hand-maintained list so newly added
  # examples and aliases are covered automatically (a previous hardcoded list
  # silently missed every alias added after it was written).
  for (nm in ferx_example()) {
    ex <- ferx_example(nm)
    expect_true(
      is.character(ex$data) && nzchar(ex$data) && file.exists(ex$data),
      label = paste0("$data for '", nm, "' must exist on disk")
    )
  }
})
test_that("emax_timecourse bundles a compartment-free model and synthetic data", {
  ex <- ferx_example("emax_timecourse")

  structural <- ferx:::.ferx_extract_blocks(ex$model)$structural_model
  expect_true(any(grepl("^EFF\\s*=", structural)))
  expect_true(any(grepl("^y\\s*=", structural)))
  expect_false(any(grepl("^(pk|ode|ode_template)\\b", structural)))
  capture.output(s <- ferx_model_inspect(ex$model))
  expect_equal(s$model_type, "compartment-free")

  dat <- utils::read.csv(ex$data)
  expect_identical(names(dat), c("ID", "TIME", "DV", "MDV"))
  expect_equal(nrow(dat), 240L)
  expect_equal(length(unique(dat$ID)), 30L)
  expect_false("AMT" %in% names(dat))
})
test_that("warfarin_derived $data points to warfarin.csv (alias check)", {
  ex <- ferx_example("warfarin_derived")
  expect_true(grepl("warfarin\\.csv$", ex$data))
})
