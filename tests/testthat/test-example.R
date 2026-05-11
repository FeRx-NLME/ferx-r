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
