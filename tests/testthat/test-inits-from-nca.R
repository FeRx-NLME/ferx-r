test_that("ferx_inits_from_nca returns ferx_inits with expected shape", {
  ex <- ferx_example("warfarin")
  inits <- ferx_inits_from_nca(ex$model, ex$data)

  expect_s3_class(inits, "ferx_inits")
  expect_named(inits, c("theta", "omega", "method", "warnings"))
  expect_true(is.numeric(inits$theta))
  expect_false(is.null(names(inits$theta)))
  expect_true(length(inits$theta) > 0L)
  expect_true(is.matrix(inits$omega))
  expect_identical(inits$method, "nca_sweep")
})

test_that("ferx_inits_from_nca respects explicit method", {
  ex <- ferx_example("warfarin")

  inits_nca <- ferx_inits_from_nca(ex$model, ex$data, method = "nca")
  expect_identical(inits_nca$method, "nca")

  inits_ebe <- ferx_inits_from_nca(ex$model, ex$data, method = "nca_ebe")
  expect_true(inits_ebe$method %in% c("nca_ebe", "nca_sweep")) # may fall back for ODE
})

test_that("ferx_inits_from_nca rejects invalid method", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_inits_from_nca(ex$model, ex$data, method = "bogus"),
    "'arg'"
  )
})

test_that("print.ferx_inits emits theta and method", {
  ex <- ferx_example("warfarin")
  inits <- ferx_inits_from_nca(ex$model, ex$data, method = "nca")
  out <- capture.output(print(inits))
  expect_true(any(grepl("ferx_inits", out)))
  expect_true(any(grepl("Theta:", out)))
  expect_true(any(grepl("method = nca", out)))
})

test_that("ferx_fit inits_from_nca arg validation", {
  expect_error(
    ferx_fit(ferx_example("warfarin")$model,
             ferx_example("warfarin")$data,
             inits_from_nca = NA),
    "TRUE/FALSE"
  )
  expect_error(
    ferx_fit(ferx_example("warfarin")$model,
             ferx_example("warfarin")$data,
             inits_from_nca = "bogus"),
    "should be one of"
  )
})
