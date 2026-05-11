test_that("ferx_predict returns a data frame", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data)
  expect_s3_class(pred, "data.frame")
})

test_that("ferx_predict has required columns: ID, TIME, PRED", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data)
  expect_true(all(c("ID", "TIME", "PRED") %in% names(pred)))
})

test_that("ferx_predict nrow matches observation rows in data", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data)
  dat  <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(pred), n_obs)
})

test_that("PRED is finite numeric with no NAs", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data)
  expect_true(is.numeric(pred$PRED))
  expect_true(all(is.finite(pred$PRED)))
})

test_that("PRED from ferx_predict(fit=) matches sdtab$PRED from ferx_fit()", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data, fit = warfarin_fit())
  sdtab_pred <- warfarin_fit()$sdtab$PRED
  expect_equal(pred$PRED, sdtab_pred, tolerance = 1e-6)
})

test_that("predict-from-fit returns data frame with required columns", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data, fit = warfarin_fit())
  expect_s3_class(pred, "data.frame")
  expect_true(all(c("ID", "TIME", "PRED") %in% names(pred)))
})

test_that("predict-from-fit PRED values are finite", {
  ex   <- ferx_example("warfarin")
  pred <- ferx_predict(ex$model, ex$data, fit = warfarin_fit())
  expect_true(all(is.finite(pred$PRED)))
})

test_that("predict-from-fit produces different PRED than population predict", {
  ex        <- ferx_example("warfarin")
  pred_pop  <- ferx_predict(ex$model, ex$data)
  pred_fit  <- ferx_predict(ex$model, ex$data, fit = warfarin_fit())
  expect_false(identical(pred_pop$PRED, pred_fit$PRED))
})

test_that("ferx_predict errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_predict("no_such_model.ferx", ex$data))
})

test_that("ferx_predict errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_predict(ex$model, "no_such_data.csv"))
})
