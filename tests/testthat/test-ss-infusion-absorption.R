# Steady-state and infusion dosing into a built-in absorption compartment
# (ferx-core #719 gaps 1 & 2). Both bundled example datasets carry the NONMEM
# 7.6.0 (ADVAN2 TRANS2, S2 = V) population prediction in their DV column, so
# ferx_predict() at the (NONMEM-matched) typical values must reproduce it - a
# cross-engine anchor. Both combinations were previously rejected at parse time.

test_that("ss_absorption: SS dosing into a built-in absorption compartment predicts", {
  ex   <- ferx_example("ss_absorption")
  pred <- ferx_predict(ex$model, ex$data)
  expect_s3_class(pred, "data.frame")
  expect_true(all(c("ID", "TIME", "PRED") %in% names(pred)))
  expect_true(all(is.finite(pred$PRED)))
})

test_that("ss_absorption: ferx_predict matches the NONMEM anchor (< 1e-4)", {
  ex   <- ferx_example("ss_absorption")
  dat  <- read.csv(ex$data)
  obs  <- dat[dat$EVID == 0, c("TIME", "DV")]
  obs$DV <- as.numeric(obs$DV)  # DV is character: dose rows use "." as a placeholder
  pred <- ferx_predict(ex$model, ex$data)
  cmp  <- merge(obs, pred[, c("TIME", "PRED")], by = "TIME")
  expect_equal(nrow(cmp), nrow(obs))
  # DV on the observation rows is the NONMEM PRED; ferx must match to < 1e-4.
  # Assert the per-point MAX relative error (matches the documented claim; a
  # mean-based tolerance could hide a single diluted outlier).
  expect_lt(max(abs(cmp$PRED - cmp$DV) / cmp$DV), 1e-4)
})

test_that("infusion_absorption: infusion into a built-in absorption compartment predicts", {
  ex   <- ferx_example("infusion_absorption")
  pred <- ferx_predict(ex$model, ex$data)
  expect_s3_class(pred, "data.frame")
  expect_true(all(c("ID", "TIME", "PRED") %in% names(pred)))
  expect_true(all(is.finite(pred$PRED)))
})

test_that("infusion_absorption: ferx_predict matches the NONMEM anchor (< 1e-4)", {
  ex   <- ferx_example("infusion_absorption")
  dat  <- read.csv(ex$data)
  obs  <- dat[dat$EVID == 0, c("TIME", "DV")]
  obs$DV <- as.numeric(obs$DV)  # DV is character: dose rows use "." as a placeholder
  pred <- ferx_predict(ex$model, ex$data)
  cmp  <- merge(obs, pred[, c("TIME", "PRED")], by = "TIME")
  expect_equal(nrow(cmp), nrow(obs))
  expect_lt(max(abs(cmp$PRED - cmp$DV) / cmp$DV), 1e-4)
})

test_that("both new absorption examples are discovered with data on disk", {
  av <- ferx_example()
  expect_true(all(c("ss_absorption", "infusion_absorption") %in% av))
  for (nm in c("ss_absorption", "infusion_absorption")) {
    ex <- ferx_example(nm)
    expect_true(file.exists(ex$model))
    expect_true(file.exists(ex$data))
  }
})
