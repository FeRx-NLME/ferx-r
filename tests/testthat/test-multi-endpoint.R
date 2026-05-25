# Multi-endpoint (per-CMT) residual error models — ferx-core #14.
#
# The emax_pkpd example declares a per-CMT [error_model] (proportional on the
# plasma endpoint CMT=2, additive on the PD endpoint CMT=3). These tests pin
# the pre-fit model-structure surface; they do not run the (heavy ODE) fit.

test_that("emax_pkpd example ships with model and data files", {
  ex <- ferx_example("emax_pkpd")
  expect_true(file.exists(ex$model))
  expect_true(file.exists(ex$data))
})

test_that("ferx_model_inspect reports the per-CMT residual structure", {
  ex <- ferx_example("emax_pkpd")
  out <- capture.output(s <- ferx_model_inspect(ex$model))

  # Residual label enumerates each endpoint's error model, ordered by CMT.
  expect_equal(s$residual, "per-CMT (CMT2=proportional, CMT3=additive)")
  expect_true(any(grepl("per-CMT", out)))
})

test_that("single-endpoint models still report a bare residual type", {
  ex <- ferx_example("warfarin")
  s <- ferx_model_inspect(ex$model)
  expect_equal(s$residual, "proportional")
})
