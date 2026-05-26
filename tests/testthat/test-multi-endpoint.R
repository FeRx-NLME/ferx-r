# Multi-endpoint (per-CMT) residual error models - ferx-core #14.
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

test_that("unrecognised per-CMT error type warns and reports 'unknown' (no NA label)", {
  # Pre-fit reparser path: an unknown error type on a CMT= line must warn and
  # fall back to "unknown" rather than emitting a label like "CMT2=NA".
  f <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(1.0)",
    "  omega ETA_CL ~ 0.1",
    "  sigma S1 ~ 0.1",
    "  sigma S2 ~ 0.1",
    "[individual_parameters]",
    "  CL = TVCL",
    "[structural_model]",
    "  ode(states=[central, effect])",
    "[error_model]",
    "  CMT=2: DV ~ bogus(S1)",
    "  CMT=3: DV ~ additive(S2)"
  ), f)
  on.exit(unlink(f), add = TRUE)

  expect_warning(
    s <- suppressMessages(ferx_model_inspect(f)),
    "Unrecognised per-CMT"
  )
  expect_equal(s$residual, "unknown")
  expect_false(grepl("NA", s$residual, fixed = TRUE))
})
