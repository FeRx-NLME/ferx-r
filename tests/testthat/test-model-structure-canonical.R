# Pins the contract that ferx_fit() returns `$model_structure` sourced from
# the Rust engine (built from the parsed CompiledModel) rather than re-parsed
# in R from the .ferx file. See ferx-core#49.
#
# Uses the shared `warfarin_fit()` helper (helper-warfarin-fit.R) so the
# real FOCEI fit is reused across test files.

test_that("fit$model_structure has the documented shape", {
  fit <- warfarin_fit()
  ms <- fit$model_structure
  expect_type(ms, "list")
  expect_named(ms, c("theta_names", "model_type", "iiv", "iov", "residual"))
})

test_that("fit$model_structure reflects what ferx-core actually parsed (warfarin = 1-cpt oral, proportional)", {
  fit <- warfarin_fit()
  ms <- fit$model_structure

  # Warfarin declares TVCL, TVV, TVKA in [parameters].
  expect_setequal(ms$theta_names, c("TVCL", "TVV", "TVKA"))

  # Structural form is unambiguous (one_cpt_oral) — model_type is set, not NULL.
  expect_equal(ms$model_type, "1-cpt oral")

  # IIV: one omega per PK parameter.
  expect_setequal(ms$iiv, c("ETA_CL", "ETA_V", "ETA_KA"))

  # No IOV in the warfarin example.
  expect_identical(ms$iov, character(0))

  # Proportional residual error model.
  expect_equal(ms$residual, "proportional")
})

test_that("ferx_model_inspect(fit) reads from the Rust-supplied structure post-fit", {
  fit <- warfarin_fit()
  out <- capture.output(s <- ferx_model_inspect(fit))
  # Returned invisibly — must equal the attached structure.
  expect_identical(s, fit$model_structure)
  expect_true(any(grepl("1-cpt oral", out)))
  expect_true(any(grepl("proportional", out)))
})
