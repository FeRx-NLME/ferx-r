# Unit tests for .ferx_format_structural(): builds the one-line structural
# summary from a model_structure list. Pure string formatting, no fit needed.

test_that(".ferx_format_structural combines type and theta names", {
  ms <- list(model_type = "1-cpt oral", theta_names = c("CL", "V", "KA"))
  expect_identical(ferx:::.ferx_format_structural(ms), "1-cpt oral  (CL, V, KA)")
})

test_that(".ferx_format_structural uses type alone when no thetas", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = "ODE", theta_names = character(0))),
    "ODE"
  )
})

test_that(".ferx_format_structural uses theta names alone when no type", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = NULL, theta_names = c("CL", "V"))),
    "CL, V"
  )
})

test_that(".ferx_format_structural falls back to 'unknown' when nothing is known", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = NULL, theta_names = character(0))),
    "unknown"
  )
})
