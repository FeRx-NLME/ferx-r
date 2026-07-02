# ferx_model_section
# ---------------------------------------------------------------------------



test_that("ferx_model(template=) created file is parseable by ferx_model_inspect()", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model(template = tmpl, path = path, edit = FALSE)
    result <- ferx_model_inspect(path)
    expect_false(is.null(result$model_type),
                 info = paste("template", tmpl, "model_type is NULL"))
  }
})
test_that("ferx_model_inspect() returns correct theta_names for warfarin", {
  ex <- ferx_example("warfarin")
  result <- ferx_model_inspect(ex$model)
  expect_equal(result$theta_names, c("TVCL", "TVV", "TVKA"))
})
test_that("ferx_model_inspect() returns correct model_type for warfarin", {
  ex <- ferx_example("warfarin")
  result <- ferx_model_inspect(ex$model)
  expect_equal(result$model_type, "1-cpt oral")
})
test_that("ferx_model_inspect() returns correct iiv for warfarin", {
  ex <- ferx_example("warfarin")
  result <- ferx_model_inspect(ex$model)
  expect_equal(result$iiv, c("ETA_CL", "ETA_V", "ETA_KA"))
})
test_that("ferx_model_inspect() returns iov = character(0) when no kappa lines", {
  ex <- ferx_example("warfarin")
  result <- ferx_model_inspect(ex$model)
  expect_equal(result$iov, character(0))
})
test_that("ferx_model_inspect() returns residual = 'proportional' for warfarin", {
  ex <- ferx_example("warfarin")
  result <- ferx_model_inspect(ex$model)
  expect_equal(result$residual, "proportional")
})
test_that("ferx_model_inspect() returns residual = 'additive'", {
  path <- write_test_model(list(
    parameters    = c("  theta TVCL(1.0, 0.001, 100.0)", "  sigma ADD_ERR ~ 0.01"),
    structural_model = "  pk one_cpt_oral(cl=1, v=10, ka=1)",
    error_model   = "  DV ~ additive(ADD_ERR)"
  ))
  on.exit(unlink(path))

  result <- ferx_model_inspect(path)
  expect_equal(result$residual, "additive")
})
test_that("ferx_model_inspect() returns residual = 'combined'", {
  path <- write_test_model(list(
    parameters    = c("  theta TVCL(1.0, 0.001, 100.0)", "  sigma ERR ~ 0.01"),
    structural_model = "  pk one_cpt_oral(cl=1, v=10, ka=1)",
    error_model   = "  DV ~ combined(ERR)"
  ))
  on.exit(unlink(path))

  result <- ferx_model_inspect(path)
  expect_equal(result$residual, "combined")
})
test_that("ferx_model_inspect() returns residual = 'unknown' and warns for unrecognised type", {
  path <- write_test_model(list(
    parameters    = "  theta TVCL(1.0, 0.001, 100.0)",
    structural_model = "  pk one_cpt_oral(cl=1, v=10, ka=1)",
    error_model   = "  DV ~ custom_error(ERR)"
  ))
  on.exit(unlink(path))

  expect_warning(result <- ferx_model_inspect(path), regexp = "Unrecognised")
  expect_equal(result$residual, "unknown")
})
test_that("ferx_model_inspect() returns list invisibly", {
  ex <- ferx_example("warfarin")
  result <- withVisible(ferx_model_inspect(ex$model))
  expect_false(result$visible)
  expect_type(result$value, "list")
})
test_that("ferx_model_inspect() prints filename in header", {
  ex <- ferx_example("warfarin")
  out <- capture.output(ferx_model_inspect(ex$model))
  expect_true(any(grepl(basename(ex$model), out, fixed = TRUE)))
})
test_that("ferx_model_inspect() errors on missing file", {
  expect_error(
    ferx_model_inspect(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})
test_that("ferx_model_inspect() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_inspect(path), regexp = "\\.ferx")
})
test_that("ferx_model_inspect() dispatches on ferx_fit object", {
  ms <- list(
    theta_names = c("TVCL", "TVV"),
    model_type  = "1-cpt oral",
    iiv         = c("ETA_CL", "ETA_V"),
    iov         = character(0),
    residual    = "proportional"
  )
  fit <- make_ferx_fit_stub(model_structure = ms, model_name = "mymodel")
  result <- ferx_model_inspect(fit)

  expect_equal(result$theta_names, c("TVCL", "TVV"))
  expect_equal(result$model_type,  "1-cpt oral")
  expect_equal(result$iiv,         c("ETA_CL", "ETA_V"))
  expect_equal(result$residual,    "proportional")
})
test_that("ferx_model_inspect() errors on ferx_fit with NULL model_structure", {
  fit <- make_ferx_fit_stub(model_structure = NULL)
  expect_error(ferx_model_inspect(fit), regexp = "No model_structure")
})
test_that("ferx_model_inspect() detects ODE model type", {
  path <- write_test_model(list(
    parameters       = "  theta TVPARAM(1.0, 0.001, 100.0)",
    structural_model = "  ode(obs_cmt=central, states=[depot, central])",
    error_model      = "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  result <- ferx_model_inspect(path)
  expect_equal(result$model_type, "ODE")
})
test_that("ferx_model_inspect() detects 2-cpt IV model type", {
  path <- write_test_model(list(
    parameters       = "  theta TVCL(5.0, 0.1, 100.0)",
    structural_model = "  pk two_cpt_iv(cl=CL, v1=V1, q=Q, v2=V2)",
    error_model      = "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  result <- ferx_model_inspect(path)
  expect_equal(result$model_type, "2-cpt IV")
})

# ---- header from test-model-inspect.R ----
write_ferx <- function(content, path = tempfile(fileext = ".ferx")) {
  writeLines(content, path)
  path
}

# --- .ferx_extract_blocks ---



test_that("ferx_model_inspect errors on missing file", {
  expect_error(ferx::ferx_model_inspect("/nonexistent/path.ferx"), "File not found")
})
test_that("ferx_model_inspect errors on non-.ferx extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("dummy", path)
  on.exit(unlink(path))
  expect_error(ferx::ferx_model_inspect(path), "\\.ferx")
})
test_that("ferx_model_inspect prints structure and returns list invisibly", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  theta TVV(10.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  out <- capture.output(s <- ferx::ferx_model_inspect(path))
  expect_type(s, "list")
  expect_named(s, c("theta_names", "model_type", "iiv", "iov", "residual"))
  expect_true(any(grepl("Structural", out)))
  expect_true(any(grepl("1-cpt oral", out)))
  expect_true(any(grepl("proportional", out)))
})
test_that("ferx_model_inspect accepts a ferx_fit object", {
  ms <- list(
    theta_names = c("TVCL", "TVV"),
    model_type  = "1-cpt oral",
    iiv         = c("ETA_CL", "ETA_V"),
    iov         = character(0),
    residual    = "proportional"
  )
  fit <- structure(
    list(model_structure = ms, model_name = "warfarin"),
    class = "ferx_fit"
  )

  out <- capture.output(s <- ferx::ferx_model_inspect(fit))
  expect_identical(s, ms)
  expect_true(any(grepl("warfarin", out)))
  expect_true(any(grepl("1-cpt oral", out)))
  expect_true(any(grepl("ETA_CL", out)))
})
test_that("ferx_model_inspect errors when ferx_fit has no model_structure", {
  fit <- structure(list(model_structure = NULL), class = "ferx_fit")
  expect_error(ferx::ferx_model_inspect(fit), "No model_structure")
})

# ---- header from test-model-structure-canonical.R ----
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

# ---- header from test-multi-endpoint.R ----
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
