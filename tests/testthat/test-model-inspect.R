write_ferx <- function(content, path = tempfile(fileext = ".ferx")) {
  writeLines(content, path)
  path
}

# --- .ferx_extract_blocks ---

test_that(".ferx_extract_blocks returns named list of block contents", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(1.0)",
    "  omega ETA_CL ~ 0.09",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)"
  ))
  on.exit(unlink(path))

  b <- ferx:::.ferx_extract_blocks(path)
  expect_type(b, "list")
  expect_true("parameters" %in% names(b))
  expect_true("structural_model" %in% names(b))
  expect_true(any(grepl("theta TVCL", b$parameters)))
  expect_true(any(grepl("pk one_cpt_oral", b$structural_model)))
})

test_that(".ferx_extract_blocks strips comments and blank lines", {
  path <- write_ferx(c(
    "[parameters]",
    "  # this is a comment",
    "  theta TVCL(1.0)  # inline comment",
    "",
    "  theta TVV(10.0)"
  ))
  on.exit(unlink(path))

  b <- ferx:::.ferx_extract_blocks(path)
  expect_equal(length(b$parameters), 2L)
  expect_false(any(grepl("#", b$parameters)))
})

test_that(".ferx_extract_blocks handles an empty file", {
  path <- write_ferx(character(0))
  on.exit(unlink(path))

  b <- ferx:::.ferx_extract_blocks(path)
  expect_equal(length(b), 0L)
})

# --- .ferx_model_type ---

test_that(".ferx_model_type detects all analytical model labels (matches Rust)", {
  # Labels mirror pk_model_type_label() in src/rust/src/lib.rs — keep in sync.
  # ferx-core #176 collapsed `*_iv_bolus` / `*_infusion` into `*_iv`; the
  # bolus-vs-infusion route is read per dose from the RATE column.
  expect_equal(ferx:::.ferx_model_type("pk one_cpt_oral(cl=CL, v=V, ka=KA)"),                                  "1-cpt oral")
  expect_equal(ferx:::.ferx_model_type("pk one_cpt_iv(cl=CL, v=V)"),                                           "1-cpt IV")
  expect_equal(ferx:::.ferx_model_type("pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2)"),                           "2-cpt oral")
  expect_equal(ferx:::.ferx_model_type("pk two_cpt_iv(cl=CL, v1=V1, q=Q, v2=V2)"),                             "2-cpt IV")
  expect_equal(ferx:::.ferx_model_type("pk three_cpt_oral(cl=CL, v1=V1, q2=Q2, v2=V2, q3=Q3, v3=V3, ka=KA)"),  "3-cpt oral")
  expect_equal(ferx:::.ferx_model_type("pk three_cpt_iv(cl=CL, v1=V1, q2=Q2, v2=V2, q3=Q3, v3=V3)"),           "3-cpt IV")
})

test_that(".ferx_model_type accepts the long `*_compartment_*` aliases", {
  expect_equal(ferx:::.ferx_model_type("pk one_compartment_oral(cl=CL, v=V, ka=KA)"),       "1-cpt oral")
  expect_equal(ferx:::.ferx_model_type("pk two_compartment_iv(cl=CL, v1=V1, q=Q, v2=V2)"),  "2-cpt IV")
})

test_that(".ferx_model_type still labels retired `*_iv_bolus` / `*_infusion` strings", {
  # These spellings were retired by ferx-core #176 and would fail at fit
  # time, but the R-side label helper accepts them so `ferx_model_inspect()`
  # on a stale .ferx file still reports something intelligible (the engine
  # will surface the migration error when the user tries to fit).
  expect_equal(ferx:::.ferx_model_type("pk one_cpt_iv_bolus(cl=CL, v=V)"),                      "1-cpt IV bolus")
  expect_equal(ferx:::.ferx_model_type("pk one_cpt_infusion(cl=CL, v=V)"),                      "1-cpt IV infusion")
  expect_equal(ferx:::.ferx_model_type("pk three_cpt_iv_bolus(cl=CL, v1=V1, q2=Q2, v2=V2, q3=Q3, v3=V3)"), "3-cpt IV bolus")
})

test_that(".ferx_model_type detects ODE at any position in lines", {
  expect_equal(ferx:::.ferx_model_type("ode(obs_cmt=central, states=[depot, central])"), "ODE")
  # ode() not at start of collapsed string when a preceding line exists
  expect_equal(
    ferx:::.ferx_model_type(c("some_prefix_line", "ode(obs_cmt=central, states=[c])")),
    "ODE"
  )
})

test_that(".ferx_model_type returns NULL for unrecognised content", {
  expect_null(ferx:::.ferx_model_type("TVCL * exp(ETA_CL)"))
  expect_null(ferx:::.ferx_model_type(character(0)))
})

# --- .ferx_parse_structure ---

test_that(".ferx_parse_structure parses a standard 1-cpt oral model", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2, 0.001, 10.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    "  theta TVKA(1.5, 0.01, 50.0)",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V  ~ 0.04",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$theta_names, c("TVCL", "TVV", "TVKA"))
  expect_equal(s$model_type,  "1-cpt oral")
  expect_equal(s$iiv,         c("ETA_CL", "ETA_V"))
  expect_equal(s$iov,         character(0))
  expect_equal(s$residual,    "proportional")
})

test_that(".ferx_parse_structure parses an ODE model", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVVMAX(4.0)",
    "  theta TVKM(6.0)",
    "  omega ETA_VMAX ~ 0.15",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  ode(obs_cmt=central, states=[depot, central])",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$model_type,  "ODE")
  expect_equal(s$theta_names, c("TVVMAX", "TVKM"))
  expect_equal(s$iiv,         "ETA_VMAX")
})

test_that(".ferx_parse_structure returns empty vectors when IIV/IOV absent", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ additive(PROP_ERR)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$iiv,      character(0))
  expect_equal(s$iov,      character(0))
  expect_equal(s$residual, "additive")
})

test_that(".ferx_parse_structure detects combined error", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ combined(PROP_ERR, ADD_ERR)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$residual, "combined")
})

test_that(".ferx_parse_structure detects LTBS via log(DV) ~ additive (case 2)", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  sigma ADD_LOG ~ 0.1",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  log(DV) ~ additive(ADD_LOG)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$residual, "additive (log-transformed)")
})

test_that(".ferx_parse_structure detects LTBS via DV ~ log_additive (case 1)", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  sigma ADD_LOG ~ 0.1",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ log_additive(ADD_LOG)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  # `log_additive` must not be misread as plain `additive`.
  expect_equal(s$residual, "additive (log-transformed)")
})

test_that(".ferx_parse_structure warns and returns 'unknown' for unrecognised error type", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "  sigma PROP_ERR ~ 0.02",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ custom_error(PROP_ERR)"
  ))
  on.exit(unlink(path))

  expect_warning(
    s <- ferx:::.ferx_parse_structure(path),
    regexp = "Unrecognised"
  )
  expect_equal(s$residual, "unknown")
})

test_that(".ferx_parse_structure returns 'unknown' when error_model block is absent", {
  path <- write_ferx(c(
    "[parameters]",
    "  theta TVCL(0.2)",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)"
  ))
  on.exit(unlink(path))

  s <- ferx:::.ferx_parse_structure(path)
  expect_equal(s$residual, "unknown")
})

# --- ferx_model_inspect ---

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
