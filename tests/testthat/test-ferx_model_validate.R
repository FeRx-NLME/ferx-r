
# ---- header from test-model-validate.R ----
# Helper: write a minimal valid warfarin-like .ferx file
write_ferx <- function(sections, path = tempfile(fileext = ".ferx")) {
  writeLines(sections, path)
  path
}

VALID_WARFARIN_SECTIONS <- c(
  "[parameters]",
  "  theta TVCL(0.134, 0.001, 10.0)",
  "  omega ETA_CL ~ 0.07",
  "  sigma PROP_ERR ~ 0.01",
  "",
  "[individual_parameters]",
  "  CL = TVCL * exp(ETA_CL)",
  "",
  "[structural_model]",
  "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
  "",
  "[error_model]",
  "  DV ~ proportional(PROP_ERR)"
)



test_that("valid warfarin model prints 'Validating:' header and returns invisibly", {
  ex  <- ferx_example("warfarin")
  out <- capture.output(res <- withVisible(ferx_model_validate(ex$model)))
  expect_true(any(grepl("Validating:", out)))
  expect_false(res$visible)
})
test_that("valid model returns ok = TRUE", {
  ex  <- ferx_example("warfarin")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
  expect_s3_class(res$diagnostics, "data.frame")
  expect_named(res$diagnostics,
               c("severity", "code", "message", "block", "line", "suggestion"))
})
test_that("valid model: $missing_sections is empty (no missing sections reported)", {
  # ferx_model_validate returns a scalar boolean, not a list.
  # We verify indirectly: output contains no [MISSING] markers.
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_model_validate(ex$model))
  expect_false(any(grepl("\\[MISSING\\]", out)))
})
test_that("valid model: $unknown_sections is empty (no unknown section markers)", {
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_model_validate(ex$model))
  expect_false(any(grepl("\\[unknown section\\]", out)))
})
test_that("model missing a required section returns FALSE", {
  path <- write_ferx(grep("^\\[error_model\\]", VALID_WARFARIN_SECTIONS,
                          value = TRUE, invert = TRUE))
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
})
test_that("model missing a required section reports the missing section name", {
  path <- write_ferx(grep("^\\[error_model\\]|^  DV ~", VALID_WARFARIN_SECTIONS,
                          value = TRUE, invert = TRUE))
  on.exit(unlink(path))
  out <- capture.output(ferx_model_validate(path))
  expect_true(any(grepl("error_model", out)))
})
test_that("model with unknown section reports it", {
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[foo]", "  bar = 1")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  out <- capture.output(ferx_model_validate(path))
  expect_true(any(grepl("unknown section|foo", out)))
})
test_that("model with unknown section is not ok", {
  # Regression for ferx-core #1040: `unknown` was printed as `[unknown section]`
  # but left out of the returned status, so `res$ok` was TRUE for a model
  # carrying a block the engine ignores.
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[foo]", "  bar = 1")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  res <- suppressWarnings(capture.output(r <- ferx_model_validate(path)))
  expect_false(isTRUE(r$ok))
})
test_that("the engine rejects an unknown section with E_UNKNOWN_BLOCK", {
  # The status above must not rest on the R-side list alone: the engine itself
  # refuses the header, so `ok` stays FALSE even if the R check were removed.
  lines <- c(VALID_WARFARIN_SECTIONS, "", "[fit_option]", "  method = focei")
  path  <- write_ferx(lines)
  on.exit(unlink(path))
  invisible(capture.output(res <- ferx_model_validate(path)))
  expect_false(isTRUE(res$ok))
  expect_true("E_UNKNOWN_BLOCK" %in% res$diagnostics$code)
  d <- res$diagnostics[res$diagnostics$code == "E_UNKNOWN_BLOCK", , drop = FALSE]
  expect_identical(d$block[[1L]], "fit_option")
  expect_true(grepl("fit_options", d$suggestion[[1L]], fixed = TRUE))
})
test_that("optional sections come from the engine, not a duplicated R list", {
  # The hard-coded vector this replaces had drifted: it omitted `covariates`
  # (and event_model, simulation, ...) so a valid bundled example printed
  # `covariates [unknown section]`, and it listed `initial_values`, which the
  # engine dropped long ago (ferx-core e5e934d).
  known <- ferx:::ferx_rust_known_blocks()
  expect_true(all(c("parameters", "individual_parameters", "structural_model",
                    "error_model", "covariates", "event_model", "simulation",
                    "adaptive_dosing", "mixture", "data_selection") %in% known))
  expect_false("initial_values" %in% known)
  expect_identical(known, sort(known))
})
test_that("a bundled example with a [covariates] block validates clean", {
  ex <- ferx_example("two_cpt_oral_cov")
  out <- capture.output(res <- ferx_model_validate(ex$model))
  expect_true(isTRUE(res$ok))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)))
})
test_that("syntax error in model body returns FALSE and prints errors", {
  lines <- c(
    "[parameters]",
    "  theta !!!BAD SYNTAX!!!",
    "",
    "[individual_parameters]",
    "  CL = 1.0",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "",
    "[error_model]",
    "  DV ~ proportional(0.01)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
})
test_that("data-dependent checks fire only when data is supplied", {
  # Reference a covariate that is not in the warfarin dataset.
  lines <- c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  omega ETA_CL ~ 0.07",
    "  sigma PROP_ERR ~ 0.01",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL) * (NOT_A_COL / 70)",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  ex <- ferx_example("warfarin")

  # Without data: the missing-covariate check cannot fire.
  res_no_data <- ferx_model_validate(path)
  expect_false(
    "E_MISSING_COVARIATE" %in% res_no_data$diagnostics$code
  )

  # With data: the check fires with a stable error code.
  res_with_data <- ferx_model_validate(path, data = ex$data)
  expect_false(isTRUE(res_with_data$ok))
  expect_true("E_MISSING_COVARIATE" %in% res_with_data$diagnostics$code)
})
test_that("diagnostics carry stable codes and severity labels", {
  lines <- c(
    "[parameters]",
    "  theta !!!BAD!!!",
    "[individual_parameters]",
    "  CL = 1.0",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "[error_model]",
    "  DV ~ proportional(0.01)"
  )
  path <- write_ferx(lines)
  on.exit(unlink(path))
  res <- ferx_model_validate(path)
  expect_false(isTRUE(res$ok))
  expect_gt(nrow(res$diagnostics), 0L)
  expect_true(all(res$diagnostics$severity %in% c("error", "warning")))
  expect_true(all(grepl("^[EW]_", res$diagnostics$code)))
})
test_that("ferx_model_validate errors on missing file", {
  expect_error(
    ferx_model_validate("no_such_file.ferx"),
    regexp = "File not found",
    ignore.case = TRUE
  )
})
test_that("ferx_model_validate errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))
  expect_error(
    ferx_model_validate(path),
    regexp = "\\.ferx",
    ignore.case = TRUE
  )
})
test_that("[derived] and [output] sections are accepted as optional, not flagged as unknown", {
  path <- write_ferx(c(
    VALID_WARFARIN_SECTIONS,
    "",
    "[derived]",
    "  KE = CL / 10",
    "",
    "[output]",
    "  CL"
  ))
  on.exit(unlink(path))
  out <- capture.output(suppressWarnings(ferx_model_validate(path)))
  expect_false(any(grepl("unknown section", out, fixed = TRUE)),
               label = "[derived] and [output] must not appear as unknown sections")
  expect_true(any(grepl("derived", out, ignore.case = TRUE)),
              label = "[derived] section should be printed as present")
  expect_true(any(grepl("output", out, ignore.case = TRUE)),
              label = "[output] section should be printed as present")
})
test_that("bundled igd_inverse_gaussian example validates (igd() parses against the engine)", {
  # Guards that the bundled inverse-Gaussian example stays parseable by the
  # pinned ferx-core: validation runs the model through the Rust engine, so a
  # valid result confirms the built-in `igd(mat, cv2)` input rate is understood
  # (the data/alias are covered by test-example.R; a full fit is skipped — the
  # anchor data is mis-specified for IG, so it stalls rather than converges).
  ex  <- ferx_example("igd_inverse_gaussian")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled biphasic_igd_absorption example validates (FR1*igd + FR2*igd parses)", {
  # Guards the shared pathway-fraction multiplier (ferx-core #388): a valid result
  # confirms the engine accepts a declared-parameter fraction (`FR*igd(...)`) and
  # two input-rate terms on one compartment, and that the fit-init fraction checks
  # (0 < FR <= 1, Σ FR ≈ 1) pass for FR1 / FR2 = 1 - FR1. The matched, recovering
  # fit is exercised by the ferx-core NONMEM anchor (freijer_biphasic_ig), not here.
  ex  <- ferx_example("biphasic_igd_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled weibull_absorption example validates (weibull() parses against the engine)", {
  # Guards that the bundled Weibull example stays parseable by the pinned
  # ferx-core: a valid result confirms the built-in `weibull(td, beta)` input
  # rate is understood by the engine (data/alias covered by test-example.R; a
  # full fit is exercised by the ferx-core NONMEM anchor, not here).
  ex  <- ferx_example("weibull_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled zero_order_absorption example validates (zero_order() parses against the engine)", {
  # Guards that the bundled zero-order example stays parseable by the pinned
  # ferx-core: a valid result confirms the built-in `zero_order(dur)` input rate
  # is understood by the engine (data/alias covered by test-example.R; the value
  # path is anchored against a NONMEM ADVAN1 modeled-duration run in ferx-core).
  ex  <- ferx_example("zero_order_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled parallel_absorption example validates (FR1*first_order + FR2*first_order parses)", {
  # Guards parallel dual-first-order absorption (ferx-core #505): a valid result
  # confirms the engine accepts the built-in `first_order(ka)` input rate composed
  # as two fraction-weighted pathways on one compartment, and that the fit-init
  # fraction checks (0 < FR <= 1, sum FR ~ 1) pass for FR1 / FR2 = 1 - FR1. The
  # matched, recovering fit is exercised by the ferx-core NONMEM anchor, not here.
  ex  <- ferx_example("parallel_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled mixed_absorption example validates (zero_order + first_order compose with fractions)", {
  # Guards mixed zero-order + first-order absorption (ferx-core #505): a valid
  # result confirms a fraction-weighted `zero_order(dur)` and `first_order(ka)` on
  # one compartment parse against the pinned engine -- the per-segment zero-order
  # channel now carries the pathway fraction (anchored to NONMEM in ferx-core).
  ex  <- ferx_example("mixed_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled per_route_lag_absorption example validates (per-route lag= parses on a parallel route)", {
  # Guards per-route absorption lag (ferx-core #856): a valid result confirms the
  # pinned engine accepts the optional `lag=` argument on a fraction-weighted
  # `first_order(ka, lag=LAG2)` pathway (IR + delayed-release), a strict superset of
  # the parallel first_order composition. The matched, recovering fit is the
  # ferx-core NONMEM $DES anchor (per_route_lag.ctl), not exercised here.
  ex  <- ferx_example("per_route_lag_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
test_that("bundled sequential_absorption example validates (zero_order + first-order compose)", {
  # Guards the sequential (zero-order depot fill -> first-order ka) composition:
  # a valid result confirms zero_order(dur) into a depot plus a hand-written
  # `- KA*depot` transfer parses against the engine (its own simulated dataset is
  # covered by test-example.R; a full fit recovers the simulation truth offline).
  ex  <- ferx_example("sequential_absorption")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res$ok))
})
