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

test_that("valid model returns TRUE", {
  ex  <- ferx_example("warfarin")
  res <- ferx_model_validate(ex$model)
  expect_true(isTRUE(res))
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
  expect_false(isTRUE(res))
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
  expect_false(isTRUE(res))
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
