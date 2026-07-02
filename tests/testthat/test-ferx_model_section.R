# ferx_model_section
# ---------------------------------------------------------------------------



test_that("ferx_model_section() returns body lines of a named section", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "parameters")
  expect_equal(result, c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"))
})
test_that("ferx_model_section() returns the last section correctly", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = c("  method = focei", "  maxiter = 300")
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "fit_options")
  expect_equal(result, c("  method = focei", "  maxiter = 300"))
})
test_that("ferx_model_section() returns character(0) for an empty section", {
  path <- write_test_model(list(
    parameters  = character(0),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- suppressMessages(ferx_model_section(path, "parameters"))
  expect_equal(result, character(0))
})
test_that("ferx_model_section() errors with available names when section is missing", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  expect_error(
    ferx_model_section(path, "odes"),
    regexp = "parameters, fit_options"
  )
})
test_that("ferx_model_section() errors on non-existent file", {
  expect_error(
    ferx_model_section(file.path(tempdir(), "no_such.ferx"), "parameters"),
    regexp = "File not found"
  )
})
test_that("ferx_model_section() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_section(path, "parameters"), regexp = "\\.ferx")
})
test_that("ferx_model_section(strip = FALSE) preserves leading whitespace (default)", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)", "    omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  expect_equal(
    ferx_model_section(path, "parameters"),
    c("  theta TVCL(1.0, 0.001, 100.0)", "    omega ETA_CL ~ 0.09")
  )
})
test_that("ferx_model_section(strip = TRUE) trims only leading whitespace", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)  ", "    omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "parameters", strip = TRUE)
  # Leading whitespace removed
  expect_equal(
    result,
    c("theta TVCL(1.0, 0.001, 100.0)  ", "omega ETA_CL ~ 0.09")
  )
  # Internal spacing untouched (between tokens)
  expect_true(grepl("theta TVCL", result[1]))
  expect_true(grepl("omega ETA_CL ~ 0.09", result[2]))
  # Trailing whitespace not trimmed
  expect_true(grepl("  $", result[1]))
})
test_that("ferx_model_section() prints '# [section]' header to console", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_section(path, "parameters"))
  expect_true(any(grepl("# [parameters]", out, fixed = TRUE)))
})
test_that("ferx_model_section() returns body invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(suppressMessages(ferx_model_section(path, "parameters")))
  expect_false(result$visible)
})
test_that("ferx_model_section() preserves blank lines within section body", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)", "", "  omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "parameters")
  expect_equal(result, c("  theta TVCL(1.0, 0.001, 100.0)", "", "  omega ETA_CL ~ 0.09"))
})
test_that("ferx_model_section() parses header with surrounding whitespace", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("  [parameters]  ", "  theta TVCL(1.0, 0.001, 100.0)"), path)
  on.exit(unlink(path))

  result <- ferx_model_section(path, "parameters")
  expect_equal(result, "  theta TVCL(1.0, 0.001, 100.0)")
})
test_that("ferx_model_section(strip = TRUE) on empty section returns character(0)", {
  path <- write_test_model(list(
    parameters  = character(0),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- suppressMessages(ferx_model_section(path, "parameters", strip = TRUE))
  expect_equal(result, character(0))
})
test_that("ferx_model(template=) created file is readable by ferx_model_section()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  expect_no_error(ferx_model_section(path, "parameters"))
})
