# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Write a minimal .ferx file from a named list of section -> lines.
write_test_model <- function(sections) {
  path  <- tempfile(fileext = ".ferx")
  lines <- character(0)
  for (nm in names(sections)) {
    lines <- c(lines, paste0("[", nm, "]"), sections[[nm]])
  }
  writeLines(lines, path)
  path
}

# ---------------------------------------------------------------------------
# ferx_model_section
# ---------------------------------------------------------------------------

test_that("ferx_model_section() returns body lines of a named section", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "parameters")
  expect_equal(result, c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"))
})

test_that("ferx_model_section() returns the last section correctly", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = focei", "  maxiter = 300")
  ))
  on.exit(unlink(path))

  result <- ferx_model_section(path, "fit_options")
  expect_equal(result, c("  method = focei", "  maxiter = 300"))
})

test_that("ferx_model_section() returns character(0) for an empty section", {
  path <- write_test_model(list(
    parameters  = character(0),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  result <- suppressMessages(ferx_model_section(path, "parameters"))
  expect_equal(result, character(0))
})

test_that("ferx_model_section() errors with available names when section is missing", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce")
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

# ---------------------------------------------------------------------------
# ferx_model_set_section
# ---------------------------------------------------------------------------

test_that("ferx_model_set_section() replaces a middle section and leaves others intact", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce"),
    error_model = c("  DV ~ proportional(PROP_ERR)")
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", c("  method = focei", "  maxiter = 500"))

  expect_equal(ferx_model_section(path, "parameters"), c("  theta TVCL(1.0, 0.001, 100.0)"))
  expect_equal(ferx_model_section(path, "fit_options"), c("  method = focei", "  maxiter = 500"))
  expect_equal(ferx_model_section(path, "error_model"), c("  DV ~ proportional(PROP_ERR)"))
})

test_that("ferx_model_set_section() replaces the last section without appending garbage", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", c("  method = focei", "  maxiter = 500"))

  result <- ferx_model_section(path, "fit_options")
  expect_equal(result, c("  method = focei", "  maxiter = 500"))

  # No extra content beyond the last section
  all_lines <- readLines(path)
  expect_false(any(is.na(all_lines)))
  last_section_start <- which(all_lines == "[fit_options]")
  expect_equal(all_lines[(last_section_start + 1):length(all_lines)],
               c("  method = focei", "  maxiter = 500"))
})

test_that("ferx_model_set_section() replaces a section with empty lines", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", character(0))

  result <- ferx_model_section(path, "fit_options")
  expect_equal(result, character(0))
  # parameters section must be unaffected
  expect_equal(ferx_model_section(path, "parameters"), c("  theta TVCL(1.0, 0.001, 100.0)"))
})

test_that("ferx_model_set_section() round-trips a section via ferx_model_section()", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  original <- ferx_model_section(path, "parameters")
  modified <- sub("TVCL\\(.*\\)", "TVCL(0.5, 0.001, 10.0)", original)
  ferx_model_set_section(path, "parameters", modified)

  expect_equal(ferx_model_section(path, "parameters"), modified)
})

test_that("ferx_model_set_section() returns path invisibly", {
  path <- write_test_model(list(fit_options = c("  method = foce")))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_set_section(path, "fit_options", c("  method = focei")))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("ferx_model_set_section() errors with available names when section is missing", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)"),
    fit_options = c("  method = foce")
  ))
  on.exit(unlink(path))

  expect_error(
    ferx_model_set_section(path, "odes", c("  d/dt(central) = 0")),
    regexp = "parameters, fit_options"
  )
})

# ---------------------------------------------------------------------------
# ferx_model_new
# ---------------------------------------------------------------------------

test_that("ferx_model_new(print = TRUE) prints to console and returns NULL invisibly", {
  result <- withVisible(ferx_model_new(print = TRUE))
  expect_null(result$value)
  expect_false(result$visible)
})

test_that("ferx_model_new(print = TRUE) works for each template", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    expect_no_error(ferx_model_new(template = tmpl, print = TRUE))
  }
})

test_that("ferx_model_new() creates a file containing the required sections", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model_new(path, edit = FALSE)

  expect_true(file.exists(path))
  content <- readLines(path)
  expect_true(any(grepl("^\\[parameters\\]", content)))
  expect_true(any(grepl("^\\[fit_options\\]", content)))
})

test_that("ferx_model_new() errors when file exists and overwrite = FALSE", {
  path <- tempfile(fileext = ".ferx")
  writeLines("existing", path)
  on.exit(unlink(path))

  expect_error(ferx_model_new(path, edit = FALSE), regexp = "already exists")
})

test_that("ferx_model_new() overwrites when overwrite = TRUE", {
  path <- tempfile(fileext = ".ferx")
  writeLines("old content", path)
  on.exit(unlink(path))

  expect_no_error(ferx_model_new(path, overwrite = TRUE, edit = FALSE))
  expect_true(any(grepl("\\[parameters\\]", readLines(path))))
})

test_that("ferx_model_new() errors on unknown template", {
  expect_error(
    ferx_model_new(print = TRUE, template = "3cpt_magical"),
    regexp = "Unknown template"
  )
})

test_that("ferx_model_new() errors when path is NULL and print = FALSE", {
  expect_error(ferx_model_new(print = FALSE), regexp = "'path' is required")
})

test_that("ferx_model_new() errors when path has wrong extension", {
  expect_error(
    ferx_model_new(tempfile(fileext = ".txt"), edit = FALSE),
    regexp = "\\.ferx"
  )
})

test_that("ferx_model_new() returns path invisibly", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- withVisible(ferx_model_new(path, edit = FALSE))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

# ---------------------------------------------------------------------------
# ferx_model_edit — overwrite guard (error path only; success opens an editor)
# ---------------------------------------------------------------------------

test_that("ferx_model_edit() errors when dest file exists and overwrite = FALSE", {
  ex   <- ferx_example("warfarin")
  dest <- tempfile()
  dir.create(dest)
  on.exit(unlink(dest, recursive = TRUE))

  file.copy(ex$model, file.path(dest, basename(ex$model)))

  expect_error(
    ferx_model_edit(ex$model, dest = dest),
    regexp = "already exists"
  )
})
