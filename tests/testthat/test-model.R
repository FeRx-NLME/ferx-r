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

# Build a minimal ferx_fit stub for ferx_model_inspect() dispatch tests.
make_ferx_fit_stub <- function(model_structure = NULL, model_name = "test_model") {
  obj <- list(model_structure = model_structure, model_name = model_name)
  class(obj) <- "ferx_fit"
  obj
}

null_editor <- function(...) invisible(NULL)

# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# ferx_model_set_section
# ---------------------------------------------------------------------------

test_that("ferx_model_set_section() replaces a middle section and leaves others intact", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce",
    error_model = "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", c("  method = focei", "  maxiter = 500"))

  expect_equal(ferx_model_section(path, "parameters"), "  theta TVCL(1.0, 0.001, 100.0)")
  expect_equal(ferx_model_section(path, "fit_options"), c("  method = focei", "  maxiter = 500"))
  expect_equal(ferx_model_section(path, "error_model"), "  DV ~ proportional(PROP_ERR)")
})

test_that("ferx_model_set_section() replaces the last section without appending garbage", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
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
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", character(0))

  result <- ferx_model_section(path, "fit_options")
  expect_equal(result, character(0))
  # parameters section must be unaffected
  expect_equal(ferx_model_section(path, "parameters"), "  theta TVCL(1.0, 0.001, 100.0)")
})

test_that("ferx_model_set_section() round-trips a section via ferx_model_section()", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  original <- ferx_model_section(path, "parameters")
  modified <- sub("TVCL\\(.*\\)", "TVCL(0.5, 0.001, 10.0)", original)
  ferx_model_set_section(path, "parameters", modified)

  expect_equal(ferx_model_section(path, "parameters"), modified)
})

test_that("ferx_model_set_section() returns path invisibly", {
  path <- write_test_model(list(fit_options = "  method = foce"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_set_section(path, "fit_options", "  method = focei"))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("ferx_model_set_section() errors with available names when section is missing", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  expect_error(
    ferx_model_set_section(path, "odes", "  d/dt(central) = 0"),
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

test_that("ferx_model_edit() rejects malformed save_as before opening the editor", {
  # Use a user-owned (non-package) file so we exercise the save_as validation
  # path without triggering the in-package copy step. The error must fire
  # before utils::file.edit() is reached.
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  expect_error(
    ferx_model_edit(path, save_as = c("a.ferx", "b.ferx")),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
  expect_error(
    ferx_model_edit(path, save_as = 42),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
  expect_error(
    ferx_model_edit(path, save_as = NA_character_),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
})

# ---------------------------------------------------------------------------
# ferx_model_show
# ---------------------------------------------------------------------------

test_that("ferx_model_show() prints '# model: <filename>' header", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl(paste0("# model: ", basename(path)), out, fixed = TRUE)))
})

test_that("ferx_model_show() prints all file lines", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = c("  method = foce", "  maxiter = 300")
  ))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl("[parameters]",             out, fixed = TRUE)))
  expect_true(any(grepl("theta TVCL",               out, fixed = TRUE)))
  expect_true(any(grepl("[fit_options]",             out, fixed = TRUE)))
  expect_true(any(grepl("method = foce",             out, fixed = TRUE)))
  expect_true(any(grepl("maxiter = 300",             out, fixed = TRUE)))
})

test_that("ferx_model_show() returns path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_show(path))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("ferx_model_show() errors on missing file", {
  expect_error(
    ferx_model_show(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})

test_that("ferx_model_show() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_show(path), regexp = "\\.ferx")
})

test_that("ferx_model_show() handles an empty .ferx file without error", {
  path <- tempfile(fileext = ".ferx")
  writeLines(character(0), path)
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl("# model:", out, fixed = TRUE)))
})

test_that("ferx_model_show() handles a file containing only comment lines", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("# this is a comment", "# another comment"), path)
  on.exit(unlink(path))

  expect_no_error(out <- capture.output(ferx_model_show(path)))
  expect_true(any(grepl("# this is a comment", out, fixed = TRUE)))
})

test_that("ferx_model_show() pipes: ferx_model_new() |> ferx_model_show()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model_new(path, edit = FALSE) |> ferx_model_show()
  expect_equal(result, path)
})

# ---------------------------------------------------------------------------
# ferx_model_section — additional coverage
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# ferx_model_set_section — additional coverage
# ---------------------------------------------------------------------------

test_that("ferx_model_set_section() errors on missing file", {
  expect_error(
    ferx_model_set_section(file.path(tempdir(), "no_such.ferx"), "parameters", character(0)),
    regexp = "File not found"
  )
})

test_that("ferx_model_set_section() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_set_section(path, "parameters", character(0)), regexp = "\\.ferx")
})

test_that("ferx_model_set_section() preserves blank separator lines between sections", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)", "", "[fit_options]", "  method = foce"), path)
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", "  method = focei")

  all_lines <- readLines(path)
  expect_true("" %in% all_lines)
})

test_that("ferx_model_set_section() preserves comment lines in other sections", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "# a comment", "  theta TVCL(1.0, 0.001, 100.0)", "[fit_options]", "  method = foce"), path)
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", "  method = focei")

  all_lines <- readLines(path)
  expect_true("# a comment" %in% all_lines)
})

test_that("ferx_model_set_section() leaves file unchanged when section not found", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  before <- readLines(path)
  expect_error(ferx_model_set_section(path, "odes", "  d/dt(central) = 0"))
  expect_equal(readLines(path), before)
})

test_that("ferx_model_set_section() correctly replaces 1 line with 3 lines", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  new_lines <- c("  method = focei", "  maxiter = 500", "  covariance = false")
  ferx_model_set_section(path, "fit_options", new_lines)

  expect_equal(ferx_model_section(path, "fit_options"), new_lines)
})

test_that("ferx_model_set_section() pipe chain: two sequential set_section calls", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  path |>
    ferx_model_set_section("parameters", "  theta TVCL(2.0, 0.001, 100.0)") |>
    ferx_model_set_section("fit_options", "  method = focei")

  expect_equal(ferx_model_section(path, "parameters"), "  theta TVCL(2.0, 0.001, 100.0)")
  expect_equal(ferx_model_section(path, "fit_options"), "  method = focei")
})

test_that("ferx_model_set_section() pipe chain: new |> set_section |> show", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model_new(path, edit = FALSE) |>
      ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
      ferx_model_show()
  )

  expect_true(any(grepl("maxiter = 999", out, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# ferx_model_new — additional coverage
# ---------------------------------------------------------------------------

test_that("ferx_model_new() each template contains all required sections", {
  required <- c("parameters", "individual_parameters", "structural_model",
                "error_model", "fit_options")
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model_new(path, template = tmpl, edit = FALSE)
    content <- readLines(path)
    for (sec in required) {
      expect_true(
        any(grepl(paste0("^\\[", sec, "\\]"), content)),
        info = paste("template", tmpl, "missing section", sec)
      )
    }
  }
})

test_that("ferx_model_new() templates do not emit the deprecated [initial_values] block", {
  # Initial values now live inline in [parameters] — re-emitting them in a
  # separate block would be redundant boilerplate that the parser silently
  # ignores. Guards against regressing the issue-16 fix.
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model_new(path, template = tmpl, edit = FALSE)
    content <- readLines(path)
    expect_false(
      any(grepl("^\\[initial_values\\]", content)),
      info = paste("template", tmpl, "still emits [initial_values]")
    )
  }
})

test_that("ferx_model_new() print=TRUE output matches what would be written to file", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model_new(path, template = tmpl, edit = FALSE)
    file_lines <- readLines(path)

    printed <- capture.output(ferx_model_new(template = tmpl, print = TRUE))
    # capture.output adds a trailing "" for the final cat("\n"); drop it
    printed <- printed[nchar(printed) > 0 | seq_along(printed) < length(printed)]

    expect_equal(printed[printed != ""], file_lines[file_lines != ""],
                 info = paste("template", tmpl, "print vs file mismatch"))
  }
})

test_that("ferx_model_new() created file is readable by ferx_model_section()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model_new(path, edit = FALSE)
  expect_no_error(ferx_model_section(path, "parameters"))
})

test_that("ferx_model_new() created file is parseable by ferx_model_inspect()", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model_new(path, template = tmpl, edit = FALSE)
    result <- ferx_model_inspect(path)
    expect_false(is.null(result$model_type),
                 info = paste("template", tmpl, "model_type is NULL"))
  }
})

test_that("ferx_model_new() pipe chain: new |> set_section", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  path |>
    ferx_model_new(edit = FALSE) |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999"))

  expect_equal(ferx_model_section(path, "fit_options"), c("  method = focei", "  maxiter = 999"))
})

# ---------------------------------------------------------------------------
# ferx_model_inspect
# ---------------------------------------------------------------------------

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

test_that("ferx_model_inspect() detects 2-cpt IV bolus model type", {
  path <- write_test_model(list(
    parameters       = "  theta TVCL(5.0, 0.1, 100.0)",
    structural_model = "  pk two_cpt_iv_bolus(cl=CL, v1=V1, q=Q, v2=V2)",
    error_model      = "  DV ~ proportional(PROP_ERR)"
  ))
  on.exit(unlink(path))

  result <- ferx_model_inspect(path)
  expect_equal(result$model_type, "2-cpt IV bolus")
})

# ---------------------------------------------------------------------------
# ferx_model_edit — additional coverage
# ---------------------------------------------------------------------------

test_that("ferx_model_edit() save_as = FALSE collapses to NULL (no copy, returns path)", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_edit(path, save_as = FALSE, .editor = null_editor))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("ferx_model_edit() save_as = NULL returns path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_edit(path, save_as = NULL, .editor = null_editor))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("ferx_model_edit() save_as = 'dest' copies file to that path", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  ferx_model_edit(path, save_as = dest, .editor = null_editor)

  expect_true(file.exists(dest))
  expect_equal(readLines(dest), readLines(path))
})

test_that("ferx_model_edit() save_as returns the dest path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  result <- withVisible(ferx_model_edit(path, save_as = dest, .editor = null_editor))
  expect_equal(result$value, dest)
  expect_false(result$visible)
})

test_that("ferx_model_edit() save_as + overwrite = FALSE errors when dest exists", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  writeLines("existing", dest)
  on.exit({ unlink(path); unlink(dest) })

  expect_error(
    ferx_model_edit(path, save_as = dest, overwrite = FALSE, .editor = null_editor),
    regexp = "already exists"
  )
})

test_that("ferx_model_edit() save_as + overwrite = TRUE silently overwrites existing dest", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  writeLines("old content", dest)
  on.exit({ unlink(path); unlink(dest) })

  expect_no_error(
    ferx_model_edit(path, save_as = dest, overwrite = TRUE, .editor = null_editor)
  )
  expect_true(file.exists(dest))
})

test_that("ferx_model_edit() copies pkg file to dest before editing", {
  ex   <- ferx_example("warfarin")
  dest <- tempfile()
  dir.create(dest)
  on.exit(unlink(dest, recursive = TRUE))

  expect_message(
    ferx_model_edit(ex$model, dest = dest, .editor = null_editor),
    regexp = "Copied to"
  )
  expect_true(file.exists(file.path(dest, basename(ex$model))))
})

test_that("ferx_model_edit() post-edit copy: .editor can modify the file before save_as", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  modifying_editor <- function(p, ...) {
    lines <- readLines(p)
    writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
  }

  ferx_model_edit(path, save_as = dest, .editor = modifying_editor)

  dest_lines <- readLines(dest)
  expect_true(any(grepl("TVV", dest_lines)))
})

test_that("ferx_model_edit() errors on missing file", {
  expect_error(
    ferx_model_edit(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})
