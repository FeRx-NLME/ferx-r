
# ---- header from test-model.R ----
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
    ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
      ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
      ferx_model_show()
  )

  expect_true(any(grepl("maxiter = 999", out, fixed = TRUE)))
})

# ---- header from test-model-pipe.R ----
# Tests for the ferx_model S3 pipe object and pipe-friendly wrappers added in #47:
#   ferx_model()       — constructor
#   print.ferx_model() — console summary
#   ferx_set_section() — pipe-friendly section replacement
#   ferx_get_section() — pipe-friendly section display
#   ferx_fit()         — ferx_model dispatch (inline, see #52)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_pipe_test_model <- function() {
  path <- tempfile(fileext = ".ferx")
  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  path
}

modifying_editor <- function(p, ...) {
  lines <- readLines(p)
  writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
}

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------



test_that("ferx_set_section() on a plain path delegates to ferx_model_set_section()", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_set_section(path, "fit_options", c("  method = focei", "  maxiter = 42"))
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 42", opts)))
})
test_that("pipe chain: ferx_model(template=) |> ferx_model_set_section() |> ferx_model_inspect()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_inspect()

  expect_equal(result$model_type, "1-cpt oral")
})
