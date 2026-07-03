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

  expect_equal(ferx_model_get_section(path, "parameters"), "  theta TVCL(1.0, 0.001, 100.0)")
  expect_equal(ferx_model_get_section(path, "fit_options"), c("  method = focei", "  maxiter = 500"))
  expect_equal(ferx_model_get_section(path, "error_model"), "  DV ~ proportional(PROP_ERR)")
})
test_that("ferx_model_set_section() replaces the last section without appending garbage", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  ferx_model_set_section(path, "fit_options", c("  method = focei", "  maxiter = 500"))

  result <- ferx_model_get_section(path, "fit_options")
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

  result <- ferx_model_get_section(path, "fit_options")
  expect_equal(result, character(0))
  # parameters section must be unaffected
  expect_equal(ferx_model_get_section(path, "parameters"), "  theta TVCL(1.0, 0.001, 100.0)")
})
test_that("ferx_model_set_section() round-trips a section via ferx_model_get_section()", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  original <- ferx_model_get_section(path, "parameters")
  modified <- sub("TVCL\\(.*\\)", "TVCL(0.5, 0.001, 10.0)", original)
  ferx_model_set_section(path, "parameters", modified)

  expect_equal(ferx_model_get_section(path, "parameters"), modified)
})
test_that("ferx_model_set_section() on a plain path returns the path invisibly", {
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

  expect_equal(ferx_model_get_section(path, "fit_options"), new_lines)
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

  expect_equal(ferx_model_get_section(path, "parameters"), "  theta TVCL(2.0, 0.001, 100.0)")
  expect_equal(ferx_model_get_section(path, "fit_options"), "  method = focei")
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
test_that("pipe chain: ferx_model(template=) |> ferx_model_set_section() |> ferx_model_inspect()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_inspect()

  expect_equal(result$model_type, "1-cpt oral")
})

# ---------------------------------------------------------------------------
# ferx_model object input (pipe + copy-on-write)
# ---------------------------------------------------------------------------

test_that("ferx_model_set_section() on ferx_model returns a ferx_model", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- ferx_model_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  expect_s3_class(result, "ferx_model")
})
test_that("ferx_model_set_section() on ferx_model returns the same object (same $model path)", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(model = path)
  result <- ferx_model_set_section(m, "fit_options", "  method = focei")
  expect_equal(result$model, m$model)
})
test_that("ferx_model_set_section() on ferx_model writes the change to disk", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_model_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  opts <- ferx_model_get_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 999", opts)))
})
test_that("ferx_model_set_section() returns ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(
    ferx_model_set_section(ferx_model(model = path), "fit_options", "  method = foce")
  )
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})
test_that("pipe chain: double ferx_model_set_section() applies both changes", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(model = path) |>
    ferx_model_set_section("fit_options",  c("  method = focei", "  maxiter = 999")) |>
    ferx_model_set_section("error_model",  "  DV ~ additive(ADD_ERR)")

  expect_true(any(grepl("maxiter = 999", ferx_model_get_section(path, "fit_options"))))
  expect_true(any(grepl("additive",      ferx_model_get_section(path, "error_model"))))
})
test_that("pipe chain: ferx_model() |> ferx_model_set_section() |> print() shows changed content", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model(model = path) |>
      ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 123")) |>
      print()
  )
  # print.ferx_model shows structural summary, not raw file lines,
  # so verify the object identity by checking the file on disk instead
  opts <- ferx_model_get_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 123", opts)))
})
test_that("ferx_model_set_section() on a ferx_model pointing at a package file does NOT modify the bundled file", {
  ex <- ferx_example("warfarin")
  before <- readLines(ex$model, warn = FALSE)
  suppressMessages(
    ferx_model(ex$model, data = ex$data) |>
      ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 7"))
  )
  after <- readLines(ex$model, warn = FALSE)
  expect_identical(before, after)
})
test_that("ferx_model_set_section() on a ferx_model pointing at a package file redirects $model to tempdir", {
  ex <- ferx_example("warfarin")
  m  <- suppressMessages(
    ferx_model_set_section(
      ferx_model(ex$model, data = ex$data),
      "fit_options", c("  method = focei", "  maxiter = 7")
    )
  )
  expect_true(startsWith(normalizePath(m$model), normalizePath(tempdir())))
  expect_true(file.exists(m$model))
  opts <- ferx_model_get_section(m$model, "fit_options")
  expect_true(any(grepl("maxiter = 7", opts)))
})
test_that("ferx_model_set_section() on a user-owned ferx_model edits in place (no copy)", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m <- ferx_model_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 12")
  )
  expect_equal(m$model, path)
  opts <- ferx_model_get_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 12", opts)))
})
test_that("ferx_model_set_section() emits a copy-on-write message for package files", {
  ex <- ferx_example("warfarin")
  expect_message(
    ferx_model_set_section(
      ferx_model(ex$model, data = ex$data),
      "fit_options", "  method = focei"
    ),
    regexp = "copying read-only package model"
  )
})
