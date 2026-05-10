# Tests for the ferx_model S3 pipe object (issue #47).
#
# Every test in this file is guarded by skip_if_no_ferx_model(), which checks
# that ferx_model() is exported. The tests are inert until #47 is merged and
# the package is reinstalled — at that point they run automatically with no
# manual changes needed.

skip_if_no_ferx_model <- function() {
  skip_if_not(
    existsMethod("ferx_model") || exists("ferx_model", envir = asNamespace("ferx"), inherits = FALSE),
    "ferx_model() not yet implemented (waiting for #47)"
  )
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_pipe_test_model <- function() {
  path <- tempfile(fileext = ".ferx")
  ferx_model_new(path, edit = FALSE)
  path
}

null_editor <- function(...) invisible(NULL)

modifying_editor <- function(p, ...) {
  lines <- readLines(p)
  writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
}

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------

test_that("ferx_model() returns an object of class 'ferx_model'", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_s3_class(ferx_model(path), "ferx_model")
})

test_that("ferx_model()$path equals the input path", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_equal(ferx_model(path)$path, normalizePath(path))
})

test_that("ferx_model() errors on missing file", {
  skip_if_no_ferx_model()
  expect_error(
    ferx_model(file.path(tempdir(), "no_such.ferx")),
    regexp = "not found|does not exist"
  )
})

test_that("ferx_model() errors on wrong extension", {
  skip_if_no_ferx_model()
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))
  expect_error(ferx_model(path), regexp = "\\.ferx")
})

# ---------------------------------------------------------------------------
# Block 2 — print.ferx_model
# ---------------------------------------------------------------------------

test_that("print.ferx_model() produces the same output as ferx_model_show()", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_equal(
    capture.output(print(ferx_model(path))),
    capture.output(ferx_model_show(path))
  )
})

test_that("print.ferx_model() returns the ferx_model object invisibly", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(print(ferx_model(path)))
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

# ---------------------------------------------------------------------------
# Block 3 — ferx_model_set_section dispatch on ferx_model
# ---------------------------------------------------------------------------

test_that("ferx_model_set_section() on ferx_model returns a ferx_model", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- ferx_model_set_section(
    ferx_model(path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  expect_s3_class(result, "ferx_model")
})

test_that("ferx_model_set_section() on ferx_model returns same path", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m <- ferx_model(path)
  result <- ferx_model_set_section(m, "fit_options", c("  method = focei"))
  expect_equal(result$path, m$path)
})

test_that("ferx_model_set_section() on ferx_model writes change to disk", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_model_set_section(
    ferx_model(path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 999", opts)))
})

test_that("ferx_model_set_section() returns ferx_model invisibly", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(
    ferx_model_set_section(ferx_model(path), "fit_options", c("  method = foce"))
  )
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

test_that("pipe chain: double ferx_model_set_section() applies both changes", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(path) |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_set_section("error_model", c("  DV ~ additive(ADD_ERR)"))

  expect_true(any(grepl("maxiter = 999", ferx_model_section(path, "fit_options"))))
  expect_true(any(grepl("additive",      ferx_model_section(path, "error_model"))))
})

# ---------------------------------------------------------------------------
# Block 4 — ferx_model_edit dispatch on ferx_model
# ---------------------------------------------------------------------------

test_that("ferx_model_edit() on ferx_model returns a ferx_model invisibly", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(
    ferx_model_edit(ferx_model(path), .editor = null_editor)
  )
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

test_that("ferx_model_edit() save_as via ferx_model captures post-edit content", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  ferx_model_edit(ferx_model(path), save_as = dest, .editor = modifying_editor)

  expect_true(file.exists(dest))
  expect_true(any(grepl("TVV", readLines(dest))))
})

# ---------------------------------------------------------------------------
# Block 5 — Full pipe chains (golden paths)
# ---------------------------------------------------------------------------

test_that("pipe chain: ferx_model_new |> ferx_model_set_section |> ferx_model_inspect", {
  skip_if_no_ferx_model()
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model_new(path, edit = FALSE) |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_inspect()

  expect_equal(result$method %||% "focei", "focei")
})

test_that("pipe chain: ferx_model_new |> ferx_model_edit |> ferx_model_show", {
  skip_if_no_ferx_model()
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model_new(path, edit = FALSE) |>
      ferx_model_edit(.editor = modifying_editor) |>
      ferx_model_show()
  )
  expect_true(any(grepl("TVV", out)))
})

test_that("pipe chain: ferx_model() |> ferx_model_set_section |> ferx_model_show", {
  skip_if_no_ferx_model()
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model(path) |>
      ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 123")) |>
      ferx_model_show()
  )
  expect_true(any(grepl("maxiter = 123", out)))
})

# ---------------------------------------------------------------------------
# Block 6 — ferx_fit dispatch on ferx_model (if #47 includes it)
# ---------------------------------------------------------------------------

test_that("ferx_fit() dispatches on ferx_model and returns a ferx_fit", {
  skip_if_no_ferx_model()
  skip("enable once ferx_fit.ferx_model dispatch is confirmed in #47")
  ex     <- ferx_example("warfarin")
  result <- ferx_model(ex$model) |> ferx_fit(data = ex$data)
  expect_s3_class(result, "ferx_fit")
})
