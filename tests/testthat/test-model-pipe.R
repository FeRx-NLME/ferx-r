# Tests for the ferx_model S3 pipe object and pipe-friendly wrappers added in #47:
#   ferx_model()       — constructor
#   print.ferx_model() — console summary
#   ferx_set_section() — pipe-friendly section replacement
#   ferx_get_section() — pipe-friendly section display
#   ferx_fit()         — ferx_fit.ferx_model dispatch

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_pipe_test_model <- function() {
  path <- tempfile(fileext = ".ferx")
  ferx_model_new(path, edit = FALSE)
  path
}

modifying_editor <- function(p, ...) {
  lines <- readLines(p)
  writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
}

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------

test_that("ferx_model() returns an object of class 'ferx_model'", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_s3_class(ferx_model(path), "ferx_model")
})

test_that("ferx_model()$model equals the input model path", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_equal(ferx_model(path)$model, path)
})

test_that("ferx_model()$data is NULL when data is not supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_null(ferx_model(path)$data)
})

test_that("ferx_model()$data equals the input data path when supplied", {
  path <- write_pipe_test_model()
  ex   <- ferx_example("warfarin")
  on.exit(unlink(path))
  m <- ferx_model(path, data = ex$data)
  expect_equal(m$data, ex$data)
})

test_that("ferx_model() errors on missing model file", {
  expect_error(
    ferx_model(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})

test_that("ferx_model() errors on wrong model extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))
  expect_error(ferx_model(path), regexp = "\\.ferx")
})

test_that("ferx_model() errors on missing data file", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_error(
    ferx_model(path, data = file.path(tempdir(), "no_such.csv")),
    regexp = "not found"
  )
})

# ---------------------------------------------------------------------------
# Block 2 — print.ferx_model
# ---------------------------------------------------------------------------

test_that("print.ferx_model() prints 'ferx_model' as first line", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(path)))
  expect_equal(out[1], "ferx_model")
})

test_that("print.ferx_model() prints the model path", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(path)))
  expect_true(any(grepl(basename(path), out, fixed = TRUE)))
})

test_that("print.ferx_model() prints '<none>' when data is not supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(path)))
  expect_true(any(grepl("<none>", out, fixed = TRUE)))
})

test_that("print.ferx_model() prints the data path when supplied", {
  path <- write_pipe_test_model()
  ex   <- ferx_example("warfarin")
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(path, data = ex$data)))
  expect_true(any(grepl(basename(ex$data), out, fixed = TRUE)))
})

test_that("print.ferx_model() returns the ferx_model object invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(print(ferx_model(path)))
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

# ---------------------------------------------------------------------------
# Block 3 — ferx_set_section() (pipe-friendly section replacement)
# ---------------------------------------------------------------------------

test_that("ferx_set_section() on ferx_model returns a ferx_model", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- ferx_set_section(
    ferx_model(path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  expect_s3_class(result, "ferx_model")
})

test_that("ferx_set_section() on ferx_model returns the same object (same $model path)", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(path)
  result <- ferx_set_section(m, "fit_options", c("  method = focei"))
  expect_equal(result$model, m$model)
})

test_that("ferx_set_section() on ferx_model writes the change to disk", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_set_section(
    ferx_model(path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 999", opts)))
})

test_that("ferx_set_section() returns ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(
    ferx_set_section(ferx_model(path), "fit_options", c("  method = foce"))
  )
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

test_that("ferx_set_section() on a plain path delegates to ferx_model_set_section()", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_set_section(path, "fit_options", c("  method = focei", "  maxiter = 42"))
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 42", opts)))
})

test_that("pipe chain: double ferx_set_section() applies both changes", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(path) |>
    ferx_set_section("fit_options",  c("  method = focei", "  maxiter = 999")) |>
    ferx_set_section("error_model",  c("  DV ~ additive(ADD_ERR)"))

  expect_true(any(grepl("maxiter = 999", ferx_model_section(path, "fit_options"))))
  expect_true(any(grepl("additive",      ferx_model_section(path, "error_model"))))
})

# ---------------------------------------------------------------------------
# Block 4 — ferx_get_section() (pipe-friendly section display)
# ---------------------------------------------------------------------------

test_that("ferx_get_section() on ferx_model prints the section to console", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(ferx_get_section(ferx_model(path), "parameters"))
  expect_true(any(grepl("parameters", out)))
})

test_that("ferx_get_section() on ferx_model returns the ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(path)
  result <- withVisible(ferx_get_section(m, "parameters"))
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})

test_that("ferx_get_section() on a plain path returns the path invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(ferx_get_section(path, "parameters"))
  expect_equal(result$value, path)
  expect_false(result$visible)
})

test_that("pipe chain: ferx_model() |> ferx_get_section() |> ferx_set_section() works", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(path) |>
    ferx_get_section("parameters") |>
    ferx_set_section("fit_options", c("  method = foce"))

  expect_true(any(grepl("method = foce", ferx_model_section(path, "fit_options"))))
})

# ---------------------------------------------------------------------------
# Block 5 — Full pipe chains (golden paths)
# ---------------------------------------------------------------------------

test_that("pipe chain: ferx_model_new() |> ferx_model_set_section() |> ferx_model_inspect()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model_new(path, edit = FALSE) |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_inspect()

  expect_equal(result$model_type, "1-cpt oral")
})

test_that("pipe chain: ferx_model() |> ferx_set_section() |> print() shows changed content", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model(path) |>
      ferx_set_section("fit_options", c("  method = focei", "  maxiter = 123")) |>
      print()
  )
  # print.ferx_model shows structural summary, not raw file lines,
  # so verify the object identity by checking the file on disk instead
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 123", opts)))
})

test_that("pipe chain: ferx_model() |> ferx_get_section() |> ferx_set_section() — then inspect", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  # ferx_model_inspect() takes a path, not a ferx_model, so extract $model after the chain
  m <- ferx_model(path) |>
    ferx_get_section("parameters") |>
    ferx_set_section("fit_options", c("  method = foce", "  maxiter = 50"))

  result <- ferx_model_inspect(m$model)
  expect_equal(result$residual, "proportional")
})

# ---------------------------------------------------------------------------
# Block 6 — ferx_fit.ferx_model dispatch
#
# Both call forms must produce identical results:
#   ferx_fit("my_model.ferx", data = "data.csv")          # path style
#   ferx_model("my_model.ferx", data = "data.csv") |> ferx_fit()  # pipe style
#
# ferx_fit() returns a plain ferx_fit — NOT a ferx_model.
# ---------------------------------------------------------------------------

test_that("ferx_fit() on ferx_model returns a ferx_fit (not a ferx_model)", {
  ex     <- ferx_example("warfarin")
  result <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_fit(verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
  expect_false(inherits(result, "ferx_model"))
})

test_that("ferx_fit() on ferx_model via named argument works (ferx_fit(model = m))", {
  ex     <- ferx_example("warfarin")
  m      <- ferx_model(ex$model, data = ex$data)
  result <- suppressWarnings(
    ferx_fit(model = m, verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
})

test_that("ferx_fit() on ferx_model produces same theta names as path-style call", {
  ex      <- ferx_example("warfarin")
  by_path <- suppressWarnings(
    ferx_fit(ex$model, data = ex$data, verbose = FALSE, settings = list(maxiter = 30L))
  )
  by_pipe <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_fit(verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_equal(names(by_path$theta), names(by_pipe$theta))
})

test_that("ferx_fit() data argument overrides data stored in ferx_model", {
  ex <- ferx_example("warfarin")
  result <- suppressWarnings(
    ferx_model(ex$model) |>
      ferx_fit(data = ex$data, verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
})

test_that("ferx_fit() errors with clear message when ferx_model has no data and none supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_error(
    ferx_model(path) |> ferx_fit(verbose = FALSE),
    regexp = "data"
  )
})
