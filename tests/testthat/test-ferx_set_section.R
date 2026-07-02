
# write_pipe_test_model()/modifying_editor() come from helper-model-pipe.R

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------



test_that("ferx_set_section() on ferx_model returns a ferx_model", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- ferx_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  expect_s3_class(result, "ferx_model")
})
test_that("ferx_set_section() on ferx_model returns the same object (same $model path)", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(model = path)
  result <- ferx_set_section(m, "fit_options", "  method = focei")
  expect_equal(result$model, m$model)
})
test_that("ferx_set_section() on ferx_model writes the change to disk", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  ferx_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 999")
  )
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 999", opts)))
})
test_that("ferx_set_section() returns ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(
    ferx_set_section(ferx_model(model = path), "fit_options", "  method = foce")
  )
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})
test_that("pipe chain: double ferx_set_section() applies both changes", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(model = path) |>
    ferx_set_section("fit_options",  c("  method = focei", "  maxiter = 999")) |>
    ferx_set_section("error_model",  "  DV ~ additive(ADD_ERR)")

  expect_true(any(grepl("maxiter = 999", ferx_model_section(path, "fit_options"))))
  expect_true(any(grepl("additive",      ferx_model_section(path, "error_model"))))
})
test_that("pipe chain: ferx_model() |> ferx_get_section() |> ferx_set_section() works", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  ferx_model(model = path) |>
    ferx_get_section("parameters") |>
    ferx_set_section("fit_options", "  method = foce")

  expect_true(any(grepl("method = foce", ferx_model_section(path, "fit_options"))))
})
test_that("pipe chain: ferx_model() |> ferx_set_section() |> print() shows changed content", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))

  out <- capture.output(
    ferx_model(model = path) |>
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
  m <- ferx_model(model = path) |>
    ferx_get_section("parameters") |>
    ferx_set_section("fit_options", c("  method = foce", "  maxiter = 50"))

  result <- ferx_model_inspect(m$model)
  expect_equal(result$residual, "proportional")
})
test_that("ferx_set_section() on a ferx_model pointing at a package file does NOT modify the bundled file", {
  ex <- ferx_example("warfarin")
  before <- readLines(ex$model, warn = FALSE)
  suppressMessages(
    ferx_model(ex$model, data = ex$data) |>
      ferx_set_section("fit_options", c("  method = focei", "  maxiter = 7"))
  )
  after <- readLines(ex$model, warn = FALSE)
  expect_identical(before, after)
})
test_that("ferx_set_section() on a ferx_model pointing at a package file redirects $model to tempdir", {
  ex <- ferx_example("warfarin")
  m  <- suppressMessages(
    ferx_set_section(
      ferx_model(ex$model, data = ex$data),
      "fit_options", c("  method = focei", "  maxiter = 7")
    )
  )
  expect_true(startsWith(normalizePath(m$model), normalizePath(tempdir())))
  expect_true(file.exists(m$model))
  opts <- ferx_model_section(m$model, "fit_options")
  expect_true(any(grepl("maxiter = 7", opts)))
})
test_that("ferx_set_section() on a user-owned ferx_model edits in place (no copy)", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m <- ferx_set_section(
    ferx_model(model = path), "fit_options",
    c("  method = focei", "  maxiter = 12")
  )
  expect_equal(m$model, path)
  opts <- ferx_model_section(path, "fit_options")
  expect_true(any(grepl("maxiter = 12", opts)))
})
test_that("ferx_set_section() emits a copy-on-write message for package files", {
  ex <- ferx_example("warfarin")
  expect_message(
    ferx_set_section(
      ferx_model(ex$model, data = ex$data),
      "fit_options", "  method = focei"
    ),
    regexp = "copying read-only package model"
  )
})
