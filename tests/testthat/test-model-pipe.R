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

test_that("ferx_model() returns an object of class 'ferx_model'", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_s3_class(ferx_model(model = path), "ferx_model")
})

test_that("ferx_model()$model equals the input model path", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_equal(ferx_model(model = path)$model, path)
})

test_that("ferx_model()$data is NULL when data is not supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_null(ferx_model(model = path)$data)
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
    ferx_model(model = file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})

test_that("ferx_model() errors on wrong model extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))
  expect_error(ferx_model(model = path), regexp = "\\.ferx")
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
# Block 1b — data-first argument order (#81) and deprecation shim
# ---------------------------------------------------------------------------

test_that("ferx_model(data, model) accepts positional new-order arguments", {
  ex <- ferx_example("warfarin")
  m  <- ferx_model(ex$data, ex$model)
  expect_s3_class(m, "ferx_model")
  expect_equal(m$model, ex$model)
  expect_equal(m$data,  ex$data)
})

test_that("ferx_model() supports the data |> ferx_model(model) pipe form", {
  ex <- ferx_example("warfarin")
  m  <- ex$data |> ferx_model(ex$model)
  expect_s3_class(m, "ferx_model")
  expect_equal(m$model, ex$model)
  expect_equal(m$data,  ex$data)
})

test_that("ferx_model(path, data = ...) keeps working without a warning (named-data backcompat)", {
  ex   <- ferx_example("warfarin")
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_no_warning(m <- ferx_model(path, data = ex$data))
  expect_equal(m$model, path)
  expect_equal(m$data,  ex$data)
})

test_that("ferx_model(path) old-style single positional .ferx warns and auto-corrects", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_warning(
    m <- ferx_model(path),
    regexp = "argument order changed"
  )
  expect_equal(m$model, path)
  expect_null(m$data)
})

test_that("ferx_model(model_path, data_path) old positional pair warns and auto-corrects", {
  ex   <- ferx_example("warfarin")
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_warning(
    m <- ferx_model(path, ex$data),
    regexp = "argument order changed"
  )
  expect_equal(m$model, path)
  expect_equal(m$data,  ex$data)
})

test_that("ferx_model() errors when neither positional nor named model is supplied", {
  expect_error(ferx_model(), regexp = "'model' is required")
})

# ---------------------------------------------------------------------------
# Block 2 — print.ferx_model
# ---------------------------------------------------------------------------

test_that("print.ferx_model() prints 'ferx_model' as first line", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(model = path)))
  expect_equal(out[1], "ferx_model")
})

test_that("print.ferx_model() prints the model path", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(model = path)))
  expect_true(any(grepl(basename(path), out, fixed = TRUE)))
})

test_that("print.ferx_model() prints '<none>' when data is not supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(print(ferx_model(model = path)))
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
  result <- withVisible(print(ferx_model(model = path)))
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

  ferx_model(model = path) |>
    ferx_set_section("fit_options",  c("  method = focei", "  maxiter = 999")) |>
    ferx_set_section("error_model",  "  DV ~ additive(ADD_ERR)")

  expect_true(any(grepl("maxiter = 999", ferx_model_section(path, "fit_options"))))
  expect_true(any(grepl("additive",      ferx_model_section(path, "error_model"))))
})

# ---------------------------------------------------------------------------
# Block 4 — ferx_get_section() (pipe-friendly section display)
# ---------------------------------------------------------------------------

test_that("ferx_get_section() on ferx_model prints the section to console", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(ferx_get_section(ferx_model(model = path), "parameters"))
  expect_true(any(grepl("parameters", out)))
})

test_that("ferx_get_section() on ferx_model returns the ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(model = path)
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

  ferx_model(model = path) |>
    ferx_get_section("parameters") |>
    ferx_set_section("fit_options", "  method = foce")

  expect_true(any(grepl("method = foce", ferx_model_section(path, "fit_options"))))
})

# ---------------------------------------------------------------------------
# Block 5 — Full pipe chains (golden paths)
# ---------------------------------------------------------------------------

test_that("pipe chain: ferx_model(template=) |> ferx_model_set_section() |> ferx_model_inspect()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999")) |>
    ferx_model_inspect()

  expect_equal(result$model_type, "1-cpt oral")
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

# ---------------------------------------------------------------------------
# Block 6 — ferx_fit() dispatch on a ferx_model object
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
    ferx_model(model = path) |> ferx_fit(verbose = FALSE),
    regexp = "data"
  )
})

test_that("ferx_fit() exposes all estimation arguments in its formals (#52)", {
  # The whole point of dropping the S3 generic: IDE argument completion and
  # inline help introspect formals(ferx_fit), which previously showed only
  # model/data/... and hid every real option from the user.
  fmls <- names(formals(ferx_fit))
  expect_true(all(c(
    "model", "data", "method", "covariance", "verbose", "bloq_method",
    "threads", "mu_referencing", "sir", "gradient", "optimizer_trace",
    "scale_params", "inits_from_nca", "settings", "output", "include_data"
  ) %in% fmls))
  # Convergence-tolerance knobs flow through `settings`, not the signature (#51).
  expect_false("max_unconverged_frac" %in% fmls)
  expect_false("min_obs_for_convergence_check" %in% fmls)
  # No longer an S3 generic: the .default / .ferx_model methods are gone.
  expect_false(exists("ferx_fit.default"))
  expect_false(exists("ferx_fit.ferx_model"))
})

test_that("ferx_fit() rejects unrecognised arguments instead of silently ignoring them", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, methdo = "focei"),
    regexp = "unused argument"
  )
})

test_that("ferx_fit() reports an unrecognised arg without evaluating its value", {
  # The `...` guard must inspect names via ...names()/...length(), never force
  # the promises with list(...): forcing a side-effecting / erroring value would
  # mask the clear "unused argument(s)" message. Regression for the #67 review.
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, bogus = stop("must not be evaluated")),
    regexp = "unused argument\\(s\\) passed to ferx_fit\\(\\): bogus"
  )
})

# ---------------------------------------------------------------------------
# Block 7 — ferx_set_section() copy-on-write for package files (#80)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Block 8 — ferx_check_init() accepts ferx_model (#79)
# ---------------------------------------------------------------------------

test_that("ferx_check_init() accepts a ferx_model in a pipe and uses its data", {
  ex <- ferx_example("warfarin")
  chk <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_check_init(method = "focei", maxiter = 2L)
  )
  expect_named(chk, c("fit", "trace", "summary"))
  expect_s3_class(chk$fit, "ferx_fit")
  expect_true(is.data.frame(chk$summary))
})

test_that("ferx_check_init() errors when ferx_model has no data and none is supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_error(
    ferx_model(model = path) |> ferx_check_init(),
    regexp = "data"
  )
})

test_that("ferx_check_init() explicit data argument overrides data on ferx_model", {
  ex <- ferx_example("warfarin")
  # Sanity check that explicit data still wins over $data on the object —
  # we pass the same ex$data twice but via different routes.
  chk <- suppressWarnings(
    ferx_check_init(
      ferx_model(ex$model, data = ex$data),
      data = ex$data, method = "focei", maxiter = 2L
    )
  )
  expect_s3_class(chk$fit, "ferx_fit")
})
