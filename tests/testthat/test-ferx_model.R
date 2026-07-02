# ferx_model_section
# ---------------------------------------------------------------------------



test_that("ferx_model(print = TRUE) prints to console and returns NULL invisibly", {
  result <- withVisible(ferx_model(print = TRUE))
  expect_null(result$value)
  expect_false(result$visible)
})
test_that("ferx_model(print = TRUE) works for each template", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    expect_no_error(ferx_model(template = tmpl, print = TRUE))
  }
})
test_that("ferx_model(template=) creates a file containing the required sections", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)

  expect_true(file.exists(path))
  content <- readLines(path)
  expect_true(any(grepl("^\\[parameters\\]", content)))
  expect_true(any(grepl("^\\[fit_options\\]", content)))
})
test_that("ferx_model(template=) errors when file exists and overwrite = FALSE", {
  path <- tempfile(fileext = ".ferx")
  writeLines("existing", path)
  on.exit(unlink(path))

  expect_error(
    ferx_model(template = "1cpt_oral", path = path, edit = FALSE),
    regexp = "already exists"
  )
})
test_that("ferx_model(template=) overwrites when overwrite = TRUE", {
  path <- tempfile(fileext = ".ferx")
  writeLines("old content", path)
  on.exit(unlink(path))

  expect_no_error(
    ferx_model(template = "1cpt_oral", path = path, overwrite = TRUE, edit = FALSE)
  )
  expect_true(any(grepl("\\[parameters\\]", readLines(path))))
})
test_that("ferx_model(template=) errors on unknown template", {
  expect_error(
    ferx_model(print = TRUE, template = "3cpt_magical"),
    regexp = "Unknown template"
  )
})
test_that("ferx_model(template=) errors when path is NULL and print = FALSE", {
  expect_error(
    ferx_model(template = "1cpt_oral", print = FALSE),
    regexp = "'path' is required"
  )
})
test_that("ferx_model(template=) errors when path has wrong extension", {
  expect_error(
    ferx_model(template = "1cpt_oral", path = tempfile(fileext = ".txt"),
               edit = FALSE),
    regexp = "\\.ferx"
  )
})
test_that("ferx_model(template=) returns a ferx_model object pointing at the file", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  m <- ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  expect_s3_class(m, "ferx_model")
  expect_equal(m$model, path)
  expect_null(m$data)
})
test_that("ferx_model(template=) each template contains all required sections", {
  required <- c("parameters", "individual_parameters", "structural_model",
                "error_model", "fit_options")
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model(template = tmpl, path = path, edit = FALSE)
    content <- readLines(path)
    for (sec in required) {
      expect_true(
        any(grepl(paste0("^\\[", sec, "\\]"), content)),
        info = paste("template", tmpl, "missing section", sec)
      )
    }
  }
})
test_that("ferx_model(template=) templates do not emit the deprecated [initial_values] block", {
  # Initial values now live inline in [parameters] — re-emitting them in a
  # separate block would be redundant boilerplate that the parser silently
  # ignores. Guards against regressing the issue-16 fix.
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model(template = tmpl, path = path, edit = FALSE)
    content <- readLines(path)
    expect_false(
      any(grepl("^\\[initial_values\\]", content)),
      info = paste("template", tmpl, "still emits [initial_values]")
    )
  }
})
test_that("ferx_model(template=) print=TRUE output matches what would be written to file", {
  for (tmpl in c("1cpt_oral", "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode")) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    ferx_model(template = tmpl, path = path, edit = FALSE)
    file_lines <- readLines(path)

    printed <- capture.output(ferx_model(template = tmpl, print = TRUE))
    # capture.output adds a trailing "" for the final cat("\n"); drop it
    printed <- printed[nchar(printed) > 0 | seq_along(printed) < length(printed)]

    expect_equal(printed[printed != ""], file_lines[file_lines != ""],
                 info = paste("template", tmpl, "print vs file mismatch"))
  }
})
test_that("ferx_model(template=) pipe chain: scaffold |> set_section", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
    ferx_model_set_section("fit_options", c("  method = focei", "  maxiter = 999"))

  expect_equal(ferx_model_section(path, "fit_options"), c("  method = focei", "  maxiter = 999"))
})

# write_pipe_test_model()/modifying_editor() come from helper-model-pipe.R

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
