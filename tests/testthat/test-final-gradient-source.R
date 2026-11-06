# `fit$final_gradient_source` (#684, ferx-core #997): what kind of gradient
# `fit$final_gradient` is, on a fresh fit and through a `.fitrx` bundle.

gradient_fit <- local({
  fits <- list()
  function(optimizer = "auto", gradient = "auto", report = TRUE) {
    key <- paste(c(optimizer, gradient, report), collapse = "/")
    if (is.null(fits[[key]])) {
      ex <- ferx_example("warfarin")
      settings <- list(maxiter = 30L, optimizer = optimizer,
                       report_final_gradient = report)
      fits[[key]] <<- suppressWarnings(ferx_fit(
        ex$model, ex$data, method = "focei", verbose = FALSE,
        covariance = FALSE, gradient = gradient, settings = settings
      ))
    }
    fits[[key]]
  }
})

test_that("each kind of gradient is labelled with its source", {
  skip_on_cran()
  # A gradient-based optimizer reports its own.
  fit <- gradient_fit()
  expect_identical(fit$final_gradient_source, "optimizer")
  expect_true(length(fit$final_gradient) > 0L)

  # BOBYQA has none, so one is computed after the fit: exact where the model
  # can use the analytic outer gradient, a central difference where it cannot.
  fit <- gradient_fit("bobyqa")
  expect_identical(fit$final_gradient_source, "analytic")
  expect_true(length(fit$final_gradient) > 0L)

  fit <- gradient_fit("bobyqa", gradient = "fd")
  expect_identical(fit$final_gradient_source, "finite_difference")
  expect_true(length(fit$final_gradient) > 0L)
})

test_that("no gradient means no source", {
  skip_on_cran()
  fit <- gradient_fit("bobyqa", report = FALSE)
  expect_null(fit$final_gradient)
  expect_null(fit$final_gradient_source)
})

test_that("the gradient and its source survive a .fitrx round trip", {
  skip_on_cran()
  for (fit in list(gradient_fit(), gradient_fit("bobyqa"),
                   gradient_fit("bobyqa", gradient = "fd"))) {
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    loaded <- ferx_load_fit(path)
    expect_identical(loaded$final_gradient_source, fit$final_gradient_source)
    expect_equal(loaded$final_gradient, fit$final_gradient, tolerance = 1e-12)
  }

  fit <- gradient_fit("bobyqa", report = FALSE)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_null(loaded$final_gradient)
  expect_null(loaded$final_gradient_source)
})

test_that("the bundle holds the gradient and its source in the engine's keys", {
  skip_on_cran()
  fit <- gradient_fit("bobyqa", gradient = "fd")
  wire <- fitrx_fit_json(fit$bundle)
  expect_identical(wire$final_gradient_source, "finite_difference")
  expect_equal(as.numeric(unlist(wire$final_gradient)), fit$final_gradient,
               tolerance = 1e-12)
})

test_that("a one-component gradient reads back from its array", {
  skip_on_cran()
  fit <- fitrx_edit_fit(gradient_fit(), function(w) {
    w$final_gradient <- list(0.25)
    w
  })
  expect_type(fitrx_fit_json(fit$bundle)$final_gradient, "list")
  expect_equal(fit$final_gradient, 0.25)
})

test_that("a non-finite gradient does not cost the bundle its engine load", {
  skip_on_cran()
  # JSON has no NaN: it would be written as a null the engine's `Vec<f64>`
  # refuses. The gradient goes, its label with it, and the bundle still loads.
  fit <- gradient_fit("bobyqa", gradient = "fd")
  expect_false(is.null(fit[["final_gradient"]]))
  bytes <- ferx:::ferx_rust_test_nan_gradient_bundle(fit$bundle)
  wire <- fitrx_fit_json(bytes)
  expect_null(wire$final_gradient)
  expect_null(wire$final_gradient_source)
  loaded <- ferx:::.ferx_fit_with_bundle(fit, bytes)
  expect_null(loaded[["final_gradient"]])
  expect_null(loaded$final_gradient_source)
})

test_that("a bundle saved before the fields existed loads with neither", {
  skip_on_cran()
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(gradient_fit(), path)
  old <- fitrx_edit_bundle(path, function(wire) {
    wire$final_gradient <- NULL
    wire$final_gradient_source <- NULL
    wire
  })
  on.exit(unlink(old), add = TRUE)
  loaded <- ferx_load_fit(old)
  expect_null(loaded$final_gradient)
  expect_null(loaded$final_gradient_source)
})
