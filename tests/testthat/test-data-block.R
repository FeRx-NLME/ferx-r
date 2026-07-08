# Tests for the model-file `[data]` block fallback (#254): fit/model/simulate/
# predict resolve the dataset from the model file when `data` is omitted.

# Helper: write a tiny NONMEM CSV + a one-cpt IV model into a fresh directory.
# `with_block` controls whether the model declares a `[data]` block pointing at
# the sibling CSV (relative path, so the engine resolves it against the model
# file's directory).
make_model <- function(dir, with_block = TRUE, csv_name = "mydata.csv") {
  csv <- file.path(dir, csv_name)
  writeLines(
    c(
      "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV",
      "1,0,.,1,100,1,0,1",
      "1,0.25,3.3287,0,.,1,0,0",
      "1,0.5,2.9472,0,.,1,0,0",
      "1,1,2.8394,0,.,1,0,0",
      "1,2,2.2970,0,.,1,0,0",
      "2,0,.,1,100,1,0,1",
      "2,0.25,3.1,0,.,1,0,0",
      "2,0.5,2.7,0,.,1,0,0",
      "2,1,2.5,0,.,1,0,0",
      "2,2,2.0,0,.,1,0,0"
    ),
    csv
  )

  body <- c(
    "[parameters]",
    "  theta TVCL(4.0, 0.01, 100.0)",
    "  theta TVV(30.0, 1.0, 500.0)",
    "",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V  ~ 0.04",
    "",
    "  sigma PROP_ERR ~ 0.04 (sd)",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV  * exp(ETA_V)",
    "",
    "[structural_model]",
    "  pk one_cpt_iv(cl=CL, v=V)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)",
    "",
    "[fit_options]",
    "  method     = foce",
    "  maxiter    = 5",
    "  covariance = false"
  )
  if (with_block) {
    body <- c(body, "", "[data]", paste0("  path = ", csv_name))
  }
  model <- file.path(dir, "model.ferx")
  writeLines(body, model)
  list(model = model, data = csv)
}

test_that(".ferx_model_data_path() resolves a declared [data] block", {
  dir <- file.path(tempdir(), "ferx_data_block_yes")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- make_model(dir, with_block = TRUE)

  got <- .ferx_model_data_path(paths$model)
  expect_false(is.null(got))
  expect_true(file.exists(got))
  expect_identical(normalizePath(got), normalizePath(paths$data))
})

test_that(".ferx_model_data_path() returns NULL without a [data] block", {
  dir <- file.path(tempdir(), "ferx_data_block_no")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- make_model(dir, with_block = FALSE)

  expect_null(.ferx_model_data_path(paths$model))
})

test_that("ferx_model()$data is populated from the [data] block", {
  dir <- file.path(tempdir(), "ferx_data_block_model")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- make_model(dir, with_block = TRUE)

  m <- ferx_model(model = paths$model)
  expect_false(is.null(m$data))
  expect_identical(normalizePath(m$data), normalizePath(paths$data))
})

test_that("ferx_fit() without data errors when no [data] block is declared", {
  dir <- file.path(tempdir(), "ferx_data_block_fit_err")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- make_model(dir, with_block = FALSE)

  expect_error(ferx_fit(paths$model), "No data supplied")
})

test_that("ferx_fit() picks up the [data] block dataset end-to-end", {
  skip_on_cran()
  dir <- file.path(tempdir(), "ferx_data_block_fit_ok")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- make_model(dir, with_block = TRUE)

  fit <- ferx_fit(paths$model, method = "foce", covariance = FALSE)
  expect_s3_class(fit, "ferx_fit")
})
