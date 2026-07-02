
# ---- header from test-frem.R ----
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Write a warfarin CSV with synthetic WT and AGE covariates.
write_warfarin_with_covariates <- function(dir) {
  ex <- ferx_example("warfarin")
  dat <- read.csv(ex$data)

  # Per-subject covariates (10 subjects).
  cov_tbl <- data.frame(
    ID  = 1:10,
    WT  = c(70, 80, 65, 90, 55, 75, 85, 60, 72, 68),
    AGE = c(35, 45, 28, 55, 22, 40, 50, 30, 38, 33)
  )
  dat <- merge(dat, cov_tbl, by = "ID", sort = FALSE)
  path <- file.path(dir, "warfarin_cov.csv")
  write.csv(dat, path, row.names = FALSE, quote = FALSE)
  path
}

# Write the base warfarin model. FREM takes its covariates from the model's
# [covariates] block, so the base model declares WT and AGE there.
write_warfarin_model <- function(dir) {
  path <- file.path(dir, "warfarin.ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.2, 0.001, 10.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    "  theta TVKA(1.5, 0.01, 50.0)",
    "",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V  ~ 0.04",
    "  omega ETA_KA ~ 0.30",
    "",
    "  sigma PROP_ERR ~ 0.02",
    "",
    "[covariates]",
    "  WT  continuous",
    "  AGE continuous",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV  * exp(ETA_V)",
    "  KA = TVKA * exp(ETA_KA)",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)",
    "",
    "[fit_options]",
    "  method = focei"
  ), path)
  path
}

# Same model with no [covariates] block — FREM has nothing to transform and
# should error.
write_warfarin_model_no_cov_block <- function(dir) {
  path <- file.path(dir, "warfarin_nocov.ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.2, 0.001, 10.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    "  theta TVKA(1.5, 0.01, 50.0)",
    "",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V  ~ 0.04",
    "  omega ETA_KA ~ 0.30",
    "",
    "  sigma PROP_ERR ~ 0.02",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV  * exp(ETA_V)",
    "  KA = TVKA * exp(ETA_KA)",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)",
    "",
    "[fit_options]",
    "  method = focei"
  ), path)
  path
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("ferx_model_to_frem() errors on missing model path", {
  expect_error(ferx_model_to_frem(model = "", data = "x.csv", covariates = "WT"),
               "non-empty path")
})
test_that("ferx_model_to_frem() errors on non-existent model file", {
  expect_error(ferx_model_to_frem(model = "/no/such/file.ferx",
                            data = "/no/such/data.csv",
                            covariates = "WT"),
               "Model file not found")
})
test_that("ferx_model_to_frem() errors on missing data argument", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)
  on.exit(unlink(model_path))

  expect_error(ferx_model_to_frem(model = model_path, covariates = "WT"),
               "data")
})
test_that("ferx_model_to_frem() errors when the model has no [covariates] block", {
  tmp <- tempdir()
  model_path <- write_warfarin_model_no_cov_block(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)
  on.exit(unlink(c(model_path, data_path)))

  expect_error(
    ferx_model_to_frem(model = model_path, data = data_path),
    "\\[covariates\\] block"
  )
})
test_that("ferx_model_to_frem() errors when the covariates filter names an undeclared covariate", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)        # declares WT, AGE
  data_path  <- write_warfarin_with_covariates(tmp)
  on.exit(unlink(c(model_path, data_path)))

  expect_error(
    ferx_model_to_frem(model = model_path, data = data_path, covariates = "SEX"),
    "not declared"
  )
})
test_that("ferx_model_to_frem() errors when `fit` is not a ferx_fit object", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)
  on.exit(unlink(c(model_path, data_path)))

  expect_error(
    ferx_model_to_frem(model = model_path, data = data_path, fit = list(theta = 1)),
    "ferx_fit"
  )
})
# ---------------------------------------------------------------------------
# End-to-end FREM transformation
# ---------------------------------------------------------------------------
test_that("ferx_model_to_frem() returns a ferx_model referencing the generated files", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  result <- ferx_model_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  # It is a regular ferx_model pointing at the generated model + data files.
  expect_s3_class(result, "ferx_model")
  expect_true(file.exists(result$model))
  expect_true(file.exists(result$data))
})
test_that("ferx_model_to_frem() uses all declared covariates when `covariates` omitted", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)   # declares WT, AGE
  data_path  <- write_warfarin_with_covariates(tmp)

  # `covariates` omitted → all covariates from the model's [covariates] block.
  result <- ferx_model_to_frem(model = model_path, data = data_path, output_dir = tmp)

  # Both WT (FREMTYPE 100) and AGE (FREMTYPE 200) pseudo-obs are emitted.
  frem_data <- read.csv(result$data)
  expect_equal(sort(unique(frem_data$FREMTYPE)), c(0, 100, 200))
})
test_that("ferx_model_to_frem() `covariates` filters to a subset of the block", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)   # declares WT, AGE
  data_path  <- write_warfarin_with_covariates(tmp)

  # Filter to WT only — AGE is declared but excluded from the FREM model.
  result <- ferx_model_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = "WT",
    output_dir = tmp
  )

  # Only WT (FREMTYPE 100) pseudo-obs; AGE (200) is excluded.
  frem_data <- read.csv(result$data)
  expect_equal(sort(unique(frem_data$FREMTYPE)), c(0, 100))
})
test_that("FREM dataset has correct row count", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  result <- ferx_model_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  frem_data <- read.csv(result$data)

  # Original: 10 subjects x (1 dose + 11 obs) = 120 rows.
  # FREM adds: 10 subjects x 2 covariates = 20 pseudo-obs.
  expect_equal(nrow(frem_data), 140L)

  # FREMTYPE column exists with correct values.
  expect_true("FREMTYPE" %in% names(frem_data))
  expect_equal(sum(frem_data$FREMTYPE == 100), 10L)  # WT pseudo-obs
  expect_equal(sum(frem_data$FREMTYPE == 200), 10L)  # AGE pseudo-obs
  expect_equal(sum(frem_data$FREMTYPE == 0),  120L)  # original rows
})
test_that("FREM model can be fitted", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  result <- ferx_model_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  # Pass the ferx_model returned by ferx_model_to_frem() straight to ferx_fit().
  fit <- ferx_fit(
    result,
    method     = "focei",
    verbose    = FALSE,
    covariance = FALSE,
    settings   = list(maxiter = 3L)
  )

  expect_s3_class(fit, "ferx_fit")
  expect_true(is.finite(fit$ofv))

  # Omega should be 5x5.
  expect_equal(nrow(fit$omega), 5L)
  expect_equal(ncol(fit$omega), 5L)

  # All omega diagonals positive.
  expect_true(all(diag(fit$omega) > 0))
})

test_that("ferx_model_to_frem() `fit` seeds the generated model's PK theta/omega inits (#239)", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  base_fit <- ferx_fit(
    model_path, data_path,
    method     = "focei",
    verbose    = FALSE,
    covariance = FALSE,
    settings   = list(maxiter = 3L)
  )
  expect_s3_class(base_fit, "ferx_fit")

  result <- ferx_model_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp,
    fit        = base_fit
  )

  # The generated FREM model's declared theta inits are the fitted values,
  # not the base model file's original declared inits (0.2, 10.0, 1.5).
  frem_lines <- readLines(result$model)
  init_of <- function(name) {
    line <- grep(sprintf("^\\s*theta %s\\(", name), frem_lines, value = TRUE)
    expect_length(line, 1L)
    as.numeric(sub(sprintf("^\\s*theta %s\\(([^,]+),.*$", name), "\\1", line))
  }
  expect_equal(init_of("TVCL"), unname(base_fit$theta[["TVCL"]]))
  expect_equal(init_of("TVV"),  unname(base_fit$theta[["TVV"]]))
  expect_equal(init_of("TVKA"), unname(base_fit$theta[["TVKA"]]))

  # It still fits.
  fit <- ferx_fit(
    result,
    method     = "focei",
    verbose    = FALSE,
    covariance = FALSE,
    settings   = list(maxiter = 3L)
  )
  expect_s3_class(fit, "ferx_fit")
  expect_true(is.finite(fit$ofv))
})
