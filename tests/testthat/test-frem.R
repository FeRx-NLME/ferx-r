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

# Write the base warfarin model.
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

test_that("ferx_to_frem() errors on missing model path", {
  expect_error(ferx_to_frem(model = "", data = "x.csv", covariates = "WT"),
               "non-empty path")
})

test_that("ferx_to_frem() errors on non-existent model file", {
  expect_error(ferx_to_frem(model = "/no/such/file.ferx",
                            data = "/no/such/data.csv",
                            covariates = "WT"),
               "Model file not found")
})

test_that("ferx_to_frem() errors on missing data argument", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)
  on.exit(unlink(model_path))

  expect_error(ferx_to_frem(model = model_path, covariates = "WT"),
               "data")
})

test_that("ferx_to_frem() errors on empty covariates", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)
  on.exit(unlink(c(model_path, data_path)))

  expect_error(ferx_to_frem(model = model_path, data = data_path,
                            covariates = character(0)),
               "covariates")
})

test_that("ferx_to_frem() errors when categorical is not subset of covariates", {
  tmp <- tempdir()
  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)
  on.exit(unlink(c(model_path, data_path)))

  expect_error(
    ferx_to_frem(model = model_path, data = data_path,
                 covariates = c("WT", "AGE"),
                 categorical = "SEX"),
    "not in `covariates`"
  )
})

# ---------------------------------------------------------------------------
# End-to-end FREM transformation
# ---------------------------------------------------------------------------

test_that("ferx_to_frem() returns ferx_frem with correct metadata", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  result <- ferx_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  # Class.
  expect_s3_class(result, "ferx_frem")

  # Output files exist.
  expect_true(file.exists(result$model_path))
  expect_true(file.exists(result$data_path))

  # Total etas: 3 PK + 2 covariate = 5.
  expect_equal(result$n_total_etas, 5L)

  # FREMTYPE map.
  expect_equal(length(result$fremtype_map), 2L)
  expect_equal(as.integer(result$fremtype_map[["WT"]]),  100L)
  expect_equal(as.integer(result$fremtype_map[["AGE"]]), 200L)

  # Covariate means: WT = 72.0, AGE = 37.6.
  expect_equal(unname(result$covariate_means[["WT"]]),  72.0, tolerance = 0.01)
  expect_equal(unname(result$covariate_means[["AGE"]]), 37.6, tolerance = 0.01)

  # Covariate variances: WT ≈ 111.56, AGE ≈ 99.38.
  expect_equal(unname(result$covariate_variances[["WT"]]),  111.56, tolerance = 0.5)
  expect_equal(unname(result$covariate_variances[["AGE"]]),  99.38, tolerance = 0.5)
})

test_that("FREM dataset has correct row count", {
  tmp <- tempfile("frem_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))

  model_path <- write_warfarin_model(tmp)
  data_path  <- write_warfarin_with_covariates(tmp)

  result <- ferx_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  frem_data <- read.csv(result$data_path)

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

  result <- ferx_to_frem(
    model      = model_path,
    data       = data_path,
    covariates = c("WT", "AGE"),
    output_dir = tmp
  )

  fit <- ferx_fit(
    result$model_path,
    result$data_path,
    method   = "focei",
    verbose  = FALSE,
    settings = list(maxiter = 3L, covariance = FALSE)
  )

  expect_s3_class(fit, "ferx_fit")
  expect_true(is.finite(fit$ofv))

  # Omega should be 5x5.
  expect_equal(nrow(fit$omega), 5L)
  expect_equal(ncol(fit$omega), 5L)

  # All omega diagonals positive.
  expect_true(all(diag(fit$omega) > 0))
})
