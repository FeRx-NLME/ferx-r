# ferx_save_fit() / ferx_load_fit() round-trip and CLI-flag tests.
#
# These exercise the .fitrx bundle format. The cross-tool test (Rust CLI
# writes, R reads) lives in ferx-core; here we only verify the R-side
# writer + reader are inverses, and that the `output` argument of
# ferx_fit() saves a valid bundle.

test_that("ferx_save_fit + ferx_load_fit round-trip a real fit", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_s3_class(loaded, "ferx_fit")
  expect_equal(loaded$method, fit$method)
  expect_equal(loaded$converged, fit$converged)
  expect_equal(loaded$ofv, fit$ofv, tolerance = 1e-12)
  expect_equal(loaded$aic, fit$aic, tolerance = 1e-12)
  expect_equal(loaded$bic, fit$bic, tolerance = 1e-12)
  expect_equal(unname(loaded$theta), unname(fit$theta), tolerance = 1e-12)
  expect_equal(names(loaded$theta), names(fit$theta))
  expect_equal(loaded$omega, fit$omega, tolerance = 1e-12)
  expect_equal(unname(loaded$sigma), unname(fit$sigma), tolerance = 1e-12)
  expect_equal(loaded$n_subjects, fit$n_subjects)
  expect_equal(loaded$n_obs, fit$n_obs)
  expect_equal(loaded$shrinkage_eta, fit$shrinkage_eta, tolerance = 1e-12)
  expect_equal(loaded$model_name, fit$model_name)
  expect_equal(loaded$ferx_version, fit$ferx_version)
})

test_that("covariance step results survive round-trip", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "no covariance step in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_true(is.matrix(loaded$cov_matrix))
  expect_equal(dim(loaded$cov_matrix), dim(fit$cov_matrix))
  expect_equal(loaded$cov_matrix, fit$cov_matrix, tolerance = 1e-12,
               ignore_attr = "dimnames")
  expect_equal(unname(loaded$se_theta), unname(fit$se_theta), tolerance = 1e-12)
})

test_that("predictions survive round-trip (per-row equality)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_true(is.data.frame(loaded$sdtab))
  expect_equal(nrow(loaded$sdtab), nrow(fit$sdtab))
  numeric_cols <- intersect(
    c("TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES", "EBE_OFV"),
    intersect(names(loaded$sdtab), names(fit$sdtab))
  )
  for (col in numeric_cols) {
    expect_equal(
      as.numeric(loaded$sdtab[[col]]),
      as.numeric(fit$sdtab[[col]]),
      tolerance = 1e-9,
      info = sprintf("column %s", col)
    )
  }
})

test_that("ferx_fit(output = ...) writes a valid bundle", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  fit <- ferx_fit(
    ex$model, ex$data,
    method = "focei", verbose = FALSE, covariance = FALSE,
    settings = list(maxiter = 10L),
    output = path
  )
  expect_true(file.exists(path))

  loaded <- ferx_load_fit(path)
  expect_s3_class(loaded, "ferx_fit")
  expect_equal(loaded$ofv, fit$ofv, tolerance = 1e-12)
})

test_that("include_data = TRUE bundles the input CSV", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  fit <- ferx_fit(
    ex$model, ex$data,
    method = "focei", verbose = FALSE, covariance = FALSE,
    settings = list(maxiter = 5L),
    output = path, include_data = TRUE
  )
  expect_true(file.exists(path))

  entries <- utils::unzip(path, list = TRUE)$Name
  expect_true("data.csv" %in% entries)

  loaded <- ferx_load_fit(path)
  expect_true(!is.na(loaded$data_path))
  expect_true(file.exists(loaded$data_path))
})

test_that("manifest reports format_version 1", {
  skip_on_cran()
  fit <- warfarin_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)

  staging <- tempfile("fitrx_inspect_")
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  dir.create(staging)
  utils::unzip(path, exdir = staging)
  manifest <- jsonlite::read_json(file.path(staging, "manifest.json"),
                                  simplifyVector = TRUE)
  expect_equal(manifest$format_version, "1")
  expect_true("fit.json" %in% manifest$entries)
})

test_that("bad path errors cleanly", {
  expect_error(ferx_load_fit("does-not-exist.fitrx"), "does not exist")
  bogus <- tempfile(fileext = ".fitrx")
  writeLines("not a zip", bogus)
  on.exit(unlink(bogus), add = TRUE)
  expect_error(ferx_load_fit(bogus))
})

test_that("ferx_save_fit rejects non-ferx_fit inputs", {
  expect_error(ferx_save_fit(list(), tempfile()), "ferx_fit")
})
