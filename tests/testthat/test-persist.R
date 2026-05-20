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

test_that("save -> load -> save again works (warfarin)", {
  # End-to-end check: load a fit from a bundle and save it again. Warfarin's
  # numeric-string IDs come back from read.csv as integers, so this path
  # exercises only the numeric branch of the predictions writer.
  skip_on_cran()
  fit <- warfarin_fit()
  first <- tempfile(fileext = ".fitrx")
  second <- tempfile(fileext = ".fitrx")
  on.exit(unlink(c(first, second)), add = TRUE)

  ferx_save_fit(fit, first)
  reloaded <- ferx_load_fit(first)
  ferx_save_fit(reloaded, second)
  twice <- ferx_load_fit(second)
  expect_equal(twice$ofv, fit$ofv, tolerance = 1e-12)
})

test_that(".fitrx_write_predictions_csv handles character sdtab$ID", {
  # Regression for the Copilot review on PR #21: the writer used to assume
  # sdtab$ID was numeric and called max(sdtab$ID, ...), which errors on
  # character input. This synthesises a fit where IDs are alphanumeric (as
  # they would be after loading a study with non-numeric subject IDs from a
  # .fitrx) and verifies that saving and reloading succeed.
  fake <- structure(
    list(
      sdtab = data.frame(
        ID = c("PT001", "PT001", "PT002"),
        TIME = c(0.5, 1.0, 0.5),
        DV = c(5, 6, 4),
        PRED = c(5, 6, 4), IPRED = c(5, 6, 4),
        CWRES = c(0, 0, 0), IWRES = c(0, 0, 0),
        EBE_OFV = c(-1, -1, -2), N_OBS = c(2L, 2L, 1L),
        stringsAsFactors = FALSE
      ),
      ebe_etas = data.frame(
        ID = c("PT001", "PT002"),
        ETA_CL = c(0.1, -0.1),
        stringsAsFactors = FALSE
      ),
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 1L, n_iterations = 1L,
      method = "FOCEI", method_chain = "FOCEI",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "Enzyme AD",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  expect_silent(ferx_save_fit(fake, path))
  loaded <- ferx_load_fit(path)
  expect_equal(as.character(loaded$sdtab$ID), c("PT001", "PT001", "PT002"))
  expect_equal(loaded$n_subjects, 2L)
})

test_that("init_as_sd flags survive ferx_save_fit / ferx_load_fit round-trip", {
  skip_if_not(nzchar(Sys.getenv("FERX_RUN_REAL_FIT", "")),
              "Set FERX_RUN_REAL_FIT=1 to run real-fit persist tests")
  fit <- warfarin_fit()
  # Inject flags so we can round-trip them without a real (sd)-annotated fit
  fit$omega_init_as_sd <- c(TRUE, FALSE)
  fit$sigma_init_as_sd <- TRUE
  fit$kappa_init_as_sd <- logical(0L)

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$omega_init_as_sd, c(TRUE, FALSE))
  expect_equal(loaded$sigma_init_as_sd, TRUE)
  expect_equal(loaded$kappa_init_as_sd, logical(0L))
})

test_that("ferx_load_fit on old .fitrx without init_as_sd produces empty logicals", {
  # Build a fake wire list that omits init_as_sd (pre-PR#57 bundle)
  wire_omega <- list(
    names = list("ETA_CL"),
    matrix = list(list(0.09)),
    fixed = list(FALSE),
    log_transformed = list(TRUE),
    shrinkage = list(0.1),
    se = NULL,
    param_corr = NULL
    # init_as_sd intentionally absent
  )
  wire_sigma <- list(
    names = list("PROP_ERR"),
    estimates = list(0.01),
    fixed = list(FALSE),
    types = list("proportional"),
    se = NULL
    # init_as_sd intentionally absent
  )
  wire <- list(
    method = "focei", method_chain = list("focei"),
    converged = TRUE, ofv = -100, aic = -90, bic = -80,
    n_obs = 10L, n_subjects = 2L, n_parameters = 2L, n_iterations = 5L,
    interaction = TRUE, wall_time_secs = 1.0,
    gradient_method_inner = "", gradient_method_outer = "",
    covariance_status = "not_requested",
    omega = wire_omega, sigma = wire_sigma,
    theta = list(names = list("TVCL"), estimates = list(1.0),
                 fixed = list(FALSE), transform = list("log"), se = NULL)
  )
  result <- ferx:::.fitrx_wire_to_fit(wire)
  expect_equal(result$omega_init_as_sd, logical(0L))
  expect_equal(result$sigma_init_as_sd, logical(0L))
  expect_equal(result$kappa_init_as_sd, logical(0L))
})
