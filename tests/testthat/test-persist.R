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

test_that("LTBS residual label survives a .fitrx round-trip", {
  skip_on_cran()
  # A fit of a log-transform-both-sides model must report its residual as
  # "additive (log-transformed)" (set by the Rust glue's residual_label from
  # CompiledModel.log_transform), and that label must survive save/load via the
  # r_extras-persisted model_structure.
  md <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.2, 0.001, 10.0)",
    "  theta TVV(10.0, 0.1, 500.0)",
    "  theta TVKA(1.5, 0.01, 50.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma ADD_LOG ~ 0.15 (sd)",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "  KA = TVKA",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  log(DV) ~ additive(ADD_LOG)"
  ), md)
  dd <- tempfile(fileext = ".csv")
  writeLines(c(
    "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV",
    "1,0,.,1,100,1,0,1",
    "1,1,8.5,0,.,1,0,0",
    "1,4,5.1,0,.,1,0,0",
    "1,12,1.9,0,.,1,0,0",
    "2,0,.,1,120,1,0,1",
    "2,1,9.8,0,.,1,0,0",
    "2,4,6.0,0,.,1,0,0",
    "2,12,2.2,0,.,1,0,0"
  ), dd)
  on.exit(unlink(c(md, dd)), add = TRUE)

  fit <- ferx_fit(md, data = dd, verbose = FALSE,
                  settings = list(maxiter = 30L))
  expect_equal(fit$model_structure$residual, "additive (log-transformed)")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_equal(loaded$model_structure$residual, "additive (log-transformed)")
})

test_that(".fitrx_build_iov_wire uses empty list when shrinkage_kappa_by_occ is absent", {
  fake_iov_fit <- list(
    omega_iov              = matrix(0.09, 1, 1),
    kappa_names            = "KAPPA_CL",
    kappa_fixed            = FALSE,
    se_kappa               = NULL,
    shrinkage_kappa        = 0.15,
    shrinkage_kappa_by_occ = NULL,
    omega_iov_param_corr   = NULL
  )
  wire <- ferx:::.fitrx_build_iov_wire(fake_iov_fit)
  expect_equal(wire$shrinkage_kappa_by_occ, list())
})

test_that("shrinkage_kappa_by_occ survives a ferx_save/ferx_load round-trip", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  skip_if(is.null(fit$shrinkage_kappa_by_occ),
          "IOV fit did not produce shrinkage_kappa_by_occ - skip round-trip")
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_equal(colnames(loaded$shrinkage_kappa_by_occ),
               colnames(fit$shrinkage_kappa_by_occ))
  expect_equal(nrow(loaded$shrinkage_kappa_by_occ),
               nrow(fit$shrinkage_kappa_by_occ))
  expect_equal(loaded$shrinkage_kappa_by_occ$occ,
               fit$shrinkage_kappa_by_occ$occ)
})

test_that("covariate_names and input_columns survive ferx_save/ferx_load round-trip", {
  # Regression guard: covariate_names was not serialised before PR#126, so it
  # silently vanished after a save/load cycle. ferx_fit() populates both fields
  # from the engine result; the test injects them into a fake fit so the check
  # runs fast without a real fit.
  fake <- structure(
    list(
      sdtab = data.frame(
        ID = 1L, TIME = 1.0, DV = 1.0,
        PRED = 1.0, IPRED = 1.0,
        CWRES = 0.0, IWRES = 0.0,
        EBE_OFV = -1.0, N_OBS = 1L
      ),
      ebe_etas      = data.frame(ID = 1L, ETA_CL = 0.0),
      theta         = c(TVCL = 1.0),
      omega         = matrix(0.04, 1L, 1L,
                             dimnames = list("ETA_CL", "ETA_CL")),
      eta_names     = "ETA_CL",
      sigma         = c(prop = 0.05),
      sigma_names   = "prop",
      sigma_types   = "proportional",
      ofv = 0, aic = 2, bic = 4,
      n_obs = 1L, n_subjects = 1L,
      n_parameters = 1L, n_iterations = 1L,
      method = "FOCEI", method_chain = "FOCEI",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "FD", gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_,
      covariate_names = c("WT", "AGE"),
      input_columns  = c("ID", "TIME", "DV", "EVID", "AMT", "WT", "AGE")
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fake, path)
  loaded <- ferx_load_fit(path)
  expect_equal(loaded$covariate_names, c("WT", "AGE"))
  expect_equal(loaded$input_columns,
               c("ID", "TIME", "DV", "EVID", "AMT", "WT", "AGE"))
})

