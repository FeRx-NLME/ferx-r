
# ---- header from test-bayes.R ----
# Tier-1 tests for the Bayes (method = "bayes") R-side surface: the R-hat
# convergence flag, the print/summary rendering of `$bayes`, and the persist
# wire builder. These exercise the wrapper logic with a mocked posterior so
# they run in milliseconds — the real MCMC engine round-trip is covered by the
# (slow) integration tests in test-fit.R / test-persist.R.

# A minimal but complete `$bayes` posterior summary, shaped exactly as
# `fit_result_to_list()` (Rust) and `ferx_load_fit()` (persist) produce it.
fake_bayes <- function(max_rhat = 1.004) {
  list(
    n_chains          = 2,
    n_warmup          = 300,
    n_draws_per_chain = 600,
    n_divergent       = 0,
    max_rhat          = max_rhat,
    param_names       = c("TVCL", "TVV", "OMEGA_CL", "SIGMA_PROP"),
    mean              = c(0.133, 8.0, 0.09, 0.10),
    sd                = c(0.01, 0.4, 0.02, 0.01),
    q025              = c(0.115, 7.2, 0.05, 0.08),
    median            = c(0.133, 8.0, 0.09, 0.10),
    q975              = c(0.151, 8.8, 0.13, 0.12),
    rhat              = c(1.001, 1.002, 1.004, 1.000),
    ess_bulk          = c(820, 760, 540, 910),
    ess_tail          = c(800, 740, 520, 900),
    mcse              = c(3e-4, 0.01, 8e-4, 3e-4)
  )
}

# .ferx_bayes_rhat_flag ------------------------------------------------------



# print.ferx_fit — BAYES posterior block -------------------------------------




# summary.ferx_fit / print.ferx_summary --------------------------------------



# .fitrx_build_bayes_wire ----------------------------------------------------




test_that(".fitrx_build_bayes_wire returns NULL when there is no posterior", {
  expect_null(ferx:::.fitrx_build_bayes_wire(list(bayes = NULL)))
  expect_null(ferx:::.fitrx_build_bayes_wire(list()))
})
test_that(".fitrx_build_bayes_wire serialises every posterior field", {
  wire <- ferx:::.fitrx_build_bayes_wire(list(bayes = fake_bayes()))
  expect_false(is.null(wire))
  # Scalars survive as plain numerics; param_names become a JSON-friendly list.
  expect_equal(wire$n_chains, 2)
  expect_equal(wire$n_draws_per_chain, 600)
  expect_equal(wire$max_rhat, 1.004)
  expect_equal(unlist(wire$param_names), c("TVCL", "TVV", "OMEGA_CL", "SIGMA_PROP"))
  # Per-parameter vectors round-trip through the unwrap helpers unchanged.
  expect_equal(ferx:::.fitrx_unwrap_opt_num_vec(wire$mean), fake_bayes()$mean)
  expect_equal(ferx:::.fitrx_unwrap_opt_num_vec(wire$ess_bulk), fake_bayes()$ess_bulk)
})

# ---- header from test-persist.R ----
# ferx_save_fit() / ferx_load_fit() round-trip and CLI-flag tests.
#
# These exercise the .fitrx bundle format. The cross-tool test (Rust CLI
# writes, R reads) lives in ferx-core; here we only verify the R-side
# writer + reader are inverses, and that the `output` argument of
# ferx_fit() saves a valid bundle.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.fitrx_wire_to_fit      <- getFromNamespace(".fitrx_wire_to_fit",      "ferx")
.fitrx_build_iov_wire   <- getFromNamespace(".fitrx_build_iov_wire",   "ferx")





















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
test_that("bare (extension-less) output path defaults to .fitrx and round-trips (#268)", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  dir <- tempfile("fitrx_ext_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  bare <- file.path(dir, "run1_base_diag")
  ferx_save_fit(fit, bare)

  # The archive is written with the conventional `.fitrx` extension, not a
  # bare name and not Info-ZIP's `.zip`.
  expect_true(file.exists(paste0(bare, ".fitrx")))
  expect_false(file.exists(bare))
  expect_false(file.exists(paste0(bare, ".zip")))

  # ... and loads back from the same bare path the caller saved to.
  loaded <- ferx_load_fit(bare)
  expect_s3_class(loaded, "ferx_fit")
  expect_equal(unname(loaded$theta), unname(fit$theta), tolerance = 1e-12)

  # An explicit extension is honoured verbatim (not forced to .fitrx).
  custom <- file.path(dir, "run1.bundle")
  ferx_save_fit(fit, custom)
  expect_true(file.exists(custom))
  expect_s3_class(ferx_load_fit(custom), "ferx_fit")
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
test_that("the model-selection surface survives round-trip", {
  skip_on_cran()
  fit <- warfarin_fit_cov()

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  # `bic_inputs` and `left_init` travel at wire level (they are ferx-core
  # FitWire fields); the three derived verdicts travel under r_extras.
  expect_equal(loaded$bic_inputs, fit$bic_inputs)
  expect_equal(loaded$left_init, fit$left_init)
  expect_equal(loaded$stalled_at_init, fit$stalled_at_init)
  expect_equal(loaded$estimate_near_boundary, fit$estimate_near_boundary)
  expect_equal(loaded$max_abs_correlation, fit$max_abs_correlation,
               tolerance = 1e-12)

  for (ty in c("mixed", "fixed", "iiv", "random")) {
    expect_equal(ferx_bic(loaded, ty), ferx_bic(fit, ty), tolerance = 1e-12)
  }
  # The condition number and the packed-Omega layout travel at wire level too.
  # Reading only the wire spelling on save wrote a null for every bundle this
  # package writes, so the condition-number gate landed in `skipped` (not
  # `failures`) on every reload - compare both halves of the verdict.
  expect_equal(
    loaded$cov_condition_number %||% loaded$condition_number,
    fit$condition_number %||% fit$cov_condition_number,
    tolerance = 1e-8
  )
  expect_equal(loaded$omega_is_diagonal, fit$omega_is_diagonal)
  expect_equal(loaded$kappa_is_diagonal, fit$kappa_is_diagonal)

  # The engine's stalled_at_init() checks init/estimate shapes before reading
  # `left_init`, so the inits have to survive for the verdict to mean anything.
  expect_equal(loaded$theta_init, fit$theta_init, tolerance = 1e-12)
  expect_equal(loaded$sigma_init, fit$sigma_init, tolerance = 1e-12)
  expect_equal(loaded$omega_init_dim, as.integer(fit$omega_init_dim))
  expect_equal(loaded$omega_init, as.numeric(fit$omega_init), tolerance = 1e-12)

  strict_loaded <- check_strictness(loaded, max_condition_number = 1e10)
  strict_fresh  <- check_strictness(fit, max_condition_number = 1e10)
  expect_equal(strict_loaded$failures, strict_fresh$failures)
  expect_equal(strict_loaded$skipped, strict_fresh$skipped)
  expect_equal(strict_loaded$passed, strict_fresh$passed)
})

test_that("tri-state verdicts stay NA across round-trip when unrecorded", {
  skip_on_cran()
  fit <- warfarin_fit_cov()

  # An optimizer that recorded no verdict leaves these NA. NA serialises to
  # JSON null, and assigning NULL back would drop the key entirely - so
  # `is.na(loaded$stalled_at_init)` would error on a zero-length vector
  # instead of reporting "no verdict".
  fit$stalled_at_init <- NA
  fit$estimate_near_boundary <- NA
  fit$max_abs_correlation <- NA_real_

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_true(is.na(loaded$stalled_at_init))
  expect_true(is.na(loaded$estimate_near_boundary))
  expect_true(is.na(loaded$max_abs_correlation))
  expect_length(loaded$stalled_at_init, 1L)
  expect_length(loaded$max_abs_correlation, 1L)
})

test_that("bayes posterior summary survives round-trip", {
  skip_on_cran()
  fit <- warfarin_bayes_fit()
  skip_if(is.null(fit$bayes), "no bayes result in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$bayes))
  b0 <- fit$bayes
  b1 <- loaded$bayes
  expect_equal(b1$n_chains, b0$n_chains)
  expect_equal(b1$n_warmup, b0$n_warmup)
  expect_equal(b1$n_draws_per_chain, b0$n_draws_per_chain)
  expect_equal(b1$max_rhat, b0$max_rhat, tolerance = 1e-12)
  expect_equal(b1$param_names, b0$param_names)
  expect_equal(b1$mean, b0$mean, tolerance = 1e-12)
  expect_equal(b1$q975, b0$q975, tolerance = 1e-12)
  expect_equal(b1$ess_bulk, b0$ess_bulk, tolerance = 1e-9)
})
test_that("SAEM conditional-distribution results survive round-trip (real fit)", {
  skip_on_cran()
  fit <- warfarin_saem_conddist_fit()
  skip_if(is.null(fit$cond_dist), "no cond_dist result in this fit")

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$cond_dist))
  expect_equal(loaded$cond_dist$data$ID, fit$cond_dist$data$ID)
  expect_equal(loaded$cond_dist$data$ETA, fit$cond_dist$data$ETA)
  expect_equal(loaded$cond_dist$data$COND_MEAN, fit$cond_dist$data$COND_MEAN,
               tolerance = 1e-9)
  expect_equal(loaded$cond_dist$data$COND_SD, fit$cond_dist$data$COND_SD,
               tolerance = 1e-9)
  expect_equal(loaded$cond_dist$data$COND_MODE, fit$cond_dist$data$COND_MODE,
               tolerance = 1e-9)
})
test_that(".fitrx_has_cond_dist / .fitrx_write_conddist_csv (synthetic)", {
  no_cd <- structure(list(cond_dist = NULL), class = "ferx_fit")
  expect_false(ferx:::.fitrx_has_cond_dist(no_cd))

  cd_data <- data.frame(
    ID = c("1", "1", "2", "2"),
    ETA = c("ETA_CL", "ETA_V", "ETA_CL", "ETA_V"),
    COND_MEAN = c(0.05, -0.02, -0.03, 0.01),
    COND_SD = c(0.12, 0.15, 0.11, 0.14),
    COND_MODE = c(0.04, -0.01, -0.02, 0.02),
    stringsAsFactors = FALSE
  )
  with_cd <- structure(
    list(cond_dist = list(data = cd_data, shrinkage = c(0.1, 0.05),
                           nsamp = 20L, burnin = 5L)),
    class = "ferx_fit"
  )
  expect_true(ferx:::.fitrx_has_cond_dist(with_cd))

  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_conddist_csv(with_cd, path)
  written <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(written$ID, c(1L, 1L, 2L, 2L))
  expect_equal(written$ETA, cd_data$ETA)
  expect_equal(written$COND_MEAN, cd_data$COND_MEAN, tolerance = 1e-12)
})
test_that("conddist.csv survives a ferx_save_fit / ferx_load_fit round-trip (synthetic)", {
  cd_data <- data.frame(
    ID = c("1", "1", "2", "2"),
    ETA = c("ETA_CL", "ETA_V", "ETA_CL", "ETA_V"),
    COND_MEAN = c(0.05, -0.02, -0.03, 0.01),
    COND_SD = c(0.12, 0.15, 0.11, 0.14),
    COND_MODE = c(0.04, -0.01, -0.02, 0.02),
    stringsAsFactors = FALSE
  )
  fake <- structure(
    list(
      theta = c(TVCL = 1.0, TVV = 10.0),
      omega = matrix(c(0.04, 0, 0, 0.09), 2L, 2L,
                     dimnames = list(c("ETA_CL", "ETA_V"), c("ETA_CL", "ETA_V"))),
      eta_names = c("ETA_CL", "ETA_V"),
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      cond_dist = list(data = cd_data, shrinkage = c(0.1, 0.05),
                        nsamp = 20L, burnin = 5L),
      # Numeric-looking IDs: ebe_etas$ID and cond_dist$data$ID must come back
      # the same type after a round-trip (#248 review).
      ebe_etas = data.frame(ID = c("1", "2"), ETA_CL = c(0.04, -0.02),
                             ETA_V = c(-0.01, 0.02), stringsAsFactors = FALSE),
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 2L, n_iterations = 1L,
      method = "SAEM", method_chain = "SAEM",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = c(0, 0), shrinkage_eps = 0,
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

  ferx_save_fit(fake, path)
  entries <- utils::unzip(path, list = TRUE)$Name
  expect_true("conddist.csv" %in% entries)

  loaded <- ferx_load_fit(path)
  expect_equal(loaded$cond_dist$data$ID, cd_data$ID)
  expect_equal(loaded$cond_dist$data$ETA, cd_data$ETA)
  expect_equal(loaded$cond_dist$data$COND_MEAN, cd_data$COND_MEAN, tolerance = 1e-12)
  # nsamp/burnin aren't recoverable from conddist.csv alone (#244); NA (not
  # 0L) so "unknown after reload" isn't confused with "zero draws retained".
  expect_equal(loaded$cond_dist$nsamp, NA_integer_)
  expect_equal(loaded$cond_dist$burnin, NA_integer_)
  # Shrinkage is recomputed from cond_mean + omega, matching ferx-core's own
  # conddist.csv reader formula.
  expect_equal(
    unname(loaded$cond_dist$shrinkage),
    c(
      1 - stats::sd(cd_data$COND_MEAN[cd_data$ETA == "ETA_CL"]) / sqrt(0.04),
      1 - stats::sd(cd_data$COND_MEAN[cd_data$ETA == "ETA_V"]) / sqrt(0.09)
    ),
    tolerance = 1e-9
  )
  # ebe_etas$ID and cond_dist$data$ID must have the same type after a
  # round-trip, so joins between the two tables don't silently break for
  # numeric-looking subject IDs (#248 review).
  expect_type(loaded$ebe_etas$ID, "character")
  expect_type(loaded$cond_dist$data$ID, "character")
})
test_that("recomputed shrinkage matches ferx-core's near-singular-omega floor (synthetic)", {
  # ferx-core's own conddist.csv reader (src/io/fitrx.rs) treats
  # Omega_jj < SAEM_OMEGA_DIAG_FLOOR (1e-6) as NaN, not just Omega_jj <= 0;
  # the R reimplementation must match or shrinkage silently diverges for a
  # near-singular (but positive) omega diagonal.
  cd_data <- data.frame(
    ID = c("1", "2", "3"),
    ETA = c("ETA_CL", "ETA_CL", "ETA_CL"),
    COND_MEAN = c(0.05, -0.03, 0.01),
    COND_SD = c(0.12, 0.11, 0.10),
    COND_MODE = c(0.04, -0.02, 0.02),
    stringsAsFactors = FALSE
  )
  omega <- matrix(1e-8, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL"))
  shrink <- ferx:::.fitrx_cond_dist_shrinkage(cd_data, omega, "ETA_CL")
  expect_true(is.na(shrink))
})
test_that("a fit with no cond_dist has no conddist.csv entry and loads NULL", {
  fake <- structure(
    list(
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      cond_dist = NULL,
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

  ferx_save_fit(fake, path)
  entries <- utils::unzip(path, list = TRUE)$Name
  expect_false("conddist.csv" %in% entries)

  loaded <- ferx_load_fit(path)
  expect_null(loaded$cond_dist)
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
test_that("init_as_sd flags survive a ferx_save_fit / ferx_load_fit round-trip (synthetic)", {
  # CI-runnable regression for the .fitrx writer dropping the SD-init flags.
  # The flags were read on load (omega/sigma/iov wire blocks) but never
  # written on save, so they silently vanished on round-trip. A synthetic
  # fit (no real estimation) keeps this in the default CI test pass; the
  # gated test below additionally exercises the path on a real fit.
  fake <- structure(
    list(
      theta = c(TVCL = 1.0, TVV = 10.0),
      omega = matrix(c(0.04, 0, 0, 0.09), 2L, 2L,
                     dimnames = list(c("ETA_CL", "ETA_V"), c("ETA_CL", "ETA_V"))),
      eta_names = c("ETA_CL", "ETA_V"),
      omega_init_as_sd = c(TRUE, FALSE),
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      sigma_init_as_sd = TRUE,
      # Minimal IOV so the iov wire (and kappa_init_as_sd) is written.
      omega_iov = matrix(0.02, 1L, 1L,
                         dimnames = list("KAPPA_CL", "KAPPA_CL")),
      kappa_names = "KAPPA_CL",
      kappa_fixed = FALSE,
      kappa_init_as_sd = TRUE,
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 2L, n_iterations = 1L,
      method = "FOCEI", method_chain = "FOCEI",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = c(0, 0), shrinkage_eps = 0,
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

  ferx_save_fit(fake, path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$omega_init_as_sd, c(TRUE, FALSE))
  expect_equal(loaded$sigma_init_as_sd, TRUE)
  expect_equal(loaded$kappa_init_as_sd, TRUE)
})
test_that("fit$trace survives a ferx_save_fit / ferx_load_fit round-trip (synthetic)", {
  # fit$trace_path is a temp file and typically won't exist on reload; the
  # bundled trace.csv is what actually needs to round-trip (see #228).
  trace_path <- write_fake_trace(n = 3L, method = "gn")
  on.exit(unlink(trace_path), add = TRUE)
  fake <- structure(
    list(
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 1L, n_iterations = 1L,
      method = "GN", method_chain = "GN",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "Enzyme AD",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_,
      trace_path = trace_path,
      trace = ferx_trace(trace_path)
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fake, path)
  loaded <- ferx_load_fit(path)

  expect_s3_class(loaded$trace, "data.frame")
  expect_equal(nrow(loaded$trace), 3L)
  expect_equal(loaded$trace$iter, fake$trace$iter)
  expect_equal(loaded$trace$ofv, fake$trace$ofv, tolerance = 1e-12)
  # grad_norm/step_norm are set by write_fake_trace(); lm_lambda etc. are NA
  # in a "gn" trace's non-applicable columns and must stay NA (not NaN or 0)
  # after the JSON-free CSV round-trip.
  expect_equal(loaded$trace$grad_norm, fake$trace$grad_norm, tolerance = 1e-12)
  expect_true(all(is.na(loaded$trace$lm_lambda)))
})
test_that("fit$trace bundles even as a header-only (zero-row) data.frame", {
  # A valid zero-iteration trace must still round-trip -- gating on
  # `nrow() > 0L` would silently drop it (see PR #246 review).
  trace_path <- write_fake_trace(n = 3L, method = "gn")
  on.exit(unlink(trace_path), add = TRUE)
  empty_trace <- ferx_trace(trace_path)[0, , drop = FALSE]
  fake <- structure(
    list(
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 1L, n_iterations = 1L,
      method = "GN", method_chain = "GN",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "Enzyme AD",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_,
      trace_path = trace_path,
      trace = empty_trace
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fake, path)
  expect_true("trace.csv" %in% unzip(path, list = TRUE)$Name)
  loaded <- ferx_load_fit(path)
  expect_s3_class(loaded$trace, "data.frame")
  expect_equal(nrow(loaded$trace), 0L)
})
test_that("fit$trace is absent from the bundle (and NULL on load) when no trace was collected", {
  fake <- structure(
    list(
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

  ferx_save_fit(fake, path)
  expect_setequal(unzip(path, list = TRUE)$Name, c(
    "manifest.json", "fit.json", "ebes.csv", "predictions.csv",
    "model.ferx", "warnings.txt"
  ))
  loaded <- ferx_load_fit(path)
  expect_null(loaded$trace)
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
  wire <- .fitrx_build_iov_wire(fake_iov_fit)
  expect_equal(wire$shrinkage_kappa_by_occ, list())
})
test_that(".fitrx_build_iov_wire omits the #1031 weight keys for an unweighted model", {
  # Regression: `list(kappa_weights = NULL)` keeps the *name*, so write_json()
  # emitted `"kappa_weights": null`. ferx-core declares the field as
  # `#[serde(default)] Vec<Option<String>>`, and serde(default) only fires for
  # a missing key - an explicit null fails with "invalid type: null, expected
  # a sequence". The keys must be absent, not null.
  fake_iov_fit <- list(
    omega_iov       = matrix(0.09, 1, 1),
    kappa_names     = "KAPPA_CL",
    kappa_fixed     = FALSE,
    shrinkage_kappa = 0.15,
    kappa_weights   = NULL
  )
  wire <- .fitrx_build_iov_wire(fake_iov_fit)
  expect_false("kappa_weights" %in% names(wire))
  expect_false("kappa_weight_typical" %in% names(wire))
  json <- as.character(jsonlite::toJSON(wire, auto_unbox = TRUE,
                                        null = "null", na = "null"))
  expect_false(grepl("kappa_weight", json, fixed = TRUE))
})
test_that(".fitrx_build_iov_wire writes the #1031 weight keys as null-holed arrays", {
  fake_iov_fit <- list(
    omega_iov            = matrix(0.09, 2, 2),
    kappa_names          = c("KAPPA_CL", "KAPPA_V"),
    kappa_fixed          = c(FALSE, FALSE),
    shrinkage_kappa      = c(0.15, 0.2),
    kappa_weights        = c(NA_character_, "NARM"),
    kappa_weight_typical = c(NA_real_, 400)
  )
  wire <- .fitrx_build_iov_wire(fake_iov_fit)
  expect_equal(wire$kappa_weights, list(NA_character_, "NARM"))
  expect_equal(wire$kappa_weight_typical, list(NA_real_, 400))
  # `na = "null"` is what the writer uses; the holes must survive as JSON null
  # so `Vec<Option<..>>` keeps every entry lined up with its kappa.
  expect_equal(as.character(jsonlite::toJSON(wire$kappa_weights,
                                             auto_unbox = TRUE, na = "null")),
               '[null,"NARM"]')
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

# ---- header from test-persist-helpers.R ----
# Unit tests for the pure-R serialization helpers in persist.R. These convert
# between the in-memory ferx_fit fields and the cross-language .fitrx wire
# schema, and need no model fit to exercise — so they live in the fast tier.

# ---------------------------------------------------------------------------
# Method name <-> token round-trip
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Covariance status <-> token
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Error model inference from sigma types
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Matrix <-> wire (row-major)
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Confidence-interval <-> wire
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Optional scalar wrappers (NULL/NA/empty collapse to NULL)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Named SE unwrap + subject string IDs
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-subject OFV / N_OBS extraction from sdtab
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# EBEs CSV writer — header-only path when there are no random effects
# ---------------------------------------------------------------------------


test_that(".fitrx_method_to_token normalises labels and defaults to focei", {
  expect_identical(ferx:::.fitrx_method_to_token(NULL), "focei")
  expect_identical(ferx:::.fitrx_method_to_token(""), "focei")
  expect_identical(ferx:::.fitrx_method_to_token("FOCEI"), "focei")   # case-insensitive
  expect_identical(ferx:::.fitrx_method_to_token("FOCE"), "foce")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce_gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn-hybrid"), "foce_gn_hybrid")
  expect_identical(ferx:::.fitrx_method_to_token("SAEM"), "saem")
  expect_identical(ferx:::.fitrx_method_to_token("Mystery"), "mystery") # unknown -> tolower
})
test_that("method label/token round-trips through the canonical display names", {
  for (lbl in c("FOCE", "FOCEI", "FOCE-GN", "FOCE-GN-Hybrid", "SAEM")) {
    expect_identical(ferx:::.fitrx_method_label(ferx:::.fitrx_method_to_token(lbl)), lbl)
  }
})
test_that(".fitrx_covariance_status_to_token normalises and defaults", {
  expect_identical(ferx:::.fitrx_covariance_status_to_token(NULL), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token(""), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("Computed"), "computed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("FAILED"), "failed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("NotRequested"), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("sir_fallback"), "sir_fallback")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("SirFallback"), "sir_fallback")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("weird"), "weird")
})
test_that(".fitrx_error_model_from_sigma_types classifies sigma vectors", {
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(NULL), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(character(0)), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("proportional", "Proportional")), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive")), "additive")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive", "proportional")), "combined")
})
test_that(".fitrx_matrix_to_wire emits row-major data and dims", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE) # [[1,2,3],[4,5,6]]
  w <- ferx:::.fitrx_matrix_to_wire(m)
  expect_identical(w$rows, 2L)
  expect_identical(w$cols, 3L)
  expect_identical(w$data, c(1, 2, 3, 4, 5, 6)) # row-major
})
test_that(".fitrx_matrix_to_wire returns NULL for NULL / non-matrix input", {
  expect_null(ferx:::.fitrx_matrix_to_wire(NULL))
  expect_null(ferx:::.fitrx_matrix_to_wire(c(1, 2, 3)))
})
test_that("matrix round-trips through wire form unchanged", {
  m <- matrix(c(0.1, 0.0, 0.0, 0.2), nrow = 2)
  expect_equal(ferx:::.fitrx_matrix_from_wire(ferx:::.fitrx_matrix_to_wire(m)), m)
})
test_that(".fitrx_ci_to_wire accepts a 2-column matrix", {
  ci <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE) # rows (1,2) (3,4)
  w  <- ferx:::.fitrx_ci_to_wire(ci)
  expect_length(w, 2L)
  expect_identical(w[[1L]], c(1, 2))
  expect_identical(w[[2L]], c(3, 4))
})
test_that(".fitrx_ci_to_wire accepts the legacy flat layout and rejects odd length", {
  expect_equal(ferx:::.fitrx_ci_to_wire(c(1, 2, 3, 4)), list(c(1, 2), c(3, 4)))
  expect_null(ferx:::.fitrx_ci_to_wire(c(1, 2, 3)))   # odd length
  expect_null(ferx:::.fitrx_ci_to_wire(NULL))
})
test_that("CI round-trips matrix -> wire -> matrix", {
  ci <- matrix(c(0.5, 1.5, 2.5, 3.5), ncol = 2, byrow = TRUE)
  back <- ferx:::.fitrx_unwrap_ci(ferx:::.fitrx_ci_to_wire(ci))
  expect_equal(back, ci)
  expect_null(ferx:::.fitrx_unwrap_ci(NULL))
})
test_that(".fitrx_subject_string_ids extracts character IDs or NULL", {
  expect_null(ferx:::.fitrx_subject_string_ids(list()))
  expect_null(ferx:::.fitrx_subject_string_ids(list(ebe_etas = data.frame(x = 1))))
  fit <- list(ebe_etas = data.frame(ID = c(1, 2, 10), eta = c(0, 0, 0)))
  expect_identical(ferx:::.fitrx_subject_string_ids(fit), c("1", "2", "10"))
})
test_that(".fitrx_per_subject_ofv_nobs pulls one row per subject from sdtab", {
  expect_null(ferx:::.fitrx_per_subject_ofv_nobs(list()))
  expect_null(ferx:::.fitrx_per_subject_ofv_nobs(list(sdtab = data.frame(ID = 1)))) # missing cols
  fit <- list(sdtab = data.frame(
    ID = c("A", "A", "B"),
    EBE_OFV = c(1.5, 1.5, 2.5),
    N_OBS = c(3L, 3L, 4L)
  ))
  res <- ferx:::.fitrx_per_subject_ofv_nobs(fit)
  expect_identical(res$id, c("A", "B"))
  expect_identical(res$ofv, c(1.5, 2.5))
  expect_identical(res$n_obs, c(3L, 4L))
})
test_that(".fitrx_write_ebes_csv writes a stable header when ebes are absent", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_ebes_csv(list(ebe_etas = NULL), path)
  expect_identical(readLines(path), "ID,ofv_contribution,n_obs")
})

# ---- header from test-sde.R ----
# Tests for SDE (diffusion / EKF) surface on the R side.
# All tests use make_fake_fit() from helper-trace.R — no real fit needed.

# uses_sde flag in print.ferx_fit ----------------------------------------





# uses_sde flag in summary / print.ferx_summary ---------------------------






# DIFF_* theta parameters flow through naturally --------------------------


# ferx_estimates() with DIFF_* theta -------------------------------------


# ferx_save_fit / ferx_load_fit round-trip --------------------------------




test_that("uses_sde survives a ferx_save_fit / ferx_load_fit round-trip (TRUE)", {
  fit <- make_fake_fit(
    uses_sde = TRUE,
    theta    = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    omega    = matrix(0.09, 1, 1),
    sigma    = 1.0,
    sigma_names = "ADD"
  )
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp))
  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)
  expect_true(isTRUE(fit2$uses_sde))
})
test_that("uses_sde survives a ferx_save_fit / ferx_load_fit round-trip (FALSE)", {
  fit <- make_fake_fit(
    uses_sde = FALSE,
    omega    = matrix(0.09, 1, 1),
    sigma    = 1.0,
    sigma_names = "ADD"
  )
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp))
  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)
  expect_false(isTRUE(fit2$uses_sde))
})

# ---- header from test-covtab.R ----
# Covariate table (`fit$covtab`) + types (`fit$covariate_types`), surfaced from
# a `[covariates]` block. See ferx-core #182/#184 for the underlying DSL + table.

# Cached fit on the covariate example so the multiple tests below share one
# estimation (mirrors helper-warfarin-fit.R's memoised pattern).
cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("two_cpt_oral_cov")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "focei", verbose = FALSE,
        covariance = FALSE, settings = list(maxiter = 5L)
      )
    }
    fit
  }
})








test_that("covtab and covariate_types survive a .fitrx save/load round-trip", {
  fit <- cov_fit()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)

  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)

  expect_s3_class(fit2$covtab, "data.frame")
  expect_identical(names(fit2$covtab), names(fit$covtab))
  expect_equal(nrow(fit2$covtab), nrow(fit$covtab))
  expect_type(fit2$covtab$ID, "character")
  expect_identical(fit2$covariate_types, fit$covariate_types)
})

# ---- header from test-data-selection.R ----
# Tests for [data_selection] / ferx_apply_selection() / ferx_fit(ignore=) support.

# Local alias avoids the ::: operator.
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# .sel_attr()/warfarin_sel_fit() come from helper-selection.R

# ---------------------------------------------------------------------------
# 1. Plain warfarin fit has NULL exclusions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 2. [data_selection] model populates exclusions
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 3. ferx_fit(data, ignore = ...) produces exclusions matching model-file path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4. ferx_apply_selection() pure-R preview
# ---------------------------------------------------------------------------







# ---------------------------------------------------------------------------
# 5. ferx_data passed directly to ferx_fit()
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6. print.ferx_fit() shows DATA SELECTION block
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 7. ferx_runlog() shows exclusion line
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 8. persist round-trip preserves exclusions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 9. print.ferx_data
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 10. ferx_apply_selection(excluded = TRUE) on a ferx_data
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 11. accept= and ignore_ids= coverage for ferx_apply_selection()
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 12. ferx_apply_selection() data.frame input
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 13. ferx_apply_selection(fit, excluded = TRUE) replays fired conditions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 14. ferx_apply_selection(fit, excluded = TRUE) NULL exclusions path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 15. ferx_fit() type validation stops for ignore/accept/ignore_ids
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 16. ferx_fit() accept= and ignore_ids= settings_parts code paths
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 17. print.ferx_fit fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 18. ferx_runlog fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 19. .cmp() operator coverage (pure R, via ferx_apply_selection() with inline df)
# ---------------------------------------------------------------------------

# .sel_test_df() comes from helper-selection.R





# ---------------------------------------------------------------------------
# 20. .parse_filter_expr() coverage: character rhs, &&, !found
# ---------------------------------------------------------------------------





test_that("ferx_save/ferx_load round-trip preserves exclusions", {
  skip_on_cran()
  fit  <- warfarin_sel_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_false(is.null(loaded$exclusions))
  expect_equal(loaded$exclusions$n_obs_excluded,   fit$exclusions$n_obs_excluded)
  expect_equal(loaded$exclusions$n_records_total,  fit$exclusions$n_records_total)
  expect_equal(loaded$exclusions$fired_ignore,     fit$exclusions$fired_ignore)
  expect_equal(loaded$exclusions$excluded_subject_ids, fit$exclusions$excluded_subject_ids)
})

test_that("sample-size-weighted IOV survives a ferx_save_fit / ferx_load_fit round-trip", {
  # ferx-core #1031: `kappa K ~ g2 weight = NARM`. The .fitrx `iov` block must
  # carry both the weight source text and the typical (median) arm size, and
  # both must come back as JSON *arrays* - a single weighted kappa is the
  # common MBMA case, and auto_unbox would otherwise write a bare scalar that
  # the engine's Vec<Option<..>> reader rejects.
  fake <- structure(
    list(
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L,
                     dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      omega_iov = matrix(2.0, 1L, 1L,
                         dimnames = list("KAPPA_CL", "KAPPA_CL")),
      kappa_names = "KAPPA_CL",
      kappa_fixed = FALSE,
      kappa_init_as_sd = TRUE,
      kappa_weights = c(KAPPA_CL = "NARM"),
      kappa_weight_typical = c(KAPPA_CL = 400),
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 2L, n_iterations = 1L,
      method = "FOCE", method_chain = "FOCE",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "analytic",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fake, path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$kappa_weights, c(KAPPA_CL = "NARM"))
  expect_equal(loaded$kappa_weight_typical, c(KAPPA_CL = 400))
})
test_that("an unweighted IOV fit writes no kappa weight fields", {
  # A pre-#1031 bundle carries neither field, and an unweighted model must
  # keep producing exactly that bundle.
  fake <- structure(
    list(
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05),
      sigma_names = "prop",
      sigma_types = "proportional",
      omega_iov = matrix(0.02, 1L, 1L,
                         dimnames = list("KAPPA_CL", "KAPPA_CL")),
      kappa_names = "KAPPA_CL",
      kappa_fixed = FALSE,
      kappa_init_as_sd = FALSE,
      ofv = 0, aic = 2, bic = 4,
      n_obs = 3L, n_subjects = 2L, n_parameters = 2L, n_iterations = 1L,
      method = "FOCE", method_chain = "FOCE",
      converged = TRUE,
      warnings = character(),
      shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "analytic",
      gradient_method_outer = "N/A",
      covariance_status = "NotRequested",
      model_source = "model fake\n",
      data_path = NA_character_
    ),
    class = "ferx_fit"
  )
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)

  ferx_save_fit(fake, path)
  loaded <- ferx_load_fit(path)

  expect_null(loaded$kappa_weights)
  expect_null(loaded$kappa_weight_typical)
})
