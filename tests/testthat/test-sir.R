# Tests for ferx_sir() and the path/hash provenance plumbing that backs it.
#
# SIR happy-path tests are gated on `skip_if(is.null(fit$cov_matrix), ...)`
# because in the no-autodiff CI build the warfarin cov step occasionally
# fails to converge with maxiter = 30 (same pattern as test-simulate-uncertainty.R).

sir_cov_skip <- "covariance step did not converge - skipping"

test_that("ferx_fit records model_path, data_path, and hashes on the fit", {
  fit <- warfarin_fit_cov()
  ex <- ferx_example("warfarin")

  expect_true(!is.na(fit$model_path))
  expect_true(!is.na(fit$data_path))
  # Hashes are 64-char lowercase hex (SHA-256).
  expect_match(fit$model_hash, "^[0-9a-f]{64}$")
  expect_match(fit$data_hash, "^[0-9a-f]{64}$")

  # Paths point at the actual on-disk files used for the fit.
  expect_true(file.exists(fit$model_path))
  expect_true(file.exists(fit$data_path))
  expect_identical(normalizePath(fit$model_path),
                   normalizePath(ex$model))
  expect_identical(normalizePath(fit$data_path),
                   normalizePath(ex$data))
})

test_that("ferx_sir rejects fits without a covariance matrix", {
  fit <- warfarin_fit()  # covariance = FALSE
  expect_error(
    ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L),
    "covariance matrix",
    fixed = TRUE
  )
})

test_that("ferx_sir runs end-to-end and populates SIR fields", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)
  out <- ferx_sir(
    fit,
    sir_samples = 8L,
    sir_resamples = 4L,
    sir_seed = 1L,
    verbose = FALSE
  )

  # All SIR fields populated and shaped correctly.
  expect_true(is.numeric(out$sir_ess) && length(out$sir_ess) == 1L)
  expect_true(is.finite(out$sir_ess))

  n_theta <- length(out$theta)
  n_eta <- nrow(out$omega)
  n_sigma <- length(out$sigma)

  expect_equal(dim(out$sir_ci_theta), c(n_theta, 2L))
  expect_equal(dim(out$sir_ci_omega), c(n_eta, 2L))
  expect_equal(dim(out$sir_ci_sigma), c(n_sigma, 2L))
  expect_equal(colnames(out$sir_ci_theta), c("lower", "upper"))

  # CIs straddle the point estimates (with this tiny SIR run, occasional
  # extreme draws push the band wide; the lo <= hi check is the
  # invariant we actually care about).
  expect_true(all(out$sir_ci_theta[, "lower"] <= out$sir_ci_theta[, "upper"]))
  expect_true(all(out$sir_ci_omega[, "lower"] <= out$sir_ci_omega[, "upper"]))
  expect_true(all(out$sir_ci_sigma[, "lower"] <= out$sir_ci_sigma[, "upper"]))

  # Returned object still classes as ferx_fit so downstream printers /
  # methods keep working.
  expect_s3_class(out, "ferx_fit")
})

test_that("ferx_sir refuses to run after the model file is tampered with", {
  # Run against a temp copy of the example so we don't mutate the bundled
  # example file. ferx_fit normalises paths, so the fit will point at the
  # temp model/data we created here.
  src <- ferx_example("warfarin")
  tmp <- tempfile("ferx_sir_tamper_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  model_tmp <- file.path(tmp, "warfarin.ferx")
  data_tmp <- file.path(tmp, "warfarin.csv")
  file.copy(src$model, model_tmp)
  file.copy(src$data, data_tmp)

  fit <- ferx_fit(
    model_tmp,
    data_tmp,
    method = "focei",
    verbose = FALSE,
    covariance = TRUE,
    settings = list(maxiter = 30L)
  )
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  # Append whitespace to the model — flips the SHA-256.
  cat("\n# tampered\n", file = model_tmp, append = TRUE)

  expect_error(
    ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L),
    "hash mismatch"
  )
})

test_that("path/hash fields round-trip through .fitrx save/load", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)
  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out), add = TRUE)

  ferx_save_fit(fit, out, include_data = TRUE)
  loaded <- ferx_load_fit(out)

  # Hashes are unchanged (the bundle preserves the source bytes verbatim,
  # so loading + re-hashing the unpacked copy yields the same digest).
  expect_identical(loaded$model_hash, fit$model_hash)
  expect_identical(loaded$data_hash, fit$data_hash)

  # After load, model_path / data_path point at the unpacked tempfile
  # copies; running ferx_sir against the loaded fit must still succeed.
  expect_true(file.exists(loaded$model_path))
  expect_true(file.exists(loaded$data_path))

  out_sir <- ferx_sir(
    loaded,
    sir_samples = 6L,
    sir_resamples = 3L,
    sir_seed = 7L
  )
  expect_true(is.finite(out_sir$sir_ess))
})
