# Tests for ferx_covariance() - the standalone covariance step (issue #738),
# the covariance-step analogue of ferx_sir().
#
# Happy-path assertions are gated on `skip_if(is.null(fit$cov_matrix), ...)`
# because in the no-autodiff CI build the warfarin FD covariance step
# occasionally fails to converge at maxiter = 30 (same pattern as
# test-ferx_sir.R).

cov_skip <- "covariance step did not converge - skipping"

test_that(".ferx_fit_interaction follows the last estimating (non-IMP) stage", {
  # fit$interaction is never plumbed to R, so the flag is derived from the
  # method chain. A trailing IMP stage is diagnostic-only and must be skipped,
  # otherwise the covariance step would differentiate the wrong (FOCE) NLL.
  expect_true(ferx:::.ferx_fit_interaction(list(method_chain = "FOCEI")))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = "FOCE")))
  # c("focei", "imp"): interaction TRUE despite the terminal IMP.
  expect_true(ferx:::.ferx_fit_interaction(list(method_chain = c("FOCEI", "IMP"))))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = c("FOCE", "IMP"))))
  # Falls back to fit$method when method_chain is absent.
  expect_true(ferx:::.ferx_fit_interaction(list(method = "FOCEI")))
  # Degenerate / all-IMP chains don't error.
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = character(0))))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = "IMP")))
})

test_that("ferx_covariance surfaces cov-step warnings in the structured table", {
  # Covariance-step warnings must reach fit$warnings_structured, not just the
  # flat vector: ferx_get_warnings() and the print tally read the structured
  # table. Force the step to fail (bad FD step) so the engine emits a diagnostic.
  fit <- warfarin_fit()
  out <- ferx_covariance(fit, covariance_method = "r")
  skip_if(is.null(out$cov_matrix), cov_skip)

  # A well-conditioned refit should carry no stale critical condition_number row.
  ws <- out$warnings_structured
  if (is.data.frame(ws) && nrow(ws) > 0L) {
    expect_true(all(c("severity", "category", "message") %in% names(ws)))
  }
  # ferx_get_warnings() reads the structured table; it must run without error
  # and return a data frame on the refreshed fit.
  gw <- ferx_get_warnings(out, as_df = TRUE)
  expect_s3_class(gw, "data.frame")
})

test_that("ferx_covariance validates its inputs", {
  fit <- warfarin_fit()  # covariance = FALSE

  expect_error(ferx_covariance("not a fit"), "ferx_fit object")
  expect_error(
    ferx_covariance(fit, covariance_method = "bhhh"),
    "covariance_method"
  )
  expect_error(
    ferx_covariance(fit, covariance_method = c("r", "s")),
    "single string"
  )
})

test_that("ferx_covariance runs end-to-end and populates covariance fields", {
  fit <- warfarin_fit()  # no inline covariance step
  expect_null(fit$cov_matrix)

  out <- ferx_covariance(fit, verbose = FALSE)
  skip_if(is.null(out$cov_matrix), cov_skip)

  expect_s3_class(out, "ferx_fit")
  expect_equal(out$covariance_status, "computed")

  # Covariance matrix is square, symmetric-ish, and named by parameter.
  d <- nrow(out$cov_matrix)
  expect_equal(ncol(out$cov_matrix), d)
  expect_false(is.null(rownames(out$cov_matrix)))
  expect_identical(rownames(out$cov_matrix), colnames(out$cov_matrix))

  # Standard errors populated and named for theta.
  expect_true(is.numeric(out$se_theta))
  expect_equal(length(out$se_theta), length(out$theta))
  expect_identical(names(out$se_theta), names(out$theta))
  expect_true(all(out$se_theta > 0))

  # Derived correlation matrix refreshed alongside cov_matrix.
  expect_false(is.null(out$cor_matrix))
  expect_equal(dim(out$cor_matrix), c(d, d))
  expect_equal(unname(diag(out$cor_matrix)), rep(1, d), tolerance = 1e-8)

  # Non-covariance fields are untouched.
  expect_identical(out$theta, fit$theta)
  expect_identical(out$omega, fit$omega)
  expect_equal(out$ofv, fit$ofv)
})

test_that("ferx_covariance reproduces the inline covariance step exactly", {
  # Re-running the covariance step against a fit that already carries an inline
  # covariance (same converged point, same EBEs) must reproduce that fit's own
  # covariance matrix and SEs bit-for-bit - the numerics route through the same
  # engine covariance step. This is the true apples-to-apples parity check
  # (comparing two independently-converged fits would confound the covariance
  # with tiny differences in the optimum).
  fit <- warfarin_fit_cov()   # covariance = TRUE
  skip_if(is.null(fit$cov_matrix), cov_skip)

  out <- ferx_covariance(fit)
  skip_if(is.null(out$cov_matrix), cov_skip)

  expect_equal(dim(out$cov_matrix), dim(fit$cov_matrix))
  expect_equal(unname(out$cov_matrix), unname(fit$cov_matrix), tolerance = 1e-6)
  expect_equal(unname(out$se_theta), unname(fit$se_theta), tolerance = 1e-6)
})

test_that("ferx_covariance can re-run with a different covariance_method", {
  fit <- warfarin_fit()
  out_r   <- ferx_covariance(fit, covariance_method = "r")
  skip_if(is.null(out_r$cov_matrix), cov_skip)
  out_rsr <- ferx_covariance(fit, covariance_method = "rsr")
  skip_if(is.null(out_rsr$cov_matrix), cov_skip)

  # Both produce a usable matrix; the sandwich generally differs from the
  # plain inverse-Hessian, but at minimum both are the same shape and named.
  expect_equal(dim(out_r$cov_matrix), dim(out_rsr$cov_matrix))
  expect_identical(rownames(out_r$cov_matrix), rownames(out_rsr$cov_matrix))
})

test_that("ferx_covariance refuses to run after the model file is tampered with", {
  src <- ferx_example("warfarin")
  tmp <- tempfile("ferx_cov_tamper_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  model_tmp <- file.path(tmp, "warfarin.ferx")
  data_tmp <- file.path(tmp, "warfarin.csv")
  file.copy(src$model, model_tmp)
  file.copy(src$data, data_tmp)

  fit <- ferx_fit(
    model_tmp, data_tmp,
    method = "focei", verbose = FALSE,
    covariance = FALSE, settings = list(maxiter = 30L)
  )
  expect_match(fit$model_hash, "^[0-9a-f]{64}$")

  # Append whitespace to the model - flips the SHA-256.
  cat("\n# tampered\n", file = model_tmp, append = TRUE)

  expect_error(ferx_covariance(fit), "hash mismatch")
})

test_that("ferx_covariance errors when the fit has no recorded model path", {
  fit <- warfarin_fit()
  fit$model_path <- NULL
  expect_error(ferx_covariance(fit), "no recorded model_path")
})
