# Fixed-effects-only (naive-pooled) fits - no random effects at all.
#
# ferx-core #989 made an empty Omega legal on a continuous endpoint; ferx-r #290
# is the R-side half, where ferx_covariance() and ferx_sir() rejected any fit
# whose ebe_etas table had no ETA columns. At n_eta = 0 `build_ebe_etas()`
# returns NULL outright, so the guard those functions actually tripped was the
# "fit$ebe_etas is empty" one, not the "no ETA columns" one #290 names.
#
# Values are anchored against a NONMEM 7.6.0 run of the same model written as
# `$OMEGA 0 FIX` on a degenerate ETA (METHOD=0 SIGDIGITS=4). Simply *omitting*
# `$OMEGA` is not equivalent: NM-TRAN then infers a single-subject analysis,
# drops the ID grouping, and aborts on the per-subject TIME restarts.

# NONMEM reference (.ext final row / -1000000001 SE row).
NM_OFV     <- -269.63700440
NM_TVCL    <- 4.84070
NM_TVV     <- 52.8324
NM_SE_TVCL <- 0.166298   # RSR sandwich - NONMEM's $COVARIANCE default
NM_SE_TVV  <- 1.76021

pooled_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("one_cpt_iv_pooled")
      # No `covariance` argument: the model file already sets `covariance = true`,
      # and passing FALSE here conflicts with it and emits an override warning.
      fit <<- ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                       settings = list(maxiter = 300L))
    }
    fit
  }
})

test_that("the bundled pooled example is complete and declares no random effects", {
  ex <- ferx_example("one_cpt_iv_pooled")
  expect_true(file.exists(ex$model))
  expect_true(file.exists(ex$data))
  # The whole point of the example: no `omega` line and no `exp(ETA_*)` term.
  # Strip comments first - the header explains what was removed and says
  # "exp(ETA_*)" in prose, which a naive grep over the raw file matches.
  src  <- readLines(ex$model, warn = FALSE)
  code <- sub("#.*$", "", src)
  expect_false(any(grepl("^\\s*omega\\s", code)))
  expect_false(any(grepl("exp\\s*\\(\\s*ETA", code)))
  # ...but sigma is still mandatory for a continuous endpoint.
  expect_true(any(grepl("^\\s*sigma\\s", code)))
})

test_that("a model with no random effects fits and reports a 0x0 Omega", {
  fit <- pooled_fit()

  # 0x0, not NULL and not a phantom 1x1 - `nrow()` must be 0 so the eta-metadata
  # guards downstream see a number rather than NA (ferx-r #271).
  expect_true(is.matrix(fit$omega))
  expect_identical(nrow(fit$omega), 0L)
  expect_identical(ncol(fit$omega), 0L)

  # Nothing downstream may invent a random effect.
  expect_null(fit$ebe_etas)
  expect_length(fit$shrinkage_eta, 0L)
  expect_true(is.finite(fit$ofv))

  # Sigma is still doing all the work.
  expect_length(fit$sigma, 1L)
  expect_true(fit$sigma[[1]] > 0)
})

test_that("the pooled fit reproduces the NONMEM naive-pooled run", {
  fit <- pooled_fit()
  skip_if(!isTRUE(fit$converged), "pooled fit did not converge - skipping anchor")

  expect_equal(fit$ofv, NM_OFV, tolerance = 1e-4)
  expect_equal(unname(fit$theta[["TVCL"]]), NM_TVCL, tolerance = 1e-3)
  expect_equal(unname(fit$theta[["TVV"]]),  NM_TVV,  tolerance = 1e-3)
})

test_that("with no random effects PRED equals IPRED and CWRES equals IWRES", {
  fit <- pooled_fit()
  sd <- fit$sdtab
  expect_true(nrow(sd) > 0L)
  # Exactly, not approximately: any non-empty conditional step would move IPRED
  # off PRED. This is the cheapest end-to-end check that FOCEI really did
  # collapse rather than quietly carrying a phantom eta.
  expect_identical(sd$PRED, sd$IPRED)
  expect_identical(sd$CWRES, sd$IWRES)
  # No ETA columns in the diagnostic table either.
  expect_length(grep("^ETA", names(sd)), 0L)
})

test_that("ferx_covariance() accepts a fit with no random effects (#290)", {
  # The regression. Before #290 this stopped with
  # "fit$ebe_etas is empty; cannot warm-start the inner loop".
  fit <- pooled_fit()
  out <- ferx_covariance(fit, covariance_method = "rsr")
  skip_if(is.null(out$cov_matrix), "covariance step did not converge - skipping")

  expect_s3_class(out, "ferx_fit")
  expect_true(all(is.finite(out$se_theta)))
  expect_true(all(out$se_theta > 0))
  # No Omega means no Omega SEs - empty or absent, never a phantom entry.
  expect_true(is.null(out$se_omega) || length(out$se_omega) == 0L)

  # RSR is NONMEM's $COVARIANCE default. On a naive-pooled model this is not a
  # cosmetic choice: the model ignores within-subject correlation by
  # construction, so the sandwich runs about twice the naive inverse-Hessian.
  expect_equal(unname(out$se_theta[["TVCL"]]), NM_SE_TVCL, tolerance = 1e-2)
  expect_equal(unname(out$se_theta[["TVV"]]),  NM_SE_TVV,  tolerance = 1e-2)
})

test_that("ferx_covariance() default 'r' differs from 'rsr' on a pooled fit", {
  # Guards the anchor above against being satisfied by the wrong estimator: if
  # `covariance_method` were ignored, both calls would agree.
  fit <- pooled_fit()
  r   <- ferx_covariance(fit, covariance_method = "r")
  rsr <- ferx_covariance(fit, covariance_method = "rsr")
  skip_if(is.null(r$cov_matrix) || is.null(rsr$cov_matrix),
          "covariance step did not converge - skipping")
  expect_true(rsr$se_theta[["TVCL"]] > 1.5 * r$se_theta[["TVCL"]])
})

test_that("ferx_sir() accepts a fit with no random effects (#290)", {
  fit <- pooled_fit()
  cov <- ferx_covariance(fit, covariance_method = "r")
  skip_if(is.null(cov$cov_matrix), "covariance step did not converge - skipping")
  out <- ferx_sir(cov, sir_samples = 200L, sir_resamples = 50L, sir_seed = 1L)
  expect_s3_class(out, "ferx_fit")
})

test_that("print() omits the OMEGA section entirely at n_eta = 0", {
  # An empty header over an empty body reads as an estimation that failed rather
  # than one that was never requested; ferx-core's console output suppresses the
  # same section.
  fit <- pooled_fit()
  txt <- paste(capture.output(print(fit)), collapse = "\n")
  expect_false(grepl("OMEGA", txt, fixed = TRUE))
  # The sections that do apply must still be there.
  expect_true(grepl("SIGMA", txt, fixed = TRUE))
})

test_that("the estimates table reports theta and sigma but no omega rows", {
  # There is no exported `ferx_estimates()`; the table is computed internally by
  # `.ferx_compute_estimates()` and attached to the fit as `$estimates`. Its
  # column is `param`, not `parameter`.
  fit <- pooled_fit()
  est <- fit$estimates
  expect_true(is.data.frame(est))
  expect_true(all(c("TVCL", "TVV") %in% est$param))
  expect_length(grep("^(ETA|OMEGA)", est$param), 0L)
})
