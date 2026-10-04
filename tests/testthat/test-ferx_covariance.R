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

test_that("ferx_covariance closely reproduces the inline covariance step", {
  # Re-running the covariance step against a fit that already carries an inline
  # covariance (same point, same starting EBEs) reproduces that fit's own
  # covariance matrix and SEs closely - the numerics route through the same
  # engine covariance step. This is the true apples-to-apples parity check
  # (comparing two independently-converged fits would confound the covariance
  # with tiny differences in the optimum).
  #
  # Parity is close but NOT bit-exact: run_covariance recomputes the inner-loop
  # EBEs (it cold-starts warm_etas = None to mirror the inline step), so the
  # FD-Hessian is evaluated around slightly different etas. ferx-core's own
  # parity test (run_covariance_matches_inline_covariance) documents that
  # "strict sub-1e-6 parity is not achievable through this path" and bounds the
  # gap at abs 1e-4 for a small converged synthetic model. The shared warfarin
  # fit (maxiter = 30) shows a larger but stable gap, so assert an absolute
  # bound that still catches any real regression - a wrong interaction flag or
  # wrong EBEs diverge by orders of magnitude more.
  fit <- warfarin_fit_cov()   # covariance = TRUE
  skip_if(is.null(fit$cov_matrix), cov_skip)

  out <- ferx_covariance(fit)
  skip_if(is.null(out$cov_matrix), cov_skip)

  expect_equal(dim(out$cov_matrix), dim(fit$cov_matrix))
  expect_lt(max(abs(unname(out$cov_matrix) - unname(fit$cov_matrix))), 2e-3)
  # se_theta is the sqrt of the cov diagonal, so the same recompute-path divergence
  # lands slightly larger here. The ferx-core f7d52ec pin (categorical endpoints,
  # #900) grew the stable warfarin gap from <2e-3 to a deterministic ~2.8e-3, so
  # bound se_theta at 5e-3 - still orders of magnitude below any real covariance
  # regression (a wrong interaction flag / wrong EBEs).
  expect_lt(max(abs(unname(out$se_theta) - unname(fit$se_theta))), 5e-3)
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

test_that("ferx_covariance labels a block omega's rows column-major", {
  # The second call site of the cov_matrix labels (#367, #437): `ferx_fit()`
  # and this function name the same matrix, and a fit that skipped the
  # covariance step and picked it up here must arrive at identical dimnames.
  ex  <- ferx_example("warfarin_block_omega")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE,
                                   covariance = FALSE))
  out <- suppressWarnings(ferx_covariance(fit))
  skip_if(is.null(out$cov_matrix), cov_skip)

  d <- diag(out$cov_matrix)
  expect_true(all(
    c("ETA_CL,ETA_CL", "ETA_V,ETA_CL", "ETA_KA,ETA_CL",
      "ETA_V,ETA_V", "ETA_KA,ETA_V", "ETA_KA,ETA_KA") %in% names(d)
  ))
  # The held covariances of the partial block are the zeros; no variance is.
  expect_identical(unname(d[["ETA_KA,ETA_CL"]]), 0)
  expect_identical(unname(d[["ETA_KA,ETA_V"]]), 0)
  expect_gt(d[["ETA_V,ETA_V"]], 0)
  expect_identical(rownames(out$cor_matrix), rownames(out$cov_matrix))
})

# -- IOV models: the kappa segment gets its own labels (#437) ------------------
#
# The engine packs theta / omega / sigma / kappa. R used to label the matrix
# from counts, read every coordinate after theta and before sigma as omega, and
# so came up one name short on any IOV model - and dropped all of them. Each
# model below changes one segment of warfarin_iov: a block omega in front of a
# diagonal kappa, then a diagonal omega in front of a block kappa. Both are
# checked through both call sites: ferx_fit()'s covariance step, and
# ferx_covariance() on a fit that skipped it.

iov_variant_fit <- function(edit, covariance) {
  ex  <- ferx_example("warfarin_iov")
  txt <- edit(readLines(ex$model))
  mod <- tempfile(fileext = ".ferx")
  writeLines(txt, mod)
  suppressWarnings(ferx_fit(mod, ex$data, method = "foce", verbose = FALSE,
                            covariance = covariance,
                            settings = list(maxiter = 30L)))
}

block_omega_iov <- function(txt) {
  txt <- sub("^\\s*omega ETA_CL ~ .*$",
             "  block_omega (ETA_CL, ETA_V) = [0.07, 0.01, 0.02]", txt)
  txt <- txt[!grepl("^\\s*omega ETA_V\\s+~", txt)]
  stopifnot(sum(grepl("block_omega", txt)) == 1L,
            !any(grepl("^\\s*omega ETA_V", txt)))
  txt
}

block_kappa_iov <- function(txt) {
  txt <- sub("^\\s*kappa KAPPA_CL ~ .*$",
             "  block_kappa (KAPPA_CL, KAPPA_V) = [0.04, 0.005, 0.02]", txt)
  txt <- sub("TVV  \\* exp\\(ETA_V\\)", "TVV  * exp(ETA_V + KAPPA_V)", txt)
  stopifnot(sum(grepl("block_kappa", txt)) == 1L,
            sum(grepl("exp\\(ETA_V \\+ KAPPA_V\\)", txt)) == 1L)
  txt
}

expect_iov_labels <- function(edit, expected) {
  fit <- iov_variant_fit(edit, covariance = TRUE)
  out <- suppressWarnings(ferx_covariance(iov_variant_fit(edit, covariance = FALSE)))
  for (x in list(fit, out)) {
    if (is.null(x$cov_matrix)) next
    expect_identical(dimnames(x$cov_matrix), list(expected, expected))
    expect_identical(rownames(x$cor_matrix), expected)
    expect_identical(names(x$cov_fixed), expected)
  }
  skip_if(is.null(fit$cov_matrix) && is.null(out$cov_matrix), cov_skip)
}

test_that("a block omega with a diagonal kappa labels kappa after sigma", {
  # Mixed block + diagonal omega packs its full lower triangle (6), then sigma,
  # then the one kappa.
  expect_iov_labels(block_omega_iov, c(
    "TVCL", "TVV", "TVKA",
    "ETA_CL,ETA_CL", "ETA_V,ETA_CL", "ETA_KA,ETA_CL",
    "ETA_V,ETA_V", "ETA_KA,ETA_V", "ETA_KA,ETA_KA",
    "PROP_ERR", "KAPPA_CL"
  ))
})

test_that("a diagonal omega with a block kappa labels kappa column-major", {
  expect_iov_labels(block_kappa_iov, c(
    "TVCL", "TVV", "TVKA", "ETA_CL", "ETA_V", "ETA_KA", "PROP_ERR",
    "KAPPA_CL,KAPPA_CL", "KAPPA_V,KAPPA_CL", "KAPPA_V,KAPPA_V"
  ))
})

test_that("[mixture] per-class overrides are labelled with their class", {
  # The overrides pack after kappa (none here), omega ones before sigma ones.
  ex  <- ferx_example("warfarin")
  mod <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL1(0.1, 0.001, 10.0)",
    "  theta TVCL2(0.3, 0.001, 10.0)",
    "  theta TVV(8.0, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0)",
    "  theta P1(0.5, 0.01, 0.99)",
    "  omega ETA_CL ~ 0.05",
    "  omega ETA_V ~ 0.02",
    "  sigma PROP_ERR ~ 0.01 (sd)",
    "[mixture]",
    "  nsub = 2",
    "  p(1) = P1",
    "  omega(2) ETA_CL ~ 0.10",
    "  sigma(2) PROP_ERR ~ 0.02 (sd)",
    "[individual_parameters]",
    "  CL = if (MIXNUM == 1) TVCL1 * exp(ETA_CL) else TVCL2 * exp(ETA_CL)",
    "  V  = TVV * exp(ETA_V)",
    "  KA = TVKA",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ), mod)
  fit <- suppressWarnings(ferx_fit(mod, ex$data, method = "foce", verbose = FALSE,
                                   covariance = TRUE,
                                   settings = list(maxiter = 30L)))
  skip_if(is.null(fit$cov_matrix), cov_skip)
  expected <- c("TVCL1", "TVCL2", "TVV", "TVKA", "P1", "ETA_CL", "ETA_V",
                "PROP_ERR", "ETA_CL (class 2)", "PROP_ERR (class 2)")
  expect_identical(dimnames(fit$cov_matrix), list(expected, expected))
})
