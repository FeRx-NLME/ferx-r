# Shared cached FOCEI fit on the bundled warfarin example.
#
# testthat sources `helper-*.R` once at the start of `test_dir()`, so the
# closure below is materialised on first call and reused by every test in
# every file that calls `warfarin_fit()`. This avoids re-running a real
# estimation for each test file that needs a live `ferx_fit` object.
#
# `maxiter = 30L` is enough for the warfarin example to settle to a sensible
# point — tests using this helper check output shape, not convergence
# quality, so a low cap keeps CI fast.
warfarin_fit_cov <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex  <- ferx_example("warfarin")
      fit <<- ferx_fit(ex$model, ex$data,
                       method = "focei", verbose = FALSE,
                       covariance = TRUE, settings = list(maxiter = 30L))
    }
    fit
  }
})

warfarin_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "focei",
        verbose = FALSE,
        covariance = FALSE,
        settings = list(maxiter = 30L)
      )
    }
    fit
  }
})
