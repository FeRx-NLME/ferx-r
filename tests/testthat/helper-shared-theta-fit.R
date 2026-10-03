# Cached SAEM and IMPMAP fits on the bundled two_cpt_oral_cov example, whose
# covariate mu-referenced typical values TVCL and TVV1 share THETA_WT. Since
# ferx-core #1620 those fits take one joint M-step for the shared theta and
# carry an info note naming it; test-shared-theta-mstep.R pins both from R.
#
# The model's [fit_options] say `method = focei` and `covariance = true`, so
# the call-time overrides raise two R conflict warnings. Only those are
# muffled here; any other R warning still reaches the test.
two_cpt_oral_cov_fit <- function(method, settings) {
  ex <- ferx_example("two_cpt_oral_cov")
  withCallingHandlers(
    ferx_fit(
      ex$model, ex$data,
      method = method, verbose = FALSE, covariance = FALSE,
      settings = settings
    ),
    warning = function(w) {
      if (grepl("overrides it with", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Default SAEM schedule: about 0.2 s, and the estimate tests rely on it.
two_cpt_oral_cov_saem_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- two_cpt_oral_cov_fit("saem", list(seed = 42L))
    }
    fit
  }
})

# Three IMPMAP iterations (under 0.1 s) already carry the note. The default
# schedule takes ~40 s, and from R its estimates barely move with #1620
# (THETA_CRCL 0.552 before, 0.564 after), so only the note is tested on IMPMAP.
two_cpt_oral_cov_impmap_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- two_cpt_oral_cov_fit(
        "impmap", list(impmap_iterations = 3L, impmap_seed = 42L)
      )
    }
    fit
  }
})
