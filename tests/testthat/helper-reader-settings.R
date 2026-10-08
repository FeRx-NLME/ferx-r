# Shared fixtures for the "post-hoc steps use the fit's own data selection"
# tests (#462 / #416): ferx_covariance(), ferx_sir() and ferx_calc_npde() on a
# fit made with a record-only `ignore =` / `accept =` / `ignore_ids =`.
#
# The oracle is the fit's own in-fit step: the quantity is "the same rows as
# the fit", and the in-fit covariance / SIR is that by definition. Compared at
# 1e-8 relative. Before #462 the standalone step re-read the data with the
# model file's selection only and landed 2.11e-05 (`ignore`) / 2.09e-04
# (`accept`) away, silently; with the fit's selection it lands within ~4e-11.

# Two observation records (IDs 7 and 8), no dose and no whole subject.
rs_ignore <- "EVID == 0 && DV < 1.0"

# Cached warfarin fits keyed by their selection arguments; covariance and
# in-fit SIR on, so every comparison has its in-fit oracle.
rs_fit <- local({
  cache <- list()
  function(..., sir = FALSE) {
    args <- list(...)
    key <- paste(deparse(c(args, sir = sir)), collapse = "")
    if (is.null(cache[[key]])) {
      ex <- ferx_example("warfarin")
      settings <- list(maxiter = 30L)
      if (sir) {
        settings <- c(settings, list(sir_samples = 200L, sir_resamples = 100L,
                                     sir_seed = 7L))
      }
      cache[[key]] <<- suppressWarnings(do.call(ferx_fit, c(
        list(ex$model, ex$data, method = "focei", verbose = FALSE,
             covariance = TRUE, sir = sir, settings = settings),
        args
      )))
    }
    cache[[key]]
  }
})

# max |a - b| / |b|: the "standalone equals in-fit" measure.
rs_rel <- function(a, b) max(abs(unname(a) - unname(b)) / abs(unname(b)))

# A fit as one made before #462 recorded its selection: neither field.
rs_legacy <- function(fit) {
  fit$reader_settings <- NULL
  fit$population_fingerprint <- NULL
  fit
}

# The three post-hoc steps on `fit`, each returning its value or the condition
# it raised. The hash-missing warning a hand-edited fit triggers is muffled.
rs_steps <- function(fit) {
  grab <- function(expr) {
    tryCatch(suppressWarnings(expr), error = function(e) e)
  }
  list(
    covariance = grab(ferx_covariance(fit)),
    sir = grab(ferx_sir(fit, sir_samples = 50L, sir_resamples = 20L,
                        sir_seed = 1L)),
    npde = grab(ferx_calc_npde(fit, nsim = 20L, seed = 5L))
  )
}
