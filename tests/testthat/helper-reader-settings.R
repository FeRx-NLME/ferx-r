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

# warfarin_iov fitted with its occasion column passed only through `settings =`
# (`VISIT`, a copy of the file's `OCC`): the reader settings record `VISIT`, the
# model file says `OCC`.
rs_iov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_iov")
      data <- tempfile(fileext = ".csv")
      rows <- utils::read.csv(ex$data)
      rows$VISIT <- rows$OCC
      utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
      fit <<- suppressWarnings(ferx_fit(
        ex$model, data, method = "foce", verbose = FALSE, covariance = FALSE,
        settings = list(maxiter = 2L, iov_column = "VISIT")
      ))
    }
    fit
  }
})

# warfarin with a second dose per subject after the last observation, dropped
# by a record-only `ignore =`: the dose rows are not in the sdtab, so npde's
# row/ID/TIME alignment cannot see the re-read put them back (#526 review 1).
rs_late_dose_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      rows <- utils::read.csv(ex$data)
      late <- rows[rows$EVID == 1, ]
      late$TIME <- 500
      rows <- rbind(rows, late)
      rows <- rows[order(rows$ID, rows$TIME), ]
      data <- tempfile(fileext = ".csv")
      utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
      fit <<- suppressWarnings(ferx_fit(
        ex$model, data, method = "focei", verbose = FALSE, covariance = FALSE,
        ignore = "EVID == 1 && TIME > 400", settings = list(maxiter = 5L)
      ))
    }
    fit
  }
})

# A copy of `fit`'s model file with a comment appended: the fit's hash no longer
# matches it. `fit$model_path` is pointed at the copy.
rs_edited_model <- function(fit, env = parent.frame()) {
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(readLines(fit$model_path), "# edited after the fit"), model)
  fit$model_path <- model
  fit
}
