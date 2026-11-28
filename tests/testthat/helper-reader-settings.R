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
  fitrx_edit_fit(fit, function(w) {
    w$reader_settings <- NULL
    w$population_fingerprint <- NULL
    w
  })
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

# A fit that records a population fingerprint but no reader settings: the shape
# ferx-core's in-memory `fit()` writes (#1941), and the only shape besides a
# hand edit that still reaches the legacy guard (#530).
rs_fingerprint_only <- function(fit) {
  fitrx_edit_fit(fit, function(w) {
    w$reader_settings <- NULL
    w
  })
}

# warfarin with a model file stating `accept = TIME < 100`, fitted with a
# record-only `ignore = "TIME > 110"`: every row the call's clause matches the
# file's clause drops anyway, so the population is the file's, but the reader
# checks `ignore` before `accept` and credits the TIME = 120 rows to the call's
# clause (#530 gap 2).
rs_accept_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      model <- tempfile(fileext = ".ferx")
      writeLines(c(readLines(ex$model), "", "[data_selection]",
                   "  accept = TIME < 100"), model)
      fit <<- suppressWarnings(ferx_fit(
        model, ex$data, method = "focei", verbose = FALSE, covariance = TRUE,
        ignore = "TIME > 110", settings = list(maxiter = 30L)
      ))
    }
    fit
  }
})

# --- #793: from-fit readers' `selection = "fit"` ------------------------------

# warfarin_iov (doses at TIME 0 and 120) fitted with its TIME 120 doses dropped
# by a record-only `ignore =`: a selection that removes doses and no
# observation, so "all" and "fit" read the same number of rows and differ only
# in the dose the fit never saw (#793 plan, M3).
rs_mid_dose_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_iov")
      fit <<- suppressWarnings(ferx_fit(
        ex$model, ex$data, verbose = FALSE, covariance = FALSE,
        ignore = "EVID == 1 && TIME > 0", settings = list(maxiter = 5L)
      ))
    }
    fit
  }
})

# The first five subject IDs of a data frame, in file order.
rs_first_ids <- function(rows) utils::head(unique(rows$ID), 5L)

# A TTE example fitted without its first five subjects, dropped by a
# record-only `ignore_ids =`. `data` (a frame) replaces the example's file;
# only the example's own file is cached (per example), so two frames never
# share a fit. The fit's own reading has the other subjects only.
rs_tte_subjects_fit <- local({
  cache <- list()
  fit_it <- function(model, data, rows) {
    suppressWarnings(ferx_fit(model, data, verbose = FALSE, covariance = FALSE,
                              ignore_ids = rs_first_ids(rows)))
  }
  function(name = "tte_exponential", data = NULL) {
    ex <- ferx_example(name)
    if (!is.null(data)) return(fit_it(ex$model, data, data))
    if (is.null(cache[[name]])) {
      cache[[name]] <<- fit_it(ex$model, ex$data, utils::read.csv(ex$data))
    }
    cache[[name]]
  }
})

# `path` (a CSV, or a frame written to one first) without the records `drop`
# selects (a logical over its rows), as a new CSV. The kept lines are the
# file's own bytes, so every number reads back exactly as the fit read it.
rs_csv_without <- function(path, drop, env = parent.frame()) {
  if (is.data.frame(path)) {
    frame <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
    utils::write.csv(path, frame, row.names = FALSE, quote = FALSE, na = ".")
    path <- frame
  }
  lines <- readLines(path)
  rows <- utils::read.csv(path, na.strings = ".")
  stopifnot(length(drop) == nrow(rows), length(lines) == nrow(rows) + 1L)
  out <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
  writeLines(lines[c(TRUE, !drop)], out)
  out
}

# Every replicate of `out` (grouped by DRAW and SIM) has `fit$sdtab`'s rows:
# its ID and TIME columns, in order.
expect_sdtab_rows <- function(out, fit) {
  groups <- split(out[c("ID", "TIME")],
                  interaction(out$DRAW, out$SIM, drop = TRUE))
  testthat::expect_gt(length(groups), 0L)
  for (g in groups) {
    testthat::expect_identical(g$ID, as.character(fit$sdtab$ID))
    testthat::expect_identical(g$TIME, as.numeric(fit$sdtab$TIME))
  }
}
