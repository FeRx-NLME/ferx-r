#' Non-parametric case bootstrap
#'
#' Resample whole subjects with replacement, refit the same model to each
#' replicate, and summarise the spread of the estimates into bias, standard
#' errors and confidence intervals. Feature target: PsN's `bootstrap`.
#'
#' The bootstrap runs in the engine (ferx-core's `ferx-tools` crate); this
#' function is a thin wrapper. The percentile estimator, the exclusion filters,
#' the parameter naming and the `sample_size` syntax all live there, so an R
#' session and the `ferx bootstrap` CLI produce the same numbers from the same
#' seed.
#'
#' ## What the bootstrap buys over the covariance step
#'
#' ferx's standard errors are the asymptotic `R^-1` matrix, symmetric by
#' construction. The bootstrap distribution assumes neither normality nor
#' symmetry: it survives a failed or non-positive-definite covariance step, and
#' it shows asymmetry in poorly identified parameters instead of averaging it
#' away. Both intervals are returned side by side (`ci_lower` / `ci_upper` for
#' the percentile interval, `ci_lower_normal` / `ci_upper_normal` for the
#' normal approximation built from the bootstrap standard error) precisely so
#' the disagreement between them is visible.
#'
#' ## Exclusions
#'
#' A replicate whose minimization terminated, or whose estimate sits on a
#' boundary, is not a draw from the sampling distribution of a converged fit,
#' so by default it is excluded from the statistics (PsN's defaults). The
#' replicate is still fitted and still appears in `$raw`; only the summary
#' drops it. `$diagnostics` reports how many were dropped and why. To change
#' the criteria after a run has finished, use [ferx_bootstrap_summarize()]
#' instead of refitting - but that needs `directory` to have been set.
#'
#' ## Cost
#'
#' This is `samples` (+ 1) complete fits. At PsN's rule-of-thumb 200 samples a
#' model that takes a minute to fit takes over three hours, so `threads` is the
#' argument that matters most. `progress` (on by default in an interactive
#' session) draws a bar as the replicates come in, so a long run says where it
#' is. Note that the run is **not** interruptible with Ctrl-C once the
#' replicates are underway - the bar advancing is not a Ctrl-C handler, and a
#' Ctrl-C pressed while the bar is up is consumed rather than queued, so it does
#' not abort the run when the fits finish either. Stopping a run means killing
#' the session, so set `directory` before starting a long one and re-run with
#' `resume = TRUE`, which refits only the replicates that directory is missing.
#'
#' ## Resuming an interrupted run
#'
#' `resume = TRUE` continues the run in `directory`, refitting only the sample
#' indices its `raw_results.csv` does not already carry. It is sound because a
#' replicate's draw is a pure function of `(seed, index)`: a reused replicate is
#' bit-for-bit the one an uninterrupted run would have produced, so a resumed
#' run's artefacts are identical to that run's. The base fit is reused too and
#' is deliberately *not* refitted - a second fit can land on a slightly
#' different optimum, and with `update_inits = TRUE` that would start the new
#' replicates from a different point than the ones already on disk.
#'
#' The engine refuses to resume from a directory that belongs to a different
#' run: the model and data hashes, the parameter names and the settings that
#' shape the replicates are recorded alongside them and are checked before any
#' is reused. Those settings are `samples`, `seed`, `sample_size`,
#' `stratify_on`, `run_base_model`, `update_inits`, `keep_covariance` and
#' `dofv` - all of them pinned, so a resumed run cannot *extend* an earlier
#' one. Calling again with a larger `samples` is refused rather than topping
#' the run up; a bigger bootstrap means a fresh `directory`.
#'
#' A replicate whose fit *errored* is carried forward as a failure rather than
#' refitted, matching PsN - a fit that failed usually fails again.
#' `retry_failed = TRUE` (which needs `resume = TRUE`) refits those instead, and
#' is for a failure that was transient - an out-of-memory kill, a full disk -
#' rather than one of the model.
#'
#' @param model Path to a `.ferx` model file, or a `ferx_model` object from
#'   [ferx_model()].
#' @param data Path to a NONMEM-format CSV. `NULL` (the default) uses the
#'   model's `[data]` block.
#' @param samples Number of bootstrap datasets. Default 200, PsN's default and
#'   its stated rule of thumb for standard errors.
#' @param seed Master seed. Each replicate's draw is derived from this and its
#'   own index, so a run is reproducible and independent of `threads`.
#' @param threads Replicates to fit concurrently. `NULL` (the default) uses the
#'   engine default.
#' @param stratify_on Name of a data column defining the resampling strata, or
#'   `NULL` (the default) for an unstratified run. Stratifying on study or arm
#'   keeps each replicate's composition close to the original dataset's.
#' @param sample_size How many subjects each replicate draws. `NULL` (the
#'   default) draws as many as the original dataset has - the only case where a
#'   stratified proportional split is guaranteed to sum back to the request. A
#'   single number sets one total for the replicate; a **named** numeric vector
#'   (`c("1001" = 12, "1002" = 24)`) sets an explicit count per stratum and is
#'   only meaningful together with `stratify_on`. Every count must be a whole
#'   number >= 1. It is PsN's `-sample_size=1001=>12,1002=>24`, spelled the R
#'   way.
#' @param update_inits Start each replicate from the base fit's final estimates
#'   (PsN default, `TRUE`). Requires `run_base_model = TRUE`.
#' @param run_base_model Fit the original dataset first. Required for
#'   `update_inits`, for the bias column, and for the normal-approximation
#'   intervals. Default `TRUE`.
#' @param keep_covariance Run the covariance step for each replicate. Off by
#'   default: it is most of the cost, and the bootstrap standard error comes
#'   from the spread of the estimates, not from any per-replicate `R^-1`.
#'   Turning it on adds `se_<parameter>` columns to `$raw` and is what the two
#'   covariance-step exclusion filters read.
#' @param dofv Compute the change in objective function value: evaluate each
#'   replicate's parameter vector on the *original* dataset with no estimation
#'   (NONMEM `MAXEVAL=0`) and subtract the original fit's OFV. Roughly doubles
#'   the run time. Default `FALSE`.
#' @param skip_minimization_terminated Exclude replicates whose minimization
#'   terminated. Default `TRUE` (PsN's default).
#' @param skip_estimate_near_boundary Exclude replicates whose estimate sits on
#'   a boundary. Default `TRUE` (PsN's default).
#' @param skip_covariance_step_terminated,skip_with_covstep_warnings Exclude
#'   replicates on the covariance step's outcome. Both default `FALSE`, and
#'   both require `keep_covariance = TRUE` - asking for one without it is an
#'   error rather than a silently dropped filter.
#' @param ci Two-sided confidence level, in percent. Default 95.
#' @param directory Where to write the CSV artefacts (`raw_results.csv`,
#'   `bootstrap_results.csv`, `bootstrap_diagnostics.csv`, the individual and
#'   key files, and `delta_ofv.csv` when `dofv = TRUE`). `NULL`, the default,
#'   computes everything in memory and writes nothing - an R user has the data
#'   frames in hand and rarely wants eight files appearing in the working
#'   directory. The CLI defaults the other way. Set it if you want to be able
#'   to call [ferx_bootstrap_summarize()] later.
#' @param resume Continue an interrupted run in `directory` instead of starting
#'   a fresh one, refitting only the replicates that directory does not already
#'   hold. Default `FALSE`. Needs `directory`, and needs the arguments that
#'   shape the replicates to match the ones the directory was written with -
#'   `samples`, `seed`, `sample_size`, `stratify_on`, `run_base_model`,
#'   `update_inits`, `keep_covariance` and `dofv`, plus the model and the data
#'   themselves. A mismatch is an error, not a fresh run: resuming with a
#'   larger `samples` does not extend the earlier run.
#' @param retry_failed Refit the replicates a resumed run finds recorded as
#'   *failed*, instead of carrying the failure forward. Default `FALSE` (PsN's
#'   default). Needs `resume = TRUE`.
#' @param progress Show a progress bar while the replicates fit. Default
#'   `interactive()`, so a script or a knitted document prints nothing. Uses
#'   the cli package when it is installed and [utils::txtProgressBar()]
#'   otherwise. The bar is drawn by this R session while the fits run in the
#'   engine, and reports what the run will actually do: a `resume = TRUE` run
#'   whose `directory` already holds most of the replicates counts only the
#'   ones still to fit.
#' @param verbose Print a one-line run header to stderr. Default `FALSE`.
#'
#' @return An object of class `ferx_bootstrap`, a list with:
#'   \describe{
#'     \item{parameters}{Data frame, one row per parameter: `parameter`,
#'       `original`, `mean`, `bias`, `standard_error`, `median`, `ci_lower`,
#'       `ci_upper`, `ci_lower_normal`, `ci_upper_normal`. This is
#'       `bootstrap_results.csv`.}
#'     \item{raw}{Data frame, one row per fit with the original dataset first
#'       (`sample = 0`), carrying the estimates *and* the termination
#'       diagnostics the exclusion filters read. This is `raw_results.csv`, and
#'       it is what makes a custom R-side filter or re-summarisation possible.}
#'     \item{diagnostics}{Data frame of `statistic` / `value`: the run counts,
#'       the exclusion tallies (`excluded: ...`) and the diagnostic means
#'       (`mean: ...`).}
#'     \item{delta_ofv}{Data frame of `sample` / `delta_ofv` when
#'       `dofv = TRUE`, otherwise `NULL`.}
#'     \item{parameter_names, subject_ids, n_completed, n_included,
#'       chi_square_df, confidence_level, samples, seed, model, data,
#'       directory}{Run metadata. `chi_square_df` is the number of estimated
#'       parameters, the reference degrees of freedom for the `dofv` panel of
#'       [plot.ferx_bootstrap()].}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' bs <- ferx_bootstrap(ex$model, ex$data, samples = 200, seed = 12345,
#'                      threads = 8)
#' bs
#' bs$parameters
#' plot(bs)
#'
#' # Stratified, with artefacts on disk so the run can be re-summarised later
#' bs <- ferx_bootstrap(ex$model, ex$data, samples = 200,
#'                      stratify_on = "SEX", directory = "warfarin-bootstrap")
#' ferx_bootstrap_summarize("warfarin-bootstrap",
#'                          skip_estimate_near_boundary = FALSE)
#'
#' # Continue that run after it was interrupted: only the replicates missing
#' # from the directory are fitted again.
#' bs <- ferx_bootstrap(ex$model, ex$data, samples = 200,
#'                      stratify_on = "SEX", directory = "warfarin-bootstrap",
#'                      resume = TRUE)
#' }
#'
#' @seealso [ferx_bootstrap_summarize()] to change the exclusion criteria on a
#'   finished run, [ferx_sir()] and [ferx_covariance()] for the parametric
#'   uncertainty routes.
#' @family fitting
#' @export
ferx_bootstrap <- function(model,
                           data = NULL,
                           samples = 200,
                           seed = 1,
                           threads = NULL,
                           stratify_on = NULL,
                           sample_size = NULL,
                           update_inits = TRUE,
                           run_base_model = TRUE,
                           keep_covariance = FALSE,
                           dofv = FALSE,
                           skip_minimization_terminated = TRUE,
                           skip_estimate_near_boundary = TRUE,
                           skip_covariance_step_terminated = FALSE,
                           skip_with_covstep_warnings = FALSE,
                           ci = 95,
                           directory = NULL,
                           resume = FALSE,
                           retry_failed = FALSE,
                           progress = interactive(),
                           verbose = FALSE) {
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model <- model$model
  }
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("`model` must be a path to a .ferx file or a ferx_model object.")
  }
  if (!file.exists(model)) {
    stop("Model file not found: ", model)
  }
  if (!is.null(data)) {
    if (!is.character(data) || length(data) != 1L || is.na(data)) {
      stop("`data` must be a single path to a NONMEM-format CSV, or NULL.")
    }
    if (!file.exists(data)) stop("Data file not found: ", data)
  }
  .ferx_check_count(samples, "samples", min = 1L)
  .ferx_check_count(seed, "seed", min = 0L)
  if (!is.null(threads)) .ferx_check_count(threads, "threads", min = 1L)
  if (!is.numeric(ci) || length(ci) != 1L || is.na(ci) || ci <= 0 || ci >= 100) {
    stop("`ci` must be a confidence level in (0, 100), e.g. 95.")
  }
  for (nm in c("update_inits", "run_base_model", "keep_covariance", "dofv",
               "skip_minimization_terminated", "skip_estimate_near_boundary",
               "skip_covariance_step_terminated", "skip_with_covstep_warnings",
               "resume", "retry_failed", "progress", "verbose")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("`", nm, "` must be TRUE or FALSE.")
    }
  }
  if (update_inits && !run_base_model) {
    stop("`update_inits = TRUE` needs `run_base_model = TRUE`: the replicates ",
         "start from the base fit's final estimates.")
  }
  # The engine checks these too, and its messages are good ones - but they name
  # the CLI flags. Checking here as well means the R user is told about
  # `directory =` and `resume =` before paying for a model compile.
  if (resume && is.null(directory)) {
    stop("`resume = TRUE` continues a run in a directory, so it needs ",
         "`directory` naming one an earlier run wrote.")
  }
  if (resume && !dir.exists(directory)) {
    stop("`resume = TRUE` needs a directory an earlier run wrote, and there ",
         "is nothing at `", directory, "`. Drop `resume` to start a fresh run ",
         "there, or check the path.")
  }
  if (retry_failed && !resume) {
    stop("`retry_failed = TRUE` refits replicates a previous run recorded as ",
         "failed, so it only means anything with `resume = TRUE`.")
  }
  if (!is.null(stratify_on) &&
      (!is.character(stratify_on) || length(stratify_on) != 1L ||
       is.na(stratify_on) || !nzchar(stratify_on))) {
    stop("`stratify_on` must be a single column name, or NULL.")
  }
  ss <- .ferx_bootstrap_sample_size(sample_size)

  # An unstratified run with per-stratum counts is a user error the engine
  # cannot diagnose: it sees only the assembled string and one nameless stratum.
  if (length(ss$keys) > 0L && is.null(stratify_on)) {
    stop("`sample_size` with per-stratum counts needs `stratify_on`.")
  }

  dir_arg <- if (is.null(directory)) "" else normalizePath(directory, mustWork = FALSE)

  res <- ferx_rust_bootstrap(
    normalizePath(model),
    if (is.null(data)) "" else normalizePath(data),
    as.integer(samples),
    as.numeric(seed),
    ss$keys,
    ss$values,
    stratify_on %||% "",
    update_inits,
    run_base_model,
    keep_covariance,
    if (is.null(threads)) 0L else as.integer(threads),
    skip_minimization_terminated,
    skip_estimate_near_boundary,
    skip_covariance_step_terminated,
    skip_with_covstep_warnings,
    dofv,
    dir_arg,
    as.numeric(ci),
    resume,
    retry_failed,
    verbose,
    # NULL is what tells the engine not to open a reporting thread at all.
    if (isTRUE(progress)) .ferx_bootstrap_progress_handler() else NULL
  )
  if (length(res) == 0L) {
    stop("ferx_bootstrap: backend returned no result.")
  }

  res$raw <- .ferx_blank_error_to_na(res$raw)
  res$samples <- as.integer(samples)
  res$seed <- as.numeric(seed)
  res$model <- normalizePath(model, mustWork = FALSE)
  res$data <- if (is.null(res$data_path) || !nzchar(res$data_path)) {
    NA_character_
  } else {
    normalizePath(res$data_path, mustWork = FALSE)
  }
  res$data_path <- NULL
  # Re-normalise now that the engine has created it: before the run the
  # path may not exist yet, and `normalizePath` leaves an unresolved path
  # (a macOS /var vs /private/var symlink, say) untouched.
  res$directory <- if (nzchar(dir_arg)) {
    normalizePath(dir_arg, mustWork = FALSE)
  } else {
    NULL
  }
  class(res) <- c("ferx_bootstrap", "list")
  res
}

#' Re-summarise a finished bootstrap run
#'
#' Recompute the statistics of a bootstrap run from its `raw_results.csv` under
#' a different set of exclusion criteria - PsN's `-summarize`. Refits nothing:
#' the estimates are already on disk, and so are the termination diagnostics
#' the filters read, so relaxing (say) `skip_minimization_terminated` is a
#' re-read rather than a re-run.
#'
#' This is the recovery path for a run where too many replicates were filtered
#' out to resolve a percentile interval. It rewrites `bootstrap_results.csv`
#' and `bootstrap_diagnostics.csv` in `directory` in place; the per-replicate
#' files the original run wrote are still valid and are left alone.
#'
#' Needs a run that was given a `directory`. [ferx_bootstrap()] with the default
#' `directory = NULL` computes in memory and writes nothing, so there is nothing
#' to re-summarise.
#'
#' @param directory Path to a finished bootstrap run directory - one containing
#'   `raw_results.csv`.
#' @inheritParams ferx_bootstrap
#'
#' @return An object of class `ferx_bootstrap` with `$parameters`,
#'   `$diagnostics`, `$raw` (read back from `raw_results.csv`), `$samples`,
#'   `$n_completed`, `$n_included`, `$confidence_level` and `$directory`.
#'   `$delta_ofv` is populated when the original run wrote a `delta_ofv.csv`.
#'   `$model`, `$data` and `$seed` are *not* set: a finished run directory does
#'   not record the model path, the data path or the seed, so the printed
#'   header omits those lines.
#'
#' @examples
#' \dontrun{
#' bs <- ferx_bootstrap(ferx_example("warfarin")$model, samples = 200,
#'                      directory = "warfarin-bootstrap")
#'
#' # Keep the boundary cases after all, without refitting
#' relaxed <- ferx_bootstrap_summarize("warfarin-bootstrap",
#'                                     skip_estimate_near_boundary = FALSE)
#' relaxed$parameters
#' }
#'
#' @seealso [ferx_bootstrap()]
#' @family fitting
#' @export
ferx_bootstrap_summarize <- function(directory,
                                     skip_minimization_terminated = TRUE,
                                     skip_estimate_near_boundary = TRUE,
                                     skip_covariance_step_terminated = FALSE,
                                     skip_with_covstep_warnings = FALSE,
                                     ci = 95) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory)) {
    stop("`directory` must be a single path to a bootstrap run directory.")
  }
  if (!is.numeric(ci) || length(ci) != 1L || is.na(ci) || ci <= 0 || ci >= 100) {
    stop("`ci` must be a confidence level in (0, 100), e.g. 95.")
  }
  for (nm in c("skip_minimization_terminated", "skip_estimate_near_boundary",
               "skip_covariance_step_terminated", "skip_with_covstep_warnings")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("`", nm, "` must be TRUE or FALSE.")
    }
  }
  if (!dir.exists(directory)) {
    stop("Bootstrap directory not found: ", directory)
  }
  raw_path <- file.path(directory, "raw_results.csv")
  if (!file.exists(raw_path)) {
    stop("`", directory, "` is not a finished bootstrap run: no raw_results.csv. ",
         "Re-run ferx_bootstrap() with `directory` set.")
  }

  dir_abs <- normalizePath(directory)
  res <- ferx_rust_bootstrap_summarize(
    dir_abs,
    skip_minimization_terminated,
    skip_estimate_near_boundary,
    skip_covariance_step_terminated,
    skip_with_covstep_warnings,
    as.numeric(ci)
  )
  if (length(res) == 0L) {
    stop("ferx_bootstrap_summarize: backend returned no result.")
  }

  # The engine's `resummarize` returns the statistics only - it never re-reads
  # the estimates into memory for the caller. Reading the same file back here
  # is a convenience, not a second implementation: nothing is recomputed from
  # it, it is handed over as `$raw`.
  res$raw <- .ferx_blank_error_to_na(
    .ferx_raw_csv_types(
      read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
    )
  )
  dofv_path <- file.path(dir_abs, "delta_ofv.csv")
  res$delta_ofv <- if (file.exists(dofv_path)) {
    read.csv(dofv_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    NULL
  }
  # Everything else a `ferx_bootstrap()` result carries is a fact about the
  # invocation (model, data, seed) that the run directory does not record;
  # `samples` is the exception - the engine wrote it to the diagnostics.
  res$samples <- as.integer(.ferx_diag_value(res$diagnostics, "samples_requested"))
  res$directory <- dir_abs
  class(res) <- c("ferx_bootstrap", "list")
  res
}

#' Print a bootstrap result
#'
#' @param x A `ferx_bootstrap` object.
#' @param digits Significant digits for the parameter table. Default 4.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @export
print.ferx_bootstrap <- function(x, digits = 4, ...) {
  cat("ferx bootstrap\n")
  if (!is.null(x$model)) cat("  Model:  ", x$model, "\n", sep = "")
  if (!is.null(x$data) && !is.na(x$data)) cat("  Data:   ", x$data, "\n", sep = "")
  requested <- .ferx_diag_value(x$diagnostics, "samples_requested")
  cat(sprintf("  Samples: %s completed of %s requested; %s included\n",
              format(x$n_completed), format(requested), format(x$n_included)))

  excl <- x$diagnostics[startsWith(as.character(x$diagnostics$statistic), "excluded: "), , drop = FALSE]
  for (i in seq_len(nrow(excl))) {
    cat(sprintf("    excluded %s: %s\n",
                format(excl$value[i]),
                sub("^excluded: ", "", excl$statistic[i])))
  }
  if (!is.null(x$directory)) cat("  Wrote:  ", x$directory, "\n", sep = "")

  cat(sprintf("\n%g%% confidence intervals (percentile / normal approximation):\n",
              x$confidence_level))
  tbl <- x$parameters
  num <- vapply(tbl, is.numeric, logical(1))
  tbl[num] <- lapply(tbl[num], signif, digits = digits)
  print(tbl, row.names = FALSE)

  if (nrow(tbl) > 0L && all(is.na(tbl$ci_lower))) {
    cat(sprintf(paste0("\nNote: %s samples cannot resolve a %g%% percentile ",
                       "interval; only the normal-approximation interval is ",
                       "shown. PsN's rule of thumb is 200 samples.\n"),
                format(x$n_included), x$confidence_level))
  }
  invisible(x)
}

#' Plot a bootstrap parameter distribution
#'
#' One histogram per parameter of the bootstrap estimates - the panel PsN
#' generates through xpose4's `boot.hist` - with the original estimate drawn as
#' a solid line and the confidence interval as dashed lines. When the run used
#' `dofv = TRUE`, a final panel shows the distribution of the change in
#' objective function value against its chi-square reference density, on
#' `chi_square_df` degrees of freedom.
#'
#' The histograms show every replicate that produced an estimate, including any
#' the exclusion filters dropped from the statistics; the lines drawn over them
#' come from `$parameters` and so reflect the included set only. The counts in
#' `$diagnostics` say how far the two differ. (Deciding which replicate is
#' excluded is the engine's job - see [ferx_bootstrap_summarize()] to change
#' the criteria.)
#'
#' @param x A `ferx_bootstrap` object.
#' @param y Ignored. Present only for consistency with the `plot(x, y, ...)`
#'   generic; must be `NULL` if supplied.
#' @param parameters Character vector selecting which parameters to plot.
#'   `NULL` (the default) plots all of them.
#' @param breaks Passed to [hist()]. Default `"FD"`.
#' @param ... Ignored.
#'
#' @return Invisibly, the character vector of parameters plotted.
#'
#' @examples
#' \dontrun{
#' bs <- ferx_bootstrap(ferx_example("warfarin")$model, samples = 200)
#' plot(bs)
#' plot(bs, parameters = c("TVCL", "TVV"))
#' }
#'
#' @family diagnostics
#' @export
plot.ferx_bootstrap <- function(x, y = NULL, parameters = NULL,
                                breaks = "FD", ...) {
  if (!is.null(y)) {
    stop("plot.ferx_bootstrap() does not use `y`; did you mean `parameters = `?")
  }
  pars <- x$parameters$parameter
  if (!is.null(parameters)) {
    unknown <- setdiff(parameters, pars)
    if (length(unknown) > 0L) {
      stop("Unknown parameter(s): ", paste(unknown, collapse = ", "))
    }
    pars <- parameters
  }
  # `sample == 0` is the fit on the original dataset, not a draw.
  reps <- x$raw[x$raw$sample > 0L, , drop = FALSE]
  n_panel <- length(pars) + as.integer(!is.null(x$delta_ofv))
  if (n_panel == 0L) {
    stop("Nothing to plot: the bootstrap summary has no parameters.")
  }

  n_col <- min(3L, n_panel)
  n_row <- ceiling(n_panel / n_col)
  op <- par(mfrow = c(n_row, n_col), mar = c(4, 4, 2.5, 1))
  on.exit(par(op), add = TRUE)

  for (p in pars) {
    row <- x$parameters[match(p, x$parameters$parameter), ]
    vals <- if (p %in% names(reps)) reps[[p]] else numeric(0)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
      plot.new()
      title(main = p, sub = "no estimates")
      next
    }
    hist(vals, breaks = .ferx_hist_breaks(vals, breaks), col = "grey85",
         border = "white", main = p, xlab = "estimate")
    if (!is.na(row$original)) {
      abline(v = row$original, lwd = 2)
    }
    lo <- if (is.na(row$ci_lower)) row$ci_lower_normal else row$ci_lower
    hi <- if (is.na(row$ci_upper)) row$ci_upper_normal else row$ci_upper
    if (!is.na(lo)) abline(v = lo, lty = 2)
    if (!is.na(hi)) abline(v = hi, lty = 2)
  }

  if (!is.null(x$delta_ofv)) {
    d <- x$delta_ofv$delta_ofv
    d <- d[is.finite(d)]
    if (length(d) == 0L) {
      plot.new()
      title(main = "delta OFV", sub = "no values")
    } else {
      hist(d, breaks = .ferx_hist_breaks(d, breaks), col = "grey85",
           border = "white", freq = FALSE, main = "delta OFV",
           xlab = "delta OFV")
      df <- x$chi_square_df
      if (!is.null(df) && !is.na(df) && df > 0) {
        grid <- seq(max(1e-8, min(d)), max(d), length.out = 200)
        lines(grid, dchisq(grid, df = df), lwd = 2)
      }
    }
  }
  invisible(pars)
}

# -- internals --------------------------------------------------------------

# The progress handler the engine calls while a run proceeds.
#
# The engine reports one event per completed fit and this draws it. Events are
# coalesced on the way over - the engine redraws on a timer, so a fast run can
# skip straight from "started" to a replicate count - so every stage opens the
# bar it needs if it is not already open, rather than assuming the previous
# stage arrived. Nothing here can throw: an error in a progress bar must not
# lose a run that may have taken hours, and the caller drops it silently, so a
# broken bar would otherwise fail invisibly and repeatedly.
#
# `envir` is the frame the bar belongs to: cli ties a progress bar's lifetime to
# a calling frame, and the handler's own frame is gone before the next event.
.ferx_bootstrap_progress_handler <- function(envir = parent.frame()) {
  # The bars are drawn from the engine's poll loop, long after this factory has
  # returned, so a lazy `parent.frame()` resolves against a call stack that no
  # longer holds the fit - cli then falls back to the global environment, whose
  # bar nothing closes when a run errors out.
  force(envir)
  use_cli <- requireNamespace("cli", quietly = TRUE)
  state <- new.env(parent = emptyenv())
  state$kind <- NULL      # "base", "replicate" or "dofv"
  state$bar <- NULL
  state$total <- NA_integer_
  state$pos <- 0L
  state$replicates <- NA_integer_
  state$announced <- FALSE

  labels <- c(base = "Fitting the base model",
              replicate = "Bootstrap replicates",
              dofv = "delta-OFV evaluations")

  close_bar <- function() {
    if (is.null(state$bar)) return(invisible(NULL))
    if (use_cli) {
      cli::cli_progress_done(id = state$bar)
    } else {
      close(state$bar)
    }
    state$bar <- NULL
    state$kind <- NULL
    state$total <- NA_integer_
    state$pos <- 0L
  }

  # Open the bar for `kind` unless the one on screen already has that kind and
  # that total. The total is half the key: a coalesced "started" leaves the
  # replicate bar opened against NA, and reopening on the first counted event is
  # what corrects it. `total` is NA for the base fit - a single fit of unknown
  # length, so cli draws a spinner and the fallback prints one line.
  open_bar <- function(kind, total) {
    total <- if (is.null(total)) NA_integer_ else as.integer(total)
    if (identical(state$kind, kind) && identical(state$total, total)) {
      return(invisible(NULL))
    }
    close_bar()
    state$kind <- kind
    state$total <- total
    state$pos <- 0L
    state$bar <- if (use_cli) {
      cli::cli_progress_bar(labels[[kind]], total = total, .envir = envir)
    } else if (is.na(total) || total <= 0L) {
      # Nothing to count against: the base fit is one fit of unknown length, and
      # a fully resumed run has no replicate left to fit. txtProgressBar() stops
      # on `max <= min`, so say it in a line instead of drawing a bar.
      message(labels[[kind]], " ...")
      NULL
    } else {
      utils::txtProgressBar(min = 0, max = total, style = 3)
    }
    invisible(NULL)
  }

  # Step the bar for `kind`, if that is the one on screen, and never backwards.
  # `completed` is handed out by a `fetch_add` inside the engine's `par_iter`
  # and reported from whichever worker finished the fit, so 5 can reach the slot
  # this session polls before 4 does. Setting it back skews cli's ETA, which
  # reads the rate, and in the text bar redraws a shorter bar.
  advance <- function(kind, completed) {
    if (is.null(state$bar) || !identical(state$kind, kind)) {
      return(invisible(NULL))
    }
    if (length(completed) != 1L || is.na(completed) || completed <= state$pos) {
      return(invisible(NULL))
    }
    state$pos <- completed
    if (use_cli) {
      cli::cli_progress_update(id = state$bar, set = completed, .envir = envir)
    } else {
      utils::setTxtProgressBar(state$bar, completed)
    }
    invisible(NULL)
  }

  # Advance a spinner that has no count to set. cli only animates a bar on an
  # update, and the engine redelivers the latest event on its timer, so the base
  # fit's repeated "started" is what keeps the spinner moving.
  tick_bar <- function() {
    if (is.null(state$bar) || !use_cli) return(invisible(NULL))
    cli::cli_progress_update(id = state$bar, .envir = envir)
    invisible(NULL)
  }

  function(stage, completed, total, base_fit) {
    tryCatch({
      if (identical(stage, "started")) {
        state$replicates <- total
        # `completed` carries the reused count on this one event: replicates a
        # `directory` already holds, which this run will not fit again.
        if (completed > 0L && !state$announced) {
          state$announced <- TRUE
          message("Reusing ", completed, " replicates already on disk")
        }
        if (isTRUE(base_fit)) {
          open_bar("base", NA_integer_)
          tick_bar()
        } else {
          open_bar("replicate", total)
        }
      } else if (identical(stage, "base_done")) {
        open_bar("replicate", state$replicates)
      } else if (identical(stage, "replicate")) {
        # The delta-OFV pass starts only once every replicate is in, so a
        # replicate event arriving after it is a late one under the same
        # out-of-order delivery: step the bar while it is up, never put it back.
        if (!identical(state$kind, "dofv")) open_bar("replicate", total)
        advance("replicate", completed)
      } else if (identical(stage, "dofv")) {
        # Any delta-OFV event retires the replicate bar, not only the one
        # carrying 1: which evaluation reports first is not knowable, so keying
        # the swap on the count would draw the early ones onto the replicate bar
        # and then restart the delta-OFV bar under them.
        open_bar("dofv", total)
        advance("dofv", completed)
      } else if (identical(stage, "finished")) {
        close_bar()
      }
    }, error = function(e) NULL)
    invisible(NULL)
  }
}


# Split an R `sample_size` into the (names, values) pair the Rust glue takes.
# The PsN `1001=>12,1002=>24` string is assembled and parsed on the Rust side,
# so this only classifies the shape.
.ferx_bootstrap_sample_size <- function(sample_size) {
  if (is.null(sample_size)) {
    return(list(keys = character(0), values = numeric(0)))
  }
  if (!is.numeric(sample_size) || length(sample_size) == 0L ||
      any(is.na(sample_size))) {
    stop("`sample_size` must be a number, or a named numeric vector of ",
         "per-stratum counts, or NULL.")
  }
  # A stratum count of 0 is accepted all the way down to `SampleSize::Total(0)`,
  # where it becomes replicates drawing no subjects at all - a fit error far
  # from its cause.
  if (any(sample_size != trunc(sample_size)) || any(sample_size < 1)) {
    stop("`sample_size` counts must be whole numbers >= 1.")
  }
  nms <- names(sample_size)
  if (is.null(nms)) {
    if (length(sample_size) != 1L) {
      stop("`sample_size`: an unnamed vector must be a single number; name the ",
           "elements (c(`1001` = 12, `1002` = 24)) for per-stratum counts.")
    }
    return(list(keys = character(0), values = as.numeric(sample_size)))
  }
  if (any(!nzchar(nms))) {
    stop("`sample_size`: either every element is named with its stratum, or none is.")
  }
  list(keys = as.character(nms), values = as.numeric(unname(sample_size)))
}

# `breaks = "FD"` and the other named rules derive a bin width from the spread
# of the data: with a single value the IQR is 0, the width comes out NA, and
# hist() aborts on `if (h > 0)`. One bin is the only honest picture of one
# estimate anyway.
.ferx_hist_breaks <- function(vals, breaks) {
  if (length(vals) < 2L) 1L else breaks
}

# raw_results.csv is untyped text: the engine writes the four status flags as
# 0/1 and leaves `error` blank for a replicate that succeeded, so read.csv
# hands back integer flags - and, when no replicate errored, an all-blank
# `error` column type-converted to logical NA. The in-memory path builds the
# same frame with logical flags and a character `error`; make the two agree.
.ferx_raw_csv_types <- function(raw) {
  if (is.null(raw)) return(raw)
  for (nm in c("minimization_successful", "estimate_near_boundary",
               "covariance_step_successful", "covariance_step_warnings")) {
    if (nm %in% names(raw)) raw[[nm]] <- as.logical(raw[[nm]])
  }
  if ("error" %in% names(raw)) raw$error <- as.character(raw$error)
  raw
}

# The engine writes "" for a replicate that did not error, matching
# raw_results.csv; NA_character_ is what an R user expects to test for.
.ferx_blank_error_to_na <- function(raw) {
  if (is.null(raw) || !"error" %in% names(raw)) return(raw)
  raw$error[!nzchar(raw$error)] <- NA_character_
  raw
}

.ferx_check_count <- function(value, name, min) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value != trunc(value) || value < min) {
    stop("`", name, "` must be a single whole number >= ", min, ".")
  }
  invisible(TRUE)
}

.ferx_diag_value <- function(diagnostics, statistic) {
  if (is.null(diagnostics)) return(NA_real_)
  i <- match(statistic, as.character(diagnostics$statistic))
  if (is.na(i)) NA_real_ else diagnostics$value[i]
}
