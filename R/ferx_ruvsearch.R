#' Residual-error model search
#'
#' Runs a residual-error search - Pharmpy's \code{ruvsearch} - over a model's
#' \code{[error_model]} block. Each iteration adds one residual-error feature
#' to the model the last iteration kept, fits every candidate, and accepts the
#' largest improvement the likelihood-ratio test calls significant at
#' \code{p_value}. Up to \code{max_iter} features are stacked this way.
#'
#' There is no search space to state: the candidates are the four families
#' Pharmpy tests, and \code{skip} is what narrows them.
#'   \describe{
#'     \item{\code{IIV_on_RUV}}{A log-normal per-subject scale on the residual
#'       SD (\code{ETA_RUV}). Needs an estimation method with eta-epsilon
#'       interaction (\code{focei}, \code{imp}, \code{impmap}, \code{saem}); on
#'       a \code{foce} model it is not tested, and the run says so in
#'       \code{$notes}.}
#'     \item{\code{power}}{The proportional loading raised to an estimated
#'       exponent, \code{power(sigma, P)}.}
#'     \item{\code{combined}}{An additive component beside the proportional
#'       one.}
#'     \item{\code{time_varying}}{Every sigma scaled by a theta for records
#'       below a time-after-dose cutoff - \code{groups - 1} candidates, at the
#'       \code{i / groups} quantiles of time after dose.}
#'   }
#' \code{power} and \code{combined} are tested together or not at all, as in
#' Pharmpy, and accepting one of them retires both.
#'
#' The search always starts from a plain proportional error model: when the
#' input is not one, its proportional counterpart is fitted first and becomes
#' the parent of iteration one (\code{base_model_id} says which). The final
#' comparison then checks the accepted stack against the input as well, so a
#' search can return the model it started from.
#'
#' A candidate has to pass the strictness gate before it can be accepted, and
#' the step table always carries the termination status and the gate verdict
#' beside the p-value - a candidate that stalled at its initial estimates
#' carries a dOFV that says nothing about the error model it was testing.
#'
#' @section Two entry forms:
#' Pass either \code{config} - a \code{.ferxsearch} file, which is the
#' reproducible artifact - or the inline arguments. The inline form is rendered
#' into the same configuration and handed to the engine's own loader, so the
#' two cannot disagree.
#'
#' The run arguments (\code{threads}, \code{retries}, \code{resume},
#' \code{directory}, \code{progress}) say how to run a search rather than what
#' to search, and are available to both forms.
#'
#' A \code{[space]} section, or a \code{[rank]} asking for a BIC or a
#' \code{cutoff}, is refused by name when the file is read: ruvsearch selects
#' on the likelihood-ratio test at \code{p_value}, not on a ranking criterion.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the search runs on, so like \code{model} it
#'   cannot be given beside \code{config} - the file's own \code{data} key
#'   says which dataset that file searches.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param groups Time-after-dose bins the time-varying candidates cut at:
#'   \code{groups - 1} candidates, at the \code{i / groups} quantiles.
#'   \code{NULL} keeps the engine default (4).
#' @param p_value The likelihood-ratio level a candidate must reach.
#'   \code{NULL} keeps the engine default (0.001).
#' @param skip Families never tested, any of \code{"IIV_on_RUV"},
#'   \code{"power"}, \code{"combined"}, \code{"time_varying"}. \code{NULL}
#'   tests them all.
#' @param max_iter Iterations, 1 to 3 as in Pharmpy. \code{NULL} keeps the
#'   engine default (3).
#' @param cwres_prescreen Screen the candidates on the parent's CWRES first and
#'   refit only the winner - Pharmpy's own path, and much cheaper than fitting
#'   every candidate to the data. \code{NULL} keeps the engine default
#'   (\code{FALSE}: every candidate is fitted).
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the per-iteration journals, \code{steps.csv},
#'   \code{models/<id>.ferx} and \code{final.ferx} are written. \code{NULL}
#'   keeps the run in memory, which also means it cannot be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}.
#' @param progress Print the engine's iteration progress to the console.
#'
#' @return An object of class \code{ferx_ruvsearch}:
#'   \describe{
#'     \item{steps}{The step table, one row per fitted model, in the engine's
#'       own column order: \code{iteration}, \code{candidate},
#'       \code{feature}, \code{screened}, \code{parent_ofv}, \code{ofv},
#'       \code{dofv}, \code{df}, \code{p_value}, \code{alpha},
#'       \code{significant}, \code{cwres_dofv}, \code{selected},
#'       \code{converged}, \code{passed}, \code{failures}, \code{seconds} -
#'       followed by \code{family} (the feature's family, which is what
#'       \code{skip} names and what a selection retires) and \code{note} (why
#'       a row could not be compared, on its own; \code{failures} carries it
#'       too when the gate had nothing to say). A row with \code{screened =
#'       TRUE} is a CWRES pre-screen fit: its \code{ofv} is on the CWRES scale
#'       rather than the data's, and \code{cwres_dofv} is its improvement over
#'       the CWRES base - only the feature the screen picked was refitted to
#'       the data.}
#'     \item{fit}{The final model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{model_text}{Every fitted model's text, named by candidate id, so a
#'       form the search rejected can be read or refitted without re-running
#'       it. \code{final_model_path} is the \code{final.ferx} written beside
#'       the run.}
#'     \item{final_features, final_families}{The features the final model
#'       carries beyond the base, in the order they were accepted. Empty when
#'       the search returned the input or the proportional base.}
#'     \item{input_model, input_ofv}{The model file's own error model, fitted
#'       as given.}
#'     \item{base_model_id, base_ofv}{The parent of iteration one:
#'       \code{"input"}, or \code{"base"} when a proportional base had to be
#'       derived.}
#'     \item{final_model, final_model_id, final_ofv}{The winning model's text,
#'       its id, and its objective function value.}
#'     \item{options}{The search as the engine read it: \code{groups},
#'       \code{p_value}, \code{skip}, \code{max_iter}, \code{cwres_prescreen}
#'       and \code{cutoff} (the \code{df = 1} chi-square cutoff at
#'       \code{p_value}).}
#'     \item{n_iterations}{Iterations actually run.}
#'     \item{candidates}{The runner's candidate table when the run wrote one -
#'       see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once (a family not
#'       tested and why, a reversion), and whether it was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("one_cpt_iv")
#'
#' # Inline: no file to author for a one-off
#' res <- ferx_ruvsearch(
#'   model     = ex$model,
#'   data      = ex$data,
#'   p_value   = 0.05,
#'   max_iter  = 2,
#'   directory = "ruvsearch-run-1"
#' )
#' res
#' res$steps[, c("iteration", "feature", "family", "dofv", "p_value", "passed")]
#' summary(res)
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res2 <- ferx_ruvsearch(config = ex$search, directory = "ruvsearch-run-2")
#' }
#'
#' @seealso \code{\link{ferx_search_config}}, \code{\link{ferx_search_results}},
#'   \code{\link{ferx_covsearch}}, \code{\link{ferx_modelsearch}}
#' @family search
#' @export
ferx_ruvsearch <- function(model = NULL,
                           data = NULL,
                           config = NULL,
                           groups = NULL,
                           p_value = NULL,
                           skip = NULL,
                           max_iter = NULL,
                           cwres_prescreen = NULL,
                           threads = NULL,
                           retries = NULL,
                           directory = NULL,
                           resume = FALSE,
                           progress = interactive()) {
  what <- "ferx_ruvsearch"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      groups = groups,
      p_value = p_value,
      skip = skip,
      max_iter = max_iter,
      cwres_prescreen = cwres_prescreen
    ),
    what
  )

  paths <- if (nzchar(config_path)) {
    list(model = "", data = "")
  } else {
    .ferx_search_model_data(model, data, what)
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_ruvsearch(
    config_path     = config_path,
    model_path      = paths$model,
    data_path       = paths$data,
    groups          = .ferx_search_count(groups, "groups", what, below = 0L),
    p_value         = .ferx_search_scalar(p_value, "p_value", what, positive = TRUE),
    skip            = .ferx_ruv_skip(skip, what),
    max_iter        = .ferx_search_count(max_iter, "max_iter", what, below = 0L),
    cwres_prescreen = .ferx_search_flag3(cwres_prescreen, "cwres_prescreen", what),
    threads         = .ferx_search_count(threads, "threads", what, below = 0L),
    retries         = .ferx_search_count(retries, "retries", what, below = -1L, min = 0L),
    resume          = .ferx_search_bool(resume, "resume", what),
    directory       = dir_arg,
    progress        = .ferx_search_bool(progress, "progress", what)
  )

  # The engine's 17 columns in the engine's order, then `family` and `note` -
  # the feature's family (what `skip` names and what a selection retires) and
  # the reason a row could not be compared, on its own.
  steps <- data.frame(
    iteration   = as.integer(raw$iteration),
    candidate   = as.character(raw$candidate),
    feature     = .ferx_search_chr(raw$feature),
    screened    = as.logical(raw$screened),
    parent_ofv  = .ferx_search_num(raw$parent_ofv),
    ofv         = .ferx_search_num(raw$ofv),
    dofv        = .ferx_search_num(raw$dofv),
    df          = as.integer(.ferx_search_num(raw$df)),
    p_value     = .ferx_search_num(raw$p_value),
    alpha       = .ferx_search_num(raw$alpha),
    significant = .ferx_search_lgl(raw$significant),
    cwres_dofv  = .ferx_search_num(raw$cwres_dofv),
    selected    = as.logical(raw$selected),
    converged   = .ferx_search_lgl(raw$converged),
    passed      = as.logical(raw$passed),
    failures    = .ferx_search_chr(raw$failures),
    seconds     = .ferx_search_num(raw$seconds),
    family      = .ferx_search_chr(raw$family),
    note        = .ferx_search_chr(raw$note),
    stringsAsFactors = FALSE
  )

  model_text <- stats::setNames(as.character(raw$model_text),
                                as.character(raw$model_id))

  # The winning model as a file: the run's own `final.ferx` when it wrote one
  # (the engine seeds the final estimates into it), otherwise a temporary copy
  # so the fit below has a model file to name and read.
  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-ruvsearch-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path, data = as.character(raw$data))
  }

  result <- list(
    steps            = steps,
    fit              = fit,
    model_text       = model_text,
    options          = list(
      groups          = as.integer(raw$opt_groups),
      p_value         = as.numeric(raw$opt_p_value),
      skip            = as.character(raw$opt_skip),
      max_iter        = as.integer(raw$opt_max_iter),
      cwres_prescreen = isTRUE(raw$opt_cwres_prescreen),
      cutoff          = as.numeric(raw$opt_cutoff)
    ),
    input_model      = as.character(raw$input_model),
    input_ofv        = as.numeric(raw$input_ofv),
    base_model_id    = as.character(raw$base_id),
    base_ofv         = as.numeric(raw$base_ofv),
    final_model      = as.character(raw$final_model),
    final_model_path = final_path,
    final_model_id   = as.character(raw$final_id),
    final_ofv        = as.numeric(raw$final_ofv),
    final_features   = as.character(raw$final_features),
    final_families   = as.character(raw$final_families),
    n_iterations     = as.integer(raw$n_iterations),
    summary_text     = as.character(raw$summary),
    candidates       = .ferx_search_candidates(
      dir_arg,
      .ferx_ruv_step_dirs(steps, as.character(raw$base_id))
    ),
    model            = as.character(raw$model),
    data             = as.character(raw$data),
    directory        = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config           = if (nzchar(config_path)) config_path else NA_character_,
    notes            = as.character(raw$notes),
    cancelled        = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_ruvsearch", "ferx_search_result")
  result
}

# The runner directories this run wrote, in the order the engine wrote them:
# the input, the proportional base when one had to be derived, and per
# iteration the CWRES pre-screen (`screen-N`) and the fits to the data
# (`iteration-N`). Taken from the rows the run itself produced, so a shorter
# run in a directory an earlier one used cannot inherit its candidate tables.
.ferx_ruv_step_dirs <- function(steps, base_model_id) {
  dirs <- "input"
  if (!identical(base_model_id, "input")) dirs <- c(dirs, "base")
  for (i in sort(unique(steps$iteration[steps$iteration > 0L]))) {
    rows <- steps$screened[steps$iteration == i]
    if (any(rows, na.rm = TRUE)) dirs <- c(dirs, sprintf("screen-%d", i))
    if (any(!rows, na.rm = TRUE)) dirs <- c(dirs, sprintf("iteration-%d", i))
  }
  dirs
}

# The `skip` families, as the engine spells them. A family is a closed set, so
# a misspelling is caught here rather than becoming a TOML key the loader
# refuses several layers down.
.ferx_ruv_skip <- function(skip, what) {
  if (is.null(skip)) return(character(0))
  if (!is.character(skip) || anyNA(skip)) {
    stop(sprintf("%s: `skip` must be a character vector of family names or NULL", what))
  }
  known <- c("IIV_on_RUV", "power", "combined", "time_varying")
  hit <- match(tolower(skip), tolower(known))
  if (anyNA(hit)) {
    stop(sprintf("%s: `skip` does not name a residual-error family: %s. Use any of %s",
                 what,
                 paste(skip[is.na(hit)], collapse = ", "),
                 paste(known, collapse = ", ")))
  }
  unique(known[hit])
}

#' @param x A \code{ferx_ruvsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_ruvsearch
#' @export
print.ferx_ruvsearch <- function(x, digits = 4, ...) {
  cat("ferx residual-error search (likelihood-ratio test at p = ",
      format(x$options$p_value), ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  # Fixed decimals rather than significant digits: a residual-error search is
  # decided on differences of a few OFV units, which four significant digits
  # would round away on a four-figure objective function value.
  ofv3 <- function(v) formatC(v, format = "f", digits = 3)
  cat(sprintf("  OFV:   %s (input) -> %s (%s, dOFV %s)\n",
              ofv3(x$input_ofv), ofv3(x$final_ofv), x$final_model_id,
              ofv3(x$input_ofv - x$final_ofv)))
  if (!identical(x$base_model_id, "input")) {
    cat(sprintf("  Base:  %s - the proportional model the search started from (OFV %s)\n",
                x$base_model_id, ofv3(x$base_ofv)))
  }
  cat(sprintf("  %d iteration%s, %d model%s fitted\n",
              x$n_iterations, if (x$n_iterations == 1L) "" else "s",
              nrow(x$steps), if (nrow(x$steps) == 1L) "" else "s"))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the rows below are what it reached.\n")
  }

  cat("\nCandidates by iteration:\n")
  tested <- x$steps[x$steps$iteration > 0L, , drop = FALSE]
  # A pre-screen row's `ofv` is on the CWRES scale, not the data's, so the
  # column that says which rows those are has to be on the table whenever the
  # run has any - two objective functions in one column, unlabelled, would read
  # as one.
  cols <- c("iteration", "feature", "family", "ofv", "dofv", "p_value",
            "significant", "selected", "converged", "passed")
  if (any(tested$screened, na.rm = TRUE)) {
    cols <- append(cols, c("screened", "cwres_dofv"), after = 3L)
  }
  .ferx_search_print_table(tested[, cols, drop = FALSE], digits)

  screened <- sum(x$steps$screened, na.rm = TRUE)
  if (screened > 0L) {
    cat(sprintf(paste0("\n%d row%s a CWRES pre-screen fit: its `ofv` is on the CWRES scale, ",
                       "not the data's,\nand `cwres_dofv` is its improvement over the CWRES ",
                       "base. Only the winner was refitted.\n"),
                screened, if (screened == 1L) " is" else "s are"))
  }
  # Screened rows are excluded from the count: a fit to the parent's CWRES is
  # not a candidate the gate rejected, and folding the two together would
  # report a search as riddled with failures whenever the pre-screen ran.
  gated <- sum(!x$steps$passed & !x$steps$screened, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d candidate%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }

  cat("\nSelected residual-error model:\n")
  if (length(x$final_features) == 0L) {
    cat("  (none - the ", x$final_model_id,
        " model's error model was not improved on)\n", sep = "")
  } else {
    for (i in seq_along(x$final_features)) {
      cat(sprintf("  %s  (%s)\n", x$final_features[i], x$final_families[i]))
    }
  }
  block <- tryCatch(.ferx_read_section(x$final_model_path, "error_model", strip = TRUE),
                    error = function(e) NULL)
  if (length(block)) {
    block <- block[nzchar(trimws(block))]
    for (line in block) cat("    ", line, "\n", sep = "")
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

#' @param object A \code{ferx_ruvsearch} object.
#' @rdname ferx_ruvsearch
#' @export
summary.ferx_ruvsearch <- function(object, digits = 4, ...) {
  print(object, digits = digits)

  cat("\nForms not selected, and why:\n")
  rejected <- object$steps[object$steps$iteration > 0L &
                             !(!is.na(object$steps$selected) & object$steps$selected), ,
                           drop = FALSE]
  if (nrow(rejected) == 0L) {
    cat("  (none - every candidate tested was accepted)\n")
  } else {
    rejected$reason <- ifelse(
      rejected$screened,
      # A pre-screen row was never a candidate to accept: it says which
      # feature the screen would refit, on the CWRES scale.
      sprintf("CWRES pre-screen only (CWRES dOFV %s)",
              format(signif(rejected$cwres_dofv, digits))),
      ifelse(!is.na(rejected$failures), rejected$failures,
             ifelse(!is.na(rejected$significant) & !rejected$significant,
                    sprintf("p = %s, above the %s cutoff",
                            format(signif(rejected$p_value, digits)),
                            format(object$options$p_value)),
                    "not the largest significant improvement of its iteration"))
    )
    .ferx_search_print_table(
      rejected[, c("iteration", "feature", "family", "dofv", "p_value",
                   "converged", "passed", "reason"), drop = FALSE],
      digits
    )
  }

  cat("\nEvery model fitted:\n")
  .ferx_search_print_table(
    object$steps[, c("iteration", "candidate", "feature", "screened", "ofv",
                     "dofv", "df", "p_value", "significant", "cwres_dofv",
                     "selected", "converged", "passed", "failures"),
                 drop = FALSE],
    digits
  )
  invisible(object)
}
