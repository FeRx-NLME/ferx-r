#' Stepwise covariate modelling (SCM)
#'
#' Runs a stepwise covariate search - PsN's \code{scm}, Pharmpy's
#' \code{covsearch} - over a model and a candidate set of covariate effects.
#' Each forward step fits every remaining effect on its own, keeps the largest
#' OFV drop that is significant by the likelihood-ratio test, and repeats; the
#' backward phase then removes the cheapest effect whose removal is not
#' significant.
#'
#' Every candidate is fitted with perturbed restarts (\code{retries}) and has to
#' pass the strictness gate before it can win, because a candidate that stalled
#' at its initial estimates carries an OFV that says nothing about the model.
#' The step table therefore always shows the termination status and the gate
#' verdict beside the OFV, and a candidate the gate excluded is a row with its
#' reason rather than an absence.
#'
#' @section Two entry forms:
#' Pass either \code{config} - a \code{.ferxsearch} file, which is the
#' reproducible artifact - or the inline arguments. The inline form is rendered
#' into the same configuration and handed to the engine's own loader, so the two
#' cannot disagree; \code{search_space} is MFL text, quoted verbatim, which is
#' also what makes a space portable to and from Pharmpy.
#'
#' The run arguments (\code{threads}, \code{retries}, \code{resume},
#' \code{directory}, \code{progress}) say how to run a search rather than what
#' to search, and are available to both forms.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the search runs on, so like \code{model} it
#'   cannot be given beside \code{config} - the file's own \code{data} key
#'   says which dataset that file searches.
#' @param search_space MFL text naming the candidate effects, e.g.
#'   \code{"COVARIATE?(@IIV, @CONTINUOUS, [pow,lin])"}. A \code{COVARIATE(...)}
#'   statement without \code{?} forces that effect into the base model before
#'   the search starts.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param algorithm \code{"scm-forward-then-backward"} or \code{"scm-forward"}.
#'   \code{NULL} keeps the engine default.
#' @param p_forward,p_backward Significance levels of the forward and backward
#'   tests (engine defaults 0.01 and 0.001).
#' @param max_steps Cap on the number of steps; \code{NULL} for unlimited.
#' @param adaptive_scope_reduction Stash effects that fail badly and retest them
#'   once at the end (SCM+). \code{NULL} keeps the engine default.
#' @param rank \code{[rank] type} for the strictness / ranking criterion, e.g.
#'   \code{"bic"}. \code{NULL} keeps the tool default (covsearch selects on the
#'   likelihood-ratio test).
#' @param cutoff \code{[rank] cutoff}. \code{NULL} keeps the default.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the per-step journals, \code{steps.csv} and
#'   \code{final.ferx} are written. \code{NULL} keeps the run in memory, which
#'   also means it cannot be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}.
#' @param progress Print the engine's step progress to the console.
#'
#' @return An object of class \code{ferx_covsearch}:
#'   \describe{
#'     \item{steps}{The step table, one row per candidate of every step, in the
#'       engine's own column order: \code{step}, \code{phase},
#'       \code{candidate}, \code{parameter}, \code{covariate}, \code{form},
#'       \code{parent_ofv}, \code{ofv}, \code{dofv}, \code{df},
#'       \code{p_value}, \code{alpha}, \code{significant}, \code{selected},
#'       \code{converged}, \code{passed}, \code{failures}.}
#'     \item{included}{The final relation set: \code{parameter},
#'       \code{covariate}, \code{form} and \code{origin} (the base model, a
#'       forced statement, or the forward step that added it).}
#'     \item{fit}{The final model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{base_model, final_model}{The model text the search started and
#'       ended on; \code{final_model_path} is the \code{final.ferx} written
#'       beside the run.}
#'     \item{base_ofv, final_ofv, final_step}{The two objective function values
#'       and the step that produced the winner (\code{0} is the base model).}
#'     \item{candidates}{The runner's candidate table when the run wrote one -
#'       see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once, and whether it
#'       was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("two_cpt_oral_cov")
#'
#' # Inline: no file to author for a one-off
#' res <- ferx_covsearch(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "COVARIATE?(@IIV, @CONTINUOUS, [pow,lin])",
#'   p_forward    = 0.01,
#'   p_backward   = 0.001,
#'   directory    = "covsearch-run-1"
#' )
#' res
#' res$steps[res$steps$selected, ]
#' summary(res)
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res2 <- ferx_covsearch(config = ex$search, directory = "covsearch-run-2")
#' }
#'
#' @seealso \code{\link{ferx_search_config}}, \code{\link{ferx_search_space}},
#'   \code{\link{ferx_search_results}}, \code{\link{ferx_allometry}}
#' @family search
#' @export
ferx_covsearch <- function(model = NULL,
                           data = NULL,
                           search_space = NULL,
                           config = NULL,
                           algorithm = NULL,
                           p_forward = NULL,
                           p_backward = NULL,
                           max_steps = NULL,
                           adaptive_scope_reduction = NULL,
                           rank = NULL,
                           cutoff = NULL,
                           threads = NULL,
                           retries = NULL,
                           directory = NULL,
                           resume = FALSE,
                           progress = interactive()) {
  what <- "ferx_covsearch"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      search_space = search_space,
      algorithm = algorithm,
      p_forward = p_forward,
      p_backward = p_backward,
      max_steps = max_steps,
      adaptive_scope_reduction = adaptive_scope_reduction,
      rank = rank,
      cutoff = cutoff
    ),
    what
  )

  if (nzchar(config_path)) {
    paths <- list(model = "", data = "")
    mfl <- ""
  } else {
    paths <- .ferx_search_model_data(model, data, what)
    mfl <- .ferx_search_space_text(search_space, what)
  }

  if (!is.null(algorithm)) {
    algorithm <- match.arg(algorithm, c("scm-forward-then-backward", "scm-forward"))
  }
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_covsearch(
    config_path = config_path,
    model_path  = paths$model,
    data_path   = paths$data,
    mfl         = mfl,
    algorithm   = algorithm %||% "",
    p_forward   = .ferx_search_scalar(p_forward, "p_forward", what, positive = TRUE),
    p_backward  = .ferx_search_scalar(p_backward, "p_backward", what, positive = TRUE),
    max_steps   = .ferx_search_count(max_steps, "max_steps", what, below = 0L),
    adaptive    = .ferx_search_flag3(adaptive_scope_reduction,
                                     "adaptive_scope_reduction", what),
    rank        = rank %||% "",
    rank_cutoff = .ferx_search_scalar(cutoff, "cutoff", what),
    threads     = .ferx_search_count(threads, "threads", what, below = 0L),
    retries     = .ferx_search_count(retries, "retries", what, below = -1L, min = 0L),
    resume      = .ferx_search_bool(resume, "resume", what),
    directory   = dir_arg,
    progress    = .ferx_search_bool(progress, "progress", what)
  )

  steps <- data.frame(
    step        = as.integer(raw$step),
    phase       = as.character(raw$phase),
    candidate   = as.character(raw$candidate),
    parameter   = as.character(raw$parameter),
    covariate   = as.character(raw$covariate),
    form        = as.character(raw$form),
    parent_ofv  = .ferx_search_num(raw$parent_ofv),
    ofv         = .ferx_search_num(raw$ofv),
    dofv        = .ferx_search_num(raw$dofv),
    df          = as.integer(.ferx_search_num(raw$df)),
    p_value     = .ferx_search_num(raw$p_value),
    alpha       = .ferx_search_num(raw$alpha),
    significant = .ferx_search_lgl(raw$significant),
    selected    = as.logical(raw$selected),
    converged   = .ferx_search_lgl(raw$converged),
    passed      = as.logical(raw$passed),
    failures    = .ferx_search_chr(raw$failures),
    stringsAsFactors = FALSE
  )

  included <- data.frame(
    parameter = as.character(raw$included_parameter),
    covariate = as.character(raw$included_covariate),
    form      = as.character(raw$included_form),
    origin    = as.character(raw$included_origin),
    stringsAsFactors = FALSE
  )

  # The winning model as a file: the run's own `final.ferx` when it wrote one
  # (the engine seeds the final estimates into it), otherwise a temporary copy
  # so the fit below has a model file to name and read.
  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-covsearch-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path, data = as.character(raw$data))
  }

  result <- list(
    steps            = steps,
    included         = included,
    fit              = fit,
    base_model       = as.character(raw$base_model),
    final_model      = as.character(raw$final_model),
    final_model_path = final_path,
    base_ofv         = as.numeric(raw$base_ofv),
    final_ofv        = as.numeric(raw$final_ofv),
    final_step       = as.integer(raw$final_step),
    algorithm        = as.character(raw$algorithm),
    candidates       = .ferx_search_candidates(dir_arg),
    model            = as.character(raw$model),
    data             = as.character(raw$data),
    directory        = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config           = if (nzchar(config_path)) config_path else NA_character_,
    notes            = as.character(raw$notes),
    cancelled        = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_covsearch", "ferx_search_result")
  result
}

#' @param x A \code{ferx_covsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_covsearch
#' @export
print.ferx_covsearch <- function(x, digits = 4, ...) {
  cat("ferx covariate search (", x$algorithm, ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  cat(sprintf("  OFV:   %s (base) -> %s (step %d)\n",
              format(signif(x$base_ofv, digits)),
              format(signif(x$final_ofv, digits)),
              x$final_step))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the steps below are what it reached.\n")
  }

  cat("\nSteps (the candidate selected at each step):\n")
  sel <- x$steps[!is.na(x$steps$selected) & x$steps$selected,
                 c("step", "phase", "parameter", "covariate", "form",
                   "dofv", "p_value", "converged", "passed"), drop = FALSE]
  .ferx_search_print_table(sel, digits)

  gated <- sum(!x$steps$passed, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d candidate%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }

  cat("\nFinal covariate relations:\n")
  if (nrow(x$included) == 0L) {
    cat("  (none - the base model was not improved on)\n")
  } else {
    for (i in seq_len(nrow(x$included))) {
      cat(sprintf("  %s ~ %s  %s  (%s)\n",
                  x$included$parameter[i], x$included$covariate[i],
                  x$included$form[i], x$included$origin[i]))
    }
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

#' @param object A \code{ferx_covsearch} object.
#' @rdname ferx_covsearch
#' @export
summary.ferx_covsearch <- function(object, digits = 4, ...) {
  print(object, digits = digits)
  cat("\nEvery candidate:\n")
  .ferx_search_print_table(
    object$steps[, c("step", "phase", "parameter", "covariate", "form",
                     "ofv", "dofv", "p_value", "significant", "selected",
                     "converged", "passed", "failures"), drop = FALSE],
    digits
  )
  invisible(object)
}
