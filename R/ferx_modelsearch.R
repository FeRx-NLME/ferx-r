#' Structural PK model search
#'
#' Runs a structural search - Pharmpy's \code{modelsearch} - over a base model
#' and a space of structural features: absorption, elimination, peripheral
#' compartments, transit compartments and lag time. Every candidate is built by
#' the engine's own model editor, fitted, gated, and ranked on the criterion
#' \code{rank} names (the mixed BIC unless the file says otherwise).
#'
#' The algorithm decides which candidates exist. \code{"exhaustive"} fits every
#' combination of the space in one layer; \code{"exhaustive_stepwise"} adds one
#' feature per layer to every model of the layer before it;
#' \code{"reduced_stepwise"} - the default - extends only the best model of
#' each feature-set group, which is what keeps a four-feature space to tens of
#' fits rather than hundreds.
#'
#' As with \code{\link{ferx_covsearch}}, a candidate has to pass the strictness
#' gate before it can be ranked, and the model table always carries the
#' termination status and the gate verdict beside the criterion - a candidate
#' that stalled at its initial estimates carries a criterion that says nothing
#' about the structure it was testing.
#'
#' @section Two entry forms:
#' Pass either \code{config} - a \code{.ferxsearch} file, which is the
#' reproducible artifact - or the inline arguments. The inline form is rendered
#' into the same configuration and handed to the engine's own loader, so the
#' two cannot disagree; \code{search_space} is MFL text, quoted verbatim, which
#' is also what makes a space portable to and from Pharmpy.
#'
#' The run arguments (\code{threads}, \code{retries}, \code{resume},
#' \code{directory}, \code{progress}) say how to run a search rather than what
#' to search, and are available to both forms.
#'
#' @section What the space may say:
#' \code{ABSORPTION}, \code{ELIMINATION}, \code{PERIPHERALS}, \code{TRANSITS}
#' and \code{LAGTIME}. A covariate or variability statement is another tool's
#' space and is refused by name before any fit starts, as is a feature the
#' engine cannot build (\code{ABSORPTION(SEQ-ZO-FO)}). Check a space with
#' \code{\link{ferx_search_space}} / \code{\link{ferx_search_coverage}} first:
#' both answer without fitting anything.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]} block.
#' @param search_space MFL text naming the structural space, e.g.
#'   \code{"ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])"}.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param algorithm \code{"reduced_stepwise"} (the default),
#'   \code{"exhaustive_stepwise"} or \code{"exhaustive"}. \code{NULL} keeps the
#'   engine default.
#' @param iiv_strategy How a candidate's new parameters are given a random
#'   effect: \code{"absorption_delay"} (the default - an eta on a new lag time
#'   or mean transit time only), \code{"add_diagonal"} (an eta on every new PK
#'   parameter) or \code{"no_add"}. \code{NULL} keeps the engine default.
#' @param rank \code{[rank] type} - the criterion candidates are ranked and
#'   gated on, e.g. \code{"bic"}, \code{"aic"}, \code{"ofv"}. \code{NULL} keeps
#'   the tool default (the mixed BIC).
#' @param cutoff \code{[rank] cutoff}: the improvement on the criterion a
#'   candidate must show over the base model to be selected. \code{NULL}
#'   selects the best model on the criterion alone.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the per-layer journals, \code{models.csv},
#'   \code{models/<id>.ferx} and \code{final.ferx} are written. \code{NULL}
#'   keeps the run in memory, which also means it cannot be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}.
#' @param progress Print the engine's layer progress to the console.
#'
#' @return An object of class \code{ferx_modelsearch}:
#'   \describe{
#'     \item{models}{The model table, one row per fitted model, in the engine's
#'       own column order: \code{id}, \code{parent}, \code{layer},
#'       \code{path}, \code{absorption}, \code{peripherals}, \code{transits},
#'       \code{lagtime}, \code{n_parameters}, \code{ofv}, \code{criterion},
#'       \code{d_criterion}, \code{rank}, \code{converged}, \code{passed},
#'       \code{failures}, \code{error}, \code{seconds}, \code{selected},
#'       \code{continued}, \code{reused} - followed by \code{structure}, the
#'       engine's one-line rendering of the same four structural columns.}
#'     \item{fit}{The winning model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{model_text}{Every candidate's model text, named by model id, so
#'       the model the table ranked second can be read or refitted without
#'       re-running the search. \code{final_model_path} is the
#'       \code{final.ferx} written beside the run.}
#'     \item{input_model, base_model_id, base_structure}{The model the search
#'       started from, and the root of the space it searched (they differ when
#'       the input lies off the space and had to be moved onto it).}
#'     \item{base_ofv, base_criterion, final_ofv, final_criterion,
#'       final_model_id}{The two objective function values, the two criterion
#'       values, and the id of the model that won.}
#'     \item{criterion, algorithm, iiv_strategy, n_layers}{What the search
#'       ranked on and how it searched.}
#'     \item{candidates}{The runner's candidate table when the run wrote one -
#'       see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once, and whether it
#'       was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Inline: no file to author for a one-off
#' res <- ferx_modelsearch(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
#'   directory    = "modelsearch-run-1"
#' )
#' res
#' res$models[res$models$passed, c("id", "structure", "criterion", "rank")]
#' summary(res)
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res2 <- ferx_modelsearch(config = ex$search, directory = "modelsearch-run-2")
#' }
#'
#' @seealso \code{\link{ferx_search_space}}, \code{\link{ferx_search_coverage}},
#'   \code{\link{ferx_covsearch}}, \code{\link{ferx_search_results}}
#' @family search
#' @export
ferx_modelsearch <- function(model = NULL,
                             data = NULL,
                             search_space = NULL,
                             config = NULL,
                             algorithm = NULL,
                             iiv_strategy = NULL,
                             rank = NULL,
                             cutoff = NULL,
                             threads = NULL,
                             retries = NULL,
                             directory = NULL,
                             resume = FALSE,
                             progress = interactive()) {
  what <- "ferx_modelsearch"
  config_path <- .ferx_search_entry_form(
    config, model,
    list(
      search_space = search_space,
      algorithm = algorithm,
      iiv_strategy = iiv_strategy,
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
    algorithm <- match.arg(
      algorithm,
      c("reduced_stepwise", "exhaustive_stepwise", "exhaustive")
    )
  }
  if (!is.null(iiv_strategy)) {
    # `fullblock` is Pharmpy's fourth strategy and the engine refuses it by
    # name (a block over the new and existing eta is a variability search's
    # move), so it stays in the choices rather than being silently unavailable.
    iiv_strategy <- match.arg(
      iiv_strategy,
      c("absorption_delay", "add_diagonal", "no_add", "fullblock")
    )
  }
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_modelsearch(
    config_path  = config_path,
    model_path   = paths$model,
    data_path    = paths$data,
    mfl          = mfl,
    algorithm    = algorithm %||% "",
    iiv_strategy = iiv_strategy %||% "",
    rank         = rank %||% "",
    rank_cutoff  = .ferx_search_scalar(cutoff, "cutoff", what),
    threads      = .ferx_search_count(threads, "threads", what, below = 0L),
    retries      = .ferx_search_count(retries, "retries", what, below = -1L, min = 0L),
    resume       = .ferx_search_bool(resume, "resume", what),
    directory    = dir_arg,
    progress     = .ferx_search_bool(progress, "progress", what)
  )

  # The engine's 21 columns in the engine's order, then `structure` - the
  # engine's own one-line rendering of the four structural columns, which the
  # table needs to be readable and R must not spell for itself.
  models <- data.frame(
    id           = as.character(raw$id),
    parent       = .ferx_search_chr(raw$parent),
    layer        = as.integer(raw$layer),
    path         = .ferx_search_chr(raw$path),
    absorption   = as.character(raw$absorption),
    peripherals  = as.integer(raw$peripherals),
    transits     = as.character(raw$transits),
    lagtime      = as.character(raw$lagtime),
    n_parameters = as.integer(.ferx_search_num(raw$n_parameters)),
    ofv          = .ferx_search_num(raw$ofv),
    criterion    = .ferx_search_num(raw$criterion),
    d_criterion  = .ferx_search_num(raw$d_criterion),
    rank         = as.integer(.ferx_search_num(raw$rank)),
    converged    = .ferx_search_lgl(raw$converged),
    passed       = as.logical(raw$passed),
    failures     = .ferx_search_chr(raw$failures),
    error        = .ferx_search_chr(raw$error),
    seconds      = .ferx_search_num(raw$seconds),
    selected     = as.logical(raw$selected),
    continued    = as.logical(raw$continued),
    reused       = as.logical(raw$reused),
    structure    = as.character(raw$structure),
    stringsAsFactors = FALSE
  )

  model_text <- stats::setNames(as.character(raw$model_text),
                                as.character(raw$model_id))

  # The winning model as a file: the run's own `final.ferx` when it wrote one
  # (the engine seeds the final estimates into it), otherwise a temporary copy
  # so the fit below has a model file to name and read.
  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-modelsearch-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path, data = as.character(raw$data))
  }

  final_row <- models[models$id == as.character(raw$final_id), , drop = FALSE]
  result <- list(
    models           = models,
    fit              = fit,
    model_text       = model_text,
    input_model      = as.character(raw$input_model),
    base_model_id    = as.character(raw$base_id),
    base_structure   = as.character(raw$base_structure),
    base_ofv         = .ferx_search_num(raw$base_ofv),
    base_criterion   = .ferx_search_num(raw$base_criterion),
    final_model      = as.character(raw$final_model),
    final_model_path = final_path,
    final_model_id   = as.character(raw$final_id),
    final_ofv        = if (nrow(final_row) == 1L) final_row$ofv else NA_real_,
    final_criterion  = .ferx_search_num(raw$final_criterion),
    criterion        = as.character(raw$criterion_label),
    algorithm        = as.character(raw$algorithm),
    iiv_strategy     = as.character(raw$iiv_strategy),
    n_layers         = as.integer(raw$n_layers),
    summary_text     = as.character(raw$summary),
    candidates       = .ferx_search_candidates(dir_arg),
    model            = as.character(raw$model),
    data             = as.character(raw$data),
    directory        = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config           = if (nzchar(config_path)) config_path else NA_character_,
    notes            = as.character(raw$notes),
    cancelled        = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_modelsearch", "ferx_search_result")
  result
}

#' @param x A \code{ferx_modelsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_modelsearch
#' @export
print.ferx_modelsearch <- function(x, digits = 4, ...) {
  cat("ferx structural model search (", x$algorithm, ", ranked on ", x$criterion,
      ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  cat("  Base:  ", x$base_model_id, " - ", x$base_structure, "\n", sep = "")
  cat(sprintf("  %s:   %s (base) -> %s (%s)\n",
              x$criterion,
              format(signif(x$base_criterion, digits)),
              format(signif(x$final_criterion, digits)),
              x$final_model_id))
  cat(sprintf("  %d model%s over %d layer%s\n",
              nrow(x$models), if (nrow(x$models) == 1L) "" else "s",
              x$n_layers, if (x$n_layers == 1L) "" else "s"))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the models below are what it reached.\n")
  }

  cat("\nRanked models (best first):\n")
  ranked <- x$models[!is.na(x$models$rank), , drop = FALSE]
  ranked <- ranked[order(ranked$rank), , drop = FALSE]
  .ferx_search_print_table(
    ranked[, c("rank", "id", "structure", "ofv", "criterion", "d_criterion",
               "converged", "passed"), drop = FALSE],
    digits
  )

  gated <- sum(!x$models$passed, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d model%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

#' @param object A \code{ferx_modelsearch} object.
#' @rdname ferx_modelsearch
#' @export
summary.ferx_modelsearch <- function(object, digits = 4, ...) {
  print(object, digits = digits)
  cat("\nEvery model:\n")
  .ferx_search_print_table(
    object$models[, c("id", "parent", "layer", "structure", "ofv", "criterion",
                      "d_criterion", "rank", "converged", "passed", "continued",
                      "reused", "failures", "error"), drop = FALSE],
    digits
  )
  invisible(object)
}
