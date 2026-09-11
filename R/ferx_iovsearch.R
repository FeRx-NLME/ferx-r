#' Inter-occasion variability search
#'
#' Runs an inter-occasion variability search - Pharmpy's \code{iovsearch} -
#' over a base model: which parameters carry a kappa, and whether the etas
#' those kappas sit beside are still worth keeping. Every candidate is built by
#' the engine's own model editor, fitted, gated, and ranked on the criterion
#' \code{rank} names (the BIC(random) unless the file says otherwise).
#'
#' @section Two stages, not one table:
#' An iovsearch is two searches in sequence. The full-IOV model - a kappa on
#' every candidate parameter - is fitted first; step 1 then removes kappas from
#' it, and step 2, taking step 1's winner as its parent, removes the etas of
#' the parameters that kept a kappa. Every row of \code{$models} carries the
#' step it was fitted in, and \code{$steps} is the engine's own per-step
#' ranking.
#'
#' @section What the base model needs:
#' The occasions are read by the engine from the base model's
#' \code{iov_column} (\code{[fit_options] iov_column = OCC}). A base without
#' one is refused by name: the dataset's occasions were never read, so no
#' candidate could be fitted. The base needs no \code{kappa} of its own - the
#' search adds them - and the column must carry at least two occasion values.
#'
#' @section Labels:
#' The engine's own \code{description} column is written in \emph{parameter}
#' names (\code{IIV([CL]+[V]);IOV([CL])}), which is Pharmpy's spelling and what
#' \code{models.csv} carries. The \code{structure}, \code{eta_labels},
#' \code{kappa_labels} and \code{kappa_block_labels} columns beside it are the
#' same structure in the model's \emph{declared} random-effect names
#' (\code{IIV([ETA_CL]+[ETA_V]);IOV([KAPPA_CL])}), read off each candidate's
#' own text. A random effect the model does not name falls back to
#' \code{OMEGA(i,i)} / \code{KAPPA<i>}, per the output label convention.
#'
#' @section Two entry forms:
#' Pass either \code{config} - a \code{.ferxsearch} file, which is the
#' reproducible artifact - or the inline arguments. The inline form is rendered
#' into the same configuration and handed to the engine's own loader, so the
#' two cannot disagree; \code{search_space} is MFL text, quoted verbatim.
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
#' @param search_space MFL text naming the parameters to try a kappa on, e.g.
#'   \code{"IOV?(*, exp)"} - the statement takes parameters and an effect, and
#'   only the exponential form is searchable. \code{NULL} takes Pharmpy's
#'   default: every parameter carrying a free eta.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param column The occasion column. \code{NULL} takes the base model's
#'   \code{iov_column}; a value that disagrees with it is an error rather than
#'   a silent override.
#' @param distribution How the added kappas are declared:
#'   \code{"same-as-iiv"} (the default - kappas blocked as their etas are),
#'   \code{"disjoint"} (one \code{kappa} line each), \code{"joint"} (one
#'   \code{block_kappa} over all of them) or \code{"explicit"} (the blocks
#'   \code{groups} names). \code{NULL} keeps the engine default.
#' @param groups The kappa blocks for \code{distribution = "explicit"}, as a
#'   list of character vectors of parameter names, e.g.
#'   \code{list(c("CL", "V"), "KA")}.
#' @param block_retries Extra starts per kappa beyond two in a candidate's
#'   largest \code{block_kappa}, on top of the run's starts. \code{NULL} keeps
#'   the engine default (2).
#' @param rank \code{[rank] type} - the criterion candidates are ranked and
#'   gated on. \code{NULL} keeps the tool default, the BIC(random), which is
#'   what \code{"bic"} means for this tool.
#' @param cutoff \code{[rank] cutoff}: the improvement on the criterion a
#'   candidate must show over its step's parent to replace it. \code{NULL}
#'   takes any improvement, which is Pharmpy's default.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the per-step journals, \code{models.csv},
#'   \code{models/<id>.ferx} and \code{final.ferx} are written. \code{NULL}
#'   keeps the run in memory, which also means it cannot be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}.
#' @param progress Print the engine's step progress to the console.
#'
#' @return An object of class \code{ferx_iovsearch}:
#'   \describe{
#'     \item{models}{The model table, one row per fitted model, in the engine's
#'       own column order: \code{id}, \code{parent}, \code{step},
#'       \code{description}, \code{etas}, \code{kappas},
#'       \code{kappa_blocks}, \code{n_parameters}, \code{ofv},
#'       \code{criterion}, \code{d_criterion}, \code{rank}, \code{converged},
#'       \code{passed}, \code{failures}, \code{error}, \code{starts},
#'       \code{seconds}, \code{selected} - followed by \code{step_kind},
#'       \code{eta_labels}, \code{kappa_labels}, \code{kappa_block_labels} and
#'       \code{structure} (the same structure in the model's declared random
#'       effect names).}
#'     \item{steps}{The engine's per-step ranking, one row per model ranked in
#'       a step: \code{step}, \code{kind}, \code{parent}, \code{id},
#'       \code{criterion}, \code{d_criterion} (positive is better),
#'       \code{rank}, \code{best}.}
#'     \item{fit}{The winning model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{model_text}{Every candidate's model text, named by model id, so a
#'       structure the search rejected can be read or refitted without
#'       re-running it. \code{final_model_path} is the \code{final.ferx}
#'       written beside the run.}
#'     \item{input_model, input_structure, input_description, input_ofv,
#'       input_criterion}{The model the search started from, as fitted.}
#'     \item{final_model_id, final_structure, final_description,
#'       final_etas, final_kappas, final_kappa_blocks, final_criterion}{The
#'       winning model, the random effects it carries and the blocks over the
#'       kappas.}
#'     \item{criterion, distribution, column, n_steps}{What the search ranked
#'       on, how the kappas were declared and which column carried the
#'       occasions.}
#'     \item{options}{The search as the engine read it: \code{distribution},
#'       \code{column}, \code{groups}, \code{block_retries}, \code{starts} and
#'       \code{cutoff}.}
#'     \item{candidates}{The runner's candidate table when the run wrote one -
#'       see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once, and whether it
#'       was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin_iov")
#'
#' # Inline: no file to author for a one-off
#' res <- ferx_iovsearch(
#'   model     = ex$model,
#'   data      = ex$data,
#'   directory = "iovsearch-run-1"
#' )
#' res
#' res$models[, c("id", "step", "structure", "criterion", "rank")]
#' summary(res)
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res2 <- ferx_iovsearch(config = ex$search, directory = "iovsearch-run-2")
#' }
#'
#' @seealso \code{\link{ferx_iivsearch}}, \code{\link{ferx_modelsearch}},
#'   \code{\link{ferx_search_space}}, \code{\link{ferx_search_results}}
#' @family search
#' @export
ferx_iovsearch <- function(model = NULL,
                           data = NULL,
                           search_space = NULL,
                           config = NULL,
                           column = NULL,
                           distribution = NULL,
                           groups = NULL,
                           block_retries = NULL,
                           rank = NULL,
                           cutoff = NULL,
                           threads = NULL,
                           retries = NULL,
                           directory = NULL,
                           resume = FALSE,
                           progress = interactive()) {
  what <- "ferx_iovsearch"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      search_space = search_space,
      column = column,
      distribution = distribution,
      groups = groups,
      block_retries = block_retries,
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
    # An IOV space is optional: without one the engine tries every parameter
    # that carries a free eta, which is Pharmpy's own default.
    mfl <- .ferx_search_space_text(search_space, what, required = FALSE)
  }

  if (!is.null(distribution)) {
    distribution <- match.arg(
      distribution,
      c("same-as-iiv", "disjoint", "joint", "explicit")
    )
  }
  if (!is.null(column) &&
        (!is.character(column) || length(column) != 1L || is.na(column))) {
    stop(what, ": `column` must be a single column name or NULL")
  }
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_iovsearch(
    config_path   = config_path,
    model_path    = paths$model,
    data_path     = paths$data,
    mfl           = mfl,
    column        = column %||% "",
    distribution  = distribution %||% "",
    groups        = .ferx_iov_groups(groups, what),
    block_retries = .ferx_search_count(block_retries, "block_retries", what,
                                       below = -1L, min = 0L),
    rank          = rank %||% "",
    rank_cutoff   = .ferx_search_scalar(cutoff, "cutoff", what),
    threads       = .ferx_search_count(threads, "threads", what, below = 0L),
    retries       = .ferx_search_count(retries, "retries", what,
                                       below = -1L, min = 0L),
    resume        = .ferx_search_bool(resume, "resume", what),
    directory     = dir_arg,
    progress      = .ferx_search_bool(progress, "progress", what)
  )

  # The engine's 19 columns in the engine's order, then the five this surface
  # adds: the step the row belongs to and the same structure in the model's
  # declared random-effect names.
  models <- data.frame(
    id                 = as.character(raw$id),
    parent             = .ferx_search_chr(raw$parent),
    step               = as.integer(raw$step),
    description        = as.character(raw$description),
    etas               = .ferx_search_chr(raw$etas),
    kappas             = .ferx_search_chr(raw$kappas),
    kappa_blocks       = .ferx_search_chr(raw$kappa_blocks),
    n_parameters       = as.integer(.ferx_search_num(raw$n_parameters)),
    ofv                = .ferx_search_num(raw$ofv),
    criterion          = .ferx_search_num(raw$criterion),
    d_criterion        = .ferx_search_num(raw$d_criterion),
    rank               = as.integer(.ferx_search_num(raw$rank)),
    converged          = .ferx_search_lgl(raw$converged),
    passed             = as.logical(raw$passed),
    failures           = .ferx_search_chr(raw$failures),
    error              = .ferx_search_chr(raw$error),
    starts             = as.integer(raw$starts),
    seconds            = .ferx_search_num(raw$seconds),
    selected           = as.logical(raw$selected),
    step_kind          = .ferx_search_chr(raw$step_kind),
    eta_labels         = .ferx_search_chr(raw$eta_labels),
    kappa_labels       = .ferx_search_chr(raw$kappa_labels),
    kappa_block_labels = .ferx_search_chr(raw$kappa_block_labels),
    structure          = as.character(raw$structure),
    stringsAsFactors = FALSE
  )

  steps <- data.frame(
    step        = as.integer(raw$s_step),
    kind        = as.character(raw$s_kind),
    parent      = as.character(raw$s_parent),
    id          = as.character(raw$s_id),
    criterion   = .ferx_search_num(raw$s_criterion),
    d_criterion = .ferx_search_num(raw$s_d_criterion),
    rank        = as.integer(.ferx_search_num(raw$s_rank)),
    best        = as.logical(raw$s_best),
    stringsAsFactors = FALSE
  )

  model_text <- stats::setNames(as.character(raw$model_text),
                                as.character(raw$model_id))

  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-iovsearch-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path, data = as.character(raw$data))
  }

  final_row <- models[models$id == as.character(raw$final_id), , drop = FALSE]
  result <- list(
    models             = models,
    steps              = steps,
    fit                = fit,
    model_text         = model_text,
    input_model        = as.character(raw$input_model),
    input_structure    = as.character(raw$input_structure),
    input_description  = as.character(raw$input_description),
    input_ofv          = .ferx_search_num(raw$input_ofv),
    input_criterion    = .ferx_search_num(raw$input_criterion),
    final_model        = as.character(raw$final_model),
    final_model_path   = final_path,
    final_model_id     = as.character(raw$final_id),
    final_structure    = as.character(raw$final_structure),
    final_description  = as.character(raw$final_description),
    final_etas         = as.character(raw$final_etas),
    final_kappas       = as.character(raw$final_kappas),
    final_kappa_blocks = .ferx_iiv_block_list(as.character(raw$final_kappa_blocks)),
    final_ofv          = if (nrow(final_row) == 1L) final_row$ofv else NA_real_,
    final_criterion    = .ferx_search_num(raw$final_criterion),
    criterion          = as.character(raw$criterion_label),
    distribution       = as.character(raw$distribution),
    column             = .ferx_search_chr(raw$column),
    n_steps            = as.integer(raw$n_steps),
    options            = list(
      distribution  = as.character(raw$distribution),
      column        = .ferx_search_chr(raw$column),
      groups        = lapply(as.character(raw$groups),
                             function(g) trimws(strsplit(g, ",", fixed = TRUE)[[1L]])),
      block_retries = as.integer(raw$block_retries),
      starts        = as.integer(raw$starts_base),
      cutoff        = .ferx_search_num(raw$opt_cutoff)
    ),
    summary_text       = as.character(raw$summary),
    # The runner directories this run wrote, in the order the engine wrote
    # them: the input, the full-IOV model, and the two removal steps. Taken
    # from the rows the run itself produced, so a shorter run in a directory an
    # earlier one used cannot inherit its candidate tables (#336 review).
    candidates         = .ferx_search_candidates(
      dir_arg,
      c("input", "iov-all",
        sprintf("step-%d", sort(unique(models$step[models$step > 0L]))))
    ),
    model              = as.character(raw$model),
    data               = as.character(raw$data),
    directory          = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config             = if (nzchar(config_path)) config_path else NA_character_,
    notes              = as.character(raw$notes),
    cancelled          = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_iovsearch", "ferx_search_result")
  result
}

# `groups` as the binding takes it: one comma-separated parameter list per
# block. A closed shape checked here rather than becoming a TOML array the
# loader refuses several layers down.
.ferx_iov_groups <- function(groups, what) {
  if (is.null(groups)) return(character(0))
  if (is.character(groups)) groups <- as.list(groups)
  if (!is.list(groups) || length(groups) == 0L) {
    stop(what, ": `groups` must be a list of character vectors of parameter ",
         "names, e.g. list(c(\"CL\", \"V\"), \"KA\")")
  }
  vapply(groups, function(g) {
    if (!is.character(g) || length(g) == 0L || anyNA(g) || !all(nzchar(g))) {
      stop(what, ": `groups` must be a list of character vectors of parameter ",
           "names, e.g. list(c(\"CL\", \"V\"), \"KA\")")
    }
    paste(trimws(g), collapse = ",")
  }, character(1))
}

#' @param x A \code{ferx_iovsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_iovsearch
#' @export
print.ferx_iovsearch <- function(x, digits = 4, ...) {
  cat("ferx inter-occasion variability search (", x$distribution,
      ", ranked on ", x$criterion, ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  cat("  Occasions: ", x$column %||% "(from the model)", "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  cat("  Input: ", x$input_structure, "\n", sep = "")
  cat(sprintf("  %s:   %s (input) -> %s (%s)\n",
              x$criterion,
              format(signif(x$input_criterion, digits)),
              format(signif(x$final_criterion, digits)),
              x$final_model_id))
  cat(sprintf("  %d model%s over %d step%s\n",
              nrow(x$models), if (nrow(x$models) == 1L) "" else "s",
              x$n_steps, if (x$n_steps == 1L) "" else "s"))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the models below are what it reached.\n")
  }

  # One table per step: removing kappas and removing etas are two decisions,
  # taken against two parents.
  for (s in unique(x$steps$step)) {
    one <- x$steps[x$steps$step == s, , drop = FALSE]
    cat(sprintf("\nStep %d (%s), parent %s:\n", s, one$kind[1L], one$parent[1L]))
    one$structure <- vapply(one$id, function(i) .ferx_iov_structure_of(x, i),
                            character(1))
    .ferx_search_print_table(
      one[, c("rank", "id", "structure", "criterion", "d_criterion", "best"),
          drop = FALSE],
      digits
    )
  }

  gated <- sum(!x$models$passed, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d model%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }

  cat("\nSelected variability structure (", x$final_model_id, "):\n", sep = "")
  cat("  IIV:  ",
      if (any(nzchar(x$final_etas))) paste(x$final_etas, collapse = ", ") else "(none)",
      "\n", sep = "")
  cat("  IOV:  ",
      if (any(nzchar(x$final_kappas))) paste(x$final_kappas, collapse = ", ") else "(none)",
      "\n", sep = "")
  pairs <- .ferx_iiv_pairs(x$final_kappa_blocks)
  if (length(pairs) > 0L) {
    cat("  IOV correlated:  ", paste(pairs, collapse = ", "), "\n", sep = "")
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

# The labelled structure of a model id, from the table the run returned.
.ferx_iov_structure_of <- function(x, id) {
  hit <- x$models$structure[x$models$id == id]
  if (length(hit) == 0L) "" else hit[1L]
}

#' @param object A \code{ferx_iovsearch} object.
#' @rdname ferx_iovsearch
#' @export
summary.ferx_iovsearch <- function(object, digits = 4, ...) {
  print(object, digits = digits)

  cat("\nStructures not selected, and why:\n")
  rejected <- object$models[!object$models$selected & object$models$step > 0L, ,
                            drop = FALSE]
  if (nrow(rejected) == 0L) {
    cat("  (none - every structure fitted was the winner of its step)\n")
  } else {
    rejected$reason <- ifelse(
      !is.na(rejected$error), rejected$error,
      ifelse(!is.na(rejected$failures), rejected$failures,
             ifelse(is.na(rejected$rank),
                    "not ranked (no usable fit)",
                    ifelse(rejected$rank == 1L,
                           "best of its step, but not the model the search ended on",
                           sprintf("ranked %d on the %s of its step",
                                   rejected$rank, object$criterion)))))
    .ferx_search_print_table(
      rejected[, c("step", "step_kind", "id", "structure", "criterion",
                   "d_criterion", "rank", "converged", "passed", "reason"),
               drop = FALSE],
      digits
    )
  }

  cat("\nEvery model:\n")
  .ferx_search_print_table(
    object$models[, c("id", "parent", "step", "structure", "n_parameters",
                      "ofv", "criterion", "d_criterion", "rank", "starts",
                      "converged", "passed", "failures", "error"),
                  drop = FALSE],
    digits
  )
  invisible(object)
}
