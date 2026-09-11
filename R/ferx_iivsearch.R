#' Variability-structure search
#'
#' Runs a variability-structure search - Pharmpy's \code{iivsearch} - over a
#' base model: which parameters carry an eta, and which of those etas are
#' correlated. Every candidate is built by the engine's own model editor,
#' fitted, gated, and ranked on the criterion \code{rank} names (the BIC(iiv)
#' unless the file says otherwise).
#'
#' @section Two stages, not one table:
#' An iivsearch is two searches in sequence. The first decides the
#' \emph{number of etas} (\code{step_kind = "no_of_etas"}); the second, taking
#' the winner of the first as its parent, decides the \emph{block structure}
#' over them (\code{step_kind = "block_structure"}); a last comparison puts the
#' winner beside the input model (\code{"compare_to_input"}). Every row of
#' \code{$models} carries the stage it was fitted in, and \code{$steps} is the
#' engine's own per-stage ranking, so a search that added an eta and then
#' blocked two of them reads as the two decisions it was.
#'
#' \code{algorithm = "simultaneous_stepwise"} is the exception: it decides the
#' block structure as it adds each eta, so it has no second stage and refuses
#' \code{correlation_algorithm}.
#'
#' A candidate whose largest block has more than two etas is fitted from more
#' starting points than a diagonal one (\code{block_retries} per eta beyond
#' two, on top of \code{retries + 1}), because a block covariance is the part
#' of an omega most likely to land in a local optimum. The starts each
#' candidate actually got are on its row.
#'
#' @section Labels:
#' The engine's own \code{description} column is written in \emph{parameter}
#' names (\code{[CL,V]+[KA]}), which is Pharmpy's spelling and what
#' \code{models.csv} carries. The \code{structure}, \code{eta_labels} and
#' \code{block_labels} columns beside it are the same structure in the model's
#' \emph{declared} random-effect names (\code{[ETA_CL,ETA_V]+[ETA_KA]}), read
#' off each candidate's own text - so a model that calls its eta something
#' other than \code{ETA_<P>} is labelled as it is written. An eta the model
#' does not name falls back to \code{OMEGA(i,i)}, per the output label
#' convention.
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
#' @section What the space may say:
#' \code{IIV} and \code{COVARIANCE} statements only. Both take two arguments:
#' \code{IIV?(@PK, exp)} makes every PK parameter's eta exploratory,
#' \code{IIV(CL, exp)} keeps CL's whatever the search decides,
#' \code{COVARIANCE?(IIV, [CL, V])} offers the two as a block and
#' \code{COVARIANCE(IIV, [CL, V])} forces it. Only the exponential form is
#' searchable, since it is the one the eta edits can reverse. A structural or
#' covariate statement is another tool's space and is refused by name before
#' any fit starts.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the search runs on, so like \code{model} it
#'   cannot be given beside \code{config} - the file's own \code{data} key
#'   says which dataset that file searches.
#' @param search_space MFL text naming the variability space, e.g.
#'   \code{"IIV?(@PK, exp); COVARIANCE?(IIV, *)"}.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param algorithm \code{"top_down_exhaustive"} (the default),
#'   \code{"bottom_up_stepwise"}, \code{"simultaneous_stepwise"} or
#'   \code{"skip"} (no eta-count stage; the block stage alone). \code{NULL}
#'   keeps the engine default.
#' @param correlation_algorithm \code{"top_down_exhaustive"} or \code{"skip"}.
#'   \code{NULL} keeps the engine's derivation: the block stage runs after any
#'   algorithm but \code{"simultaneous_stepwise"}.
#' @param as_fullblock Pharmpy's \code{as_fullblock}: a bottom-up candidate
#'   blocks every eta it carries rather than adding the new one diagonally.
#'   \code{NULL} keeps the engine default (\code{FALSE}).
#' @param block_retries Extra starts per eta beyond two in a candidate's
#'   largest block, on top of the run's starts. \code{NULL} keeps the engine
#'   default (2).
#' @param rank \code{[rank] type} - the criterion candidates are ranked and
#'   gated on. \code{NULL} keeps the tool default, the BIC(iiv), which is what
#'   \code{"bic"} means for this tool.
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
#' @return An object of class \code{ferx_iivsearch}:
#'   \describe{
#'     \item{models}{The model table, one row per fitted model, in the engine's
#'       own column order: \code{id}, \code{parent}, \code{step},
#'       \code{description}, \code{etas}, \code{blocks}, \code{n_parameters},
#'       \code{ofv}, \code{criterion}, \code{d_criterion}, \code{rank},
#'       \code{converged}, \code{passed}, \code{failures}, \code{error},
#'       \code{starts}, \code{seconds}, \code{selected} - followed by
#'       \code{step_kind} (the stage the row belongs to), \code{eta_labels},
#'       \code{block_labels} and \code{structure} (the same structure in the
#'       model's declared eta names).}
#'     \item{steps}{The engine's per-stage ranking, one row per model ranked in
#'       a step: \code{step}, \code{kind}, \code{parent}, \code{id},
#'       \code{criterion}, \code{d_criterion} (positive is better),
#'       \code{rank}, \code{best}.}
#'     \item{fit}{The winning model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{model_text}{Every candidate's model text, named by model id, so a
#'       structure the search rejected can be read or refitted without
#'       re-running it. \code{final_model_path} is the \code{final.ferx}
#'       written beside the run.}
#'     \item{input_model, input_structure, input_description}{The model the
#'       search started from, and its structure in declared eta names and in
#'       the engine's parameter names.}
#'     \item{base_model_id, base_ofv, base_criterion}{The root of the search:
#'       \code{"input"}, or \code{"base"} when one had to be derived.}
#'     \item{final_model_id, final_structure, final_description, final_etas,
#'       final_blocks, final_ofv, final_criterion}{The winning model, its
#'       structure, the etas it carries and the blocks over them.}
#'     \item{criterion, algorithm, correlation_algorithm, block_stage,
#'       n_steps}{What the search ranked on and how it searched.}
#'     \item{options}{The search as the engine read it: \code{algorithm},
#'       \code{correlation_algorithm}, \code{as_fullblock},
#'       \code{block_retries}, \code{starts} and \code{cutoff}.}
#'     \item{candidates}{The runner's candidate table when the run wrote one -
#'       see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once, and whether it
#'       was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("one_cpt_iv")
#'
#' # Inline: no file to author for a one-off
#' res <- ferx_iivsearch(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
#'   directory    = "iivsearch-run-1"
#' )
#' res
#' res$models[, c("id", "step_kind", "structure", "criterion", "rank")]
#' res$steps
#' summary(res)
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res2 <- ferx_iivsearch(config = ex$search, directory = "iivsearch-run-2")
#' }
#'
#' @seealso \code{\link{ferx_iovsearch}}, \code{\link{ferx_modelsearch}},
#'   \code{\link{ferx_search_space}}, \code{\link{ferx_search_results}}
#' @family search
#' @export
ferx_iivsearch <- function(model = NULL,
                           data = NULL,
                           search_space = NULL,
                           config = NULL,
                           algorithm = NULL,
                           correlation_algorithm = NULL,
                           as_fullblock = NULL,
                           block_retries = NULL,
                           rank = NULL,
                           cutoff = NULL,
                           threads = NULL,
                           retries = NULL,
                           directory = NULL,
                           resume = FALSE,
                           progress = interactive()) {
  what <- "ferx_iivsearch"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      search_space = search_space,
      algorithm = algorithm,
      correlation_algorithm = correlation_algorithm,
      as_fullblock = as_fullblock,
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
    mfl <- .ferx_search_space_text(search_space, what)
  }

  if (!is.null(algorithm)) {
    algorithm <- match.arg(
      algorithm,
      c("top_down_exhaustive", "bottom_up_stepwise", "simultaneous_stepwise",
        "skip")
    )
  }
  if (!is.null(correlation_algorithm)) {
    correlation_algorithm <- match.arg(
      correlation_algorithm,
      c("top_down_exhaustive", "skip")
    )
  }
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_iivsearch(
    config_path           = config_path,
    model_path            = paths$model,
    data_path             = paths$data,
    mfl                   = mfl,
    algorithm             = algorithm %||% "",
    correlation_algorithm = correlation_algorithm %||% "",
    as_fullblock          = .ferx_search_flag3(as_fullblock, "as_fullblock", what),
    block_retries         = .ferx_search_count(block_retries, "block_retries", what,
                                               below = -1L, min = 0L),
    rank                  = rank %||% "",
    rank_cutoff           = .ferx_search_scalar(cutoff, "cutoff", what),
    threads               = .ferx_search_count(threads, "threads", what, below = 0L),
    retries               = .ferx_search_count(retries, "retries", what,
                                               below = -1L, min = 0L),
    resume                = .ferx_search_bool(resume, "resume", what),
    directory             = dir_arg,
    progress              = .ferx_search_bool(progress, "progress", what)
  )

  # The engine's 18 columns in the engine's order, then the four this surface
  # adds: the stage the row belongs to (an iivsearch is two searches, and a
  # table that lost which is which would report one it did not run) and the
  # same structure in the model's declared eta names.
  models <- data.frame(
    id           = as.character(raw$id),
    parent       = .ferx_search_chr(raw$parent),
    step         = as.integer(raw$step),
    description  = as.character(raw$description),
    etas         = .ferx_search_chr(raw$etas),
    blocks       = .ferx_search_chr(raw$blocks),
    n_parameters = as.integer(.ferx_search_num(raw$n_parameters)),
    ofv          = .ferx_search_num(raw$ofv),
    criterion    = .ferx_search_num(raw$criterion),
    d_criterion  = .ferx_search_num(raw$d_criterion),
    rank         = as.integer(.ferx_search_num(raw$rank)),
    converged    = .ferx_search_lgl(raw$converged),
    passed       = as.logical(raw$passed),
    failures     = .ferx_search_chr(raw$failures),
    error        = .ferx_search_chr(raw$error),
    starts       = as.integer(raw$starts),
    seconds      = .ferx_search_num(raw$seconds),
    selected     = as.logical(raw$selected),
    step_kind    = .ferx_search_chr(raw$step_kind),
    eta_labels   = .ferx_search_chr(raw$eta_labels),
    block_labels = .ferx_search_chr(raw$block_labels),
    structure    = as.character(raw$structure),
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

  # The winning model as a file: the run's own `final.ferx` when it wrote one
  # (the engine seeds the final estimates into it), otherwise a temporary copy
  # so the fit below has a model file to name and read.
  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-iivsearch-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path, data = as.character(raw$data))
  }

  final_row <- models[models$id == as.character(raw$final_id), , drop = FALSE]
  result <- list(
    models                = models,
    steps                 = steps,
    fit                   = fit,
    model_text            = model_text,
    input_model           = as.character(raw$input_model),
    input_structure       = as.character(raw$input_structure),
    input_description     = as.character(raw$input_description),
    base_model_id         = as.character(raw$base_id),
    base_ofv              = .ferx_search_num(raw$base_ofv),
    base_criterion        = .ferx_search_num(raw$base_criterion),
    final_model           = as.character(raw$final_model),
    final_model_path      = final_path,
    final_model_id        = as.character(raw$final_id),
    final_structure       = as.character(raw$final_structure),
    final_description     = as.character(raw$final_description),
    final_etas            = as.character(raw$final_etas),
    final_blocks          = .ferx_iiv_block_list(as.character(raw$final_blocks)),
    final_ofv             = if (nrow(final_row) == 1L) final_row$ofv else NA_real_,
    final_criterion       = .ferx_search_num(raw$final_criterion),
    criterion             = as.character(raw$criterion_label),
    algorithm             = as.character(raw$algorithm),
    correlation_algorithm = .ferx_search_chr(raw$correlation_algorithm),
    block_stage           = isTRUE(raw$block_stage),
    n_steps               = as.integer(raw$n_steps),
    options               = list(
      algorithm             = as.character(raw$algorithm),
      correlation_algorithm = .ferx_search_chr(raw$correlation_algorithm),
      as_fullblock          = isTRUE(raw$as_fullblock),
      block_retries         = as.integer(raw$block_retries),
      starts                = as.integer(raw$starts_base),
      cutoff                = .ferx_search_num(raw$opt_cutoff)
    ),
    summary_text          = as.character(raw$summary),
    # The runner directories this run wrote - the root fit, the derived base
    # and one per step it actually fitted - rather than whatever the directory
    # holds: a shorter run in a directory an earlier one used would otherwise
    # inherit its candidate tables (#336 review).
    candidates            = .ferx_search_candidates(
      dir_arg,
      c("input", "base",
        sprintf("step-%d", sort(unique(models$step[models$step > 0L]))))
    ),
    model                 = as.character(raw$model),
    data                  = as.character(raw$data),
    directory             = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config                = if (nzchar(config_path)) config_path else NA_character_,
    notes                 = as.character(raw$notes),
    cancelled             = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_iivsearch", "ferx_search_result")
  result
}

# The block column as a list of blocks: `;` between blocks, `,` within one.
# A structure with no block is `character(0)`, not `""`.
.ferx_iiv_block_list <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(list())
  lapply(strsplit(x, ";", fixed = TRUE)[[1L]],
         function(b) trimws(strsplit(b, ",", fixed = TRUE)[[1L]]))
}

# The correlated pairs of a block, in the label convention's own spelling:
# `ETA_V ~ ETA_CL`, the off-diagonal form print.ferx_fit uses.
.ferx_iiv_pairs <- function(blocks) {
  out <- character(0)
  for (b in blocks) {
    if (length(b) < 2L) next
    for (i in seq_along(b)) {
      for (j in seq_len(i - 1L)) {
        out <- c(out, sprintf("%s ~ %s", b[i], b[j]))
      }
    }
  }
  out
}

#' @param x A \code{ferx_iivsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_iivsearch
#' @export
print.ferx_iivsearch <- function(x, digits = 4, ...) {
  cat("ferx variability-structure search (", x$algorithm,
      if (x$block_stage) " + block structure" else "",
      ", ranked on ", x$criterion, ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  cat("  Input: ", .ferx_iiv_show(x$input_structure), "\n", sep = "")
  if (!identical(x$base_model_id, "input")) {
    cat("  Base:  ", x$base_model_id, " (criterion ",
        format(signif(x$base_criterion, digits)), ")\n", sep = "")
  }
  cat(sprintf("  %s:   %s (%s) -> %s (%s)\n",
              x$criterion,
              format(signif(x$base_criterion, digits)), x$base_model_id,
              format(signif(x$final_criterion, digits)), x$final_model_id))
  cat(sprintf("  %d model%s over %d step%s\n",
              nrow(x$models), if (nrow(x$models) == 1L) "" else "s",
              x$n_steps, if (x$n_steps == 1L) "" else "s"))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the models below are what it reached.\n")
  }

  # One table per stage, never one table with the stage lost: the number of
  # etas and the block structure are two decisions, taken against two parents.
  for (k in unique(x$steps$kind)) {
    rows <- x$steps[x$steps$kind == k, , drop = FALSE]
    for (s in unique(rows$step)) {
      one <- rows[rows$step == s, , drop = FALSE]
      cat(sprintf("\nStep %d (%s), parent %s (%s):\n", s, k, one$parent[1L],
                  .ferx_iiv_show(.ferx_iiv_structure_of(x, one$parent[1L]))))
      one$structure <- vapply(one$id, function(i) .ferx_iiv_structure_of(x, i),
                              character(1))
      .ferx_search_print_table(
        one[, c("rank", "id", "structure", "criterion", "d_criterion", "best"),
            drop = FALSE],
        digits
      )
    }
  }

  gated <- sum(!x$models$passed, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d model%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }

  cat("\nSelected variability structure (", x$final_model_id, "):\n", sep = "")
  if (length(x$final_etas) == 0L || !any(nzchar(x$final_etas))) {
    cat("  (no eta - every parameter is fixed across subjects)\n")
  } else {
    cat("  IIV:         ", paste(x$final_etas, collapse = ", "), "\n", sep = "")
  }
  pairs <- .ferx_iiv_pairs(x$final_blocks)
  if (length(pairs) > 0L) {
    cat("  Correlated:  ", paste(pairs, collapse = ", "), "\n", sep = "")
  } else if (any(nzchar(x$final_etas))) {
    cat("  Correlated:  (none - every eta is diagonal)\n")
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

# The labelled structure of a model id, from the table the run returned.
.ferx_iiv_structure_of <- function(x, id) {
  hit <- x$models$structure[x$models$id == id]
  if (length(hit) == 0L) "" else hit[1L]
}

# An empty structure prints as what it is, not as an empty string.
.ferx_iiv_show <- function(s) {
  if (length(s) != 1L || is.na(s) || !nzchar(s)) "(no eta)" else s
}

#' @param object A \code{ferx_iivsearch} object.
#' @rdname ferx_iivsearch
#' @export
summary.ferx_iivsearch <- function(object, digits = 4, ...) {
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
    object$models[, c("id", "parent", "step", "step_kind", "structure",
                      "n_parameters", "ofv", "criterion", "d_criterion",
                      "rank", "starts", "converged", "passed", "failures",
                      "error"), drop = FALSE],
    digits
  )
  invisible(object)
}
