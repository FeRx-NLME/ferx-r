#' Automatic model development (AMD)
#'
#' Runs the whole model-development pipeline - Pharmpy's \code{amd} - over one
#' model and one search space: the structural model, the variability structure,
#' the residual-error model, inter-occasion variability, allometric scaling and
#' the covariate model, each decided by the tool that decides it
#' (\code{\link{ferx_modelsearch}}, \code{\link{ferx_iivsearch}},
#' \code{\link{ferx_ruvsearch}}, \code{\link{ferx_iovsearch}},
#' \code{\link{ferx_allometry}}, \code{\link{ferx_covsearch}}).
#'
#' Every step starts from the model the previous step selected, seeded with its
#' estimates, so a step never re-derives what the one before it decided. One
#' search space describes them all: the engine partitions it by statement kind
#' and hands each tool only the statements it can read
#' (\code{ABSORPTION} / \code{PERIPHERALS} / \code{LAGTIME} to the structural
#' step, \code{IIV} / \code{COVARIANCE(IIV, ...)} to the variability step,
#' \code{IOV} to the occasion step, \code{ALLOMETRY} to the scaling step,
#' \code{COVARIATE} to the covariate step). A step the space says nothing about
#' is skipped, and the step table says so in as many words.
#'
#' A \code{[rank]} criterion is narrowed the same way: the two steps that select
#' by a likelihood-ratio test (\code{ruvsearch}, \code{covsearch}) keep their
#' own p-values, and the steps that rank candidates get the criterion.
#'
#' @section The report is the product:
#' A pipeline that shows only its final model cannot be audited: a step whose
#' winner beat its siblings because they stalled reads exactly like a step that
#' found a real improvement. So the strictness verdict
#' (\code{passed} / \code{failures}) and the termination status
#' (\code{converged}) are columns at both levels - on every candidate of every
#' step, and on the model each step selected - beside the criterion, the
#' \code{dOFV} and the wall clock. \code{print()} shows the step table;
#' \code{summary()} adds every step's candidates.
#'
#' The same two tables are written to \code{directory} as \code{steps.csv} and
#' \code{candidates.csv}, with one subdirectory per step holding that tool's
#' own fuller record, and \code{final.ferx} at the top.
#'
#' @section Strategies:
#' \describe{
#'   \item{\code{"default"}}{structural, IIV, residual, IOV, allometry,
#'     covariates - Pharmpy's own order.}
#'   \item{\code{"reevaluation"}}{the default order, then IIV and residual
#'     again: both were decided before the model had its covariates.}
#'   \item{\code{"SIR"}}{structural, IIV, residual.}
#'   \item{\code{"SRI"}}{structural, residual, IIV.}
#'   \item{\code{"RSI"}}{residual, structural, IIV.}
#' }
#' \code{skip} leaves a named step out of whichever order is chosen.
#'
#' @section Retries:
#' \code{retries} says how many perturbed restarts every candidate of every
#' step gets on top of the fit from its exact initial estimates.
#' \code{retries_on} is the separate pass over a model a step has already
#' selected, which refits it from perturbed starts to check that its optimum is
#' the one it landed on: \code{"all_final"} (the default) runs that pass after
#' each step, \code{"final"} only on the pipeline's final model, and
#' \code{"skip"} not at all. The pass never adopts a fit the strictness gate
#' rejected, and its row is in the candidate table with \code{tool = "retries"}.
#'
#' @section Two entry forms:
#' Pass either \code{config} - a \code{.ferxsearch} file, which is the
#' reproducible artifact - or the inline arguments. The inline form is rendered
#' into the same configuration and handed to the engine's own loader, so the
#' two cannot disagree. Only the file can carry a per-tool section
#' (\code{[modelsearch]}, \code{[iivsearch]}, \code{[covsearch]}, ...); the
#' inline form runs each step at its own defaults.
#'
#' The run arguments (\code{threads}, \code{retries}, \code{resume},
#' \code{directory}, \code{progress}) say how to run a pipeline rather than
#' what to search, and are available to both forms.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the pipeline runs on, so like \code{model} it
#'   cannot be given beside \code{config}.
#' @param search_space MFL search space, as a string or a character vector of
#'   lines. Required in the inline form.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param strategy Step order: \code{"default"}, \code{"reevaluation"},
#'   \code{"SIR"}, \code{"SRI"} or \code{"RSI"}. \code{NULL} keeps the engine
#'   default (\code{"default"}).
#' @param retries_on Which selected models get the perturbed-restart pass:
#'   \code{"all_final"}, \code{"final"} or \code{"skip"}. \code{NULL} keeps the
#'   engine default (\code{"all_final"}).
#' @param skip Steps left out of the pipeline, any of \code{"structural"},
#'   \code{"iivsearch"}, \code{"residual"}, \code{"iovsearch"},
#'   \code{"allometry"}, \code{"covariates"}. A step the space says nothing
#'   about is skipped anyway; this is for one the space does describe.
#' @param rank The \code{[rank] type} the ranking steps use, e.g.
#'   \code{"bic"} or \code{"ofv"}. \code{NULL} keeps each tool's own default,
#'   which is a BIC.
#' @param cutoff The \code{[rank] cutoff}. \code{NULL} keeps the default.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the fit from the
#'   exact initials. \code{NULL} keeps the engine default.
#' @param directory Where the per-step directories, \code{steps.csv},
#'   \code{candidates.csv} and \code{final.ferx} are written. \code{NULL} runs
#'   the pipeline in a temporary directory that is removed when the call
#'   returns - the tables come back on the object either way, but nothing is
#'   left to read afterwards and the run cannot be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}.
#' @param progress Print the engine's step progress to the console.
#'
#' @return An object of class \code{ferx_amd}:
#'   \describe{
#'     \item{steps}{The step table, one row per planned step - skipped ones
#'       included - in the engine's own column order, with the pipeline
#'       position beside the step name (a \code{reevaluation} run has two rows
#'       called \code{iivsearch}): \code{index}, \code{step}, \code{tool},
#'       \code{rerun}, \code{directory}, \code{status} (\code{ran} /
#'       \code{skipped} / \code{failed}), \code{reason}, \code{criterion},
#'       \code{value_before}, \code{value_after}, \code{d_value},
#'       \code{ofv_before}, \code{ofv_after}, \code{d_ofv}, \code{candidates},
#'       \code{selected}, \code{seconds}, \code{converged}, \code{passed} and
#'       \code{notes}. \code{converged} and \code{passed} are the termination
#'       status and the strictness verdict of the model the step selected,
#'       taken from that model's row in \code{candidates}; they are \code{NA}
#'       for a step that selected no candidate of its own and carried the
#'       model it was handed forward, whose verdict is on the row of the step
#'       that produced it.}
#'     \item{candidates}{Every candidate of every step, in the engine's own
#'       column order: \code{step}, \code{tool}, \code{id}, \code{parent},
#'       \code{description}, \code{criterion}, \code{value}, \code{d_value},
#'       \code{ofv}, \code{d_ofv}, \code{rank}, \code{converged},
#'       \code{passed}, \code{failures}, \code{error}, \code{note},
#'       \code{seconds}, \code{selected}. \code{criterion} names what
#'       \code{value} is on, which is not always the step's own: a ruvsearch
#'       pre-screen candidate is fitted to the parent's CWRES, so its number is
#'       on the CWRES scale and must not be compared with a data OFV.
#'       \code{tool} is \code{"start"} for the pipeline's own first fit and
#'       \code{"retries"} for a perturbed-restart pass.}
#'     \item{tools}{One entry per step that ran - a step that failed included,
#'       carrying its reason and an empty candidate table - named by position
#'       and step. Each holds what the step decided (\code{status},
#'       \code{reason}, \code{criterion}, \code{selected}, \code{seconds}) and
#'       its own slice of \code{candidates}. When the run kept its
#'       \code{directory}, the tool's own record is nested beside that, read
#'       back from the step's directory with \code{\link{ferx_search_results}}:
#'       \code{models} or \code{steps} (whichever table that tool writes),
#'       \code{input_model_path} and \code{final_model_path}, and
#'       \code{model_paths} for every candidate model it kept. These are the
#'       engine's own tables rather than \code{ferx_modelsearch} /
#'       \code{ferx_iivsearch} objects: the pipeline returns each step's
#'       candidates adapted to one shape and leaves the tool's fuller record on
#'       disk, so a classed object would have to be rebuilt in R from those
#'       files - a second construction path, and still without the fit and the
#'       options only a direct call returns. Run the tool on its own when you
#'       want its object.}
#'     \item{fit}{The final model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when the pipeline ended with no fit in hand.}
#'     \item{input_model, input_ofv}{The model the pipeline started from, and
#'       its objective function value.}
#'     \item{final_model, final_model_path, final_ofv, d_ofv}{The model it
#'       ended on - already seeded from its own estimates - the file it was
#'       written to, its objective function value, and the change against the
#'       starting model (negative is an improvement).}
#'     \item{options}{The pipeline as the engine read it: \code{strategy},
#'       \code{retries_on} and \code{skip}.}
#'     \item{summary_text}{The engine's own report, as \code{ferx amd} prints
#'       it.}
#'     \item{notes, cancelled}{What the pipeline wants said once (a step that
#'       failed, a retries pass that could not run), and whether it was stopped
#'       early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin_amd")
#'
#' # What would run, before paying for it
#' ferx_amd_plan(model = ex$model, data = ex$data,
#'               search_space = "ABSORPTION(FO); PERIPHERALS(0..1); IIV?(@PK, exp)")
#'
#' # Reproducible: the .ferxsearch file is the artifact
#' res <- ferx_amd(config = ex$search, directory = "amd-run")
#' res
#' res$steps[, c("index", "step", "status", "d_value", "passed", "selected")]
#' summary(res)
#'
#' # Inline, with the residual step left out
#' res2 <- ferx_amd(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "ABSORPTION(FO); PERIPHERALS(0..1); IIV?(@PK, exp)",
#'   strategy     = "SIR",
#'   skip         = "residual",
#'   directory    = "amd-run-2"
#' )
#' }
#'
#' @seealso \code{\link{ferx_amd_plan}} for the pipeline without the fitting,
#'   \code{\link{ferx_search_config}}, \code{\link{ferx_modelsearch}},
#'   \code{\link{ferx_iivsearch}}, \code{\link{ferx_ruvsearch}},
#'   \code{\link{ferx_iovsearch}}, \code{\link{ferx_allometry}},
#'   \code{\link{ferx_covsearch}}
#' @family search
#' @export
ferx_amd <- function(model = NULL,
                     data = NULL,
                     search_space = NULL,
                     config = NULL,
                     strategy = NULL,
                     retries_on = NULL,
                     skip = NULL,
                     rank = NULL,
                     cutoff = NULL,
                     threads = NULL,
                     retries = NULL,
                     directory = NULL,
                     resume = FALSE,
                     progress = interactive()) {
  what <- "ferx_amd"
  args <- .ferx_amd_args(model, data, search_space, config, strategy,
                         retries_on, skip, rank, cutoff, what)

  # A run directory is not optional for the engine: AMD materialises every
  # step's input model as a file, because that is what the next tool reads.
  # `directory = NULL` therefore means a temporary one, removed on the way out
  # rather than never written.
  dir_arg <- .ferx_search_directory(directory, what)
  keep <- nzchar(dir_arg)
  if (!keep && .ferx_search_bool(resume, "resume", what)) {
    stop(what, ": `resume = TRUE` needs the `directory` of the run to resume; ",
         "without one there are no journalled fits to reuse")
  }
  if (!keep) {
    dir_arg <- tempfile(pattern = "ferx-amd-")
    dir.create(dir_arg, recursive = TRUE)
    on.exit(unlink(dir_arg, recursive = TRUE), add = TRUE)
  }

  raw <- ferx_rust_amd(
    config_path = args$config_path,
    model_path  = args$model,
    data_path   = args$data,
    mfl         = args$mfl,
    strategy    = args$strategy,
    retries_on  = args$retries_on,
    skip        = args$skip,
    rank        = args$rank,
    rank_cutoff = .ferx_search_scalar(cutoff, "cutoff", what),
    threads     = .ferx_search_count(threads, "threads", what, below = 0L),
    retries     = .ferx_search_count(retries, "retries", what, below = -1L, min = 0L),
    resume      = .ferx_search_bool(resume, "resume", what),
    directory   = dir_arg,
    progress    = .ferx_search_bool(progress, "progress", what)
  )

  # Every candidate of every step, as `CANDIDATE_COLUMNS` orders it.
  candidates <- data.frame(
    step        = as.integer(raw$c_step),
    tool        = as.character(raw$c_tool),
    id          = as.character(raw$c_id),
    parent      = .ferx_search_chr(raw$c_parent),
    description = as.character(raw$c_description),
    criterion   = as.character(raw$c_criterion),
    value       = .ferx_search_num(raw$c_value),
    d_value     = .ferx_search_num(raw$c_d_value),
    ofv         = .ferx_search_num(raw$c_ofv),
    d_ofv       = .ferx_search_num(raw$c_d_ofv),
    rank        = as.integer(.ferx_search_num(raw$c_rank)),
    converged   = .ferx_search_lgl(raw$c_converged),
    passed      = as.logical(raw$c_passed),
    failures    = .ferx_search_chr(raw$c_failures),
    error       = .ferx_search_chr(raw$c_error),
    note        = .ferx_search_chr(raw$c_note),
    seconds     = .ferx_search_num(raw$c_seconds),
    selected    = as.logical(raw$c_selected),
    stringsAsFactors = FALSE
  )

  # The step table, as `STEP_COLUMNS` orders it, with the pipeline position in
  # front. `d_value` and `d_ofv` are the same subtraction the engine's own
  # steps.csv writes, and the two gate columns are joined from the candidate
  # table rather than derived: `converged` and `passed` belong to the model the
  # step selected, which is a row of `candidates`.
  steps <- data.frame(
    index        = as.integer(raw$s_index),
    step         = as.character(raw$s_step),
    tool         = as.character(raw$s_tool),
    rerun        = as.logical(raw$s_rerun),
    directory    = as.character(raw$s_dir),
    status       = as.character(raw$s_status),
    reason       = .ferx_search_chr(raw$s_reason),
    criterion    = .ferx_search_chr(raw$s_criterion),
    value_before = .ferx_search_num(raw$s_value_before),
    value_after  = .ferx_search_num(raw$s_value_after),
    ofv_before   = .ferx_search_num(raw$s_ofv_before),
    ofv_after    = .ferx_search_num(raw$s_ofv_after),
    candidates   = as.integer(raw$s_candidates),
    selected     = .ferx_search_chr(raw$s_selected),
    seconds      = .ferx_search_num(raw$s_seconds),
    notes        = .ferx_search_chr(raw$s_notes),
    stringsAsFactors = FALSE
  )
  steps$d_value <- steps$value_after - steps$value_before
  steps$d_ofv <- steps$ofv_after - steps$ofv_before
  verdict <- .ferx_amd_step_verdict(steps$index, candidates)
  steps$converged <- verdict$converged
  steps$passed <- verdict$passed
  # The engine's `STEP_COLUMNS`, in the engine's order, with the pipeline
  # position in front and the two joined columns plus the step's notes behind.
  steps <- steps[, c("index", "step", "tool", "rerun", "directory", "status",
                     "reason", "criterion", "value_before", "value_after",
                     "d_value", "ofv_before", "ofv_after", "d_ofv",
                     "candidates", "selected", "seconds",
                     "converged", "passed", "notes")]

  # The winning model as a file: the run's own `final.ferx` when one survives
  # the call (the engine seeds the final estimates into it), otherwise a
  # temporary copy so the printed sections have a file to read.
  final_path <- file.path(dir_arg, "final.ferx")
  if (!keep || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-amd-final-", fileext = ".ferx")
    writeLines(as.character(raw$final_model), final_path)
  }

  fit <- if (is.null(raw$final_fit)) {
    NULL
  } else {
    .ferx_fit_from_raw(raw$final_fit, model = final_path,
                       data = as.character(raw$data))
  }

  result <- list(
    steps            = steps,
    candidates       = candidates,
    tools            = .ferx_amd_tools(steps, candidates,
                                       if (keep) dir_arg else NA_character_),
    fit              = fit,
    options          = list(
      strategy   = as.character(raw$strategy),
      retries_on = as.character(raw$retries_on),
      skip       = as.character(raw$skip_steps)
    ),
    input_model      = as.character(raw$input_model),
    input_ofv        = .ferx_search_num(raw$input_ofv),
    final_model      = as.character(raw$final_model),
    final_model_path = final_path,
    final_ofv        = .ferx_search_num(raw$final_ofv),
    d_ofv            = .ferx_search_num(raw$d_ofv),
    summary_text     = as.character(raw$summary),
    model            = as.character(raw$model),
    data             = as.character(raw$data),
    directory        = if (keep) dir_arg else NA_character_,
    steps_csv        = if (keep) as.character(raw$steps_csv) else NA_character_,
    candidates_csv   = if (keep) as.character(raw$candidates_csv) else NA_character_,
    config           = if (nzchar(args$config_path)) args$config_path else NA_character_,
    notes            = as.character(raw$notes),
    cancelled        = isTRUE(raw$cancelled)
  )
  # Partial matching would answer `res$fit` with `res$final_ofv` if the element
  # were dropped for being NULL, so it is kept explicitly.
  result["fit"] <- list(fit)
  class(result) <- c("ferx_amd", "ferx_search_result")
  result
}

#' The AMD pipeline as planned, without fitting anything
#'
#' Returns the steps \code{\link{ferx_amd}} would run for this space and this
#' model, in order, each either runnable or skipped with the reason - the plan
#' the engine computes before its first fit, so it is the same list the run
#' itself walks rather than a second derivation of it.
#'
#' The model and the dataset are read (a step is skipped when the starting
#' model declares no \code{iov_column}, which only the model says), but nothing
#' is fitted.
#'
#' @inheritParams ferx_amd
#'
#' @return A data frame with one row per planned step: \code{index},
#'   \code{step}, \code{tool}, \code{rerun}, \code{directory} (the name the
#'   step's files would go under) and \code{skipped} (the reason, or \code{NA}
#'   when the step will run). The options as the engine read them are attached
#'   as the \code{options} attribute.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin_amd")
#' ferx_amd_plan(config = ex$search)
#' ferx_amd_plan(model = ex$model, data = ex$data,
#'               search_space = "IIV?(@PK, exp)", strategy = "SIR")
#' }
#'
#' @seealso \code{\link{ferx_amd}}
#' @family search
#' @export
ferx_amd_plan <- function(model = NULL,
                          data = NULL,
                          search_space = NULL,
                          config = NULL,
                          strategy = NULL,
                          retries_on = NULL,
                          skip = NULL) {
  what <- "ferx_amd_plan"
  args <- .ferx_amd_args(model, data, search_space, config, strategy,
                         retries_on, skip, NULL, NULL, what)
  raw <- ferx_rust_amd_plan(
    config_path = args$config_path,
    model_path  = args$model,
    data_path   = args$data,
    mfl         = args$mfl,
    strategy    = args$strategy,
    retries_on  = args$retries_on,
    skip        = args$skip
  )
  out <- data.frame(
    index     = as.integer(raw$index),
    step      = as.character(raw$step),
    tool      = as.character(raw$tool),
    rerun     = as.logical(raw$rerun),
    directory = as.character(raw$directory),
    skipped   = .ferx_search_chr(raw$skipped),
    stringsAsFactors = FALSE
  )
  attr(out, "options") <- list(
    strategy   = as.character(raw$strategy),
    retries_on = as.character(raw$retries_on),
    skip       = as.character(raw$skip_steps),
    iov_column = .ferx_search_chr(raw$iov_column),
    model      = as.character(raw$model),
    data       = as.character(raw$data)
  )
  out
}

# -- Internals ---------------------------------------------------------------

# The step names `[amd] skip` knows, as the engine spells them.
.ferx_amd_steps <- c("structural", "iivsearch", "residual", "iovsearch",
                     "allometry", "covariates")

# Everything both entry points validate the same way: the entry form, the
# space, and the three `[amd]` keys. Kept in one place so the plan cannot
# accept a pipeline the run would refuse, or the other way round.
.ferx_amd_args <- function(model, data, search_space, config, strategy,
                           retries_on, skip, rank, cutoff, what) {
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      search_space = search_space,
      strategy     = strategy,
      retries_on   = retries_on,
      skip         = skip,
      rank         = rank,
      cutoff       = cutoff
    ),
    what
  )

  if (nzchar(config_path)) {
    paths <- list(model = "", data = "")
    mfl <- ""
  } else {
    paths <- .ferx_search_model_data(model, data, what)
    # A pipeline with no space is a residual-error search with five skipped
    # steps, which is `ferx_ruvsearch()` under another name.
    mfl <- .ferx_search_space_text(search_space, what)
  }

  if (!is.null(strategy)) {
    # `match.arg()` partial-matches, so "S" would resolve to "SIR" - harmless
    # here, since every choice is spelled out in the error when it does not.
    strategy <- match.arg(strategy,
                          c("default", "reevaluation", "SIR", "SRI", "RSI"))
  }
  if (!is.null(retries_on)) {
    retries_on <- match.arg(retries_on, c("all_final", "final", "skip"))
  }
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  list(
    config_path = config_path,
    model       = paths$model,
    data        = paths$data,
    mfl         = mfl,
    strategy    = strategy %||% "",
    retries_on  = retries_on %||% "",
    skip        = .ferx_amd_skip(skip, what),
    rank        = rank %||% ""
  )
}

# The `skip` steps, as the engine spells them. A closed set, so a misspelling
# is caught here rather than becoming a TOML key the loader refuses several
# layers down.
.ferx_amd_skip <- function(skip, what) {
  if (is.null(skip)) return(character(0))
  if (!is.character(skip) || anyNA(skip)) {
    stop(sprintf("%s: `skip` must be a character vector of step names or NULL", what))
  }
  hit <- match(tolower(skip), .ferx_amd_steps)
  if (anyNA(hit)) {
    stop(sprintf("%s: `skip` does not name a pipeline step: %s. Use any of %s",
                 what,
                 paste(skip[is.na(hit)], collapse = ", "),
                 paste(.ferx_amd_steps, collapse = ", ")))
  }
  unique(.ferx_amd_steps[hit])
}

# The termination status and the gate verdict of the model each step selected.
#
# The step table has neither of its own: a step reports what it selected, and
# whether that model converged and passed is a fact about the *fit*, which is a
# row of the candidate table. The last selected row of the step is the one -
# a retries pass that improved on the step's winner is filed under the same
# step index and marked selected in its turn, and it is the model the step
# actually handed on.
#
# A step that selected nothing - a residual search whose candidates all failed
# their likelihood-ratio test, say - has no row of its own to read, and the
# model it carried forward belongs to the step that produced it. That is NA
# here rather than the previous step's verdict copied onto a row it is not
# about.
.ferx_amd_step_verdict <- function(index, candidates) {
  converged <- rep(NA, length(index))
  passed <- rep(NA, length(index))
  for (i in seq_along(index)) {
    rows <- which(candidates$step == index[i] &
                    !is.na(candidates$selected) & candidates$selected)
    if (length(rows) == 0L) next
    row <- rows[length(rows)]
    converged[i] <- candidates$converged[row]
    passed[i] <- candidates$passed[row]
  }
  list(converged = converged, passed = passed)
}

# The per-step view: each step that ran, with its slice of the pipeline's
# candidate table and - when the run kept its directory - the tool's *own*
# tables, read back from that step's directory with `ferx_search_results()`,
# the same reader a user would call on it.
#
# These are tables rather than `ferx_modelsearch` / `ferx_iivsearch` objects
# (#356 review). The engine's `AmdResult` carries no per-tool result struct:
# each step's rows are adapted into the pipeline's own `CandidateRow` shape and
# the tool's fuller record is left in its directory. Building classed objects
# here would mean a second construction path for each one, fed from CSVs rather
# than from the binding that builds them everywhere else - the drift risk the
# shared search contract exists to avoid - and they would still be missing what
# only the binding returns (`$fit`, the options as the engine read them, every
# candidate's model text). So what is nested is what the engine actually wrote:
# the tool's `models.csv` / `steps.csv`, its `input.ferx` and `final.ferx`, and
# the models it kept.
.ferx_amd_tools <- function(steps, candidates, directory = NA_character_) {
  ran <- steps[steps$status != "skipped", , drop = FALSE]
  if (nrow(ran) == 0L) return(list())
  out <- lapply(seq_len(nrow(ran)), function(i) {
    rows <- candidates[candidates$step == ran$index[i], , drop = FALSE]
    rownames(rows) <- NULL
    entry <- list(
      index      = ran$index[i],
      step       = ran$step[i],
      tool       = ran$tool[i],
      directory  = ran$directory[i],
      status     = ran$status[i],
      reason     = ran$reason[i],
      criterion  = ran$criterion[i],
      selected   = ran$selected[i],
      seconds    = ran$seconds[i],
      candidates = rows
    )
    c(entry, .ferx_amd_tool_files(directory, ran$directory[i]))
  })
  names(out) <- sprintf("%02d-%s", ran$index, ran$step)
  out
}

# What one step's directory holds, when the run kept one: the tool's own table
# (`models.csv` for the three tools that rank models, `steps.csv` for the two
# that test them), the model it was handed, the model it selected, and the text
# of every model it kept. Absent keys for a run held in memory, which wrote
# none of them.
.ferx_amd_tool_files <- function(directory, step_dir) {
  if (is.na(directory) || !nzchar(directory)) return(list())
  dir <- file.path(directory, step_dir)
  if (!dir.exists(dir)) return(list())

  out <- list(path = dir)
  for (type in c("models", "steps")) {
    file <- file.path(dir, paste0(type, ".csv"))
    if (!file.exists(file)) next
    tab <- tryCatch(ferx_search_results(dir, type = type),
                    error = function(e) NULL)
    if (!is.null(tab)) out[[type]] <- tab
  }
  for (nm in c("input", "final")) {
    file <- file.path(dir, paste0(nm, ".ferx"))
    if (file.exists(file)) out[[paste0(nm, "_model_path")]] <- file
  }
  models <- list.files(file.path(dir, "models"), pattern = "\\.ferx$",
                       full.names = TRUE)
  if (length(models) > 0L) {
    out$model_paths <- stats::setNames(models, sub("\\.ferx$", "",
                                                   basename(models)))
  }
  out
}

# A fixed-decimal OFV. A pipeline step is decided on differences of a few OFV
# units, which four significant digits would round away on a four-figure
# objective function value.
.ferx_amd_num <- function(v, digits = 3) {
  ifelse(is.na(v), "-", formatC(v, format = "f", digits = digits))
}

#' @param x A \code{ferx_amd} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_amd
#' @export
print.ferx_amd <- function(x, digits = 4, ...) {
  cat("ferx AMD pipeline (", x$options$strategy, " strategy, retries ",
      x$options$retries_on, ")\n", sep = "")
  cat("  Model: ", x$model, "\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  cat(sprintf("  OFV:   %s (start) -> %s (final, dOFV %s)\n",
              .ferx_amd_num(x$input_ofv), .ferx_amd_num(x$final_ofv),
              .ferx_amd_num(x$d_ofv)))
  ran <- sum(x$steps$status != "skipped")
  # The steps' wall clock, not the candidates' summed seconds: candidates are
  # fitted in parallel, so their sum is thread-seconds and would report a run
  # as several times longer than it took.
  cat(sprintf("  %d of %d step%s ran, %d model%s fitted, %s s across the steps\n",
              ran, nrow(x$steps), if (nrow(x$steps) == 1L) "" else "s",
              nrow(x$candidates), if (nrow(x$candidates) == 1L) "" else "s",
              formatC(sum(x$steps$seconds, na.rm = TRUE),
                      format = "f", digits = 1)))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the pipeline finished; the steps below are what it reached.\n")
  }

  cat("\nPipeline:\n")
  .ferx_search_print_table(
    x$steps[, c("index", "step", "status", "criterion", "value_before",
                "value_after", "d_value", "d_ofv", "candidates", "converged",
                "passed", "seconds"), drop = FALSE],
    digits
  )

  # A step that did not run, or ran and failed, is the part of the table a
  # reader most needs in words: `status` alone does not say why.
  aside <- x$steps[x$steps$status != "ran", , drop = FALSE]
  if (nrow(aside) > 0L) {
    cat("\nSteps that did not run:\n")
    for (i in seq_len(nrow(aside))) {
      cat(sprintf("  %d %-12s %s: %s\n", aside$index[i], aside$step[i],
                  aside$status[i], aside$reason[i]))
    }
  }

  cat("\nSelected:\n")
  chose <- x$steps[x$steps$status == "ran", , drop = FALSE]
  if (nrow(chose) == 0L) {
    cat("  (nothing - no step ran)\n")
  } else {
    for (i in seq_len(nrow(chose))) {
      cat(sprintf("  %-12s %s\n", chose$step[i],
                  if (is.na(chose$selected[i])) "(nothing)" else chose$selected[i]))
    }
  }

  gated <- sum(!x$candidates$passed, na.rm = TRUE)
  if (gated > 0L) {
    cat(sprintf("\n%d candidate%s excluded by the strictness gate; summary() lists them.\n",
                gated, if (gated == 1L) "" else "s"))
  }
  failed <- x$steps$step[x$steps$status == "failed"]
  if (length(failed) > 0L) {
    cat("Steps that failed: ", paste(failed, collapse = ", "),
        " - the pipeline carried on from the model each was handed.\n", sep = "")
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

#' @param object A \code{ferx_amd} object.
#' @rdname ferx_amd
#' @export
summary.ferx_amd <- function(object, digits = 4, ...) {
  print(object, digits = digits)

  cols <- c("id", "description", "criterion", "value", "d_value", "d_ofv",
            "rank", "converged", "passed", "selected", "seconds")
  for (tool in object$tools) {
    cat(sprintf("\nStep %d - %s (%s), ranked on %s:\n", tool$index, tool$step,
                tool$tool, if (is.na(tool$criterion)) "-" else tool$criterion))
    .ferx_search_print_table(tool$candidates[, cols, drop = FALSE], digits)
  }

  rejected <- object$candidates[!is.na(object$candidates$passed) &
                                  !object$candidates$passed, , drop = FALSE]
  cat("\nCandidates the strictness gate excluded, and why:\n")
  if (nrow(rejected) == 0L) {
    cat("  (none - every candidate fitted passed the gate)\n")
  } else {
    .ferx_search_print_table(
      rejected[, c("step", "tool", "id", "description", "converged",
                   "failures", "error"), drop = FALSE],
      digits
    )
  }
  invisible(object)
}
