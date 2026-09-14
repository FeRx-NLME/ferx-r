#' Global model search
#'
#' Searches a model space \emph{globally} rather than one axis at a time.
#' Every structural category the space names becomes an axis of a grid with
#' that category's values as its alleles, every optional \code{COVARIATE?} pair
#' becomes an axis whose alleles are \code{none} and each of its forms, and the
#' search walks the grid: \code{"exhaustive"} fits every point, \code{"ga"} -
#' the default - runs pyDarwin's genetic algorithm over it and fits a fraction
#' of them.
#'
#' Where \code{\link{ferx_modelsearch}} decides the structure and
#' \code{\link{ferx_covsearch}} the covariate relations, each conditional on
#' the other being fixed, a global search decides them together. That is the
#' case for using it - a covariate that only pays off once a second
#' compartment is in the model is invisible to a stepwise search - and the
#' cost is the fit count: an exhaustive grid is the product of its axes, and
#' is refused outright above \code{max_models} rather than silently truncated.
#'
#' @section A different name from ferx_fit()'s global_search:
#' \code{ferx_fit(settings = list(global_search = TRUE))} is a global
#' \emph{optimizer} phase inside the estimation of one model - a wider search
#' for the parameter values of a fixed model. \code{ferx_globalsearch()} is a
#' global search over \emph{models}. The two are unrelated; the name here is
#' the engine's, the \code{.ferxsearch} file's (\code{[globalsearch]}), the
#' \code{ferx globalsearch} command's and Pharmpy/pyDarwin's, so a search is
#' portable between them.
#'
#' @section Ranking, and the charges the criterion cannot see:
#' \code{[rank] type} defaults to \code{"penalized"} here - pyDarwin's
#' fitness, the objective function plus a charge per estimated parameter and
#' per unhealthy fit - where every other tool defaults to a BIC. Any criterion
#' may be named instead.
#'
#' Whatever the criterion, the search charges three things on top of it, and
#' ranks on the sum (\code{fitness}, always finite):
#' \describe{
#'   \item{\code{charge_non_influential}}{Per gene of a genome that changed
#'     nothing in the rendered model - a covariate on a parameter the
#'     structural choice removed. A tie-break towards the simpler genotype
#'     among models that are identical in every other way.}
#'   \item{\code{charge_gate}}{A fit the \code{[strictness]} gate refused. It
#'     cannot be selected, but it still steers a genetic algorithm.}
#'   \item{\code{charge_crash}}{A candidate that produced no fit at all -
#'     it does not compile, or the fit errored. Large, and finite, so the
#'     algorithm can still rank it.}
#' }
#' The three are columns of \code{$models} beside \code{criterion} and
#' \code{fitness}, so a genome that lost to a tie-break penalty does not read
#' as though it lost on OFV. \code{$penalties} is the effective schedule the
#' run charged.
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
#' \code{ABSORPTION}, \code{ELIMINATION}, \code{PERIPHERALS}, \code{TRANSITS},
#' \code{LAGTIME} and \code{COVARIATE?}. A variability statement
#' (\code{IIV?} / \code{IOV?}) or \code{ALLOMETRY} is another tool's axis and
#' is refused by name before any fit starts, as is a feature the engine cannot
#' build. Check a space with \code{\link{ferx_search_space}} /
#' \code{\link{ferx_search_coverage}} first: both answer without fitting
#' anything.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the search runs on, so like \code{model} it
#'   cannot be given beside \code{config}.
#' @param search_space MFL text naming the grid, e.g.
#'   \code{"PERIPHERALS(0..1); LAGTIME([OFF,ON]); COVARIATE?(CL, WT, pow)"}.
#' @param config Path to a \code{.ferxsearch} file. Mutually exclusive with the
#'   arguments that state the search.
#' @param algorithm \code{"ga"} (the default) or \code{"exhaustive"}.
#'   \code{NULL} keeps the engine default. An unrecognised value is an error
#'   naming it - there is no partial matching here, because
#'   \code{algorithm = "e"} choosing an exhaustive grid over a four-axis space
#'   is not a typo worth guessing at.
#' @param iiv_strategy How a candidate's new structural parameters are given a
#'   random effect: \code{"absorption_delay"} (the default),
#'   \code{"add_diagonal"} or \code{"no_add"}. \code{NULL} keeps the engine
#'   default.
#' @param max_models The largest grid \code{"exhaustive"} will enumerate (500
#'   unless the file says otherwise). A bigger grid is an error naming the
#'   size, not a truncated search; the genetic algorithm has no such cap.
#' @param ga A named list of \code{[globalsearch.ga]} knobs -
#'   \code{population_size}, \code{generations}, \code{crossover_rate},
#'   \code{mutation_rate}, \code{gene_mutation_probability}, \code{elites},
#'   \code{tournament_size}, \code{downhill_period}, \code{niches},
#'   \code{niche_radius}, \code{niche_penalty}, \code{sharing_alpha},
#'   \code{final_downhill}, \code{seed}. The key list and its defaults come
#'   from the engine, so an unknown knob is refused by name. Only meaningful
#'   for \code{algorithm = "ga"}.
#' @param penalties A named list of \code{[rank.penalties]} charges to
#'   override - \code{theta}, \code{omega}, \code{sigma}, \code{convergence},
#'   \code{covariance}, \code{correlation}, \code{max_correlation},
#'   \code{condition_number}, \code{max_condition_number},
#'   \code{non_influential}, \code{crash}, \code{gate}. The first nine are the
#'   \code{"penalized"} criterion's own schedule; the last three are the
#'   search-level charges, which apply under any criterion.
#' @param rank \code{[rank] type} - the criterion the fitness is built on,
#'   e.g. \code{"penalized"}, \code{"bic"}, \code{"aic"}, \code{"ofv"}.
#'   \code{NULL} keeps the tool default (\code{"penalized"}).
#' @param cutoff \code{[rank] cutoff}: the margin, on the fitness scale, by
#'   which a candidate must beat the \emph{input} model to be selected - not
#'   the base model, as in the stepwise tools, because a global search has no
#'   base beyond where it started. \code{NULL} selects the best eligible
#'   candidate outright.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per candidate on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the per-batch journals, \code{models.csv},
#'   \code{generations.csv}, \code{models/<id>.ferx} and \code{final.ferx} are
#'   written. \code{NULL} keeps the run in memory, which also means it cannot
#'   be resumed.
#' @param resume Reuse the fits already journalled in \code{directory}. The
#'   genetic algorithm proposes the same genomes for the same \code{seed} on
#'   the same grid, which is what makes a resume reuse its journal. Refused
#'   before the first fit when there is nothing to resume from - no
#'   \code{directory}, or a \code{directory} no earlier run wrote - rather
#'   than quietly searching from scratch.
#' @param progress Print the engine's batch progress to the console.
#'
#' @return An object of class \code{ferx_globalsearch}:
#'   \describe{
#'     \item{models}{The model table, one row per grid point evaluated (and
#'       the input), in the engine's own column order: \code{id},
#'       \code{parent}, \code{step}, \code{genome}, \code{absorption},
#'       \code{elimination}, \code{peripherals}, \code{transits},
#'       \code{lagtime}, \code{covariates}, \code{n_parameters}, \code{ofv},
#'       \code{criterion}, \code{fitness}, \code{rank}, \code{converged},
#'       \code{passed}, \code{failures}, \code{error}, \code{seconds},
#'       \code{selected}, \code{non_influential}, \code{duplicate_of},
#'       \code{reused} - followed by \code{charge_non_influential},
#'       \code{charge_gate} and \code{charge_crash}, the decomposition of what
#'       the search added to the criterion to reach the fitness.}
#'     \item{fit}{The winning model's fit as a \code{ferx_fit}, or \code{NULL}
#'       when a degraded resume could not recover it.}
#'     \item{model_text}{Every evaluated model's text, named by model id.}
#'     \item{axes, space_size}{The grid: a named list of each axis's allele
#'       labels, and the number of points they multiply out to. A double
#'       rather than an integer, because the product of enough axes exceeds
#'       R's integer range while remaining a grid the genetic algorithm can
#'       search.}
#'     \item{generations}{The genetic algorithm's trajectory - \code{index},
#'       \code{best} (the model id), \code{best_fitness}, \code{mean_fitness},
#'       \code{polished} - or \code{NULL} for an exhaustive search.}
#'     \item{input_model, input_ofv, input_criterion, input_fitness}{The model
#'       the search started from, and what it scored.}
#'     \item{final_model, final_model_path, final_model_id, final_criterion,
#'       final_fitness}{The winning model's text, the \code{final.ferx} written
#'       beside the run, its id, and the two values it won on.}
#'     \item{criterion, algorithm, iiv_strategy, max_models, cutoff}{What the
#'       search ranked on and how it searched.}
#'     \item{penalties, ga}{The effective \code{[rank.penalties]} schedule and
#'       \code{[globalsearch.ga]} settings, as named numeric vectors.}
#'     \item{n_fitted}{Models actually fitted, as against grid points
#'       evaluated: a duplicate genome and a cached one cost no fit.}
#'     \item{candidates}{The runner's candidate table for every batch the run
#'       wrote one - see \code{\link{ferx_search_results}}.}
#'     \item{notes, cancelled}{What the search wants said once, and whether it
#'       was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("two_cpt_oral_global")
#'
#' # Inline: a small grid, enumerated exhaustively
#' res <- ferx_globalsearch(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "PERIPHERALS(0..1); COVARIATE?(CL, WT, pow)",
#'   algorithm    = "exhaustive",
#'   directory    = "globalsearch-run-1"
#' )
#' res
#' res$models[, c("id", "genome", "criterion", "charge_gate", "fitness", "rank")]
#'
#' # The genetic algorithm over the same grid, seeded so a resume can reuse
#' # the journal it writes
#' res2 <- ferx_globalsearch(
#'   model        = ex$model,
#'   data         = ex$data,
#'   search_space = "PERIPHERALS(0..1); COVARIATE?(CL, WT, pow)",
#'   algorithm    = "ga",
#'   ga           = list(population_size = 8, generations = 3, seed = 1),
#'   directory    = "globalsearch-run-2"
#' )
#' res2$generations
#'
#' # Reproducible: the bundled .ferxsearch is the artifact, and states the
#' # whole search - a four-axis grid, enumerated, ranked on penalized fitness
#' res3 <- ferx_globalsearch(config = ex$search, directory = "globalsearch-run-3")
#' summary(res3)
#' }
#'
#' @seealso \code{\link{ferx_modelsearch}} and \code{\link{ferx_covsearch}}
#'   for the stepwise searches of the same two axes,
#'   \code{\link{ferx_search_space}}, \code{\link{ferx_search_coverage}},
#'   \code{\link{ferx_search_results}}
#' @family search
#' @export
ferx_globalsearch <- function(model = NULL,
                              data = NULL,
                              search_space = NULL,
                              config = NULL,
                              algorithm = NULL,
                              iiv_strategy = NULL,
                              max_models = NULL,
                              ga = NULL,
                              penalties = NULL,
                              rank = NULL,
                              cutoff = NULL,
                              threads = NULL,
                              retries = NULL,
                              directory = NULL,
                              resume = FALSE,
                              progress = interactive()) {
  what <- "ferx_globalsearch"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      search_space = search_space,
      algorithm    = algorithm,
      iiv_strategy = iiv_strategy,
      max_models   = max_models,
      ga           = ga,
      penalties    = penalties,
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
    mfl <- .ferx_search_space_text(search_space, what)
  }

  # Exact, not partial: `algorithm = "e"` choosing an exhaustive enumeration
  # over a grid whose size the user has not looked at is not a typo worth
  # guessing at, and the engine would accept neither spelling (#364).
  algorithm <- .ferx_search_choice(algorithm, c("ga", "exhaustive"),
                                   "algorithm", what)
  # `fullblock` is Pharmpy's fourth strategy and the engine refuses it by name
  # (a block over the new and existing eta is a variability search's move), so
  # it stays in the choices rather than being silently unavailable.
  iiv_strategy <- .ferx_search_choice(
    iiv_strategy,
    c("absorption_delay", "add_diagonal", "no_add", "fullblock"),
    "iiv_strategy", what
  )
  if (!is.null(rank) && (!is.character(rank) || length(rank) != 1L || is.na(rank))) {
    stop(what, ": `rank` must be a single string or NULL")
  }

  keys <- .ferx_globalsearch_option_keys()
  ga_toml <- .ferx_search_toml_table(ga, "ga", keys$ga_name, keys$ga_kind, what)
  pen_toml <- .ferx_search_toml_table(
    penalties, "penalties", keys$penalty_name,
    rep("number", length(keys$penalty_name)), what
  )

  dir_arg <- .ferx_search_directory(directory, what, resume)
  raw <- ferx_rust_globalsearch(
    config_path    = config_path,
    model_path     = paths$model,
    data_path      = paths$data,
    mfl            = mfl,
    algorithm      = algorithm %||% "",
    iiv_strategy   = iiv_strategy %||% "",
    max_models     = .ferx_search_count(max_models, "max_models", what,
                                        below = 0L, min = 1L),
    ga_keys        = ga_toml$keys,
    ga_values      = ga_toml$values,
    penalty_keys   = pen_toml$keys,
    penalty_values = pen_toml$values,
    rank           = rank %||% "",
    rank_cutoff    = .ferx_search_scalar(cutoff, "cutoff", what),
    threads        = .ferx_search_count(threads, "threads", what, below = 0L),
    retries        = .ferx_search_count(retries, "retries", what,
                                        below = -1L, min = 0L),
    resume         = .ferx_search_bool(resume, "resume", what),
    directory      = dir_arg,
    progress       = .ferx_search_bool(progress, "progress", what)
  )

  # The engine's 24 columns in the engine's order, then the three charges -
  # computed by the glue from the same penalty schedule the search ranked on,
  # so `criterion + charges == fitness` holds row for row.
  models <- data.frame(
    id              = as.character(raw$id),
    parent          = .ferx_search_chr(raw$parent),
    step            = as.character(raw$step),
    genome          = .ferx_search_chr(raw$genome),
    absorption      = .ferx_search_chr(raw$absorption),
    elimination     = .ferx_search_chr(raw$elimination),
    peripherals     = as.integer(.ferx_search_num(raw$peripherals)),
    transits        = .ferx_search_chr(raw$transits),
    lagtime         = .ferx_search_chr(raw$lagtime),
    covariates      = .ferx_search_chr(raw$covariates),
    n_parameters    = as.integer(.ferx_search_num(raw$n_parameters)),
    ofv             = .ferx_search_num(raw$ofv),
    criterion       = .ferx_search_num(raw$criterion),
    fitness         = .ferx_search_num(raw$fitness),
    rank            = as.integer(.ferx_search_num(raw$rank)),
    converged       = .ferx_search_lgl(raw$converged),
    passed          = as.logical(raw$passed),
    failures        = .ferx_search_chr(raw$failures),
    error           = .ferx_search_chr(raw$error),
    seconds         = .ferx_search_num(raw$seconds),
    selected        = as.logical(raw$selected),
    non_influential = as.integer(raw$non_influential),
    duplicate_of    = .ferx_search_chr(raw$duplicate_of),
    reused          = as.logical(raw$reused),
    charge_non_influential = .ferx_search_num(raw$charge_non_influential),
    charge_gate            = .ferx_search_num(raw$charge_gate),
    charge_crash           = .ferx_search_num(raw$charge_crash),
    stringsAsFactors = FALSE
  )

  model_text <- stats::setNames(as.character(raw$model_text),
                                as.character(raw$model_id))

  generations <- if (length(raw$g_index) == 0L) {
    NULL
  } else {
    data.frame(
      index        = as.integer(raw$g_index),
      best         = as.character(raw$g_best),
      best_fitness = .ferx_search_num(raw$g_best_fitness),
      mean_fitness = .ferx_search_num(raw$g_mean_fitness),
      polished     = as.integer(raw$g_polished),
      stringsAsFactors = FALSE
    )
  }

  # The winning model as a file: the run's own `final.ferx` when it wrote one
  # (the engine seeds the final estimates into it), otherwise a temporary copy
  # so the fit below has a model file to name and read.
  final_path <- if (nzchar(dir_arg)) file.path(dir_arg, "final.ferx") else NA_character_
  if (is.na(final_path) || !file.exists(final_path)) {
    final_path <- tempfile(pattern = "ferx-globalsearch-final-", fileext = ".ferx")
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
    axes             = raw$axes,
    # A double, not an integer: the grid is the product of its axes, and 31
    # binary axes already exceed R's 32-bit integer range. `as.integer()`
    # would answer that with NA and a coercion warning, for a search the
    # genetic algorithm runs perfectly well (it only ever evaluates its
    # population). The engine counts it as a `usize` and the glue hands it
    # over as a double for the same reason.
    space_size       = as.numeric(raw$space_size),
    generations      = generations,
    input_model      = as.character(raw$input_model),
    input_ofv        = .ferx_search_num(raw$input_ofv),
    input_criterion  = .ferx_search_num(raw$input_criterion),
    input_fitness    = .ferx_search_num(raw$input_fitness),
    final_model      = as.character(raw$final_model),
    final_model_path = final_path,
    final_model_id   = as.character(raw$final_id),
    final_criterion  = if (nrow(final_row) == 1L) final_row$criterion else NA_real_,
    final_fitness    = .ferx_search_num(raw$final_fitness),
    criterion        = as.character(raw$criterion_label),
    algorithm        = as.character(raw$algorithm),
    iiv_strategy     = as.character(raw$iiv_strategy),
    max_models       = as.integer(raw$max_models),
    cutoff           = .ferx_search_num(raw$rank_cutoff),
    penalties        = stats::setNames(as.numeric(raw$penalty_value),
                                       as.character(raw$penalty_name)),
    ga               = stats::setNames(as.numeric(raw$ga_value),
                                       as.character(raw$ga_name)),
    n_fitted         = as.integer(raw$n_fitted),
    summary_text     = as.character(raw$summary),
    # The batch directories this run wrote, named by the run's own rows rather
    # than scanned off disk: a shorter run in a directory an earlier one used
    # would otherwise inherit its candidate tables (#336 review).
    candidates       = .ferx_search_candidates(dir_arg, unique(models$step)),
    model            = as.character(raw$model),
    data             = as.character(raw$data),
    directory        = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config           = if (nzchar(config_path)) config_path else NA_character_,
    notes            = as.character(raw$notes),
    cancelled        = isTRUE(raw$cancelled)
  )
  class(result) <- c("ferx_globalsearch", "ferx_search_result")
  result
}

# The engine's own key lists for `[globalsearch.ga]` and `[rank.penalties]`,
# fetched once per session. R needs them to refuse an unknown knob by name
# before a config file is rendered; keeping them on the engine's side means a
# knob ferx-core adds does not have to be restated here.
.ferx_globalsearch_option_keys <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      raw <- ferx_rust_globalsearch_option_keys()
      cached <<- list(
        ga_name      = as.character(raw$ga_name),
        ga_kind      = as.character(raw$ga_kind),
        penalty_name = as.character(raw$penalty_name)
      )
    }
    cached
  }
})

# A choice argument, matched exactly. `match.arg()` would accept a prefix, and
# for `algorithm` that is the difference between a genetic algorithm and an
# exhaustive enumeration of a grid nobody sized (#364).
.ferx_search_choice <- function(x, choices, name, what) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s: `%s` must be a single string or NULL", what, name))
  }
  if (!x %in% choices) {
    stop(sprintf("%s: `%s = \"%s\"` is not one of %s", what, name, x,
                 paste0("\"", choices, "\"", collapse = ", ")))
  }
  x
}

# A named list of engine keys, as the parallel key / value character vectors
# the binding renders into a TOML sub-table.
#
# `known` and `kinds` come from the engine, so a knob this package has never
# heard of is still checked against the list the loader would check it
# against - and refused here, by name, rather than surfacing three layers down
# as a serde `unknown field`.
.ferx_search_toml_table <- function(x, arg, known, kinds, what) {
  empty <- list(keys = character(0), values = character(0))
  if (is.null(x)) return(empty)
  if (!is.list(x) && !is.numeric(x) && !is.logical(x)) {
    stop(sprintf("%s: `%s` must be a named list or NULL", what, arg))
  }
  x <- as.list(x)
  if (length(x) == 0L) return(empty)
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop(sprintf("%s: every element of `%s` must be named", what, arg))
  }
  if (anyDuplicated(nms)) {
    stop(sprintf("%s: `%s` names a setting twice: %s", what, arg,
                 paste(unique(nms[duplicated(nms)]), collapse = ", ")))
  }
  unknown <- setdiff(nms, known)
  if (length(unknown)) {
    stop(sprintf("%s: `%s` has no setting called %s. The settings are: %s",
                 what, arg,
                 paste(paste0("`", unknown, "`"), collapse = ", "),
                 paste(known, collapse = ", ")))
  }

  values <- character(length(x))
  for (i in seq_along(x)) {
    v <- x[[i]]
    kind <- kinds[match(nms[i], known)]
    label <- sprintf("%s$%s", arg, nms[i])
    if (identical(kind, "flag")) {
      values[i] <- if (.ferx_search_bool(v, label, what)) "true" else "false"
      next
    }
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v)) {
      stop(sprintf("%s: `%s` must be a single finite number", what, label))
    }
    if (identical(kind, "count")) {
      if (v < 0 || v != round(v)) {
        stop(sprintf("%s: `%s` must be a whole number, zero or more", what, label))
      }
      values[i] <- sprintf("%.0f", v)
      next
    }
    # 15 significant digits and no exponent: the engine's crash charge is
    # 99999999, which `format()`'s 7-digit default would render as `1e+08` -
    # a number, but not the one the user asked for.
    values[i] <- format(v, digits = 15, scientific = FALSE)
  }
  list(keys = nms, values = values)
}

#' @param x A \code{ferx_globalsearch} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_globalsearch
#' @export
print.ferx_globalsearch <- function(x, digits = 4, ...) {
  cat("ferx global model search (", x$algorithm, ", ranked on ", x$criterion,
      ")\n", sep = "")
  cat("  Data:  ", x$data, "\n", sep = "")
  if (!is.na(x$directory)) cat("  Wrote: ", x$directory, "\n", sep = "")
  # `%s` and an explicit format: `space_size` is a double so it can exceed R's
  # integer range, and `%d` refuses a double outright.
  cat(sprintf("  Grid:  %s point%s over %d ax%s\n",
              format(x$space_size, big.mark = ",", scientific = FALSE,
                     trim = TRUE),
              if (isTRUE(x$space_size == 1)) "" else "s",
              length(x$axes), if (length(x$axes) == 1L) "is" else "es"))
  for (nm in names(x$axes)) {
    cat("    ", nm, ": ", paste(x$axes[[nm]], collapse = " | "), "\n", sep = "")
  }
  cat(sprintf("  Input: OFV %s, %s %s\n",
              format(signif(x$input_ofv, digits)),
              x$criterion,
              format(signif(x$input_criterion, digits))))
  cat(sprintf("  Fitness: %s (input) -> %s (%s)\n",
              format(signif(x$input_fitness, digits)),
              format(signif(x$final_fitness, digits)),
              x$final_model_id))
  cat(sprintf("  %d model%s fitted of %d evaluated\n",
              x$n_fitted, if (x$n_fitted == 1L) "" else "s", nrow(x$models)))
  if (isTRUE(x$cancelled)) {
    cat("  Cancelled before the search finished; the models below are what it reached.\n")
  }

  cat("\nRanked models (best first):\n")
  ranked <- x$models[!is.na(x$models$rank), , drop = FALSE]
  ranked <- ranked[order(ranked$rank), , drop = FALSE]
  .ferx_search_print_table(
    ranked[, c("rank", "id", "genome", "ofv", "criterion", "fitness",
               "converged", "passed"), drop = FALSE],
    digits
  )

  charged <- x$models[x$models$charge_gate > 0 | x$models$charge_crash > 0 |
                        x$models$charge_non_influential > 0, , drop = FALSE]
  charged <- charged[!is.na(charged$id), , drop = FALSE]
  if (nrow(charged) > 0L) {
    cat("\nModels the search charged beyond the criterion:\n")
    .ferx_search_print_table(
      charged[, c("id", "genome", "criterion", "charge_non_influential",
                  "charge_gate", "charge_crash", "fitness"), drop = FALSE],
      digits
    )
  }

  if (!is.null(x$generations)) {
    cat("\nGenerations:\n")
    .ferx_search_print_table(x$generations, digits)
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}

#' @param object A \code{ferx_globalsearch} object.
#' @rdname ferx_globalsearch
#' @export
summary.ferx_globalsearch <- function(object, digits = 4, ...) {
  print(object, digits = digits)
  cat("\nEvery model:\n")
  .ferx_search_print_table(
    object$models[, c("id", "step", "genome", "n_parameters", "ofv",
                      "criterion", "fitness", "rank", "converged", "passed",
                      "duplicate_of", "reused", "failures", "error"),
                  drop = FALSE],
    digits
  )
  cat("\nPenalty schedule (what the fitness charged):\n")
  print(object$penalties)
  invisible(object)
}
