# -- Shared internals for the model-space search tools -----------------------
#
# One place for what every search tool needs, so covsearch, allometry and the
# tools that follow (modelsearch, ruvsearch, iivsearch / iovsearch, and the AMD
# pipeline over them - ferx-r #334) cannot each grow their own answer to the
# same question. Three rules are enforced here rather than restated per tool:
#
#  1. **Two entry forms, one grammar.** A tool takes either a `.ferxsearch`
#     file or inline arguments; the inline form is rendered into a config file
#     by the *engine-facing* layer and loaded by the engine's own loader. R
#     never parses or synthesises MFL, so the two forms cannot disagree about
#     what a search space means.
#  2. **The tables are the engine's.** A step or model table is built from the
#     columns the binding returns, in the engine's own order. Nothing is
#     recomputed in R, so a run's object and its `steps.csv` are one table.
#  3. **A missing value is NA, not a sentinel.** The bindings pass `NaN` for an
#     absent number and `""` for an absent tri-state logical (a candidate that
#     never fitted has no `converged`, which is not `FALSE`).

# The glue's tri-state logical: "" is NA, "true"/"false" are themselves.
.ferx_search_lgl <- function(x) {
  x <- as.character(x)
  out <- rep(NA, length(x))
  out[x %in% c("true", "TRUE")] <- TRUE
  out[x %in% c("false", "FALSE")] <- FALSE
  out
}

# NaN is the glue's "no value" for a numeric column.
.ferx_search_num <- function(x) {
  x <- as.numeric(x)
  x[is.nan(x)] <- NA_real_
  x
}

# "" is the glue's "no value" for a character column.
.ferx_search_chr <- function(x) {
  x <- as.character(x)
  x[!nzchar(x)] <- NA_character_
  x
}

# A model argument may be a path or a `ferx_model` object, exactly as
# `ferx_bootstrap()` accepts both. Returns the two paths the bindings take,
# with "" standing for "the model file's own [data] block".
.ferx_search_model_data <- function(model, data, what) {
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model <- model$model
  }
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop(sprintf("%s: `model` must be a single path to a .ferx file, or a ferx_model object", what))
  }
  if (!file.exists(model)) {
    stop(sprintf("%s: model file not found: %s", what, model))
  }
  if (is.null(data)) {
    data <- .ferx_model_data_path(model)
  }
  if (is.null(data)) {
    stop(sprintf(
      paste0("%s: no data supplied. Pass `data`, or add a `[data]` block ",
             "(`path = ...`) to the model file."),
      what
    ))
  }
  if (!is.character(data) || length(data) != 1L || is.na(data) || !file.exists(data)) {
    stop(sprintf("%s: data file not found: %s", what, as.character(data)[1L]))
  }
  list(model = normalizePath(model), data = normalizePath(data))
}

# A single finite number, or NaN for "not stated" - the sentinel the bindings
# read as "keep the engine's default".
#
# `Inf` has to be refused rather than passed through: the bindings emit a key
# only when its value is finite, so an infinite argument would be dropped on
# the way to the config file and the search would run as though it had never
# been given - the failure mode this sentinel exists to avoid, arriving by the
# other door. NA is already refused above for the same reason.
.ferx_search_scalar <- function(x, name, what, positive = FALSE) {
  if (is.null(x)) return(NaN)
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s: `%s` must be a single number or NULL", what, name))
  }
  if (!is.finite(x)) {
    stop(sprintf("%s: `%s` must be a finite number (got %s); use NULL to keep the engine default",
                 what, name, format(x)))
  }
  if (positive && x <= 0) {
    stop(sprintf("%s: `%s` must be positive", what, name))
  }
  as.numeric(x)
}

# A count for the bindings' integer sentinels: `below` (0 for threads, -1 for
# retries and max_steps) means "not stated".
.ferx_search_count <- function(x, name, what, below = 0L, min = 1L) {
  if (is.null(x)) return(as.integer(below))
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != as.integer(x)) {
    stop(sprintf("%s: `%s` must be a single whole number or NULL", what, name))
  }
  if (x < min) {
    stop(sprintf("%s: `%s` must be at least %d", what, name, min))
  }
  as.integer(x)
}

# A tri-state logical for the bindings: -1 unset, 0 FALSE, 1 TRUE.
.ferx_search_flag3 <- function(x, name, what) {
  if (is.null(x)) return(-1L)
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s: `%s` must be TRUE, FALSE or NULL", what, name))
  }
  if (x) 1L else 0L
}

.ferx_search_bool <- function(x, name, what) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s: `%s` must be TRUE or FALSE", what, name))
  }
  x
}

# "" is the bindings' "no directory": the run stays in memory, writes nothing
# and cannot be resumed.
.ferx_search_directory <- function(directory, what) {
  if (is.null(directory)) return("")
  if (!is.character(directory) || length(directory) != 1L || is.na(directory)) {
    stop(sprintf("%s: `directory` must be a single path or NULL", what))
  }
  normalizePath(directory, mustWork = FALSE)
}

# The mutual exclusion that keeps the inline form from becoming a second way of
# saying what the file already says. `config` states the search; the run knobs
# (`threads`, `retries`, `resume`, `directory`, `progress`) say how to run it
# and stay available to both forms - a user resuming a file-driven run should
# not have to edit the file to do it.
.ferx_search_entry_form <- function(config, model, data, inline, what) {
  named <- names(inline)[!vapply(inline, is.null, logical(1))]
  if (!is.null(config)) {
    if (!is.character(config) || length(config) != 1L || is.na(config)) {
      stop(sprintf("%s: `config` must be a single path to a .ferxsearch file", what))
    }
    if (!file.exists(config)) {
      stop(sprintf("%s: config file not found: %s", what, config))
    }
    # `data` belongs with `model`, not with the run knobs: the file's `data =`
    # key names the dataset the search runs on, so accepting a second one here
    # would silently search a different dataset than the one asked for.
    if (!is.null(model) || !is.null(data) || length(named) > 0L) {
      stop(sprintf(
        paste0("%s: `config` states the whole search, so `%s` cannot be given ",
               "beside it. Edit the .ferxsearch file, or drop `config` and pass ",
               "the arguments inline."),
        what,
        paste(c(if (!is.null(model)) "model", if (!is.null(data)) "data", named),
              collapse = "`, `")
      ))
    }
    return(normalizePath(config))
  }
  if (is.null(model)) {
    stop(sprintf("%s: pass either `config` (a .ferxsearch file) or `model`", what))
  }
  ""
}

# The MFL search space, quoted verbatim. Several lines are pasted with newlines
# so a multi-line space can be given as a character vector; nothing here reads
# the grammar, which is the engine's.
.ferx_search_space_text <- function(search_space, what, required = TRUE) {
  if (is.null(search_space)) {
    if (required) {
      stop(sprintf("%s: `search_space` is required in the inline form (MFL text)", what))
    }
    return("")
  }
  if (!is.character(search_space) || anyNA(search_space)) {
    stop(sprintf("%s: `search_space` must be MFL text", what))
  }
  text <- paste(search_space, collapse = "\n")
  # `""`, `character(0)` and a vector of blank lines are the same thing to the
  # engine as no `[space]` section at all, so a tool that needs a space has to
  # refuse them here too - `required` would otherwise be a rule about the
  # argument being absent rather than about the search having a space. For AMD
  # that mattered most: an empty space plans every step but the residual one as
  # skipped, which is `ferx_ruvsearch()` wearing six rows (#356 review).
  if (required && !nzchar(trimws(text))) {
    stop(sprintf(
      "%s: `search_space` is empty; the inline form needs MFL text stating what to search",
      what
    ))
  }
  text
}

# The candidate tables of a run that wrote any.
#
# A stepwise tool gives each phase its own runner directory (`base/`,
# `forward-1/`, `backward-2/`, ...), each with the `candidates.csv` the runner
# writes, so the table a user wants is those stacked with the phase named. A
# run held in memory wrote none, and a cancelled run left a partial table -
# both are `ferx_search_results()`'s own business, so this only decides which
# directories to ask about.
#
# `runs` names those directories, in order, and every tool passes the ones its
# own result reached. Scanning the directory instead would fold in whatever an
# earlier, longer run left behind: the engine rewrites the steps it executes
# but removes nothing, so a second run with fewer steps in the same directory
# would come back carrying candidate rows for steps it never fitted - a result
# object contradicting its own step table (#336 review). A directory the
# current run did not write is simply absent, which is why each one is checked
# rather than assumed.
.ferx_search_candidates <- function(directory, runs = NULL) {
  if (!nzchar(directory) || !dir.exists(directory)) return(NULL)

  read_one <- function(dir, label) {
    tab <- tryCatch(ferx_search_results(dir), error = function(e) NULL)
    if (is.null(tab) || nrow(tab) == 0L) return(NULL)
    cbind(run = label, tab, stringsAsFactors = FALSE)
  }

  parts <- list(read_one(directory, ""))
  subs <- if (is.null(runs)) {
    list.dirs(directory, recursive = FALSE, full.names = TRUE)
  } else {
    file.path(directory, unique(runs))
  }
  for (s in subs) {
    if (!dir.exists(s)) next
    parts[[length(parts) + 1L]] <- read_one(s, basename(s))
  }
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0L) return(NULL)
  out <- do.call(rbind, parts)
  # A single-directory tool has nothing to name, so the column would be a
  # column of "".
  if (all(!nzchar(out$run))) out$run <- NULL
  rownames(out) <- NULL
  out
}

# Print helper: a data frame, rounded, without row names.
.ferx_search_print_table <- function(df, digits = 4) {
  if (is.null(df) || nrow(df) == 0L) {
    cat("  (no rows)\n")
    return(invisible(NULL))
  }
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, digits = digits)
  print(df, row.names = FALSE)
  invisible(NULL)
}

# Print helper: the notes a tool emits once (effects already in the base model,
# symbols that resolved to nothing, a journal that could not be written).
.ferx_search_print_notes <- function(notes) {
  if (length(notes) == 0L) return(invisible(NULL))
  cat("\nNotes:\n")
  for (n in notes) cat("  - ", n, "\n", sep = "")
  invisible(NULL)
}
