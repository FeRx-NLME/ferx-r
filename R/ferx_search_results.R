#' Read a search run's table
#'
#' Reads a table a search run writes into its directory - the runner's
#' \code{candidates.csv} by default, a structural search's \code{models.csv}
#' with \code{type = "models"}, or a stepwise search's \code{steps.csv} with
#' \code{type = "steps"} - and returns it typed. The column list comes from the
#' engine for every table, never a copy maintained here, so the R data frame
#' and the file on disk are the same table.
#'
#' The candidate table is one row per candidate the run was given and
#' in the order it was given them, including the ones that failed. A candidate
#' the strictness gate excluded carries \emph{why} in \code{failures} rather
#' than being absent, because a candidate missing from a report cannot be told
#' apart from one that was never generated.
#'
#' A cancelled run writes \code{candidates.partial.csv} instead, leaving any
#' complete table beside it untouched; this function reads the complete table
#' when there is one and falls back to the partial table otherwise. The column
#' list comes from the engine (never a copy maintained here), so a column the
#' engine adds shows up as a column here.
#'
#' The engine writes an empty cell for "there is no value" - a candidate that
#' never fitted has no OFV, a row that succeeded has nothing to say about
#' whether a resume would retry it. Those cells come back as \code{NA}, not
#' \code{NaN} and not \code{""}.
#'
#' @param directory Path to the run directory, or directly to the
#'   \code{candidates.csv} / \code{candidates.partial.csv} /
#'   \code{models.csv} / \code{steps.csv} file itself.
#' @param type Which table to read: \code{"candidates"} (the default) for the
#'   runner's own candidate table, written by every tool, \code{"models"} for
#'   the model table \code{\link{ferx_modelsearch}},
#'   \code{\link{ferx_iivsearch}} and \code{\link{ferx_iovsearch}} write, or
#'   \code{"steps"} for the step table \code{\link{ferx_covsearch}} and
#'   \code{\link{ferx_ruvsearch}} write. A run has more than one: the candidate
#'   table is one row per fit the runner was asked for, the model or step table
#'   one row per model the search itself decided on. Several tools write a file
#'   called \code{models.csv}, and both stepwise tools one called
#'   \code{steps.csv}; which schema a file carries is read off its own header,
#'   and reported as the \code{tool} attribute.
#' @param partial Which table to read: \code{NULL} (the default) prefers the
#'   complete table and falls back to the partial one, \code{TRUE} demands the
#'   partial table, \code{FALSE} demands the complete one. When
#'   \code{directory} names a file directly, a value that disagrees with the
#'   file named is an error rather than an ignored argument.
#'
#' @return For \code{type = "steps"}, a data frame with the engine's 17 step
#'   columns for whichever stepwise tool wrote the file - covsearch's
#'   (\code{step}, \code{phase}, \code{candidate}, \code{parameter},
#'   \code{covariate}, \code{form}, ...) or ruvsearch's (\code{iteration},
#'   \code{candidate}, \code{feature}, \code{screened}, ...) - typed the same
#'   way, and carrying a \code{tool} attribute naming the one it matched.
#'   For \code{type = "models"}, a data frame with the engine's model-table
#'   columns for whichever tool wrote the file - modelsearch's 21 (\code{id},
#'   \code{parent}, \code{layer}, \code{path}, the four structural columns,
#'   \code{criterion}, \code{rank}, ...), iivsearch's 18 (\code{id},
#'   \code{parent}, \code{step}, \code{description}, \code{etas},
#'   \code{blocks}, ...) or iovsearch's 19 (..., \code{kappas},
#'   \code{kappa_blocks}, ...) - typed the same way and carrying the same
#'   \code{tool} attribute. Otherwise a data frame with the
#'   engine's 15 candidate columns: \code{id},
#'   \code{parent}, \code{hash}, \code{features}, \code{criterion} (numeric),
#'   \code{ofv} (numeric), \code{converged} (logical), \code{passed}
#'   (logical), \code{failures}, \code{skipped}, \code{seconds} (numeric),
#'   \code{error}, \code{retryable} (logical), \code{duplicate_of} and
#'   \code{reused} (logical). Carries attributes \code{path} (the file read)
#'   and \code{partial} (whether it is a cancelled run's table).
#'
#' @examples
#' \dontrun{
#' res <- ferx_search_results("search-run-1")
#' res[res$passed, c("id", "features", "criterion")]
#' # Candidates the gate excluded, and why
#' res[!res$passed, c("id", "failures")]
#' }
#'
#' @seealso \code{\link{ferx_search_config}}, \code{\link{check_strictness}}
#' @family search
#' @export
ferx_search_results <- function(directory, partial = NULL,
                               type = c("candidates", "models", "steps")) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory)) {
    stop("'directory' must be a single path")
  }
  if (!is.null(partial) && (!is.logical(partial) || length(partial) != 1L ||
                            is.na(partial))) {
    stop("'partial' must be TRUE, FALSE, or NULL")
  }
  type <- match.arg(type)

  # The model table has no partial twin: a cancelled structural search still
  # writes one `models.csv`, holding the models it reached.
  if (type == "models") return(.ferx_read_model_table(directory, partial))
  if (type == "steps") return(.ferx_read_step_table(directory, partial))

  complete_path <- file.path(directory, "candidates.csv")
  partial_path  <- file.path(directory, "candidates.partial.csv")

  if (!dir.exists(directory) && file.exists(directory)) {
    # A file was passed directly. `partial` still means what it says: it
    # demands a table of that kind, so a mismatch is an error rather than a
    # silently ignored argument.
    path <- directory
    is_partial <- grepl("candidates\\.partial\\.csv$", path)
    if (!is.null(partial) && partial != is_partial) {
      stop("`", path, "` is a ", if (is_partial) "partial" else "complete",
           " candidate table, but partial = ", partial, " was asked for")
    }
  } else if (isTRUE(partial)) {
    if (!file.exists(partial_path)) {
      stop("No partial candidate table in ", directory,
           " (looked for candidates.partial.csv)")
    }
    path <- partial_path
    is_partial <- TRUE
  } else if (isFALSE(partial)) {
    if (!file.exists(complete_path)) {
      stop("No candidate table in ", directory, " (looked for candidates.csv)")
    }
    path <- complete_path
    is_partial <- FALSE
  } else if (file.exists(complete_path)) {
    path <- complete_path
    is_partial <- FALSE
  } else if (file.exists(partial_path)) {
    path <- partial_path
    is_partial <- TRUE
  } else {
    stop("No candidate table in ", directory,
         " (looked for candidates.csv and candidates.partial.csv)")
  }

  raw <- utils::read.csv(path, colClasses = "character", check.names = FALSE)

  # The engine owns the column list; a table missing one is a table this
  # version of ferx cannot read, and saying which column is missing beats an
  # NULL column surfacing three lines later.
  expected <- ferx_rust_search_table_columns()
  missing <- setdiff(expected, names(raw))
  if (length(missing)) {
    stop("`", path, "` is not a search candidate table - missing column",
         if (length(missing) == 1L) " " else "s ",
         paste(missing, collapse = ", "))
  }
  # Engine order first, anything the file carries beyond it after - a column a
  # newer engine wrote is worth keeping, not silently dropping.
  raw <- raw[, c(expected, setdiff(names(raw), expected)), drop = FALSE]

  num <- c("criterion", "ofv", "seconds")
  lgl <- c("converged", "passed", "retryable", "reused")
  # Every remaining column, not only the engine's own: a column a newer engine
  # wrote is kept above, so its empty cells must become NA like the rest.
  chr <- setdiff(names(raw), c(num, lgl))
  for (col in num) raw[[col]] <- .ferx_csv_num(raw[[col]])
  for (col in lgl) raw[[col]] <- .ferx_csv_lgl(raw[[col]])
  for (col in chr) raw[[col]] <- .ferx_csv_chr(raw[[col]])

  attr(raw, "path") <- path
  attr(raw, "partial") <- is_partial
  raw
}

# The engine writes an empty cell for a non-finite / absent number.
.ferx_csv_num <- function(x) {
  x <- as.character(x)
  suppressWarnings(as.numeric(ifelse(nzchar(x), x, NA_character_)))
}

# ... and for a logical with nothing to say (`converged` on a candidate that
# never fitted, `retryable` on one that did not fail).
.ferx_csv_lgl <- function(x) {
  x <- tolower(as.character(x))
  out <- rep(NA, length(x))
  out[x == "true"] <- TRUE
  out[x == "false"] <- FALSE
  out
}

.ferx_csv_chr <- function(x) {
  x <- as.character(x)
  ifelse(nzchar(x), x, NA_character_)
}

# The model table of a structural or variability search (`models.csv`), typed
# the same way as the candidate table above and against the engine's own column
# lists.
#
# Three tools write a file of that name - modelsearch's 21 structural columns,
# iivsearch's 18 and iovsearch's 19 - so the file's own header says which tool
# wrote it, and a file matching none is named as such rather than typed against
# the wrong schema. The three are told apart by a column only one of them has
# (`layer`, `blocks`, `kappas`), which is why a partial match is not enough.
#
# There is no `models.partial.csv`: a cancelled search writes the one table
# with the rows it reached, and says so on the object it returns. `partial` is
# therefore only accepted as FALSE / NULL, rather than being ignored.
.ferx_read_model_table <- function(directory, partial) {
  if (isTRUE(partial)) {
    stop("A search writes no partial model table; a cancelled run's ",
         "models.csv holds the models it reached")
  }
  path <- if (!dir.exists(directory) && file.exists(directory)) {
    directory
  } else {
    file.path(directory, "models.csv")
  }
  if (!file.exists(path)) {
    stop("No model table in ", directory, " (looked for models.csv)")
  }

  raw <- utils::read.csv(path, colClasses = "character", check.names = FALSE)
  schemas <- list(
    modelsearch = ferx_rust_modelsearch_columns(),
    iivsearch   = ferx_rust_iivsearch_columns(),
    iovsearch   = ferx_rust_iovsearch_columns()
  )
  hit <- vapply(schemas, function(cols) length(setdiff(cols, names(raw))) == 0L,
                logical(1))
  if (!any(hit)) {
    stop("`", path, "` is not a search model table - it carries none of the ",
         "structural columns (", paste(schemas$modelsearch, collapse = ", "),
         "), the variability columns (",
         paste(schemas$iivsearch, collapse = ", "),
         ") or the inter-occasion columns (",
         paste(schemas$iovsearch, collapse = ", "), ")")
  }
  tool <- names(schemas)[hit][1L]
  expected <- schemas[[tool]]
  # Engine order first, anything the file carries beyond it after - a column a
  # newer engine wrote is worth keeping, not silently dropping.
  raw <- raw[, c(expected, setdiff(names(raw), expected)), drop = FALSE]

  int <- intersect(c("layer", "step", "peripherals", "n_parameters", "rank",
                     "starts"), names(raw))
  num <- intersect(c("ofv", "criterion", "d_criterion", "seconds"), names(raw))
  lgl <- intersect(c("converged", "passed", "selected", "continued", "reused"),
                   names(raw))
  chr <- setdiff(names(raw), c(int, num, lgl))
  for (col in int) raw[[col]] <- as.integer(.ferx_csv_num(raw[[col]]))
  for (col in num) raw[[col]] <- .ferx_csv_num(raw[[col]])
  for (col in lgl) raw[[col]] <- .ferx_csv_lgl(raw[[col]])
  for (col in chr) raw[[col]] <- .ferx_csv_chr(raw[[col]])

  attr(raw, "path") <- path
  attr(raw, "partial") <- FALSE
  attr(raw, "tool") <- tool
  raw
}

# The step table of a stepwise search (`steps.csv`), typed the same way as the
# tables above and against the engine's own column lists.
#
# covsearch and ruvsearch both write a file of that name, 17 columns each and
# only `candidate` in common - so the file's own header says which tool wrote
# it, and a file matching neither is named as such rather than typed against
# the wrong schema. There is no `steps.partial.csv`: a cancelled run writes the
# one table with the rows it reached, and says so on the object it returns.
.ferx_read_step_table <- function(directory, partial) {
  if (isTRUE(partial)) {
    stop("A stepwise search writes no partial step table; a cancelled run's ",
         "steps.csv holds the steps it reached")
  }
  path <- if (!dir.exists(directory) && file.exists(directory)) {
    directory
  } else {
    file.path(directory, "steps.csv")
  }
  if (!file.exists(path)) {
    stop("No step table in ", directory, " (looked for steps.csv)")
  }

  raw <- utils::read.csv(path, colClasses = "character", check.names = FALSE)
  schemas <- list(
    covsearch = ferx_rust_covsearch_columns(),
    ruvsearch = ferx_rust_ruvsearch_columns()
  )
  hit <- vapply(schemas, function(cols) length(setdiff(cols, names(raw))) == 0L,
                logical(1))
  if (!any(hit)) {
    stop("`", path, "` is not a search step table - it carries neither the ",
         "covariate columns (", paste(schemas$covsearch, collapse = ", "),
         ") nor the residual-error columns (",
         paste(schemas$ruvsearch, collapse = ", "), ")")
  }
  tool <- names(schemas)[hit][1L]
  expected <- schemas[[tool]]
  # Engine order first, anything the file carries beyond it after - a column a
  # newer engine wrote is worth keeping, not silently dropping.
  raw <- raw[, c(expected, setdiff(names(raw), expected)), drop = FALSE]

  int <- intersect(c("step", "iteration", "df"), names(raw))
  num <- intersect(c("parent_ofv", "ofv", "dofv", "p_value", "alpha",
                     "cwres_dofv", "seconds"), names(raw))
  lgl <- intersect(c("significant", "selected", "converged", "passed",
                     "screened"), names(raw))
  chr <- setdiff(names(raw), c(int, num, lgl))
  for (col in int) raw[[col]] <- as.integer(.ferx_csv_num(raw[[col]]))
  for (col in num) raw[[col]] <- .ferx_csv_num(raw[[col]])
  for (col in lgl) raw[[col]] <- .ferx_csv_lgl(raw[[col]])
  for (col in chr) raw[[col]] <- .ferx_csv_chr(raw[[col]])

  attr(raw, "path") <- path
  attr(raw, "partial") <- FALSE
  attr(raw, "tool") <- tool
  raw
}
