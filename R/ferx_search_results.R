#' Read a search run's candidate table
#'
#' Reads the \code{candidates.csv} a search run writes into its directory and
#' returns it as a typed data frame - one row per candidate the run was given,
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
#' @param directory Path to the run directory, or directly to a
#'   \code{candidates.csv} / \code{candidates.partial.csv} file.
#' @param partial Which table to read: \code{NULL} (the default) prefers the
#'   complete table and falls back to the partial one, \code{TRUE} demands the
#'   partial table, \code{FALSE} demands the complete one. When
#'   \code{directory} names a file directly, a value that disagrees with the
#'   file named is an error rather than an ignored argument.
#'
#' @return A data frame with the engine's 15 columns: \code{id},
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
ferx_search_results <- function(directory, partial = NULL) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory)) {
    stop("'directory' must be a single path")
  }
  if (!is.null(partial) && (!is.logical(partial) || length(partial) != 1L ||
                            is.na(partial))) {
    stop("'partial' must be TRUE, FALSE, or NULL")
  }

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
