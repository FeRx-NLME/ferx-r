#' Preview data-selection filtering
#'
#' Applies IGNORE/ACCEPT record-level filters to a NONMEM dataset in pure R,
#' mirroring the semantics of the \code{[data_selection]} model block. This is
#' a \strong{preview} function: the canonical, engine-validated filtering
#' happens in Rust at fit time. Use \code{ferx_selection()} to inspect which
#' records will be excluded before committing to a fit, or to pass pre-built
#' conditions directly to \code{\link{ferx_fit}}.
#'
#' @param data Character path to a NONMEM CSV file, or a data.frame already
#'   read into memory. When a path is supplied it is recorded on the returned
#'   object so \code{\link{ferx_fit}(data = <result>)} can pass it through
#'   to Rust without writing a temporary file.
#' @param ignore Character vector of filter expressions. A record is
#'   \strong{excluded} when any expression evaluates to \code{TRUE}.
#'   Expressions use column names with standard comparison operators
#'   (\code{==}, \code{!=}, \code{<}, \code{<=}, \code{>}, \code{>=});
#'   use \code{&&} for AND within one expression.
#'   Example: \code{c("DV < 0.001", "EVID != 0")}.
#' @param accept Character vector of filter expressions. A record is
#'   \strong{kept} only when all conditions pass; excluded if any fail.
#' @param ignore_ids Numeric or character vector of subject IDs to exclude
#'   entirely.
#'
#' @return An object of class \code{"ferx_data"} (also a \code{data.frame})
#'   containing the \strong{retained} records. The following attributes carry
#'   metadata:
#'   \describe{
#'     \item{\code{source_path}}{Character path supplied as \code{data}, or
#'       \code{NULL} when \code{data} was a data.frame.}
#'     \item{\code{ignore}, \code{accept}, \code{ignore_ids}}{The filter
#'       conditions used (character / character / character vectors).}
#'     \item{\code{excluded_rows}}{A data.frame of the excluded records with
#'       an extra column \code{.exclude_reason} (first matching rule string).}
#'     \item{\code{exclusions}}{A named list mirroring \code{fit$exclusions}:
#'       \code{n_records_total}, \code{n_obs_excluded}, \code{n_dose_excluded},
#'       \code{n_other_excluded}, \code{excluded_subject_ids},
#'       \code{fired_ignore}, \code{fired_accept}.}
#'   }
#'
#' @seealso \code{\link{ferx_selection_excluded}} to retrieve excluded rows;
#'   \code{\link{ferx_fit}} which accepts a \code{ferx_data} object as
#'   \code{data}.
#' @family data selection
#' @export
ferx_selection <- function(data, ignore = NULL, accept = NULL, ignore_ids = NULL) {
  source_path <- NULL
  if (is.character(data) && length(data) == 1L) {
    source_path <- data
    df <- utils::read.csv(data, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("`data` must be a file path or a data.frame")
  }
  if (!is.null(ignore) && !is.character(ignore))
    stop("`ignore` must be a character vector")
  if (!is.null(accept) && !is.character(accept))
    stop("`accept` must be a character vector")
  if (!is.null(ignore_ids) && !is.numeric(ignore_ids) && !is.character(ignore_ids))
    stop("`ignore_ids` must be a numeric or character vector")

  ignore     <- as.character(ignore)
  accept     <- as.character(accept)
  ignore_ids <- as.character(ignore_ids)

  .get_col <- function(df, name) {
    idx <- match(name, tolower(names(df)))
    if (is.na(idx)) NULL else df[[idx]]
  }

  .eval_expr <- function(expr_str, df) {
    parts <- .parse_filter_expr(expr_str)
    if (is.null(parts)) return(rep(FALSE, nrow(df)))
    Reduce("&", lapply(parts, function(p) {
      col_vals <- .get_col(df, p$col)
      if (is.null(col_vals)) return(rep(FALSE, nrow(df)))
      .cmp(col_vals, p$op, p$rhs)
    }))
  }

  n <- nrow(df)
  include <- rep(TRUE, n)
  exclude_reason <- rep(NA_character_, n)

  id_col <- .get_col(df, "id")
  if (length(ignore_ids) > 0L && !is.null(id_col)) {
    ids_char <- as.character(id_col)
    mask <- ids_char %in% ignore_ids & include
    exclude_reason[mask & is.na(exclude_reason)] <-
      paste0("ignore_subjects: ", ids_char[mask & is.na(exclude_reason)])
    include[mask] <- FALSE
  }

  fired_ignore <- character(0)
  for (expr in ignore) {
    mask <- .eval_expr(expr, df) & include
    if (any(mask)) fired_ignore <- c(fired_ignore, paste0("ignore: ", expr))
    exclude_reason[mask & is.na(exclude_reason)] <- paste0("ignore: ", expr)
    include[mask] <- FALSE
  }

  fired_accept <- character(0)
  for (expr in accept) {
    mask <- !.eval_expr(expr, df) & include
    if (any(mask)) fired_accept <- c(fired_accept, paste0("accept: ", expr))
    exclude_reason[mask & is.na(exclude_reason)] <- paste0("accept: ", expr)
    include[mask] <- FALSE
  }

  retained <- df[include, , drop = FALSE]
  excluded_df <- df[!include, , drop = FALSE]
  excluded_df$.exclude_reason <- exclude_reason[!include]

  evid_col <- .get_col(df[!include, , drop = FALSE], "evid")
  mdv_col  <- .get_col(df[!include, , drop = FALSE], "mdv")
  excl_evid <- if (!is.null(evid_col)) as.integer(evid_col) else integer(nrow(excluded_df))
  excl_mdv  <- if (!is.null(mdv_col))  as.integer(mdv_col)  else integer(nrow(excluded_df))
  n_obs_excl   <- sum(excl_evid == 0L & excl_mdv == 0L, na.rm = TRUE)
  n_dose_excl  <- sum(excl_evid %in% c(1L, 4L), na.rm = TRUE)
  n_other_excl <- nrow(excluded_df) - n_obs_excl - n_dose_excl

  excl_id_col <- .get_col(excluded_df, "id")
  incl_id_col <- .get_col(retained, "id")
  excl_subj_ids <- character(0)
  if (!is.null(excl_id_col) && !is.null(incl_id_col)) {
    all_excl <- unique(as.character(excl_id_col))
    still_in <- unique(as.character(incl_id_col))
    excl_subj_ids <- setdiff(all_excl, still_in)
  }

  exclusions <- list(
    n_records_total      = n,
    n_obs_excluded       = n_obs_excl,
    n_dose_excluded      = n_dose_excl,
    n_other_excluded     = n_other_excl,
    excluded_subject_ids = excl_subj_ids,
    fired_ignore         = fired_ignore,
    fired_accept         = fired_accept
  )

  structure(
    retained,
    class        = c("ferx_data", "data.frame"),
    source_path  = source_path,
    ignore       = ignore,
    accept       = accept,
    ignore_ids   = ignore_ids,
    excluded_rows = excluded_df,
    exclusions   = exclusions
  )
}

#' Print a ferx_data object
#'
#' Shows a compact exclusion summary (records retained, excluded counts, fired
#' rules) followed by the first rows of the retained data.
#'
#' @param x A \code{ferx_data} object.
#' @param n Number of data rows to display. Default 6.
#' @param ... Passed to \code{print.data.frame}.
#' @export
print.ferx_data <- function(x, n = 6L, ...) {
  ex      <- attr(x, "exclusions", exact = TRUE)
  n_total <- if (!is.null(ex)) ex$n_records_total else nrow(x)
  n_ret   <- nrow(x)
  n_excl  <- if (!is.null(ex))
    (ex$n_obs_excluded %||% 0L) + (ex$n_dose_excluded %||% 0L) +
      (ex$n_other_excluded %||% 0L)
  else 0L

  cat(sprintf("<ferx_data>  %d of %d records retained", n_ret, n_total))
  if (n_excl > 0L)
    cat(sprintf("  (%d obs, %d doses, %d other excluded)",
                ex$n_obs_excluded, ex$n_dose_excluded, ex$n_other_excluded))
  cat("\n")

  if (!is.null(ex)) {
    strip_pfx <- function(x) sub("^(ignore|accept): ", "", x)
    for (cond in ex$fired_ignore)
      cat(sprintf("  Fired ignore: %s\n", strip_pfx(cond)))
    for (cond in ex$fired_accept)
      cat(sprintf("  Fired accept: %s\n", strip_pfx(cond)))
    if (length(ex$excluded_subject_ids) > 0L)
      cat(sprintf("  Subjects excluded: %s\n",
                  paste(ex$excluded_subject_ids, collapse = ", ")))
  }
  if (n_excl > 0L)
    cat("  Use ferx_selection_excluded(x) to inspect excluded records.\n")
  cat("\n")

  print(utils::head(as.data.frame(x), n = n), ...)
  if (n_ret > n)
    cat(sprintf("  ... %d more rows\n", n_ret - n))

  invisible(x)
}

#' Retrieve excluded records from a ferx_data or ferx_fit object
#'
#' @param x A \code{ferx_data} object returned by \code{\link{ferx_selection}},
#'   or a \code{ferx_fit} object.
#' @param ... Ignored.
#' @return A data.frame of excluded records. For \code{ferx_data} objects the
#'   frame includes a \code{.exclude_reason} column naming the first rule that
#'   excluded each record.
#' @family data selection
#' @export
ferx_selection_excluded <- function(x, ...) UseMethod("ferx_selection_excluded")

#' @export
ferx_selection_excluded.ferx_data <- function(x, ...) attr(x, "excluded_rows")

#' @export
ferx_selection_excluded.ferx_fit <- function(x, ...) {
  ex <- x$exclusions
  if (is.null(ex)) {
    message("No data-selection rules were active for this fit.")
    return(invisible(NULL))
  }
  dp <- x$data_path
  if (is.null(dp) || !nzchar(dp) || !file.exists(dp)) {
    message("Data path not available on fit object; cannot retrieve excluded records.")
    return(invisible(NULL))
  }
  n_excl <- ex$n_obs_excluded + ex$n_dose_excluded + ex$n_other_excluded
  if (n_excl == 0L && length(ex$excluded_subject_ids) == 0L) {
    message("No records were excluded by data-selection rules.")
    return(invisible(data.frame()))
  }
  # Replay the fired conditions through the R-side preview filter to obtain the
  # actual excluded rows with their .exclude_reason column.
  # fired_ignore / fired_accept carry a "ignore: " / "accept: " prefix (matching
  # the Rust ExclusionSummary format); strip it to recover bare expressions for
  # ferx_selection(ignore=...).
  strip_prefix <- function(x) sub("^(ignore|accept): ", "", x)
  sel <- ferx_selection(
    dp,
    ignore     = strip_prefix(ex$fired_ignore),
    accept     = strip_prefix(ex$fired_accept),
    ignore_ids = ex$excluded_subject_ids
  )
  ferx_selection_excluded(sel)
}

.get_col_sel <- function(df, name) {
  idx <- match(name, tolower(names(df)))
  if (is.na(idx)) NULL else df[[idx]]
}

.parse_filter_expr <- function(expr) {
  ops <- c("==", "!=", ">=", "<=", ">", "<")
  and_parts <- strsplit(expr, "&&")[[1L]]
  result <- vector("list", length(and_parts))
  for (i in seq_along(and_parts)) {
    part <- trimws(and_parts[[i]])
    found <- FALSE
    for (op in ops) {
      pos <- regexpr(op, part, fixed = TRUE)[1L]
      if (pos < 1L) next
      lhs <- trimws(substr(part, 1L, pos - 1L))
      rhs_raw <- trimws(substr(part, pos + nchar(op), nchar(part)))
      if (!nzchar(lhs) || !nzchar(rhs_raw)) next
      rhs <- if (grepl('^["\'].*["\']$', rhs_raw)) {
        substr(rhs_raw, 2L, nchar(rhs_raw) - 1L)
      } else if (suppressWarnings(!is.na(as.numeric(rhs_raw)))) {
        as.numeric(rhs_raw)
      } else {
        rhs_raw
      }
      result[[i]] <- list(col = tolower(lhs), op = op, rhs = rhs)
      found <- TRUE
      break
    }
    if (!found) return(NULL)
  }
  result
}

.cmp <- function(lhs, op, rhs) {
  if (is.na(rhs) || (is.numeric(rhs) && is.nan(rhs))) return(rep(FALSE, length(lhs)))
  lhs_num <- suppressWarnings(as.numeric(lhs))
  if (is.numeric(rhs)) {
    lhs_use <- lhs_num
    lhs_use[is.na(lhs_use)] <- NaN
  } else {
    lhs_use <- as.character(lhs)
  }
  switch(op,
    "==" = !is.na(lhs_use) & lhs_use == rhs,
    "!=" = !is.na(lhs_use) & lhs_use != rhs,
    "<"  = !is.na(lhs_use) & lhs_use <  rhs,
    "<=" = !is.na(lhs_use) & lhs_use <= rhs,
    ">"  = !is.na(lhs_use) & lhs_use >  rhs,
    ">=" = !is.na(lhs_use) & lhs_use >= rhs,
    rep(FALSE, length(lhs_use))
  )
}
