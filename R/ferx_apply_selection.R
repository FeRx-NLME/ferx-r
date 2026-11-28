#' Preview data-selection filtering
#'
#' Applies IGNORE/ACCEPT record-level filters to a NONMEM dataset in pure R,
#' approximating the semantics of the \code{[data_selection]} model block.
#' This is a \strong{preview} function: the canonical, engine-validated
#' filtering happens in Rust at fit time. The R-side implementation is
#' intentionally lenient: an unparseable expression, or a condition on a column
#' the data does not have, is silently treated as never matching, whereas Rust
#' will error on invalid filters. A condition that never matches does opposite
#' things under the two filters: an \code{ignore} on it excludes no record, but
#' an \code{accept} on it excludes \strong{every} record, since no record
#' passes it. Use \code{ferx_apply_selection()} to inspect which records will
#' be excluded before committing to a fit, or to pass pre-built conditions
#' directly to \code{\link{ferx_fit}}.
#'
#' @param data One of:
#'   \itemize{
#'     \item a character path to a NONMEM CSV file (recorded on the returned
#'       object so \code{\link{ferx_fit}(data = <result>)} can pass it through
#'       to Rust without writing a temporary file),
#'     \item a data.frame already read into memory,
#'     \item a \code{ferx_data} object returned by an earlier call (its cached
#'       selection is reused; combine with \code{excluded = TRUE} to pull the
#'       excluded records back out), or
#'     \item a \code{ferx_fit} object, whose fired data-selection rules are
#'       replayed against the original dataset.
#'   }
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
#' @param excluded Logical. When \code{FALSE} (default) the retained records are
#'   returned as a \code{ferx_data} object. When \code{TRUE} the function returns
#'   a plain data.frame of the \strong{excluded} records instead, with an extra
#'   \code{.exclude_reason} column naming the first rule that dropped each record
#'   (this replaces the former \code{ferx_selection_excluded()} accessor).
#' @param model Optional path to the \code{.ferx} model file (or a
#'   \code{ferx_model} object) the data will be fitted with. Its
#'   \code{[data_selection]} clauses (\code{ignore}, \code{accept},
#'   \code{ignore_subjects}) are applied first, then \code{ignore},
#'   \code{accept} and \code{ignore_ids}, so the preview selects the records
#'   \code{\link{ferx_fit}(model, data)} with the same arguments would. The
#'   clauses are still evaluated by the lenient R preview described above.
#'   When the model declares a \code{[data]} column map (e.g.
#'   \code{time = TAFD}), the expressions and IDs are evaluated against the
#'   columns under the names the engine gives them, so \code{"TIME > 24"}
#'   matches the \code{TAFD} column as it does at fit time. The returned
#'   records keep the dataset's own headers. Ignored when \code{data} is a
#'   \code{ferx_data} object; for a \code{ferx_fit} the rules the fit fired are
#'   replayed, with only the column map read from the fit's model file.
#'
#' @return When \code{excluded = FALSE}, an object of class \code{"ferx_data"}
#'   (also a \code{data.frame}) containing the \strong{retained} records. When
#'   \code{excluded = TRUE}, a plain data.frame of the excluded records. The
#'   \code{ferx_data} object carries these attributes:
#'   \describe{
#'     \item{\code{source_path}}{Character path supplied as \code{data}, or
#'       \code{NULL} when \code{data} was a data.frame.}
#'     \item{\code{source_data}}{The data.frame supplied as \code{data}, all
#'       records included, or \code{NULL} when \code{data} was a path.
#'       \code{\link{ferx_fit}} fits it with the stored conditions, as it fits
#'       \code{source_path}.}
#'     \item{\code{retained_rows}}{The retained records as selected. When the
#'       object's rows no longer match them - it was subset or edited -
#'       \code{\link{ferx_fit}} fits the rows it now holds, with a warning,
#'       instead of \code{source_path} / \code{source_data}.}
#'     \item{\code{ignore}, \code{accept}, \code{ignore_ids}}{Every filter
#'       condition the preview applied, the \code{model}'s
#'       \code{[data_selection]} clauses first (character / character /
#'       character vectors). \code{\link{ferx_fit}} sends the ones its own
#'       model does not state.}
#'     \item{\code{excluded_rows}}{A data.frame of the excluded records with
#'       an extra column \code{.exclude_reason} (first matching rule string).}
#'     \item{\code{exclusions}}{A named list mirroring \code{fit$exclusions}:
#'       \code{n_records_total}, \code{n_obs_excluded}, \code{n_dose_excluded},
#'       \code{n_other_excluded}, \code{excluded_subject_ids},
#'       \code{fired_ignore}, \code{fired_accept}.}
#'   }
#'
#' @seealso \code{\link{ferx_fit}} which accepts a \code{ferx_data} object as
#'   \code{data}. Pass \code{excluded = TRUE} to retrieve the excluded rows.
#' @family data selection
#' @export
ferx_apply_selection <- function(data, ignore = NULL, accept = NULL,
                                 ignore_ids = NULL, excluded = FALSE,
                                 model = NULL) {
  # An already-selected object caches its excluded rows on an attribute; a
  # fitted object carries its fired conditions. Handle both before the raw
  # path/data.frame compute so callers can recover the selection - or, with
  # excluded = TRUE, the excluded rows - straight from those objects.
  if (inherits(data, "ferx_data")) {
    if (isTRUE(excluded)) return(attr(data, "excluded_rows"))
    return(data)
  }
  if (inherits(data, "ferx_fit")) {
    return(.ferx_selection_from_fit(data, excluded = excluded))
  }

  source_path <- NULL
  if (is.character(data) && length(data) == 1L) {
    source_path <- data
    df <- utils::read.csv(data, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("`data` must be a file path, a data.frame, a ferx_data, or a ferx_fit object")
  }
  if (!is.null(ignore) && !is.character(ignore))
    stop("`ignore` must be a character vector")
  if (!is.null(accept) && !is.character(accept))
    stop("`accept` must be a character vector")
  if (!is.null(ignore_ids) && !is.numeric(ignore_ids) && !is.character(ignore_ids))
    stop("`ignore_ids` must be a numeric or character vector")

  if (inherits(model, "ferx_model")) model <- model$model
  if (!is.null(model) && !(is.character(model) && length(model) == 1L &&
                           file.exists(model))) {
    stop("`model` must be a path to an existing .ferx file or a ferx_model object")
  }
  # The model's own `[data_selection]` clauses go first: the engine's reader
  # evaluates the block's clauses before the call's (#794).
  block <- .ferx_model_data_selection(model)
  .ferx_select(
    df, source_path,
    ignore     = unique(c(block$ignore, as.character(ignore))),
    accept     = unique(c(block$accept, as.character(accept))),
    ignore_ids = unique(c(block$ignore_ids, as.character(ignore_ids))),
    column_map = .ferx_model_column_map(model),
    excluded   = excluded
  )
}

# The selection itself: evaluate `ignore` / `accept` / `ignore_ids` exactly as
# given on `df` (renamed by `column_map` for evaluation only) and build the
# ferx_data, or the excluded rows when `excluded = TRUE`. Reads no model file,
# so a fit's replay (.ferx_selection_from_fit) applies the fit's fired clauses
# and nothing a since-edited model file states.
.ferx_select <- function(df, source_path, ignore, accept, ignore_ids,
                         column_map, excluded = FALSE) {
  ignore     <- as.character(ignore)
  accept     <- as.character(accept)
  ignore_ids <- as.character(ignore_ids)

  # Filters are evaluated on `df_eval`, the dataset under the column names the
  # engine's reader gives it (the model's `[data]` map applied); the rows
  # returned come from `df`, with the dataset's own headers, so a ferx_data
  # passed on to ferx_fit() still matches the model's map.
  df_eval <- .ferx_apply_column_map(df, column_map)

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

  fired_ignore <- character(0)
  # Subject IDs as the file spells them: the engine matches `ignore_subjects`
  # / `ignore_ids` and reports subjects verbatim, and read.csv() would turn
  # the ID `02` into 2 (#795 review, row 3).
  ids_char <- .ferx_id_text(source_path, column_map, .get_col(df_eval, "id"))
  if (length(ignore_ids) > 0L && !is.null(ids_char)) {
    mask <- ids_char %in% ignore_ids & include
    exclude_reason[mask & is.na(exclude_reason)] <-
      paste0("ignore_subjects: ", ids_char[mask & is.na(exclude_reason)])
    include[mask] <- FALSE
    # The engine names each listed subject that removed a record among its
    # fired ignore rules (#794 audit A3).
    fired_ignore <- unique(exclude_reason[mask])
  }

  for (expr in ignore) {
    mask <- .eval_expr(expr, df_eval) & include
    if (any(mask)) fired_ignore <- c(fired_ignore, paste0("ignore: ", expr))
    exclude_reason[mask & is.na(exclude_reason)] <- paste0("ignore: ", expr)
    include[mask] <- FALSE
  }

  fired_accept <- character(0)
  for (expr in accept) {
    mask <- !.eval_expr(expr, df_eval) & include
    if (any(mask)) fired_accept <- c(fired_accept, paste0("accept: ", expr))
    exclude_reason[mask & is.na(exclude_reason)] <- paste0("accept: ", expr)
    include[mask] <- FALSE
  }

  # The engine's reader lists a rule when it first removes a record, so its
  # fired lists run in the order of each rule's first excluded record, not in
  # clause order (#794). Each rule claims exactly the records whose reason it
  # is, so its first one is the first match of that reason.
  by_first_row <- function(fired) {
    fired[order(match(fired, exclude_reason))]
  }
  fired_ignore <- by_first_row(fired_ignore)
  fired_accept <- by_first_row(fired_accept)

  retained <- df[include, , drop = FALSE]
  excluded_df <- df[!include, , drop = FALSE]
  excluded_df$.exclude_reason <- exclude_reason[!include]

  evid_col <- .get_col(df_eval[!include, , drop = FALSE], "evid")
  mdv_col  <- .get_col(df_eval[!include, , drop = FALSE], "mdv")
  excl_evid <- if (!is.null(evid_col)) as.integer(evid_col) else integer(nrow(excluded_df))
  excl_mdv  <- if (!is.null(mdv_col))  as.integer(mdv_col)  else integer(nrow(excluded_df))
  n_obs_excl   <- sum(excl_evid == 0L & excl_mdv == 0L, na.rm = TRUE)
  n_dose_excl  <- sum(excl_evid %in% c(1L, 4L), na.rm = TRUE)
  n_other_excl <- nrow(excluded_df) - n_obs_excl - n_dose_excl

  excl_subj_ids <- character(0)
  if (!is.null(ids_char)) {
    all_excl <- unique(ids_char[!include])
    still_in <- unique(ids_char[include])
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

  if (isTRUE(excluded)) return(excluded_df)

  structure(
    retained,
    class        = c("ferx_data", "data.frame"),
    source_path  = source_path,
    # The whole frame a selection was made from, so ferx_fit() can hand the
    # engine the dataset and the rules, as it does for a path (#419).
    source_data  = if (is.null(source_path)) df,
    # The rows as selected, so ferx_fit() can tell that they were since subset
    # or edited and fit what the object now holds instead (#419). Shares its
    # columns with the object until one of them is modified.
    retained_rows = .ferx_selection_rows(retained),
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
#' @return \code{x}, invisibly.
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
    cat("  Pass excluded = TRUE to ferx_apply_selection() to inspect excluded records.\n")
  cat("\n")

  print(utils::head(as.data.frame(x), n = n), ...)
  if (n_ret > n)
    cat(sprintf("  ... %d more rows\n", n_ret - n))

  invisible(x)
}

# Replay a fitted object's fired data-selection rules against the original
# dataset. Backs `ferx_apply_selection(<ferx_fit>, excluded = ...)`: returns the
# retained ferx_data (excluded = FALSE) or the excluded rows (excluded = TRUE).
.ferx_selection_from_fit <- function(fit, excluded = FALSE) {
  ex <- fit$exclusions
  if (is.null(ex)) {
    message("No data-selection rules were active for this fit.")
    return(invisible(NULL))
  }
  dp <- .ferx_fit_data_path(fit)
  if (is.null(dp) || !nzchar(dp) || !file.exists(dp)) {
    message("Data path not available on fit object; cannot retrieve excluded records.")
    return(invisible(NULL))
  }
  n_excl <- ex$n_obs_excluded + ex$n_dose_excluded + ex$n_other_excluded
  if (isTRUE(excluded) && n_excl == 0L && length(ex$excluded_subject_ids) == 0L) {
    message("No records were excluded by data-selection rules.")
    return(invisible(data.frame()))
  }
  # fired_ignore / fired_accept carry an "ignore: " / "accept: " prefix (matching
  # the Rust ExclusionSummary format); strip it to recover bare expressions.
  # fired_ignore also names each listed subject as "ignore_subjects: <ID>":
  # that is an ID, not an expression, so it goes back into ignore_ids (#794
  # audit A4) - as an `ignore` clause the engine refuses it as malformed.
  fired_ignore <- as.character(ex$fired_ignore)
  is_subject <- startsWith(fired_ignore, "ignore_subjects: ")
  df <- utils::read.csv(dp, stringsAsFactors = FALSE, check.names = FALSE)
  # The worker, not ferx_apply_selection(model = ): the replay applies the
  # rules the fit fired and nothing else, so a `[data_selection]` clause added
  # to the model file since the fit must not join them. Only the file's
  # `[data]` column map is read, when the file is still there, so the replayed
  # rules see the columns the fit's reader saw.
  .ferx_select(
    df, dp,
    ignore     = sub("^ignore: ", "", fired_ignore[!is_subject]),
    accept     = sub("^accept: ", "", as.character(ex$fired_accept)),
    ignore_ids = sub("^ignore_subjects: ", "", fired_ignore[is_subject]),
    column_map = .ferx_model_column_map(.ferx_existing_path(fit$model_path)),
    excluded   = excluded
  )
}

# A ferx_data's rows as a plain data.frame - its columns and their names,
# none of the selection attributes, and row names reset to 1..n. ferx_fit()
# compares it with `retained_rows`, and fits it when the two differ. Row names
# are not data the engine reads, and dplyr renumbers them on any verb (a
# no-op `filter(TRUE)` turns 2, 3, 4 into 1, 2, 3), so they must not count as
# an edit; a change in the rows' order or content still does.
.ferx_selection_rows <- function(x) {
  n <- .row_names_info(x, 2L)
  attributes(x) <- list(names = attr(x, "names"),
                        row.names = .set_row_names(n),
                        class = "data.frame")
  x
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
      } else if (suppressWarnings(is.finite(as.numeric(rhs_raw)))) {
        # Only a finite literal is numeric. Inf/-Inf (and NaN, which the previous
        # is.na() test already rejected) fall through to the string branch and are
        # treated as label text, matching ferx-core's filter_expr parser.
        as.numeric(rhs_raw)
      } else {
        rhs_raw
      }
      # Ordered comparisons require a numeric RHS. ferx-core rejects an ordered
      # op against a non-numeric value (e.g. `BW < abc`) at parse time; mirror
      # that here by treating the whole clause as unparseable (the lenient
      # preview then never matches it: an ignore excludes nothing, an accept
      # excludes every record) rather than silently doing a lexical
      # string comparison the Rust fit never performs.
      if (!is.numeric(rhs) && op %in% c(">=", "<=", ">", "<")) return(NULL)
      result[[i]] <- list(col = tolower(lhs), op = op, rhs = rhs)
      found <- TRUE
      break
    }
    if (!found) {
      # No operator. Accept a lone column name as the NONMEM `IGNORE=C`
      # shorthand: a bare identifier `C` expands to `C == C` (drop rows whose
      # `C` column holds the literal label "C"). Mirrors ferx-core's parser
      # (src/io/filter_expr.rs); without this the preview returns no exclusions
      # while the Rust fit drops the flagged rows. Case is preserved on the RHS
      # so it matches the raw cell, with the column name lowercased for lookup.
      if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", part)) {
        result[[i]] <- list(col = tolower(part), op = "==", rhs = part)
      } else {
        return(NULL)
      }
    }
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

# The subject-ID column as text, one entry per record of `id_col`, or NULL
# when the data has no ID column. For a CSV the column is read again with no
# type conversion, after the same `[data]` column map, so `02` stays `02` as
# the engine's reader sees it, padding trimmed as the reader trims it. A
# data.frame's IDs were parsed before the
# preview saw them, so they are taken as they are.
.ferx_id_text <- function(source_path, column_map, id_col) {
  if (is.null(id_col)) return(NULL)
  if (!is.null(source_path)) {
    raw <- utils::read.csv(source_path, colClasses = "character",
                           check.names = FALSE)
    txt <- .get_col_sel(.ferx_apply_column_map(raw, column_map), "id")
    # The engine's reader trims every field before it matches or reports an
    # ID (Rust `str::trim`, Unicode whitespace), so ` 2` is subject `2`
    # (#795 review, round 2).
    if (!is.null(txt) && length(txt) == length(id_col)) {
      return(trimws(txt, whitespace = "[\\h\\v]"))
    }
  }
  as.character(id_col)
}

# `path` when it names an existing file, else NULL (a loaded fit's model_path
# can be NA or point at a file that has since gone).
.ferx_existing_path <- function(path) {
  if (is.character(path) && length(path) == 1L && !is.na(path) &&
      nzchar(path) && file.exists(path)) path else NULL
}
