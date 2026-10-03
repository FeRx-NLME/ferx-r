#' Build and validate a summary-level (MBMA) dataset
#'
#' Turns a table of published arm-level summaries - one row per study, arm and
#' timepoint, with the arm mean, its precision and the arm size - into the
#' record layout a model-based meta-analysis (MBMA) model is fitted on, and
#' checks it on the way. Every check below catches a data-preparation mistake
#' that otherwise fits without complaint: a percentage read as a proportion, a
#' standard deviation read as a standard error, an arm whose size changes
#' between visits, a study that contributes a single arm.
#'
#' Each study becomes one subject (\code{ID}), so a random effect on \code{ID}
#' is between-study variability, and each arm becomes one occasion
#' (\code{OCC}), so a \code{kappa} is between-treatment-arm variability. The
#' study is also written to a separate \code{STUDY} column: a
#' \code{theta NAME[COL, ...]} level block reads its columns from the
#' covariates, and \code{ID} is not one.
#'
#' The result is a data.frame. The fitting functions take a CSV path, so
#' write it out first:
#' \code{utils::write.csv(d, path, row.names = FALSE)}.
#'
#' @section Checks:
#' \describe{
#'   \item{Errors}{A \code{se} (or \code{sd}) that is missing, zero,
#'     negative or infinite where \code{mean} is present, an \code{n} that is
#'     missing, non-positive or infinite, an infinite \code{mean} or one
#'     outside the range of \code{scale}, a missing \code{study} or
#'     \code{arm}, a missing or infinite \code{time}, a repeated (study, arm,
#'     time) row, a missing value in a \code{levels} column, and an
#'     \code{index} column that is not 1..N without gaps. Each message names
#'     the offending rows.}
#'   \item{Warnings}{An arm whose \code{n} changes between timepoints (the
#'     per-row \code{n} is kept in \code{NARM}), a study with a single arm,
#'     and an arm without a \code{time == 0} record.}
#'   \item{Accepted}{A negative \code{mean}: change from baseline is negative
#'     by nature.}
#' }
#' Rows whose \code{mean} is missing are dropped, with a message; if that
#' leaves no row, the call is refused.
#'
#' @section Levels for a theta level block:
#' \code{levels = c("STUDY", "TIME")} adds a \code{LEVEL_IDX} column numbering
#' the observed (STUDY, TIME) combinations 1..N, in the order the engine
#' discovers the levels of \code{theta NAME[STUDY, TIME]}: ascending by the
#' first column, then the next. The counted form
#' \code{theta NAME[N]} with \code{NAME[LEVEL_IDX]} then estimates the same
#' levels as the data-driven form with \code{contrast = none}. The counted
#' form is what lets you define the cells yourself (for example, pooled visit
#' windows, via \code{index}). The table mapping each index to its cell is
#' kept in \code{attr(d, "levels")}, labelled as the engine labels a level
#' (\code{STUDY=7,TIME=4}). The levels are built from the values
#' \code{utils::write.csv()} writes (15 significant digits), which are what
#' the engine reads, and the level columns in the result hold those values:
#' 0.1 + 0.2 and 0.3 are one level, as they are in the written file.
#'
#' @param data A data.frame with one row per study, arm and timepoint.
#' @param study,arm,time,mean,n Names of the columns in \code{data} holding the
#'   study, the arm within the study, the timepoint, the arm mean and the arm
#'   size. \code{study} and \code{arm} may be numeric, character or factor;
#'   \code{time}, \code{mean} and \code{n} must be numeric.
#' @param se,sd Name of the column holding the standard error of the arm
#'   mean, or the between-subject standard deviation. Give exactly one;
#'   \code{sd} is converted to a standard error as \code{sd / sqrt(n)}.
#' @param scale The scale \code{mean} is reported on. \code{"continuous"}
#'   (default) accepts any value; \code{"proportion"} requires values in
#'   \code{[0, 1]}; \code{"percent"} requires values in \code{[0, 100]} and
#'   divides \code{mean} and \code{se} by 100.
#' @param covariates Optional names of further columns to carry through.
#'   Character and factor columns are coded as integers 1..K (factor levels
#'   order, or sorted values), with the codes kept in
#'   \code{attr(d, "codes")}; logical columns become 0/1.
#' @param levels Optional output column names - \code{"STUDY"}, \code{"TIME"}
#'   or a covariate - to build the \code{LEVEL_IDX} level index from.
#' @param index Optional name of a column in \code{data} holding your own
#'   1-based level index for a counted \code{theta NAME[N]} block. It is
#'   checked (whole numbers, 1..N, no level left out) and carried through
#'   under its own name. Not combined with \code{levels}.
#'
#' @return A \code{ferx_mbma_data} data.frame, sorted by study, arm and time,
#'   with columns \code{ID} (the study), \code{STUDY}, \code{OCC} (the arm,
#'   1..K within each study), \code{TIME}, \code{DV} (the arm mean),
#'   \code{MDV} (0), \code{NARM} (the arm size), \code{SE}, the covariates,
#'   and \code{LEVEL_IDX} or the \code{index} column. Attributes:
#'   \describe{
#'     \item{\code{arms}}{A data.frame mapping each (\code{STUDY},
#'       \code{OCC}) to the original \code{study} and \code{arm} values.}
#'     \item{\code{codes}}{A named list of the integer codes given to each
#'       character or factor column (\code{study} included).}
#'     \item{\code{levels}}{With \code{levels}: a data.frame of
#'       \code{index}, \code{label} and one column per level column.}
#'     \item{\code{scale}}{The \code{scale} argument.}
#'   }
#'
#' @examples
#' arms <- data.frame(
#'   trial = rep(c("A", "B"), each = 4),
#'   treatment = rep(rep(c("placebo", "drug"), each = 2), 2),
#'   week = rep(c(0, 12), 4),
#'   change = c(0, -1.1, 0, -2.4, 0, -0.8, 0, -2.0),
#'   sd = c(2.1, 2.3, 2.0, 2.4, 1.9, 2.2, 2.1, 2.5),
#'   n = c(120, 120, 118, 118, 95, 95, 97, 97)
#' )
#' d <- ferx_mbma_data(arms,
#'   study = "trial", arm = "treatment", time = "week",
#'   mean = "change", sd = "sd", n = "n",
#'   levels = c("STUDY", "TIME")
#' )
#' d
#' attr(d, "levels")
#' path <- tempfile(fileext = ".csv")
#' utils::write.csv(d, path, row.names = FALSE)
#'
#' @export
ferx_mbma_data <- function(data, study, arm, time, mean, se = NULL, sd = NULL,
                           n, scale = c("continuous", "proportion", "percent"),
                           covariates = NULL, levels = NULL, index = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  scale <- match.arg(scale)
  if (is.null(se) == is.null(sd)) {
    stop("Give exactly one of `se` and `sd`.", call. = FALSE)
  }
  if (!is.null(levels) && !is.null(index)) {
    stop("Give `levels` or `index`, not both.", call. = FALSE)
  }
  roles <- list(
    study = study, arm = arm, time = time, mean = mean, se = se, sd = sd,
    n = n, index = index
  )
  for (role in names(roles)) {
    col <- roles[[role]]
    if (!is.null(col) && !(is.character(col) && length(col) == 1L && !is.na(col))) {
      stop(sprintf("`%s` must be a single column name.", role), call. = FALSE)
    }
  }
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("`covariates` must be a character vector of column names.", call. = FALSE)
  }
  wanted <- c(unlist(roles, use.names = FALSE), covariates)
  missing_cols <- setdiff(wanted, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "`data` has no column(s) %s.", .mbma_list(missing_cols, quote = TRUE)
    ), call. = FALSE)
  }
  for (role in c("time", "mean", "n", "se", "sd")) {
    col <- roles[[role]]
    if (!is.null(col) && !is.numeric(data[[col]])) {
      stop(sprintf("Column `%s` (`%s`) must be numeric.", col, role), call. = FALSE)
    }
  }

  key_na <- is.na(data[[study]]) | is.na(data[[arm]]) | !is.finite(data[[time]])
  if (any(key_na)) {
    stop(sprintf(
      "`study` and `arm` must be present and `time` finite; row(s) %s.",
      .mbma_list(which(key_na))
    ), call. = FALSE)
  }

  absent <- is.na(data[[mean]])
  if (any(absent)) {
    message(sprintf(
      "Dropped %d row(s) with a missing `mean`.", sum(absent)
    ))
    data <- data[!absent, , drop = FALSE]
  }
  if (nrow(data) == 0L) {
    stop("No row has a `mean`.", call. = FALSE)
  }
  cell <- .mbma_cells(data[[study]], data[[arm]], data[[time]])

  n_val <- data[[n]]
  bad_n <- !is.finite(n_val) | n_val <= 0
  if (any(bad_n)) {
    stop(sprintf(
      "`n` must be positive and finite on every row: %s.",
      .mbma_list(cell[bad_n])
    ), call. = FALSE)
  }

  if (is.null(se)) {
    prec <- data[[sd]]
    prec_role <- "sd"
  } else {
    prec <- data[[se]]
    prec_role <- "se"
  }
  bad_prec <- !is.finite(prec) | prec <= 0
  if (any(bad_prec)) {
    stop(sprintf(
      "`%s` must be positive and finite wherever `mean` is present: %s.",
      prec_role, .mbma_list(cell[bad_prec])
    ), call. = FALSE)
  }
  if (is.null(se)) {
    prec <- prec / sqrt(n_val)
    message("Converted `sd` to a standard error as sd / sqrt(n).")
  }

  dv <- data[[mean]]
  bad_mean <- !is.finite(dv)
  if (any(bad_mean)) {
    stop(sprintf(
      "`mean` must be finite: %s.", .mbma_list(cell[bad_mean])
    ), call. = FALSE)
  }
  if (scale == "proportion") {
    out <- dv < 0 | dv > 1
    if (any(out)) {
      hint <- if (min(dv) >= 0 && max(dv) <= 100) {
        " The values look like percentages; use scale = \"percent\"."
      } else {
        ""
      }
      stop(sprintf(
        "`mean` is outside [0, 1] for scale = \"proportion\": %s.%s",
        .mbma_list(cell[out]), hint
      ), call. = FALSE)
    }
  } else if (scale == "percent") {
    out <- dv < 0 | dv > 100
    if (any(out)) {
      stop(sprintf(
        "`mean` is outside [0, 100] for scale = \"percent\": %s.",
        .mbma_list(cell[out])
      ), call. = FALSE)
    }
    dv <- dv / 100
    prec <- prec / 100
    message("Converted `mean` and `se` from percent to proportions (divided by 100).")
  }

  dup <- duplicated(cell)
  if (any(dup)) {
    stop(sprintf(
      "Each (study, arm, time) must appear once; repeated: %s.",
      .mbma_list(unique(cell[dup]))
    ), call. = FALSE)
  }

  codes <- list()
  study_code <- .mbma_code(data[[study]])
  if (!is.null(study_code$codes)) codes$study <- study_code$codes
  study_val <- study_code$values
  arm_val <- .mbma_code(data[[arm]])$values
  arm_key <- paste(study_val, arm_val, sep = "\r")
  occ <- as.integer(stats::ave(
    arm_val, study_val,
    FUN = function(a) match(a, sort(unique(a)))
  ))
  arm_label <- sprintf("study %s arm %s", data[[study]], data[[arm]])

  n_varies <- tapply(n_val, arm_key, function(v) length(unique(v)) > 1L)
  if (any(n_varies)) {
    warning(sprintf(
      "`n` changes between timepoints within %s; each row keeps its own n in NARM.",
      .mbma_list(unique(arm_label[arm_key %in% names(n_varies)[n_varies]]))
    ), call. = FALSE)
  }
  n_arms <- tapply(arm_key, study_val, function(a) length(unique(a)))
  if (any(n_arms == 1L)) {
    single <- unique(data[[study]][study_val %in% names(n_arms)[n_arms == 1L]])
    warning(sprintf(
      "Study(ies) %s contribute a single arm, so no comparison within them informs a treatment effect.",
      .mbma_list(single)
    ), call. = FALSE)
  }
  has_base <- tapply(data[[time]] == 0, arm_key, any)
  if (!all(has_base)) {
    warning(sprintf(
      "No time = 0 record for %s. Expected for change-from-baseline data; otherwise check the time axis.",
      .mbma_list(unique(arm_label[arm_key %in% names(has_base)[!has_base]]))
    ), call. = FALSE)
  }

  out_df <- data.frame(
    ID = study_val, STUDY = study_val, OCC = occ, TIME = data[[time]],
    DV = dv, MDV = 0L, NARM = n_val, SE = prec
  )
  reserved <- c(names(out_df), "LEVEL_IDX")
  if (!is.null(index) && index %in% reserved) {
    stop(sprintf(
      "`index` column `%s` would overwrite an output column; rename it first.", index
    ), call. = FALSE)
  }
  clash <- intersect(covariates, c(reserved, index))
  if (length(clash) > 0L) {
    stop(sprintf(
      "Covariate column(s) %s would overwrite an output column; rename them first.",
      .mbma_list(clash, quote = TRUE)
    ), call. = FALSE)
  }
  for (cov in covariates) {
    coded <- .mbma_code(data[[cov]])
    out_df[[cov]] <- coded$values
    if (!is.null(coded$codes)) codes[[cov]] <- coded$codes
  }

  if (!is.null(index)) {
    out_df[[index]] <- .mbma_check_index(data[[index]], index, cell)
  }

  if (!is.null(levels)) {
    allowed <- c("STUDY", "TIME", covariates)
    if (!is.character(levels) || length(levels) == 0L ||
      !all(levels %in% allowed) || anyDuplicated(levels)) {
      stop(sprintf(
        "`levels` must name distinct output columns among %s.",
        .mbma_list(allowed, quote = TRUE)
      ), call. = FALSE)
    }
    for (col in levels) {
      bad <- !is.finite(out_df[[col]])
      if (any(bad)) {
        stop(sprintf(
          "Level column `%s` must be finite on every row: %s.",
          col, .mbma_list(cell[bad])
        ), call. = FALSE)
      }
    }
  }

  ord <- order(out_df$STUDY, out_df$OCC, out_df$TIME)
  out_df <- out_df[ord, , drop = FALSE]
  rownames(out_df) <- NULL

  level_table <- NULL
  if (!is.null(levels)) {
    built <- .mbma_level_index(out_df[levels])
    out_df[levels] <- built$values
    out_df$ID <- out_df$STUDY
    out_df$LEVEL_IDX <- built$index
    level_table <- built$table
  }

  # From the sorted frame, so a normalised STUDY (levels =) matches it.
  arms <- unique(data.frame(
    STUDY = out_df$STUDY, OCC = out_df$OCC,
    study = data[[study]][ord], arm = data[[arm]][ord]
  ))
  arms <- arms[order(arms$STUDY, arms$OCC), , drop = FALSE]
  rownames(arms) <- NULL

  structure(
    out_df,
    class = c("ferx_mbma_data", "data.frame"),
    arms = arms,
    codes = codes,
    levels = level_table,
    scale = scale
  )
}

#' Print a ferx_mbma_data object
#'
#' Reports what the checks in \code{\link{ferx_mbma_data}} cannot decide for
#' you: how many studies and arms there are, how many arms each study and
#' how many timepoints each arm contributes, and the range of the
#' between-subject SD the standard errors imply (\code{SE * sqrt(NARM)}). An
#' implied SD far from what the endpoint plausibly has usually means a
#' standard deviation was read as a standard error, or the reverse.
#'
#' @param x A \code{ferx_mbma_data} object.
#' @param n Number of records to show.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.ferx_mbma_data <- function(x, n = 6L, ...) {
  arm_key <- paste(x$STUDY, x$OCC, sep = "\r")
  per_study <- tapply(arm_key, x$STUDY, function(a) length(unique(a)))
  per_arm <- tabulate(match(arm_key, unique(arm_key)))
  implied <- x$SE * sqrt(x$NARM)
  cat(sprintf(
    "<ferx_mbma_data>  %d records, %d studies, %d arms\n",
    nrow(x), length(per_study), length(unique(arm_key))
  ))
  cat(sprintf("  arms per study:      %s\n", .mbma_range(per_study)))
  cat(sprintf("  timepoints per arm:  %s\n", .mbma_range(per_arm)))
  cat(sprintf(
    "  implied SD (SE * sqrt(NARM)): %s to %s\n",
    format(min(implied), digits = 3), format(max(implied), digits = 3)
  ))
  lv <- attr(x, "levels")
  if (!is.null(lv) && "LEVEL_IDX" %in% names(x)) {
    cols <- setdiff(names(lv), c("index", "label"))
    cat(sprintf(
      "  LEVEL_IDX: %d levels of (%s); declare theta NAME[%d]\n",
      nrow(lv), paste(cols, collapse = ", "), nrow(lv)
    ))
  }
  cat("\n")
  shown <- x[seq_len(min(n, nrow(x))), , drop = FALSE]
  class(shown) <- "data.frame"
  print(shown)
  if (nrow(x) > n) cat(sprintf("... %d more records\n", nrow(x) - n))
  invisible(x)
}

# -- helpers ------------------------------------------------------------------

# "(study A, arm drug, time 12)" for every row, used to name offending rows.
.mbma_cells <- function(study, arm, time) {
  sprintf("(study %s, arm %s, time %s)", study, arm, as.character(time))
}

# Comma-separated listing, the first 10 items and a count of the rest.
.mbma_list <- function(x, quote = FALSE) {
  x <- as.character(x)
  if (quote) x <- sprintf("`%s`", x)
  if (length(x) > 10L) {
    x <- c(x[1:10], sprintf("and %d more", length(x) - 10L))
  }
  paste(x, collapse = ", ")
}

.mbma_range <- function(x) {
  if (min(x) == max(x)) as.character(min(x)) else sprintf("%d-%d", min(x), max(x))
}

# Numeric values pass through; logical becomes 0/1; character and factor
# become 1..K (factor levels order, or sorted values), with the code table.
.mbma_code <- function(x) {
  if (is.numeric(x)) {
    return(list(values = x, codes = NULL))
  }
  if (is.logical(x)) {
    return(list(values = as.integer(x), codes = NULL))
  }
  # Radix sort: C-locale order, so the codes do not depend on the collation
  # locale of the machine that built the frame.
  lv <- if (is.factor(x)) {
    base::levels(droplevels(x))
  } else {
    sort(unique(as.character(x)), method = "radix")
  }
  values <- match(as.character(x), lv)
  list(values = values, codes = stats::setNames(seq_along(lv), lv))
}

# A user-supplied counted-form index: whole numbers 1..N, every level used.
.mbma_check_index <- function(idx, col, cell) {
  if (!is.numeric(idx)) {
    stop(sprintf("`index` column `%s` must be numeric.", col), call. = FALSE)
  }
  bad <- !is.finite(idx) | idx != round(idx) | idx < 1
  if (any(bad)) {
    stop(sprintf(
      "`index` column `%s` must hold whole numbers from 1; offending value(s) %s at %s.",
      col, .mbma_list(unique(idx[bad])), .mbma_list(cell[bad])
    ), call. = FALSE)
  }
  gaps <- setdiff(seq_len(max(idx)), idx)
  if (length(gaps) > 0L) {
    stop(sprintf(
      "`index` column `%s` skips level(s) %s of 1..%d; a counted block's levels must all appear.",
      col, .mbma_list(gaps), as.integer(max(idx))
    ), call. = FALSE)
  }
  as.integer(idx)
}

# Number the observed combinations of `cols` 1..N as the engine discovers a
# level block's levels, and label each as the engine does. The engine reads
# the CSV, so a value is what utils::write.csv() writes for it (15
# significant digits): levels are matched by that text, sorted ascending by
# its value column by column, and the columns come back holding those values.
.mbma_level_index <- function(cols) {
  text <- lapply(cols, .mbma_csv_text)
  values <- lapply(text, as.numeric)
  key <- do.call(paste, c(unname(text), sep = "\r"))
  first <- which(!duplicated(key))
  first <- first[do.call(order, unname(lapply(values, `[`, first)))]
  label <- do.call(paste, c(
    lapply(names(cols), function(nm) {
      paste0(nm, "=", .mbma_level_text(text[[nm]][first]))
    }),
    sep = ","
  ))
  cells <- as.data.frame(lapply(values, `[`, first))
  list(
    index = match(key, key[first]),
    values = values,
    table = data.frame(index = seq_along(first), label = label, cells)
  )
}

# The text utils::write.csv() writes for a numeric vector.
.mbma_csv_text <- function(v) {
  con <- textConnection("out", "w", local = TRUE)
  utils::write.table(
    data.frame(v), con, sep = ",", col.names = FALSE, row.names = FALSE
  )
  close(con)
  out
}

# A level value as the engine labels it: the shortest decimal that reads
# back as the value, never in exponent form. A decimal of 15 significant
# digits or fewer reads back exactly, so the written text is that decimal;
# only an exponent needs rewriting.
.mbma_level_text <- function(s) {
  sci <- grepl("e", s, fixed = TRUE)
  s[sci] <- vapply(
    as.numeric(s[sci]), format, character(1),
    scientific = FALSE, digits = 15, trim = TRUE
  )
  s
}
