# ETA-covariate correlation table: Pearson correlations between empirical
# Bayes estimates (ETAs) and covariates in the original dataset, identifying
# which covariates are most worth testing in a formal covariate search. Only
# columns that are constant within each subject are treated as covariates.
# Stored on the fit object as `fit$eta_cov`; NULL when the model declares no
# etas, the dataset has no usable numeric covariates, or the data file can no
# longer be read. (Formerly the exported ferx_eta_cov(fit, data), which took
# the dataset explicitly; see issue #226.)
#
# `ebe_etas` is purpose-built: ID + one column per BSV eta. `data_path` is
# re-read from disk here (rather than reusing an in-memory data.frame) since
# ferx_fit() only ever holds the CSV path, not a loaded copy of the dataset.
.ferx_compute_eta_cov <- function(ebe_etas, data_path) {
  if (is.null(ebe_etas) || !is.data.frame(ebe_etas)) return(NULL)
  data <- tryCatch(
    suppressWarnings(utils::read.csv(data_path, stringsAsFactors = FALSE,
                                      check.names = FALSE)),
    error = function(e) NULL
  )
  if (is.null(data) || !is.data.frame(data)) {
    # Only warn when a path was actually recorded but couldn't be read -
    # data_path being NA (no data bundled/resolvable) is the legitimate
    # "eta_cov not computable" case and shouldn't be noisy.
    if (!is.null(data_path) && length(data_path) == 1L && !is.na(data_path) &&
        nzchar(data_path)) {
      warning("fit$eta_cov: could not read data from '", data_path,
               "'; eta-covariate correlations were not computed.",
               call. = FALSE)
    }
    return(NULL)
  }

  # Treat every non-ID column as an eta - a "^ETA" prefix filter would
  # silently drop columns from models that don't follow the conventional
  # naming.
  ebe_id  <- if ("ID" %in% names(ebe_etas)) "ID" else names(ebe_etas)[1L]
  data_id <- if ("ID" %in% names(data))     "ID" else names(data)[1L]
  eta_cols <- setdiff(names(ebe_etas), ebe_id)
  if (length(eta_cols) == 0L) return(NULL)

  # Time-varying or non-covariate columns to skip
  SKIP <- c("TIME", "DV", "AMT", "EVID", "MDV", "CMT", "RATE",
            "II", "SS", "CENS", "LLOQ", "BLQ")

  # One row per subject already in ebe_etas
  etas <- ebe_etas[, c(ebe_id, eta_cols), drop = FALSE]

  # Numeric columns in data that could be covariates
  num_cols <- names(data)[vapply(data, is.numeric, logical(1L))]
  num_cols <- setdiff(num_cols, c(data_id, SKIP))

  if (length(num_cols) == 0L) return(NULL)

  # Keep only columns that are constant per subject (heuristic). Split on the
  # ordinal subject index, not the raw ID: a dataset that reuses a subject ID in
  # a non-contiguous block has two distinct subjects sharing that ID, and
  # splitting by ID would lump their records together (a covariate that differs
  # between the two would then be wrongly flagged non-constant, and one covariate
  # row would later be double-weighted across both subjects' ETAs).
  subj <- .ferx_subject_index(data[[data_id]])
  data_sub <- do.call(rbind, lapply(
    split(data[, c(data_id, num_cols), drop = FALSE], subj),
    function(chunk) {
      row <- chunk[1L, , drop = FALSE]
      for (col in num_cols) {
        if (length(unique(chunk[[col]])) > 1L) row[[col]] <- NA_real_
      }
      row
    }
  ))
  rownames(data_sub) <- NULL

  cov_cols <- num_cols[
    vapply(num_cols,
           function(col) sum(!is.na(data_sub[[col]])) > 0L,
           logical(1L))
  ]

  if (length(cov_cols) == 0L) return(NULL)

  # `etas` (one row per subject, subject order) and `data_sub` (one row per
  # subject index, same order) align positionally, so bind them directly - an
  # ID merge would cross-join the two subjects that share a reused ID. Fall back
  # to an ID merge only when the subject counts disagree.
  if (nrow(data_sub) == nrow(etas)) {
    merged <- cbind(etas, data_sub[, cov_cols, drop = FALSE])
  } else {
    # `data_sub` can hold two rows for one textual ID (a reused ID), so index
    # with match() rather than merge(): merge would cross-join those rows and
    # double-count the subject in every correlation.
    idx    <- match(as.character(etas[[ebe_id]]), as.character(data_sub[[data_id]]))
    merged <- cbind(etas, data_sub[idx, cov_cols, drop = FALSE])
    rownames(merged) <- NULL
  }

  rows <- vector("list", length(eta_cols) * length(cov_cols))
  k    <- 0L
  for (eta in eta_cols) {
    for (cov in cov_cols) {
      k <- k + 1L
      x  <- merged[[eta]]
      y  <- merged[[cov]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3L) {
        rows[[k]] <- data.frame(eta = eta, covariate = cov,
                                r = NA_real_, p_val = NA_real_, flag = "",
                                stringsAsFactors = FALSE)
        next
      }
      ct  <- suppressWarnings(cor.test(x[ok], y[ok]))
      r   <- as.numeric(ct$estimate)
      p   <- ct$p.value
      flg <- if (!is.na(r) && abs(r) >= 0.3) "[!]" else ""
      rows[[k]] <- data.frame(eta = eta, covariate = cov,
                              r = round(r, 3), p_val = round(p, 4),
                              flag = flg, stringsAsFactors = FALSE)
    }
  }

  result <- do.call(rbind, rows)
  result <- result[order(-abs(result$r), na.last = TRUE), ]
  rownames(result) <- NULL
  result
}
