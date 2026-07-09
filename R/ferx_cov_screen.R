#' Quick covariate screen against individual parameters and ETAs
#'
#' A fast, informal covariate screen built on the fit's covariate table
#' (\code{fit$covtab}). For every covariate it reports the association with each
#' inter-individually-varying parameter, measured two ways: against the
#' subject's individual parameter estimate (the "EBE" column, from
#' \code{fit$individual_estimates}) and against the parameter's ETA (the "ETA"
#' column, from \code{fit$ebe_etas}). Only parameters that have IIV (an ETA) are
#' screened. This is a screening aid to decide what is worth a formal covariate
#' search - it is not itself a covariate test.
#'
#' Covariates are aggregated to one value per subject first: the median for
#' continuous covariates (so time-varying covariates collapse to a typical
#' value) and the most frequent level for categorical covariates.
#'
#' The association measure depends on the covariate type (from
#' \code{fit$covariate_types}): a signed Pearson correlation for continuous
#' covariates, and the correlation ratio (eta, in \code{[0, 1]}) for categorical
#' covariates. Rows are kept only when either measure is at least
#' \code{threshold} in absolute value, and are ordered by the stronger of the
#' two.
#'
#' @param fit A \code{ferx_fit} object with a \code{covtab} (the model must
#'   declare a \code{[covariates]} block).
#' @param threshold Minimum absolute association to report. Default 0.2.
#' @return Data frame with columns \code{parameter}, \code{covariate},
#'   \code{type}, \code{ebe}, and \code{eta}, ordered by descending association
#'   strength. Returned invisibly; the table is also printed. Returns an empty
#'   data frame (with a message) when nothing clears the threshold, and
#'   \code{NULL} when the fit has no covariate table or no ETAs.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("two_cpt_oral_cov")
#' fit <- ferx_fit(ex$model, ex$data, method = "focei", covariance = FALSE)
#' ferx_cov_screen(fit)
#' }
#' @family diagnostics
#' @export
ferx_cov_screen <- function(fit, threshold = 0.2) {
  if (is.null(fit$covtab) || !is.data.frame(fit$covtab)) {
    message("No covariate table (`fit$covtab`) - the model needs a ",
            "[covariates] block.")
    return(invisible(NULL))
  }
  if (is.null(fit$ebe_etas) || !is.data.frame(fit$ebe_etas)) {
    message("No ETAs (`fit$ebe_etas`) - the model declares no ",
            "inter-individual variability.")
    return(invisible(NULL))
  }

  covtab <- fit$covtab
  id_col <- if ("ID" %in% names(covtab)) "ID" else names(covtab)[1L]

  # Declared covariate names + types. Fall back to inferring from the covtab
  # columns (everything except the structural ID/TIME/EVID columns) when
  # covariate_types is absent.
  types <- fit$covariate_types
  if (!is.null(types) && length(types)) {
    cov_names <- names(types)
  } else {
    cov_names <- setdiff(names(covtab), c(id_col, "TIME", "EVID"))
    types <- stats::setNames(
      ifelse(vapply(covtab[cov_names], is.numeric, logical(1L)),
             "continuous", "categorical"),
      cov_names
    )
  }
  cov_names <- intersect(cov_names, names(covtab))
  if (!length(cov_names)) {
    message("No declared covariates found in `fit$covtab`.")
    return(invisible(NULL))
  }

  # Aggregate each covariate to one value per subject: median (continuous) or
  # the most frequent level (categorical).
  # Aggregate on the ordinal subject index, not the raw ID: a dataset that
  # reuses a subject ID in a non-contiguous block has two distinct subjects
  # sharing that ID, and keying on ID would collapse them into one row (and
  # later double-weight that row against both subjects' ETAs).
  subj <- .ferx_subject_index(covtab[[id_col]])
  ids  <- covtab[[id_col]][!duplicated(subj)]  # one ID per subject, subject order
  keys <- as.character(seq_along(ids))
  mode_fn <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NA)
    ux <- unique(v)
    ux[which.max(tabulate(match(v, ux)))]
  }
  percov <- stats::setNames(
    data.frame(as.character(ids), stringsAsFactors = FALSE), id_col
  )
  for (cov in cov_names) {
    if (identical(types[[cov]], "categorical")) {
      agg <- tapply(covtab[[cov]], subj, mode_fn)
    } else {
      agg <- tapply(covtab[[cov]], subj,
                    function(v) stats::median(v, na.rm = TRUE))
    }
    percov[[cov]] <- agg[keys]
  }

  # Subject-level ETAs and individual parameter estimates.
  etas   <- fit$ebe_etas
  eta_id <- if ("ID" %in% names(etas)) "ID" else names(etas)[1L]
  eta_cols <- setdiff(names(etas), eta_id)
  if (!length(eta_cols)) {
    message("No ETA columns found in `fit$ebe_etas`.")
    return(invisible(NULL))
  }
  etas[[eta_id]] <- as.character(etas[[eta_id]])

  ipar     <- fit$individual_estimates
  has_ipar <- is.data.frame(ipar) && nrow(ipar) > 0L
  if (has_ipar) {
    ipar_id    <- if ("ID" %in% names(ipar)) "ID" else names(ipar)[1L]
    param_cols <- setdiff(names(ipar), ipar_id)
  } else {
    param_cols <- character(0)
  }

  # Map an ETA to its individual parameter by stripping an ETA prefix/suffix
  # (e.g. ETA_CL -> CL); NA when there is no matching parameter column.
  map_param <- function(eta) {
    cand <- sub("^ETA[_.]?", "", eta, ignore.case = TRUE)
    cand <- sub("[_.]?ETA$", "", cand, ignore.case = TRUE)
    hit  <- param_cols[toupper(param_cols) == toupper(cand)]
    if (length(hit) == 1L) hit else NA_character_
  }

  # Association: signed Pearson r (continuous) or correlation ratio (eta,
  # categorical). Needs at least 3 usable pairs and non-zero variance.
  assoc <- function(response, x, type) {
    ok <- is.finite(response) & !is.na(x)
    if (sum(ok) < 3L) return(NA_real_)
    response <- response[ok]
    x <- x[ok]
    if (identical(type, "categorical")) {
      f <- factor(x)
      if (nlevels(f) < 2L) return(NA_real_)
      grand  <- mean(response)
      ss_tot <- sum((response - grand)^2)
      if (ss_tot <= 0) return(NA_real_)
      ss_bet <- sum(vapply(split(response, f),
                           function(g) length(g) * (mean(g) - grand)^2,
                           numeric(1L)))
      sqrt(ss_bet / ss_tot)
    } else {
      x <- as.numeric(x)
      if (stats::sd(x) == 0 || stats::sd(response) == 0) return(NA_real_)
      suppressWarnings(stats::cor(response, x))
    }
  }

  # `percov`, `etas` and `ipar` are each one row per subject in subject order,
  # so bind/index positionally - an ID merge would cross-join two subjects that
  # share a reused ID. Fall back to an ID match when the counts disagree.
  if (nrow(percov) == nrow(etas)) {
    merged   <- cbind(percov, etas[, eta_cols, drop = FALSE])
    ipar_idx <- if (has_ipar) {
      if (nrow(ipar) == nrow(percov)) seq_len(nrow(percov))
      else match(merged[[id_col]], as.character(ipar[[ipar_id]]))
    } else {
      NULL
    }
  } else {
    merged   <- merge(percov, etas, by.x = id_col, by.y = eta_id)
    ipar_idx <- if (has_ipar) {
      match(merged[[id_col]], as.character(ipar[[ipar_id]]))
    } else {
      NULL
    }
  }

  rows <- vector("list", length(eta_cols) * length(cov_names))
  k    <- 0L
  for (eta in eta_cols) {
    param  <- map_param(eta)
    plabel <- if (!is.na(param)) {
      param
    } else {
      sub("^ETA[_.]?", "", eta, ignore.case = TRUE)
    }
    for (cov in cov_names) {
      k <- k + 1L
      type  <- types[[cov]]
      r_eta <- assoc(merged[[eta]], merged[[cov]], type)
      r_ebe <- if (!is.na(param) && param %in% param_cols) {
        assoc(ipar[[param]][ipar_idx], merged[[cov]], type)
      } else {
        NA_real_
      }
      rows[[k]] <- data.frame(
        parameter = plabel, covariate = cov, type = type,
        ebe = round(r_ebe, 3), eta = round(r_eta, 3),
        stringsAsFactors = FALSE
      )
    }
  }
  result   <- do.call(rbind, rows)
  strength <- pmax(abs(result$ebe), abs(result$eta), na.rm = TRUE)
  result   <- result[!is.na(strength) & strength >= threshold, , drop = FALSE]
  result   <- result[order(-pmax(abs(result$ebe), abs(result$eta),
                                 na.rm = TRUE)), , drop = FALSE]
  rownames(result) <- NULL

  if (nrow(result) == 0L) {
    message(sprintf("No covariate correlations with |r| >= %.2f.", threshold))
    return(invisible(result))
  }
  print(result)
  invisible(result)
}
