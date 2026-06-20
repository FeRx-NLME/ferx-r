.dw_label <- function(dw) {
  if (dw < 1.5) "positive autocorrelation"
  else if (dw > 2.5) "negative autocorrelation"
  else "no autocorrelation"
}

#' Structured diagnostic flags from a fit result
#'
#' Returns a named list of diagnostic data frames covering IWRES autocorrelation
#' and parameter shrinkage. Designed to be inspected programmatically or used
#' as the basis for custom plots.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return A named list with elements:
#'   \describe{
#'     \item{autocorrelation}{1-row data frame: \code{dw_statistic},
#'       \code{lag1_r}, \code{flag} (character interpretation of the DW value).
#'       \code{NULL} when no IWRES autocorrelation data is available.}
#'     \item{shrinkage}{Data frame with columns \code{param}, \code{type}
#'       (\code{"eta"} or \code{"eps"}), \code{shrinkage} (proportion 0-1),
#'       \code{shrinkage_pct} (percentage). \code{NULL} when no shrinkage data
#'       is available.}
#'   }
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' diag <- check_diagnostics(fit)
#' diag$autocorrelation
#' diag$shrinkage
#' }
#' @family diagnostics
#' @export
check_diagnostics <- function(fit) {
  # Autocorrelation
  dw  <- fit$dw_statistic
  r   <- fit$iwres_lag1_r
  autocorr <- if (!is.null(dw) && !is.na(dw)) {
    data.frame(dw_statistic = dw, lag1_r = r %||% NA_real_,
               flag = .dw_label(dw),
               stringsAsFactors = FALSE)
  } else NULL

  # Shrinkage
  rows <- list()
  if (!is.null(fit$shrinkage_eta)) {
    eta_lbls <- if (!is.null(fit$eta_names) && length(fit$eta_names) == length(fit$shrinkage_eta))
      fit$eta_names else sprintf("ETA%d", seq_along(fit$shrinkage_eta))
    for (i in seq_along(fit$shrinkage_eta)) {
      sh <- fit$shrinkage_eta[i]
      if (!is.na(sh)) {
        rows[[length(rows) + 1L]] <- data.frame(
          param = eta_lbls[i], type = "eta",
          shrinkage = sh, shrinkage_pct = sh * 100,
          stringsAsFactors = FALSE)
      }
    }
  }
  eps_vec <- fit$shrinkage_eps
  if (!is.null(eps_vec)) {
    eps_vec <- eps_vec[!is.na(eps_vec)]
    if (length(eps_vec) > 0L) {
      eps_lbls <- if (!is.null(fit$sigma_names) && length(fit$sigma_names) == length(fit$shrinkage_eps))
        fit$sigma_names[!is.na(fit$shrinkage_eps)]
      else if (length(eps_vec) == 1L)
        "EPS"
      else
        sprintf("EPS%d", seq_along(eps_vec))
      for (i in seq_along(eps_vec)) {
        rows[[length(rows) + 1L]] <- data.frame(
          param = eps_lbls[i], type = "eps",
          shrinkage = eps_vec[i], shrinkage_pct = eps_vec[i] * 100,
          stringsAsFactors = FALSE)
      }
    }
  }
  shrinkage_df <- if (length(rows) > 0L) do.call(rbind, rows) else NULL

  list(autocorrelation = autocorr, shrinkage = shrinkage_df)
}

#' Correlation matrix of estimated parameters
#'
#' Converts the parameter covariance matrix (already stored on the fit object)
#' into a correlation matrix, making off-diagonal structure immediately visible.
#' A correlation close to \eqn{\pm 1} between two parameters flags a structural
#' identifiability problem in the model.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return Named correlation matrix (invisibly). Printed to the console.
#' @seealso \code{\link{ferx_estimates}} for SEs and \%RSE.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = TRUE)
#' if (!is.null(fit$cov_matrix)) ferx_cor_matrix(fit)
#' @family diagnostics
#' @export
ferx_cor_matrix <- function(fit) {
  if (is.null(fit$cov_matrix)) {
    stop(
      "Covariance matrix not available. ",
      "Run ferx_fit() with covariance = TRUE and check fit$covariance_status."
    )
  }
  se <- sqrt(diag(fit$cov_matrix))
  if (any(se <= 0, na.rm = TRUE)) {
    warning("One or more diagonal elements are non-positive; ",
            "correlation matrix may not be meaningful.")
    se[se <= 0] <- NA_real_
  }
  cor_mat <- fit$cov_matrix / outer(se, se)
  # Clip to [-1, 1] for numerical noise on the diagonal
  cor_mat[cor_mat >  1] <-  1
  cor_mat[cor_mat < -1] <- -1
  # A correlation matrix has an exactly-unit diagonal by definition. Force it:
  # cov[i,i] / sqrt(cov[i,i])^2 can drift to 0.9999999999999999, which the `> 1`
  # clip above does not catch. Leave NA where a non-positive variance was
  # flagged above (those rows/cols are already NA and warned).
  ok <- !is.na(diag(cor_mat))
  diag(cor_mat)[ok] <- 1
  print(round(cor_mat, 3))
  invisible(cor_mat)
}

#' ETA-covariate correlation table
#'
#' Computes Pearson correlations between empirical Bayes estimates (ETAs) and
#' covariates in the original dataset. Identifies which covariates are most
#' worth testing in a formal covariate search. Only columns that are constant
#' within each subject are treated as covariates.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @param data The original dataset (data frame) passed to
#'   \code{\link{ferx_fit}}.
#' @return Data frame with columns \code{eta}, \code{covariate}, \code{r},
#'   \code{p_val}, \code{flag}, sorted by descending \code{|r|}. Returned
#'   invisibly; the full table is printed to the console.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' obs <- read.csv(ex$data)
#' ferx_eta_cov(fit, obs)
#' @family diagnostics
#' @export
ferx_eta_cov <- function(fit, data) {
  if (is.null(fit$ebe_etas) || !is.data.frame(fit$ebe_etas)) {
    stop("`fit$ebe_etas` is not available.")
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a data.frame - pass the dataset used in ferx_fit().")
  }

  # ebe_etas is purpose-built: ID + one column per BSV eta. Treat every
  # non-ID column as an eta - a "^ETA" prefix filter would silently drop
  # columns from models that don't follow the conventional naming.
  ebe_id  <- if ("ID" %in% names(fit$ebe_etas)) "ID" else names(fit$ebe_etas)[1L]
  data_id <- if ("ID" %in% names(data))         "ID" else names(data)[1L]
  eta_cols <- setdiff(names(fit$ebe_etas), ebe_id)
  if (length(eta_cols) == 0L) {
    message("No ETA columns found in fit$ebe_etas.")
    return(invisible(NULL))
  }

  # Time-varying or non-covariate columns to skip
  SKIP <- c("TIME", "DV", "AMT", "EVID", "MDV", "CMT", "RATE",
            "II", "SS", "CENS", "LLOQ", "BLQ")

  # One row per subject already in ebe_etas
  etas <- fit$ebe_etas[, c(ebe_id, eta_cols), drop = FALSE]

  # Numeric columns in data that could be covariates
  num_cols <- names(data)[vapply(data, is.numeric, logical(1L))]
  num_cols <- setdiff(num_cols, c(data_id, SKIP))

  if (length(num_cols) == 0L) {
    message("No numeric covariate columns found in data.")
    return(invisible(NULL))
  }

  # Keep only columns that are constant per subject (heuristic)
  data_sub <- do.call(rbind, lapply(
    split(data[, c(data_id, num_cols), drop = FALSE], data[[data_id]]),
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

  if (length(cov_cols) == 0L) {
    message("No constant-per-subject numeric covariates found in data.")
    return(invisible(NULL))
  }

  merged <- merge(etas, data_sub[, c(data_id, cov_cols), drop = FALSE],
                  by.x = ebe_id, by.y = data_id)

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
  print(result)
  invisible(result)
}

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
  ids <- unique(covtab[[id_col]])
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
      agg <- tapply(covtab[[cov]], covtab[[id_col]], mode_fn)
    } else {
      agg <- tapply(covtab[[cov]], covtab[[id_col]],
                    function(v) stats::median(v, na.rm = TRUE))
    }
    percov[[cov]] <- agg[as.character(ids)]
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

  merged   <- merge(percov, etas, by.x = id_col, by.y = eta_id)
  ipar_idx <- if (has_ipar) {
    match(merged[[id_col]], as.character(ipar[[ipar_id]]))
  } else {
    NULL
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

#' Tidy parameter estimates table
#'
#' Extracts all estimated parameters (theta, omega diagonal, sigma) into a
#' single tidy data frame, adding percent relative standard error (\%RSE),
#' 95\% confidence intervals, and-for log/logit-transformed thetas-natural-scale
#' back-transformed estimates and CIs.
#'
#' Omega is reported on the variance scale (matching the \code{.ferx} model
#' file convention). For block omega, only the diagonal variances are included.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return A data frame with columns \code{param}, \code{transform},
#'   \code{estimate}, \code{se}, \code{rse_pct}, \code{lower_95},
#'   \code{upper_95}, \code{estimate_natural}, \code{lower_95_natural},
#'   \code{upper_95_natural}, \code{init_as_sd}. SE-derived and
#'   natural-scale columns are \code{NA} when not applicable or when the
#'   covariance step was not run. \code{init_as_sd} is \code{TRUE} for
#'   omega, sigma, or kappa rows where the user annotated the initial
#'   value with \code{(sd)} in the model file; always \code{FALSE} for
#'   theta rows.
#' @seealso \code{\link{ferx_cor_matrix}} for parameter correlations.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' ferx_estimates(fit)
#' @family diagnostics
#' @export
ferx_estimates <- function(fit) {
  rows <- list()

  # Theta
  theta_names <- names(fit$theta)
  if (is.null(theta_names)) theta_names <- paste0("THETA", seq_along(fit$theta))
  for (i in seq_along(fit$theta)) {
    se        <- if (!is.null(fit$se_theta) && length(fit$se_theta) >= i) fit$se_theta[i] else NA_real_
    transform <- if (!is.null(fit$theta_transforms) && length(fit$theta_transforms) >= i) fit$theta_transforms[i] else "identity"
    rows[[length(rows) + 1L]] <- .ferx_est_row(theta_names[i], fit$theta[i], se, transform, FALSE)
  }

  # Omega diagonal (variance scale)
  om    <- fit$omega
  if (is.null(dim(om))) om <- matrix(om, 1L, 1L)
  n_eta <- nrow(om)
  for (i in seq_len(n_eta)) {
    pname    <- if (!is.null(fit$eta_names) && length(fit$eta_names) >= i && nzchar(fit$eta_names[i])) fit$eta_names[i] else sprintf("OMEGA(%d,%d)", i, i)
    se       <- .omega_se_at(fit$se_omega, n_eta, i, i)
    init_sd  <- !is.null(fit$omega_init_as_sd) && length(fit$omega_init_as_sd) >= i && isTRUE(fit$omega_init_as_sd[i])
    rows[[length(rows) + 1L]] <- .ferx_est_row(pname, om[i, i], se, "variance", init_sd)
  }

  # Sigma
  for (i in seq_along(fit$sigma)) {
    pname         <- if (!is.null(fit$sigma_names) && length(fit$sigma_names) >= i && nzchar(fit$sigma_names[i])) fit$sigma_names[i] else sprintf("SIGMA(%d)", i)
    se            <- if (!is.null(fit$se_sigma) && length(fit$se_sigma) >= i) fit$se_sigma[i] else NA_real_
    sig_transform <- if (!is.null(fit$sigma_types) && length(fit$sigma_types) >= i) fit$sigma_types[i] else "proportional"
    init_sd       <- !is.null(fit$sigma_init_as_sd) && length(fit$sigma_init_as_sd) >= i && isTRUE(fit$sigma_init_as_sd[i])
    rows[[length(rows) + 1L]] <- .ferx_est_row(pname, fit$sigma[i], se, sig_transform, init_sd)
  }

  # Kappa (IOV diagonal)
  if (!is.null(fit$omega_iov)) {
    m_iov  <- fit$omega_iov
    if (is.null(dim(m_iov))) m_iov <- matrix(m_iov, 1L, 1L)
    n_kap  <- nrow(m_iov)
    kap_names <- if (!is.null(fit$kappa_names) && length(fit$kappa_names) == n_kap) fit$kappa_names else paste0("KAPPA", seq_len(n_kap))
    n_se   <- length(fit$se_kappa)
    n_tri  <- n_kap * (n_kap + 1L) / 2L
    is_block_se <- (n_se == n_tri && n_kap > 1L)
    diag_se_idx <- function(j) j * n_kap - j * (j - 1L) / 2L - (n_kap - j)
    for (i in seq_len(n_kap)) {
      se_idx  <- if (is_block_se) diag_se_idx(i) else i
      se      <- if (n_se >= se_idx) fit$se_kappa[se_idx] else NA_real_
      kap_type <- if (!is.null(fit$kappa_param_types) && length(fit$kappa_param_types) >= i) fit$kappa_param_types[i] else "variance"
      init_sd  <- !is.null(fit$kappa_init_as_sd) && length(fit$kappa_init_as_sd) >= i && isTRUE(fit$kappa_init_as_sd[i])
      rows[[length(rows) + 1L]] <- .ferx_est_row(kap_names[i], m_iov[i, i], se, kap_type, init_sd)
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

.ferx_inv_logit <- function(x) 1 / (1 + exp(-x))

.ferx_est_row <- function(param, estimate, se, transform = "identity", init_as_sd = FALSE) {
  rse_pct  <- if (!is.na(se) && abs(estimate) > 1e-12) abs(se / estimate) * 100 else NA_real_

  # Asymmetric CI and natural-scale back-transform per theta type
  if (transform %in% c("identity", "variance", "proportional", "additive")) {
    lower_95          <- if (!is.na(se)) estimate - 1.96 * se else NA_real_
    upper_95          <- if (!is.na(se)) estimate + 1.96 * se else NA_real_
    estimate_natural  <- NA_real_
    lower_95_natural  <- NA_real_
    upper_95_natural  <- NA_real_
  } else if (transform == "log") {
    lower_95          <- if (!is.na(se)) estimate - 1.96 * se else NA_real_
    upper_95          <- if (!is.na(se)) estimate + 1.96 * se else NA_real_
    estimate_natural  <- if (!is.na(se)) exp(estimate) else NA_real_
    lower_95_natural  <- if (!is.na(se)) exp(estimate - 1.96 * se) else NA_real_
    upper_95_natural  <- if (!is.na(se)) exp(estimate + 1.96 * se) else NA_real_
  } else if (transform %in% c("logit", "logit_probability")) {
    # theta is on the logit scale; CI is symmetric on logit then back-transformed
    lower_95          <- if (!is.na(se)) estimate - 1.96 * se else NA_real_
    upper_95          <- if (!is.na(se)) estimate + 1.96 * se else NA_real_
    estimate_natural  <- if (!is.na(se)) .ferx_inv_logit(estimate) else NA_real_
    lower_95_natural  <- if (!is.na(se)) .ferx_inv_logit(estimate - 1.96 * se) else NA_real_
    upper_95_natural  <- if (!is.na(se)) .ferx_inv_logit(estimate + 1.96 * se) else NA_real_
  } else {
    lower_95          <- if (!is.na(se)) estimate - 1.96 * se else NA_real_
    upper_95          <- if (!is.na(se)) estimate + 1.96 * se else NA_real_
    estimate_natural  <- NA_real_
    lower_95_natural  <- NA_real_
    upper_95_natural  <- NA_real_
  }

  data.frame(param            = param,
             transform        = transform,
             estimate         = estimate,
             se               = se,
             rse_pct          = rse_pct,
             lower_95         = lower_95,
             upper_95         = upper_95,
             estimate_natural = estimate_natural,
             lower_95_natural = lower_95_natural,
             upper_95_natural = upper_95_natural,
             init_as_sd       = init_as_sd,
             stringsAsFactors = FALSE)
}

# Resolve once whether cli colour output is available. Called at the top of
# ferx_warnings() and from print.ferx_fit() so the capability check is not
# repeated inside per-row loops.
#
# We deliberately do NOT gate on `isatty(stdout())`. Inside RStudio's console
# `stdout()` is a pipe (so isatty() returns FALSE), yet RStudio renders ANSI
# escape codes correctly - and RStudio is the primary interactive audience
# for ferx. `cli::num_ansi_colors()` already performs its own RStudio /
# terminal / R.app detection and returns > 1 only when colour will render,
# so it is sufficient by itself.
.ferx_use_cli <- function() {
  tryCatch(
    requireNamespace("cli", quietly = TRUE) &&
      cli::num_ansi_colors() > 1L,
    error = function(e) FALSE
  )
}

# Apply an ANSI style to text.  use_cli must be pre-computed by .ferx_use_cli()
# and passed in so the capability check is not repeated on every call.
# Falls back to the plain string when use_cli is FALSE.
.ferx_style <- function(text, style, use_cli = .ferx_use_cli()) {
  if (!use_cli) return(text)
  switch(style,
    bold   = cli::style_bold(text),
    green  = cli::col_green(cli::style_bold(text)),
    red    = cli::col_red(cli::style_bold(text)),
    yellow = cli::col_yellow(text),
    dim    = cli::col_grey(text),
    text
  )
}

# Remediation guidance keyed by the fixed category vocabulary that ferx-core
# (and the R-side additions) emit. Extend this table when core grows a new
# category. Returns NULL for unknown categories so callers can skip printing
# rather than showing a generic placeholder.
#
# message and uses_sde are used for dw_autocorrelation (positive vs negative
# DW; SDE hint suppressed when already in use) and for covariance_step (the
# specific failure mode embedded in the message text selects targeted advice).
.ferx_warning_guidance <- function(category, message = "", uses_sde = FALSE) {
  if (category == "dw_autocorrelation") {
    if (grepl("egative", message, ignore.case = TRUE)) {
      return("Negative IWRES autocorrelation suggests over-parameterisation or a misspecified error model. Consider removing a parameter or simplifying the residual model.")
    }
    sde_hint <- if (isTRUE(uses_sde)) "" else
      " For ODE models, also consider SDE process noise ([diffusion] block)."
    return(paste0(
      "Positive IWRES autocorrelation suggests missing structural dynamics.",
      " Consider transit absorption, an extra compartment, or IOV on ka/F.",
      sde_hint
    ))
  }
  if (category == "covariance_step") {
    # Omega non-PD -- checked before general non-PD because omega messages also
    # contain "not positive definite" and "eigenvalue".
    if (grepl("omega matrix is not positive definite", message, ignore.case = TRUE)) {
      return(paste0(
        "Omega is near-singular at convergence. Consider a diagonal omega ",
        "structure, fixing a small variance to a small positive constant, or ",
        "removing the corresponding ETA from the model."
      ))
    }
    # NonPdHessian path: eigenvalue list is present in the message.
    if (grepl("not positive definite", message, ignore.case = TRUE) &&
        grepl("eigenvalue", message, ignore.case = TRUE)) {
      return(paste0(
        "Inspect the eigenvalue list in the warning: a near-zero minimum ",
        "(e.g. 1e-6) suggests a near-unidentifiable parameter; a clearly ",
        "negative minimum indicates structural non-identifiability. Consider ",
        "fixing the most collinear parameter, removing it, or switching to a ",
        "diagonal omega structure."
      ))
    }
    # Non-finite or zero Hessian diagonal: parameter name(s) listed in message.
    if (grepl("ill-conditioned entries", message, ignore.case = TRUE)) {
      return(paste0(
        "The named parameter(s) have a flat or non-finite Hessian diagonal, ",
        "meaning the parameter is not informed by the data or the objective ",
        "function overflows near convergence. Consider fixing the parameter, ",
        "tightening its bounds, or increasing fd_hessian_step ",
        "(e.g. ferx_fit(..., fd_hessian_step = 0.05))."
      ))
    }
    # Model evaluation overflow/underflow.
    if (grepl("base ofv is non-finite", message, ignore.case = TRUE)) {
      return(paste0(
        "Model evaluation overflowed or underflowed at convergence. Check for ",
        "extreme parameter values, verify that all DV values are positive for ",
        "proportional error, and consider constraining thetas to physiologically ",
        "plausible ranges."
      ))
    }
    # Regularisation path -- severity is embedded in the message.
    if (grepl("covariance step regularized", message, ignore.case = TRUE)) {
      if (grepl("severity: severe", message, ignore.case = TRUE)) {
        return(paste0(
          "Severe Hessian regularisation: standard errors are likely unreliable. ",
          "Run ferx_sir() to obtain non-parametric confidence intervals, or ",
          "simplify the model structure."
        ))
      }
      if (grepl("severity: moderate", message, ignore.case = TRUE)) {
        return(paste0(
          "Moderate Hessian regularisation: standard errors should be interpreted ",
          "with caution. Run ferx_sir() to obtain non-parametric confidence ",
          "intervals as a cross-check."
        ))
      }
      # severity: minor (or any unrecognised tier from future core versions) --
      # treat as benign; minor is the only tier ferx-core emits below moderate.
      return(paste0(
        "Minor Hessian regularisation: standard errors are likely reliable. A ",
        "small eigenvalue floor was applied; this is common on smooth OFV ",
        "surfaces and is usually benign."
      ))
    }
    # Generic fallback for older or unrecognised covariance_step messages.
    return(paste0(
      "Standard errors unavailable. Check identifiability; try a simpler ",
      "omega/sigma structure or covariance = FALSE for development."
    ))
  }
  switch(category,
    convergence        = "Optimizer did not reach convergence. Try different initial values, method = c(\"saem\", \"focei\"), or settings = list(n_starts = 4L).",
    condition_number   = "Parameters are correlated/ill-scaled. Consider fixing or removing a parameter, or reparameterising.",
    optimizer_health   = "Optimizer struggled (trust region / Hessian). Inspect the trace and consider better starting values.",
    eta_normality      = "ETA distribution may be non-normal. High shrinkage or sparse data can cause this; prefer QQ-plots for diagnosis.",
    bloq_method        = "BLOQ handling note. Set method = \"focei\" explicitly to silence, or review the M3 setup.",
    sir                = "SIR uncertainty step issue. Ensure covariance = TRUE and inspect SIR tuning (sir_samples / sir_resamples).",
    importance_sampling = "Importance-sampling ESS collapsed for some subjects. Raise imp_samples / imp_proposal_df or check EBE quality.",
    data_quality       = "Data issue detected. Review the flagged observations in the dataset.",
    omega_structure    = "Mixed parameterisation in a block omega. Check the [individual_parameters] forms for the correlated etas.",
    ebe_convergence    = "Some subjects' inner EBE search did not converge. Inspect those subjects or relax inner_tol / max_unconverged_frac.",
    gradient_fallback  = "Gradient method fell back (e.g. AD -> FD or HMC -> MH). The fit is valid; expect a longer runtime.",
    mu_referencing     = "Mu-referencing was auto-detected for the listed parameters (informational).",
    optimizer_config   = "Optimizer configuration note (informational).",
    multi_start        = "Multi-start information (informational).",
    threads            = "Thread-pool sizing note. Consider matching threads to the subject count.",
    cancelled          = "The fit was cancelled before completion.",
    unused_parameter   = "A declared parameter is never referenced in [individual_parameters] or [error_model]. Remove it from [parameters] or complete the expression that uses it.",
    NULL
  )
}

#' Structured warnings from a ferx fit
#'
#' Returns or prints the warnings produced by a fit, classified by severity
#' (\code{critical}, \code{warning}, \code{info}) and category. Severity and
#' category are assigned by the ferx-core engine (for engine warnings) or by
#' the R diagnostics layer (for condition number and ETA normality); the R
#' side never re-parses message text to guess severity.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @param as_df Logical. When \code{TRUE}, returns the raw data frame
#'   (\code{severity}, \code{category}, \code{message}, \code{source_method})
#'   instead of pretty-printing. Default \code{FALSE}.
#' @return When \code{as_df = TRUE}, a data frame. Otherwise the same data
#'   frame invisibly, after printing a grouped, colour-coded summary
#'   (critical first, then warning, then info).
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' ferx_warnings(fit)
#' ferx_warnings(fit, as_df = TRUE)
#' @family diagnostics
#' @export
ferx_warnings <- function(fit, as_df = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object")
  }
  df <- fit$warnings_structured
  if (is.null(df) || !is.data.frame(df)) {
    # Backward-compat fallback for fits saved before structured warnings:
    # surface the flat character vector as undifferentiated warnings.
    msgs <- fit$warnings %||% character(0)
    df <- data.frame(
      severity      = rep("warning", length(msgs)),
      category      = rep("general", length(msgs)),
      message       = as.character(msgs),
      source_method = rep("", length(msgs)),
      stringsAsFactors = FALSE
    )
  }
  if (isTRUE(as_df)) {
    return(df)
  }

  use_cli   <- .ferx_use_cli()
  uses_sde  <- isTRUE(fit$uses_sde)
  model_lbl <- fit$model_name %||% "fit"
  cat(sprintf("ferx fit warnings  (%s)\n", model_lbl))
  cat(strrep("-", 49), "\n", sep = "")

  if (nrow(df) == 0L) {
    cat("  No warnings.\n")
    cat(strrep("-", 49), "\n", sep = "")
    return(invisible(df))
  }

  order_lvl <- c(critical = 1L, warning = 2L, info = 3L)
  sev_key <- order_lvl[df$severity]
  sev_key[is.na(sev_key)] <- 99L
  df_ord <- df[order(sev_key), , drop = FALSE]

  label_for <- function(sev) {
    switch(sev,
      critical = .ferx_style("[CRITICAL]", "red",    use_cli),
      warning  = .ferx_style("[WARNING] ", "yellow", use_cli),
      info     = .ferx_style("[INFO]    ", "dim",    use_cli),
      sprintf("[%s]", toupper(sev))
    )
  }
  wrap_indent <- function(text, indent = "            ", width = 70L) {
    parts <- strwrap(text, width = width)
    paste0(indent, parts, collapse = "\n")
  }

  for (i in seq_len(nrow(df_ord))) {
    row <- df_ord[i, ]
    cat(sprintf("%s %s\n", label_for(row$severity), row$category))
    cat(wrap_indent(row$message), "\n", sep = "")
    guide <- .ferx_warning_guidance(row$category,
                                    message  = row$message,
                                    uses_sde = uses_sde)
    if (!is.null(guide)) {
      cat(.ferx_style(wrap_indent(guide), "dim", use_cli), "\n", sep = "")
    }
    cat("\n")
  }

  n_crit <- sum(df$severity == "critical")
  n_warn <- sum(df$severity == "warning")
  n_info <- sum(df$severity == "info")
  cat(strrep("-", 49), "\n", sep = "")
  cat(sprintf("%s   %s   %s\n",
    .ferx_style(sprintf("%d CRITICAL", n_crit), if (n_crit > 0) "red" else "dim", use_cli),
    .ferx_style(sprintf("%d WARNING",  n_warn), if (n_warn > 0) "yellow" else "dim", use_cli),
    .ferx_style(sprintf("%d INFO",     n_info), "dim", use_cli)
  ))
  invisible(df)
}
