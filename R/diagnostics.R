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
#' @export
ferx_cor_matrix <- function(fit) {
  if (is.null(fit$cov_matrix)) {
    stop(
      "Covariance matrix not available. ",
      "Run ferx_fit() with covariance = TRUE and check fit$covariance_status."
    )
  }
  d  <- nrow(fit$cov_matrix)
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
#' @export
ferx_eta_cov <- function(fit, data) {
  if (is.null(fit$sdtab) || !is.data.frame(fit$sdtab)) {
    stop("`fit$sdtab` is not available.")
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("`data` must be a data.frame — pass the dataset used in ferx_fit().")
  }

  eta_cols <- grep("^ETA", names(fit$sdtab), value = TRUE)
  if (length(eta_cols) == 0L) {
    message("No ETA columns found in fit$sdtab.")
    return(invisible(NULL))
  }

  # Time-varying or non-covariate columns to skip
  SKIP <- c("TIME", "DV", "AMT", "EVID", "MDV", "CMT", "RATE",
            "II", "SS", "CENS", "LLOQ", "BLQ")

  sdtab_id <- if ("ID" %in% names(fit$sdtab)) "ID" else names(fit$sdtab)[1L]
  data_id  <- if ("ID" %in% names(data))      "ID" else names(data)[1L]

  # One row per subject from sdtab
  etas <- fit$sdtab[!duplicated(fit$sdtab[[sdtab_id]]),
                    c(sdtab_id, eta_cols), drop = FALSE]

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
                  by.x = sdtab_id, by.y = data_id)

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

#' Tidy parameter estimates table
#'
#' Extracts all estimated parameters (theta, omega diagonal, sigma) into a
#' single tidy data frame, adding percent relative standard error (\%RSE) and
#' 95\% confidence intervals when the covariance step was run.
#'
#' Omega is reported on the variance scale (matching the \code{.ferx} model
#' file convention). For block omega, only the diagonal variances are included.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return A data frame with columns \code{param}, \code{estimate}, \code{se},
#'   \code{rse_pct}, \code{lower_95}, \code{upper_95}. SE-derived columns are
#'   \code{NA} when the covariance step was not run or failed.
#' @seealso \code{\link{ferx_cor_matrix}} for parameter correlations.
#' @export
ferx_estimates <- function(fit) {
  rows <- list()

  # Theta
  theta_names <- names(fit$theta)
  if (is.null(theta_names)) theta_names <- paste0("THETA", seq_along(fit$theta))
  for (i in seq_along(fit$theta)) {
    se <- if (!is.null(fit$se_theta) && length(fit$se_theta) >= i) fit$se_theta[i] else NA_real_
    rows[[length(rows) + 1L]] <- .ferx_est_row(theta_names[i], fit$theta[i], se)
  }

  # Omega diagonal (variance scale)
  om    <- fit$omega
  if (is.null(dim(om))) om <- matrix(om, 1L, 1L)
  n_eta <- nrow(om)
  for (i in seq_len(n_eta)) {
    pname <- sprintf("OMEGA(%d,%d)", i, i)
    se    <- if (!is.null(fit$se_omega) && length(fit$se_omega) >= i) fit$se_omega[i] else NA_real_
    rows[[length(rows) + 1L]] <- .ferx_est_row(pname, om[i, i], se)
  }

  # Sigma
  for (i in seq_along(fit$sigma)) {
    pname <- sprintf("SIGMA(%d)", i)
    se    <- if (!is.null(fit$se_sigma) && length(fit$se_sigma) >= i) fit$se_sigma[i] else NA_real_
    rows[[length(rows) + 1L]] <- .ferx_est_row(pname, fit$sigma[i], se)
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

.ferx_est_row <- function(param, estimate, se) {
  rse_pct  <- if (!is.na(se) && abs(estimate) > 1e-12) abs(se / estimate) * 100 else NA_real_
  lower_95 <- if (!is.na(se)) estimate - 1.96 * se else NA_real_
  upper_95 <- if (!is.na(se)) estimate + 1.96 * se else NA_real_
  data.frame(param    = param,
             estimate = estimate,
             se       = se,
             rse_pct  = rse_pct,
             lower_95 = lower_95,
             upper_95 = upper_95,
             stringsAsFactors = FALSE)
}
