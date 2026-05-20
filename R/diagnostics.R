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
#'
#' # Inspect autocorrelation
#' diag$autocorrelation
#'
#' # Plot IWRES vs time coloured by subject, with DW annotation
#' if (!is.null(diag$autocorrelation)) {
#'   sd <- fit$sdtab
#'   plot(sd$TIME, sd$IWRES, col = as.factor(sd$ID),
#'        main = sprintf("IWRES vs Time  |  DW = %.2f  lag-1 r = %.3f",
#'                       diag$autocorrelation$dw_statistic,
#'                       diag$autocorrelation$lag1_r),
#'        xlab = "TIME", ylab = "IWRES", pch = 16, cex = 0.7)
#'   abline(h = 0, lty = 2)
#' }
#' }
.dw_label <- function(dw) {
  if (dw < 1.5) "positive autocorrelation"
  else if (dw > 2.5) "negative autocorrelation"
  else "no autocorrelation"
}

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
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' obs <- read.csv(ex$data)
#' ferx_eta_cov(fit, obs)
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
    se       <- if (!is.null(fit$se_omega) && length(fit$se_omega) >= i) fit$se_omega[i] else NA_real_
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
