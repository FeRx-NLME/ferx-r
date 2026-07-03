# Tidy parameter estimates table: theta, omega diagonal, sigma, and (for IOV
# models) kappa diagonal, with percent relative standard error (%RSE), 95%
# confidence intervals, and (for log/logit-transformed thetas) natural-scale
# back-transformed estimates and CIs. Omega is reported on the variance scale
# (matching the .ferx model file convention); for block omega, only the
# diagonal variances are included. Stored on the fit object as
# `fit$estimates`. (Formerly the exported ferx_estimates(fit); see issue #226.)
.ferx_compute_estimates <- function(fit) {
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
