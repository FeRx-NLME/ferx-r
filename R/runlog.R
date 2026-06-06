#' NONMEM-style run log for a ferx fit
#'
#' Produces a plain-text summary modelled on the information content of a
#' NONMEM \code{.lst} file: run header, model source, data summary, initial
#' vs. final parameter table with SE and \%RSE, covariance-step diagnostics,
#' shrinkage and normality diagnostics, and (where available) the final
#' gradient vector with a convergence diagnostic.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} or
#'   \code{\link{ferx_load_fit}}.
#' @param gradient_tol Numeric scalar.  Threshold for the relative gradient
#'   (\code{|grad[i] * param[i]|}); components exceeding this trigger a
#'   convergence warning.  Default \code{0.1}.
#' @param show_iterations Logical.  When \code{TRUE} (the default) and the fit
#'   was run with \code{optimizer_trace = TRUE}, an "Iteration history" section
#'   is inserted showing OFV, delta-OFV, and method-specific convergence
#'   metrics (gradient norm for FOCE/BFGS, LM-lambda for Gauss-Newton,
#'   MH-accept-rate for SAEM) per iteration.  Runs with more than 30 iterations
#'   are truncated to the first 10 and last 10 rows.  Set to \code{FALSE} to
#'   suppress this section.
#' @param verbose Logical.  When \code{TRUE} (the default) the output is
#'   printed to the console.  When \code{FALSE} it is returned invisibly as a
#'   single character string.
#'
#' @return A single character string containing the full log (invisibly when
#'   \code{verbose = TRUE}).
#'
#' @details
#' Sections that depend on fields introduced in engine version 0.1.6
#' (\code{model_text}, \code{theta_init}, \code{omega_init}, \code{sigma_init},
#' \code{obs_time_range}, \code{final_gradient}) are silently omitted when
#' those fields are absent or \code{NULL} -- so the function works against
#' older \code{.fitrx} bundles or in-memory fits that do not carry file
#' provenance.
#'
#' The data summary shows subject and observation counts when \code{fit$n_subjects}
#' and \code{fit$n_obs} are present; observation time range when
#' \code{fit$obs_time_range} is present.  Missing counts produce a partial
#' summary rather than an error.
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' ferx_runlog(fit)
#' ferx_runlog(fit, verbose = FALSE)  # capture as string
#' }
#'
#' @export
ferx_runlog <- function(fit, gradient_tol = 0.1, show_iterations = TRUE, verbose = TRUE) {
  lines <- character(0)

  .sec <- function(title) {
    c(strrep("=", 60), toupper(title), strrep("-", 60))
  }
  .blank <- function() ""

  # -- Run header -------------------------------------------------------------
  lines <- c(lines, strrep("=", 60))
  lines <- c(lines, sprintf("  ferx run log"))
  if (!is.null(fit$model_name) && nzchar(fit$model_name)) {
    lines <- c(lines, sprintf("  Model:   %s", fit$model_name))
  }
  lines <- c(lines, sprintf("  Method:  %s", toupper(fit$method %||% "unknown")))
  lines <- c(lines, sprintf("  ferx:    %s", fit$ferx_version %||% "unknown"))
  lines <- c(lines, sprintf("  Date:    %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  lines <- c(lines, strrep("=", 60))
  lines <- c(lines, .blank())

  # -- Model file -------------------------------------------------------------
  model_text <- fit$model_text
  if (!is.null(model_text) && !is.na(model_text) && nzchar(model_text)) {
    lines <- c(lines, .sec("Model file"), model_text, .blank())
  }

  # -- Data summary -----------------------------------------------------------
  lines <- c(lines, .sec("Data summary"))
  n_subj <- fit$n_subjects
  n_obs  <- fit$n_obs
  tr     <- fit$obs_time_range

  data_name <- fit$data_name %||% NULL
  if (!is.null(data_name) && nzchar(data_name)) {
    lines <- c(lines, sprintf("  Dataset:        %s", data_name))
  }
  if (!is.null(n_subj) && !is.null(n_obs)) {
    lines <- c(lines, sprintf("  Subjects:       %d", as.integer(n_subj)))
    lines <- c(lines, sprintf("  Observations:   %d", as.integer(n_obs)))
  }
  if (!is.null(tr) && length(tr) == 2L && all(is.finite(tr))) {
    lines <- c(lines, sprintf("  Time range:     %.4g to %.4g", tr[1], tr[2]))
  }
  if (!is.null(fit$input_columns) && length(fit$input_columns) > 0L) {
    lines <- c(lines, sprintf("  Input columns:  %s", paste(fit$input_columns, collapse = "  ")))
  }
  ex <- fit$exclusions
  if (!is.null(ex)) {
    lines <- c(lines, sprintf(
      "  Data selection: %d obs, %d doses, %d other excluded (of %d records)",
      as.integer(ex$n_obs_excluded  %||% 0L),
      as.integer(ex$n_dose_excluded %||% 0L),
      as.integer(ex$n_other_excluded %||% 0L),
      as.integer(ex$n_records_total %||% 0L)
    ))
    for (cond in ex$fired_ignore)
      lines <- c(lines, sprintf("    Fired ignore:  \"%s\"", cond))
    for (cond in ex$fired_accept)
      lines <- c(lines, sprintf("    Fired accept:  \"%s\"", cond))
    if (length(ex$excluded_subject_ids) > 0L)
      lines <- c(lines, sprintf("    Excl. subjs:   [%s]",
                                paste(ex$excluded_subject_ids, collapse = ", ")))
  }
  if (is.null(n_subj) && is.null(n_obs) && is.null(tr) &&
      (is.null(fit$input_columns) || length(fit$input_columns) == 0L) &&
      is.null(ex)) {
    lines <- c(lines, "  (data summary not available for this fit)")
  }
  lines <- c(lines, .blank())

  # -- Parameter table: INITIAL vs FINAL with SE / %RSE ----------------------
  has_inits <- !is.null(fit$theta_init) && length(fit$theta_init) > 0L
  has_se    <- !is.null(fit$se_theta)   && length(fit$se_theta)   > 0L

  lines <- c(lines, .sec("Parameter estimates"))

  hdr_width <- if (has_inits) 60 else 46
  if (has_se) hdr_width <- hdr_width + 20
  if (has_inits && has_se) {
    lines <- c(lines, sprintf(
      "  %-28s  %14s  %14s  %8s  %7s",
      "PARAMETER", "INITIAL", "FINAL", "SE", "%RSE"
    ))
  } else if (has_se) {
    lines <- c(lines, sprintf(
      "  %-28s  %14s  %8s  %7s",
      "PARAMETER", "FINAL", "SE", "%RSE"
    ))
  } else if (has_inits) {
    lines <- c(lines, sprintf("  %-28s  %14s  %14s", "PARAMETER", "INITIAL", "FINAL"))
  } else {
    lines <- c(lines, sprintf("  %-28s  %14s", "PARAMETER", "FINAL"))
  }
  lines <- c(lines, sprintf("  %s", strrep("-", hdr_width)))

  # Theta rows
  thetas   <- fit$theta       %||% numeric(0)
  th_names <- names(fit$theta) %||% fit$theta_names %||%
              fit$model_structure$theta_names %||%
              paste0("THETA(", seq_along(thetas), ")")
  th_init  <- if (has_inits) fit$theta_init else NULL
  th_fixed <- fit$theta_fixed %||% rep(FALSE, length(thetas))
  se_th    <- if (has_se) fit$se_theta else NULL

  for (i in seq_along(thetas)) {
    nm    <- if (i <= length(th_names)) th_names[i] else sprintf("THETA(%d)", i)
    fixed <- isTRUE(th_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", thetas[i])
    se_str    <- ""
    rse_str   <- ""
    if (has_se && !fixed && !is.null(se_th) && i <= length(se_th)) {
      se_val  <- se_th[i]
      se_str  <- if (is.finite(se_val)) sprintf("%8.4g", se_val) else sprintf("%8s", "N/A")
      rse_val <- if (is.finite(se_val) && thetas[i] != 0) abs(se_val / thetas[i]) * 100 else NA
      rse_str <- if (is.finite(rse_val)) sprintf("%6.1f%%", rse_val) else sprintf("%7s", "N/A")
    } else if (has_se && fixed) {
      se_str  <- sprintf("%8s", "FIXED")
      rse_str <- sprintf("%7s", "FIXED")
    }
    if (has_inits && has_se) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", th_init[i])
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s  %s", nm, init_str, final_str, se_str, rse_str))
    } else if (has_se) {
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s", nm, final_str, se_str, rse_str))
    } else if (has_inits) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", th_init[i])
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else {
      lines <- c(lines, sprintf("  %-28s  %s", nm, final_str))
    }
  }

  # Omega rows: bare declared name (e.g. ETA_CL); fallback OMEGA(i,i)
  om      <- fit$omega
  if (is.null(dim(om))) om <- matrix(om, 1L, 1L)
  n_eta   <- nrow(om)
  om_init <- fit$omega_init
  if (!is.null(om_init) && is.null(dim(om_init))) {
    om_init <- matrix(om_init, n_eta, n_eta)
  }
  has_om_init <- !is.null(om_init) && all(dim(om_init) == n_eta)
  om_fixed <- fit$omega_fixed %||% rep(FALSE, n_eta)
  eta_nms  <- fit$eta_names   %||% paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
  se_om    <- fit$se_omega

  for (i in seq_len(n_eta)) {
    nm    <- if (i <= length(eta_nms)) eta_nms[i] else sprintf("OMEGA(%d,%d)", i, i)
    fixed <- isTRUE(om_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", om[i, i])
    se_str  <- ""
    rse_str <- ""
    if (has_se && !is.null(se_om) && i <= length(se_om) && !fixed) {
      se_val  <- se_om[i]
      se_str  <- if (is.finite(se_val)) sprintf("%8.4g", se_val) else sprintf("%8s", "N/A")
      rse_val <- if (is.finite(se_val) && om[i, i] != 0) abs(se_val / om[i, i]) * 100 else NA
      rse_str <- if (is.finite(rse_val)) sprintf("%6.1f%%", rse_val) else sprintf("%7s", "N/A")
    } else if (has_se && fixed) {
      se_str  <- sprintf("%8s", "FIXED")
      rse_str <- sprintf("%7s", "FIXED")
    }
    if (has_inits && has_om_init && has_se) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", om_init[i, i])
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s  %s", nm, init_str, final_str, se_str, rse_str))
    } else if (has_inits && has_om_init) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", om_init[i, i])
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else if (has_inits && has_se) {
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s  %s",
        nm, sprintf("%14s", "N/A"), final_str, se_str, rse_str))
    } else if (has_se) {
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s", nm, final_str, se_str, rse_str))
    } else if (has_inits) {
      init_str <- if (has_om_init) sprintf("%14.6g", om_init[i, i]) else sprintf("%14s", "N/A")
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else {
      lines <- c(lines, sprintf("  %-28s  %s", nm, final_str))
    }
  }

  # Sigma rows: bare declared name (e.g. EPS_PROP); fallback SIGMA(i)
  sigmas    <- fit$sigma       %||% numeric(0)
  sig_nms   <- fit$sigma_names %||% paste0("SIGMA(", seq_along(sigmas), ")")
  sig_init  <- if (has_inits && !is.null(fit$sigma_init) &&
                   length(fit$sigma_init) > 0L) fit$sigma_init else NULL
  sig_fixed <- fit$sigma_fixed %||% rep(FALSE, length(sigmas))
  se_sig    <- fit$se_sigma

  for (i in seq_along(sigmas)) {
    nm    <- if (i <= length(sig_nms)) sig_nms[i] else sprintf("SIGMA(%d)", i)
    fixed <- isTRUE(sig_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", sigmas[i])
    se_str  <- ""
    rse_str <- ""
    if (has_se && !is.null(se_sig) && i <= length(se_sig) && !fixed) {
      se_val  <- se_sig[i]
      se_str  <- if (is.finite(se_val)) sprintf("%8.4g", se_val) else sprintf("%8s", "N/A")
      rse_val <- if (is.finite(se_val) && sigmas[i] != 0) abs(se_val / sigmas[i]) * 100 else NA
      rse_str <- if (is.finite(rse_val)) sprintf("%6.1f%%", rse_val) else sprintf("%7s", "N/A")
    } else if (has_se && fixed) {
      se_str  <- sprintf("%8s", "FIXED")
      rse_str <- sprintf("%7s", "FIXED")
    }
    if (has_inits && !is.null(sig_init) && has_se) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", sig_init[i])
      lines <- c(lines, sprintf("  %-28s  %s  %s  %s  %s", nm, init_str, final_str, se_str, rse_str))
    } else if (has_inits && !is.null(sig_init)) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", sig_init[i])
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else if (has_se) {
      init_str <- if (has_inits) sprintf("%14s", "N/A") else NULL
      if (!is.null(init_str)) {
        lines <- c(lines, sprintf("  %-28s  %s  %s  %s  %s", nm, init_str, final_str, se_str, rse_str))
      } else {
        lines <- c(lines, sprintf("  %-28s  %s  %s  %s", nm, final_str, se_str, rse_str))
      }
    } else if (has_inits) {
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, sprintf("%14s", "N/A"), final_str))
    } else {
      lines <- c(lines, sprintf("  %-28s  %s", nm, final_str))
    }
  }
  lines <- c(lines, .blank())

  # -- Estimation settings ----------------------------------------------------
  has_settings <- (!is.null(fit$optimizer_label)   && !is.na(fit$optimizer_label)   && nzchar(fit$optimizer_label)) ||
                  (!is.null(fit$bloq_method_label) && !is.na(fit$bloq_method_label) && nzchar(fit$bloq_method_label)) ||
                  (!is.null(fit$inits_from_nca)    && !is.na(fit$inits_from_nca)    && nzchar(fit$inits_from_nca)) ||
                  (!is.null(fit$covariate_names)   && length(fit$covariate_names) > 0L)
  if (has_settings) {
    lines <- c(lines, .sec("Estimation settings"))
    if (!is.null(fit$optimizer_label) && !is.na(fit$optimizer_label) && nzchar(fit$optimizer_label)) {
      n_st <- fit$n_starts %||% 1L
      start_str <- if (n_st > 1L) sprintf("  Optimizer:      %s  (multi-start n=%d)", fit$optimizer_label, n_st) else
                   sprintf("  Optimizer:      %s", fit$optimizer_label)
      lines <- c(lines, start_str)
    }
    if (!is.null(fit$outer_maxiter) && fit$outer_maxiter > 0L) {
      lines <- c(lines, sprintf("  Max iterations: %d", as.integer(fit$outer_maxiter)))
    }
    has_gradient_optimizer <- !is.null(fit$gradient_method_outer) &&
                              !grepl("N/A|not_applicable", fit$gradient_method_outer %||% "", ignore.case = TRUE)
    if (has_gradient_optimizer && !is.null(fit$outer_gtol) && is.finite(fit$outer_gtol)) {
      lines <- c(lines, sprintf("  Gradient tol:   %.4g", fit$outer_gtol))
    }
    if (!is.null(fit$bloq_method_label) && !is.na(fit$bloq_method_label) && nzchar(fit$bloq_method_label)) {
      lines <- c(lines, sprintf("  BLOQ method:    %s", fit$bloq_method_label))
    }
    if (!is.null(fit$inits_from_nca) && !is.na(fit$inits_from_nca) && nzchar(fit$inits_from_nca)) {
      lines <- c(lines, sprintf("  Initial values: NCA-derived (%s)", fit$inits_from_nca))
    }
    # Seeds -- only print when non-NULL and non-NA
    seed_parts <- character(0)
    if (!is.null(fit$multi_start_seed) && !is.na(fit$multi_start_seed)) seed_parts <- c(seed_parts, sprintf("multi_start=%g", fit$multi_start_seed))
    if (!is.null(fit$saem_seed)        && !is.na(fit$saem_seed))        seed_parts <- c(seed_parts, sprintf("saem=%g",        fit$saem_seed))
    if (!is.null(fit$sir_seed_used)    && !is.na(fit$sir_seed_used))    seed_parts <- c(seed_parts, sprintf("sir=%g",         fit$sir_seed_used))
    if (!is.null(fit$is_seed)          && !is.na(fit$is_seed))          seed_parts <- c(seed_parts, sprintf("is=%g",          fit$is_seed))
    if (length(seed_parts)) {
      lines <- c(lines, sprintf("  Seeds:          %s", paste(seed_parts, collapse = "  ")))
    }
    if (!is.null(fit$covariate_names) && length(fit$covariate_names) > 0L) {
      lines <- c(lines, sprintf("  Covariates:     %s", paste(fit$covariate_names, collapse = ", ")))
    }
    lines <- c(lines, .blank())
  }

  # -- OFV / AIC / BIC -------------------------------------------------------
  lines <- c(lines, .sec("Objective function"))
  lines <- c(lines,
    sprintf("  OFV   = %.6g", fit$ofv %||% NA_real_),
    sprintf("  AIC   = %.6g", fit$aic %||% NA_real_),
    sprintf("  BIC   = %.6g", fit$bic %||% NA_real_),
    sprintf("  Convergence: %s", if (isTRUE(fit$converged)) "YES" else "NO")
  )
  lines <- c(lines, .blank())

  # -- Iteration history (requires optimizer_trace = TRUE) -------------------
  has_trace <- isTRUE(show_iterations) &&
               !is.null(fit$trace_path)   &&
               !is.na(fit$trace_path)     &&
               nzchar(fit$trace_path)     &&
               file.exists(fit$trace_path)
  if (has_trace) {
    tr <- tryCatch(ferx_trace(fit), error = function(e) NULL)
    if (!is.null(tr) && nrow(tr) > 0L) {
      lines <- c(lines, .sec("Iteration history"))
      lines <- c(lines, .runlog_iter_table(tr))
      lines <- c(lines, sprintf("  Total: %d iteration(s)", nrow(tr)))
      lines <- c(lines, .blank())
    }
  }

  # -- Covariance step --------------------------------------------------------
  has_cov <- !is.null(fit$condition_number) || !is.null(fit$eigenvalues)
  if (has_cov) {
    lines <- c(lines, .sec("Covariance step"))
    if (!is.null(fit$condition_number) && is.finite(fit$condition_number)) {
      cn <- fit$condition_number
      cn_flag <- if (cn > 1000) "  [!] condition number > 1000 -- possible ill-conditioning" else ""
      lines <- c(lines, sprintf("  Condition number: %.4g%s", cn, cn_flag))
    }
    if (!is.null(fit$eigenvalues) && length(fit$eigenvalues) > 0L) {
      ev_str <- paste(sprintf("%.4g", fit$eigenvalues), collapse = "  ")
      lines <- c(lines, sprintf("  Eigenvalues:      %s", ev_str))
    }
    lines <- c(lines, .blank())
  }

  # -- Diagnostics ------------------------------------------------------------
  has_shrink <- (!is.null(fit$shrinkage_eta) && any(!is.na(fit$shrinkage_eta))) ||
                (!is.null(fit$shrinkage_eps) && !is.na(fit$shrinkage_eps))
  has_dw     <- !is.null(fit$dw_statistic) && !is.na(fit$dw_statistic)
  has_norm   <- !is.null(fit$eta_normality) && nrow(fit$eta_normality) > 0L

  if (has_shrink || has_dw || has_norm) {
    lines <- c(lines, .sec("Diagnostics"))

    if (has_shrink) {
      if (!is.null(fit$shrinkage_eta) && any(!is.na(fit$shrinkage_eta))) {
        eta_nms_sh <- fit$eta_names %||% paste0("ETA", seq_along(fit$shrinkage_eta))
        parts <- mapply(function(nm, s) {
          if (is.na(s)) NULL else sprintf("%s=%.1f%%", nm, s * 100)
        }, eta_nms_sh, fit$shrinkage_eta, SIMPLIFY = FALSE)
        parts <- Filter(Negate(is.null), parts)
        if (length(parts)) {
          lines <- c(lines, sprintf("  ETA shrinkage:  %s", paste(parts, collapse = "  ")))
        }
      }
      if (!is.null(fit$shrinkage_eps) && !is.na(fit$shrinkage_eps)) {
        lines <- c(lines, sprintf("  EPS shrinkage:  %.1f%%", fit$shrinkage_eps * 100))
      }
    }

    if (has_dw) {
      dw_label <- if (fit$dw_statistic < 1.5) "positive autocorrelation" else
                  if (fit$dw_statistic > 2.5) "negative autocorrelation" else
                  "no autocorrelation"
      lines <- c(lines, sprintf("  DW statistic:   %.4g (%s)", fit$dw_statistic, dw_label))
    }

    if (has_norm) {
      norm <- fit$eta_normality
      norm_parts <- vapply(seq_len(nrow(norm)), function(i) {
        row <- norm[i, ]
        flag <- if (!is.na(row$p_val) && row$p_val < 0.05) " [!]" else ""
        if (!is.na(row$p_val)) {
          sprintf("%s: W=%.3f p=%.3f%s", row$eta, row$W, row$p_val, flag)
        } else {
          sprintf("%s: W=%.3f", row$eta, row$W)
        }
      }, character(1L))
      lines <- c(lines, "  ETA normality (Shapiro-Wilk):")
      for (part in norm_parts) lines <- c(lines, sprintf("    %s", part))
    }

    lines <- c(lines, .blank())
  }

  # -- Final gradient ---------------------------------------------------------
  grad <- fit$final_gradient
  if (!is.null(grad) && length(grad) > 0L) {
    lines <- c(lines, .sec("Final gradient"))

    # Relative gradient: grad[i] * |param[i]| - scale-free, same units across
    # all parameter types. Packed order: thetas, omega diagonal, sigma.
    th_vals  <- abs(fit$theta  %||% numeric(0))
    om_diag  <- if (!is.null(fit$omega) && length(fit$omega) > 0L) abs(diag(fit$omega)) else numeric(0)
    sig_vals <- abs(fit$sigma  %||% numeric(0))
    param_scale <- c(th_vals, om_diag, sig_vals)

    if (length(param_scale) == length(grad)) {
      grad_rel <- grad * param_scale
      use_rel  <- TRUE
    } else {
      grad_rel <- grad
      use_rel  <- FALSE
    }

    large <- abs(grad_rel) > gradient_tol
    if (any(large)) {
      lines <- c(lines, sprintf(
        "  WARNING: %d component(s) exceed |relative gradient| > %.4g -- possible non-convergence.",
        sum(large), gradient_tol
      ))
    } else {
      lines <- c(lines, sprintf(
        "  All %d relative gradient components satisfy |grad * param| <= %.4g.",
        length(grad_rel), gradient_tol
      ))
    }

    label <- if (use_rel) "  Relative gradient (grad * |param|):" else "  Gradient (raw - parameter count mismatch):"
    lines <- c(lines, label)
    chunks <- split(grad_rel, ceiling(seq_along(grad_rel) / 8))
    for (chunk in chunks) {
      lines <- c(lines, sprintf("  [ %s ]", paste(sprintf("%.4g", chunk), collapse = "  ")))
    }
    lines <- c(lines, sprintf("  max|grad * param| = %.4g", max(abs(grad_rel))))
    lines <- c(lines, .blank())
  } else {
    lines <- c(lines, .sec("Final gradient"))
    optimizer <- fit$gradient_method_outer %||% ""
    msg <- if (grepl("N/A|not_applicable", optimizer, ignore.case = TRUE) ||
                identical(optimizer, "")) {
      "  Gradient not available for this optimizer (BOBYQA / BFGS / GN / SAEM)."
    } else {
      "  No gradient recorded (fit produced by an older engine version)."
    }
    lines <- c(lines, msg, .blank())
  }

  # -- Runtime ----------------------------------------------------------------
  lines <- c(lines, .sec("Runtime"))
  if (!is.null(fit$wall_time_secs) && is.finite(fit$wall_time_secs)) {
    lines <- c(lines, sprintf("  Wall time:   %.1f seconds", fit$wall_time_secs))
  }
  gm_inner <- fit$gradient_method_inner %||% ""
  gm_outer <- fit$gradient_method_outer %||% ""
  if (nzchar(gm_inner) || nzchar(gm_outer)) {
    lines <- c(lines, sprintf("  Gradient:    %s (inner) / %s (outer)",
      if (nzchar(gm_inner)) gm_inner else "N/A",
      if (nzchar(gm_outer)) gm_outer else "N/A"
    ))
  }
  lines <- c(lines, strrep("=", 60))

  # -- Assemble and return ----------------------------------------------------
  out <- paste(lines, collapse = "\n")
  if (isTRUE(verbose)) cat(out, "\n")
  invisible(out)
}

#' Print the full optimizer iteration history for a ferx fit
#'
#' Displays (or returns) the complete per-iteration table from the optimizer
#' trace CSV.  Unlike the "Iteration history" section in
#' \code{\link{ferx_runlog}}, this function never truncates: every iteration
#' row is shown without omission.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} with
#'   \code{optimizer_trace = TRUE}, or a path to a trace CSV file.
#' @param verbose Logical.  When \code{TRUE} (the default) the table is
#'   printed.  When \code{FALSE} it is returned invisibly as a character string.
#'
#' @return A character string (invisibly when \code{verbose = TRUE}).
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "focei",
#'                 covariance = FALSE, optimizer_trace = TRUE)
#' ferx_runlog_iters(fit)
#' }
#'
#' @seealso \code{\link{ferx_runlog}}, \code{\link{ferx_trace}},
#'   \code{\link{ferx_plot_trace}}
#' @family diagnostics
#' @export
ferx_runlog_iters <- function(fit, verbose = TRUE) {
  tr <- if (is.character(fit)) {
    ferx_trace(fit)
  } else if (inherits(fit, "ferx_fit")) {
    if (is.null(fit$trace_path) || is.na(fit$trace_path) ||
        !nzchar(fit$trace_path) || !file.exists(fit$trace_path)) {
      stop("No trace file found. Re-run with optimizer_trace = TRUE.")
    }
    ferx_trace(fit)
  } else {
    stop("`fit` must be a ferx_fit object or a path to a trace CSV")
  }

  hdr  <- c(strrep("=", 60), "ITERATION HISTORY", strrep("-", 60))
  body <- .runlog_iter_table(tr, truncate = FALSE)
  foot <- c(sprintf("  Total: %d iteration(s)", nrow(tr)), strrep("=", 60))
  out  <- paste(c(hdr, body, foot), collapse = "\n")
  if (isTRUE(verbose)) cat(out, "\n")
  invisible(out)
}

# Internal: format a per-iteration trace data frame as a character vector of
# fixed-width table lines (header + rows).  Called by ferx_runlog (truncated)
# and ferx_runlog_iters (full).
.runlog_iter_table <- function(tr, truncate = TRUE,
                               trunc_total = 30L,
                               trunc_head  = 10L,
                               trunc_tail  = 10L) {
  if (is.null(tr) || nrow(tr) == 0L) return("  (no iterations recorded)")
  n <- nrow(tr)

  # Detect method category
  methods  <- unique(tr$method[!is.na(tr$method)])
  is_saem  <- any(methods == "saem")
  is_gn    <- !is_saem && any(grepl("^gn", methods))

  # OFV proxy column (for SAEM, tr$ofv already holds cond_nll as proxy)
  ofv_main <- tr$ofv

  # Delta-OFV: use ofv_delta for GN when populated, else diff()
  if (is_gn && "ofv_delta" %in% names(tr) && any(!is.na(tr$ofv_delta))) {
    dofv <- tr$ofv_delta
  } else {
    dofv <- c(NA_real_, diff(ofv_main))
  }

  # Column visibility
  show_grad  <- !is_saem && "grad_norm"      %in% names(tr) && any(!is.na(tr$grad_norm))
  show_step  <- !is_saem && "step_norm"      %in% names(tr) && any(!is.na(tr$step_norm))
  show_lm    <- is_gn    && "lm_lambda"      %in% names(tr) && any(!is.na(tr$lm_lambda))
  show_acc   <- is_gn    && "step_accepted"  %in% names(tr) && any(!is.na(tr$step_accepted))
  show_mh    <- is_saem  && "mh_accept_rate" %in% names(tr) && any(!is.na(tr$mh_accept_rate))
  show_gamma <- is_saem  && "gamma"          %in% names(tr) && any(!is.na(tr$gamma))
  show_ebe   <- "n_ebe_unconverged" %in% names(tr) &&
                any(tr$n_ebe_unconverged > 0L, na.rm = TRUE)
  show_phase <- !is_saem && "phase" %in% names(tr) &&
                any(nzchar(tr$phase) & !is.na(tr$phase))
  show_saem_phase <- is_saem && "phase" %in% names(tr) &&
                     any(nzchar(tr$phase) & !is.na(tr$phase))

  # Header
  ofv_label  <- if (is_saem) "COND_NLL"  else "OFV"
  dofv_label <- if (is_saem) "dCOND_NLL" else "dOFV"
  hdr <- sprintf("  %4s  %13s  %13s", "ITER", ofv_label, dofv_label)
  bar <- sprintf("  %4s  %13s  %13s", "----", "-------------", "-------------")

  add_col <- function(label, w) {
    hdr <<- paste0(hdr, sprintf(paste0("  %", w, "s"), label))
    bar <<- paste0(bar, sprintf(paste0("  %", w, "s"), strrep("-", w)))
  }

  if (!is_saem) {
    if (show_grad) add_col("GRAD_NORM", 11)
    if (show_step) add_col("STEP_NORM", 11)
    if (show_lm)   add_col("LM_LAMBDA", 11)
    if (show_acc)  add_col("ACC",        5)
  } else {
    if (show_saem_phase) add_col("PHASE",     10)
    if (show_gamma)      add_col("GAMMA",      9)
    if (show_mh)         add_col("MH_ACCEPT",  9)
  }
  if (show_ebe) add_col("EBE_WARN", 8)

  # Helpers
  na_w  <- function(w) sprintf(paste0("%", w, "s"), "NA")
  fmt_f <- function(x, w, d) {
    if (is.null(x) || is.na(x) || !is.finite(x)) na_w(w)
    else sprintf(paste0("%", w, ".", d, "g"), x)
  }
  fmt_d <- function(x) {
    if (is.na(x) || !is.finite(x)) sprintf("%13s", "---")
    else sprintf("%+13.6g", x)
  }

  make_row <- function(i) {
    r <- sprintf("  %4d  %13.6g  %s", as.integer(tr$iter[i]), ofv_main[i], fmt_d(dofv[i]))
    if (!is_saem) {
      if (show_grad) r <- paste0(r, fmt_f(tr$grad_norm[i],      13, 4))
      if (show_step) r <- paste0(r, fmt_f(tr$step_norm[i],      13, 4))
      if (show_lm)   r <- paste0(r, fmt_f(tr$lm_lambda[i],      13, 4))
      if (show_acc) {
        v   <- tr$step_accepted[i]
        r   <- paste0(r, if (is.na(v)) sprintf("  %5s", "NA") else
                         if (v == 1)   sprintf("  %5s", "YES")  else sprintf("  %5s", "NO"))
      }
    } else {
      if (show_saem_phase) {
        ph <- tr$phase[i]
        r  <- paste0(r, sprintf("  %10s", if (is.na(ph) || !nzchar(ph)) "" else ph))
      }
      if (show_gamma) r <- paste0(r, fmt_f(tr$gamma[i],           11, 4))
      if (show_mh)    r <- paste0(r, fmt_f(tr$mh_accept_rate[i],  11, 4))
    }
    if (show_ebe) {
      v <- tr$n_ebe_unconverged[i]
      r <- paste0(r, if (is.na(v)) sprintf("  %8s", "NA") else sprintf("  %8d", as.integer(v)))
    }
    r
  }

  # Build rows for a set of indices, inserting phase-change separators
  build_rows <- function(indices) {
    out      <- character(0)
    prev_ph  <- NULL
    for (i in indices) {
      if (show_phase) {
        ph <- tr$phase[i]
        if (is.na(ph)) ph <- ""
        if (!is.null(prev_ph) && nzchar(ph) && ph != prev_ph) {
          label <- paste0("-- phase: ", ph, " ")
          out   <- c(out, sprintf("  %s%s", label, strrep("-", max(0L, 40L - nchar(label)))))
        }
        prev_ph <- ph
      }
      out <- c(out, make_row(i))
    }
    out
  }

  # Apply truncation
  if (truncate && n > trunc_total) {
    idx_head  <- seq_len(trunc_head)
    idx_tail  <- seq(n - trunc_tail + 1L, n)
    n_omitted <- n - trunc_head - trunc_tail
    rows <- c(
      build_rows(idx_head),
      sprintf("  ... %d rows not shown (use ferx_runlog_iters() for full table) ...", n_omitted),
      build_rows(idx_tail)
    )
  } else {
    rows <- build_rows(seq_len(n))
  }

  c(hdr, bar, rows)
}
