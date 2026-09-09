# Format the structural display string from a model_structure list.
# Combines the optional type label with the population parameter (theta) names.
.ferx_format_structural <- function(ms) {
  theta_str <- if (length(ms$theta_names) > 0L)
    paste(ms$theta_names, collapse = ", ")
  else
    NULL
  if (!is.null(ms$model_type) && !is.null(theta_str))
    sprintf("%s  (%s)", ms$model_type, theta_str)
  else if (!is.null(ms$model_type))
    ms$model_type
  else if (!is.null(theta_str))
    theta_str
  else
    "unknown"
}

# Normalise the sample-size-weighted IOV vectors (ferx-core #1031) that come
# off the FFI list or a .fitrx bundle. The engine leaves both empty unless some
# kappa declares `weight = <expr>`, so an ordinary IOV model keeps NULL here.
# When present they are padded/truncated to one entry per kappa and named by
# kappa_names, with NA marking an unweighted kappa in a model where some other
# kappa is weighted. Returns list(kappa_weights = , kappa_weight_typical = ).
.ferx_name_kappa_weights <- function(weights, typical, kappa_names, n_kappa) {
  none <- list(kappa_weights = NULL, kappa_weight_typical = NULL)
  if (is.null(weights) || length(weights) == 0L || n_kappa < 1L) return(none)
  w <- as.character(weights)
  w[!nzchar(w)] <- NA_character_
  if (all(is.na(w))) return(none)
  length(w) <- n_kappa
  tv <- suppressWarnings(as.numeric(typical %||% numeric(0)))
  length(tv) <- n_kappa
  if (!is.null(kappa_names) && length(kappa_names) == n_kappa) {
    names(w) <- kappa_names
    names(tv) <- kappa_names
  }
  list(kappa_weights = w, kappa_weight_typical = tv)
}

# Format the one-line weight annotation printed under a weighted kappa's
# estimate. `var` is the *unweighted* gamma^2 the engine reports; the number a
# reader needs next to it is the effective between-occasion SD at a typical
# weight, gamma / sqrt(W). Mirrors the ferx-core CLI output. Returns NULL when
# kappa `i` carries no weight.
.ferx_format_kappa_weight <- function(fit, i, var, name) {
  w <- fit$kappa_weights
  if (is.null(w) || length(w) < i || is.na(w[[i]])) return(NULL)
  expr <- as.character(w[[i]])
  tv <- fit$kappa_weight_typical
  n <- if (!is.null(tv) && length(tv) >= i) tv[[i]] else NA_real_
  if (!is.na(n) && is.finite(n) && n > 0 && !is.na(var) && var >= 0) {
    sprintf("%22s weight = %s  ->  SD = %.4f at %s = %.4f (kappa ~ N(0, %s/%s))",
            "", expr, sqrt(var) / sqrt(n), expr, n, name, expr)
  } else {
    sprintf("%22s weight = %s (kappa ~ N(0, %s/%s))", "", expr, name, expr)
  }
}

# Kappa labels for a model_structure list, with a sample-size-weighted kappa
# (ferx-core #1031) annotated by its weight expression: `KAPPA_EMAX (weight =
# NARM)`. `iov_weights` is absent (or all NA) for every model that declares no
# weight. `model_structure` is persisted verbatim under r_extras and read back
# with `simplifyVector = FALSE`, so after a save/load round-trip both vectors
# arrive as *lists* whose NA slots are NULL holes - `is.na()` is FALSE for
# those and `as.character(NULL)` is "NULL", which would render
# `(weight = NULL)`. Flatten first.
#
# Every surface that reports IOV must go through this, so print.ferx_fit(),
# print.ferx_summary() and ferx_model_inspect() cannot disagree about whether
# a model is weighted. Returns character(0) for a model with no IOV.
.ferx_iov_labels <- function(ms) {
  iov_lbl <- as.character(.fitrx_unwrap_opt_chr_vec(ms$iov) %||% character())
  wts <- as.character(.fitrx_unwrap_opt_chr_vec(ms$iov_weights) %||% character())
  if (length(iov_lbl) > 0L && length(wts) == length(iov_lbl)) {
    has_w <- !is.na(wts) & nzchar(wts)
    iov_lbl[has_w] <- sprintf("%s (weight = %s)", iov_lbl[has_w], wts[has_w])
  }
  iov_lbl
}

# Print the four structural detail lines (Structural / IIV / IOV / Residual).
# The caller is responsible for any preceding header line.
.ferx_print_structure <- function(ms) {
  cat(sprintf("  Structural:  %s\n", .ferx_format_structural(ms)))
  cat(sprintf("  IIV:         %s\n",
    if (length(ms$iiv) > 0L) paste(ms$iiv, collapse = ", ") else "none"))
  iov_lbl <- .ferx_iov_labels(ms)
  cat(sprintf("  IOV:         %s\n",
    if (length(iov_lbl) > 0L) paste(iov_lbl, collapse = ", ") else "none"))
  cat(sprintf("  Residual:    %s\n", ms$residual))
  invisible(NULL)
}

# Assemble the structured-warning data frame for a fit result. Core supplies
# the severity/category/message triples it classified (parallel vectors on the
# raw FFI list); the R side appends only diagnostics it computes itself and
# that core never sees: the condition-number flag and the per-ETA Shapiro-Wilk
# normality flags. Returns a data frame with columns severity, category,
# message - always a data frame (zero rows when there are no warnings).
.ferx_assemble_structured_warnings <- function(raw, result) {
  df <- data.frame(
    severity      = as.character(raw$warnings_severity      %||% character(0)),
    category      = as.character(raw$warnings_category      %||% character(0)),
    message       = as.character(raw$warnings_message       %||% character(0)),
    source_method = as.character(raw$warnings_source_method %||% character(0)),
    stringsAsFactors = FALSE
  )
  extra <- list()
  # Condition number (NONMEM convention: > 1000 flags ill-conditioning).
  if (!is.null(result$condition_number) && is.finite(result$condition_number) &&
        result$condition_number > 1000) {
    extra[[length(extra) + 1L]] <- data.frame(
      severity      = "critical",
      category      = "condition_number",
      message       = sprintf(
        "High condition number (%.1f) -- parameter space may be ill-conditioned",
        result$condition_number
      ),
      source_method = "",
      stringsAsFactors = FALSE
    )
  }
  # Per-ETA Shapiro-Wilk normality flags (computed R-side from the EBEs),
  # folded into a single warning that lists every flagged ETA (ferx-core#163).
  nn_msg <- .ferx_eta_normality_warning(result$eta_normality)
  if (!is.null(nn_msg)) {
    extra[[length(extra) + 1L]] <- data.frame(
      severity      = "warning",
      category      = "eta_normality",
      message       = nn_msg,
      source_method = "",
      stringsAsFactors = FALSE
    )
  }
  if (length(extra) > 0L) df <- rbind(df, do.call(rbind, extra))
  df
}

# Fold the per-ETA Shapiro-Wilk flags into a single warning message that lists
# every ETA flagged as possibly non-normal, with its p-value. Returns NULL when
# the normality frame is missing/empty or no ETA is flagged. Firing once instead
# of once per ETA keeps the warnings panel readable (ferx-core#163). The
# explanatory hint (high shrinkage / sparse data, prefer QQ-plots) is attached
# by the eta_normality category guidance in ferx_get_warnings().
.ferx_eta_normality_warning <- function(eta_normality) {
  if (is.null(eta_normality) || !is.data.frame(eta_normality) ||
        nrow(eta_normality) == 0L) {
    return(NULL)
  }
  flagged <- eta_normality[
    !is.na(eta_normality$flag) & nzchar(eta_normality$flag) &
      !is.na(eta_normality$p_val), ,
    drop = FALSE
  ]
  if (nrow(flagged) == 0L) {
    return(NULL)
  }
  parts <- sprintf("%s (p=%.4f)", flagged$eta, flagged$p_val)
  noun <- if (nrow(flagged) == 1L) "ETA" else "ETAs"
  sprintf(
    "Shapiro-Wilk flags possible non-normal distribution for %d %s: %s",
    nrow(flagged), noun, paste(parts, collapse = ", ")
  )
}

# Internal: look up SE for omega element (i, j) from se_omega vector.
# se_omega may be diagonal-only (length n_eta) or full lower-triangle
# (length n_eta*(n_eta+1)/2, column-major).
.omega_se_at <- function(se_omega, n_eta, i, j) {
  if (is.null(se_omega)) return(NA_real_)
  # Ensure r >= c (symmetric)
  r <- max(i, j); c <- min(i, j)
  n_lt <- n_eta * (n_eta + 1L) / 2L
  if (length(se_omega) == n_lt && n_lt != n_eta) {
    # Full lower-triangle (block omega)
    col_offset <- if (c == 1L) 0L else (c - 1L) * n_eta - (c - 1L) * (c - 2L) / 2L
    idx <- col_offset + (r - c) + 1L  # 1-based
    if (idx >= 1L && idx <= length(se_omega)) se_omega[idx] else NA_real_
  } else {
    # Diagonal-only
    if (r == c && r <= length(se_omega)) se_omega[r] else NA_real_
  }
}

.dw_label <- function(dw) {
  if (dw < 1.5) "positive autocorrelation"
  else if (dw > 2.5) "negative autocorrelation"
  else "no autocorrelation"
}

# Resolve once whether cli colour output is available. Called at the top of
# ferx_get_warnings() and from print.ferx_fit() so the capability check is not
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

.ferx_inv_logit <- function(x) 1 / (1 + exp(-x))

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

# -- The raw FFI list -> `ferx_fit` conversion --------------------------------
#
# Extracted from `ferx_fit()` so there is exactly one construction path for a
# fresh fit. The search tools (`ferx_covsearch()`, `ferx_allometry()`) receive
# the engine's winning `FitResult` as the same raw list and must present it as
# the object a user already knows; re-deriving any of this beside them would be
# a second implementation to keep in step.
#
# `gradient_arg`, `settings_used` and `model_file_opts` describe the *call* that
# produced the fit. A search candidate was fitted by the engine's own options,
# so those default to "nothing stated here".
#
# @param raw The list `ferx_rust_fit()` (or a tool binding) returned.
# @param model,data Paths to the model and dataset, for the name / path /
#   source fields the engine does not carry.
# @return A `ferx_fit` object. Saving it (`output=`) stays with the caller.
.ferx_fit_from_raw <- function(raw,
                               model,
                               data,
                               gradient_arg = "",
                               settings_used = list(),
                               model_file_opts = list(),
                               fit_started_at = Sys.time()) {
  # Structure the result
  result <- raw

  # Name the theta vector
  if (length(result$theta) > 0 && length(result$theta_names) > 0) {
    names(result$theta) <- result$theta_names
  }

  # Reshape omega into a matrix. A model with no random effects (n_eta = 0 - e.g.
  # a fixed-effects `[binary_model]`) returns an empty omega; keep it a 0x0 matrix
  # rather than a bare numeric(0), so downstream `nrow(result$omega)` is 0, not
  # NULL (the latter poisons the eta-metadata guards below with NA - #271).
  if (length(result$omega) > 0 && !is.null(result$omega_dim)) {
    d <- result$omega_dim
    result$omega <- matrix(result$omega, nrow = d, ncol = d)
    eta_nms_om <- result$eta_names
    omega_dim_nms <- if (!is.null(eta_nms_om) && length(eta_nms_om) == d) eta_nms_om else paste0("OMEGA(", seq_len(d), ",", seq_len(d), ")")
    rownames(result$omega) <- colnames(result$omega) <- omega_dim_nms
  } else {
    result$omega <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  }

  # Name SE vectors
  if (length(result$se_theta) > 0 && length(result$theta_names) > 0) {
    names(result$se_theta) <- result$theta_names
  }
  if (length(result$se_theta) == 0) result$se_theta <- NULL
  if (length(result$se_omega) == 0) result$se_omega <- NULL
  if (length(result$se_sigma) == 0) result$se_sigma <- NULL

  # SIR: NaN ess => not computed; flat [lo, hi, lo, hi, ...] => (n, 2) matrix
  if (is.null(result$sir_ess) || !is.finite(result$sir_ess)) {
    result$sir_ess <- NULL
  }
  reshape_ci <- function(v, row_names = NULL) {
    if (length(v) == 0) {
      return(NULL)
    }
    m <- matrix(v,
      ncol = 2, byrow = TRUE,
      dimnames = list(row_names, c("lower", "upper"))
    )
    m
  }
  result$sir_ci_theta <- reshape_ci(result$sir_ci_theta, result$theta_names)
  n_eta <- if (is.null(dim(result$omega))) NULL else nrow(result$omega)
  eta_names <- if (!is.null(n_eta)) {
    en <- result$eta_names
    if (!is.null(en) && length(en) == n_eta) en else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
  } else NULL
  result$sir_ci_omega <- reshape_ci(result$sir_ci_omega, eta_names)
  sn <- result$sigma_names
  sig_names <- if (!is.null(sn) && length(sn) == length(result$sigma)) sn else paste0("SIGMA(", seq_along(result$sigma), ")")
  result$sir_ci_sigma <- reshape_ci(result$sir_ci_sigma, sig_names)

  # Normalize trace_path: NULL/empty means no trace was written
  tp <- result$trace_path
  if (is.null(tp) || length(tp) == 0L || !nzchar(tp[[1L]])) {
    result$trace_path <- NULL
    result$trace <- NULL
  } else {
    .ferx_state$last_trace_path <- result$trace_path
    .ferx_state$last_trace_time <- fit_started_at
    .ferx_state$last_trace_model <- model
    # Read the trace CSV into the fit object itself, not just its path, so
    # it survives a save/load round-trip and callers don't need to re-read
    # a temp file that may since have been deleted (see ferx_trace()).
    result$trace <- tryCatch(.ferx_read_trace_csv(result$trace_path),
                              error = function(e) NULL)
  }

  # Normalize shrinkage: NaN ? NA (consistent with other optional numerics)
  if (!is.null(result$shrinkage_eta)) {
    result$shrinkage_eta[!is.finite(result$shrinkage_eta)] <- NA_real_
  }
  if (!is.null(result$shrinkage_eps) && !is.finite(result$shrinkage_eps)) {
    result$shrinkage_eps <- NA_real_
  }
  if (!is.null(result$dw_statistic) && !is.finite(result$dw_statistic)) {
    result$dw_statistic <- NA_real_
  }
  if (!is.null(result$iwres_lag1_r) && !is.finite(result$iwres_lag1_r)) {
    result$iwres_lag1_r <- NA_real_
  }

  # Normalize covariance_status: missing from older binaries ? "not_requested"
  if (is.null(result$covariance_status)) {
    result$covariance_status <- "not_requested"
  }

  # Model-selection surface (ferx-core #1177). `max_abs_correlation` arrives as
  # NaN when the fit has no covariance matrix to read one off, and the two
  # tri-state verdicts arrive absent when the engine recorded none. Normalize
  # both to NA so check_strictness() meets a single "no input" shape whether
  # the fit came from here or from ferx_load_fit().
  if (!is.null(result$max_abs_correlation) &&
        !is.finite(result$max_abs_correlation)) {
    result$max_abs_correlation <- NA_real_
  }
  result$left_init <- .fitrx_unwrap_opt_lgl(result$left_init)
  result$stalled_at_init <- .fitrx_unwrap_opt_lgl(result$stalled_at_init)
  # The packed Omega / kappa layout, same tri-state shape. ferx_save_fit()
  # writes these so a reader can undo the Cholesky parameterisation before
  # reading correlations off the covariance matrix.
  result$omega_is_diagonal <- .fitrx_unwrap_opt_lgl(result$omega_is_diagonal)
  result$kappa_is_diagonal <- .fitrx_unwrap_opt_lgl(result$kappa_is_diagonal)

  # Fitted `block_sigma` residual correlations as a tidy frame, or NULL when
  # the model declares none. The FFI ships them as parallel vectors.
  result$residual_correlations <- .ferx_residual_corr_frame(result)
  for (k in c("residual_correlation_i", "residual_correlation_j",
              "residual_correlation_rho", "residual_correlation_fixed",
              "residual_correlation_names", "se_residual_correlation")) {
    result[[k]] <- NULL
  }

  # Reshape cov_matrix into a named square matrix (param ? param)
  d <- result$cov_matrix_dim %||% 0L
  if (!is.null(result$cov_matrix) && length(result$cov_matrix) > 0L && d > 0L) {
    m <- matrix(result$cov_matrix, nrow = d, ncol = d, byrow = TRUE)
    n_theta <- length(result$theta_names)
    n_eta <- result$omega_dim %||% 0L
    n_sigma <- length(result$sigma)
    # The engine packs the `block_sigma` correlations *last* (after sigma), so
    # they have to come out of the count before it can be read as omega -
    # otherwise every trailing coordinate is off by the number of rho's and the
    # sigma rows carry the wrong labels.
    rho_nms <- .ferx_residual_corr_labels(result)
    n_rho <- length(rho_nms)
    n_omega_packed <- d - n_theta - n_sigma - n_rho
    # Determine parameterisation: diagonal (n_omega_packed == n_eta) or block
    eta_nms <- if (!is.null(result$eta_names) && length(result$eta_names) == n_eta) result$eta_names else NULL
    omega_names <- if (n_omega_packed == n_eta) {
      if (!is.null(eta_nms)) eta_nms
      else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
    } else {
      # Block lower-triangle: L(i,j) for i >= j, column-major
      nm <- character(n_omega_packed)
      k <- 0L
      for (i in seq_len(n_eta)) {
        for (j in seq_len(i)) {
          k <- k + 1L
          nm[k] <- if (!is.null(eta_nms)) sprintf("%s,%s", eta_nms[i], eta_nms[j]) else sprintf("OMEGA(%d,%d)", i, j)
        }
      }
      nm
    }
    sig_nms <- if (!is.null(result$sigma_names) && length(result$sigma_names) == n_sigma) result$sigma_names else NULL
    pnames <- c(
      result$theta_names,
      if (n_omega_packed > 0L) omega_names else character(0L),
      if (n_sigma > 0L) (if (!is.null(sig_nms)) sig_nms else paste0("SIGMA(", seq_len(n_sigma), ")")) else character(0L),
      rho_nms
    )
    if (length(pnames) == d) rownames(m) <- colnames(m) <- pnames
    result$cov_matrix <- m
  } else {
    result$cov_matrix <- NULL
  }
  result$cov_matrix_dim <- NULL

  result <- .ferx_apply_cov_sentinels(result)

  # ETA normality (Shapiro-Wilk) - computed in R from per-subject EBEs.
  # Fold every flagged ETA into a single warning rather than one per ETA
  # (ferx-core#163).
  result$eta_normality <- .ferx_compute_eta_normality(result$ebe_etas)
  nn_msg <- .ferx_eta_normality_warning(result$eta_normality)
  if (!is.null(nn_msg)) {
    result$warnings <- c(result$warnings, nn_msg)
  }

  # Reshape omega_iov into a named matrix (NULL when no IOV)
  d_iov <- result$omega_iov_dim %||% 0L
  if (!is.null(result$omega_iov) && length(result$omega_iov) > 0L && d_iov > 0L) {
    m_iov <- matrix(result$omega_iov, nrow = d_iov, ncol = d_iov)
    if (length(result$kappa_names) == d_iov) {
      rownames(m_iov) <- colnames(m_iov) <- result$kappa_names
    }
    result$omega_iov <- m_iov
    if (length(result$kappa_names) > 0L && length(result$shrinkage_kappa) == length(result$kappa_names)) {
      names(result$shrinkage_kappa) <- result$kappa_names
    }
    if (length(result$se_kappa) == 0L) {
      result$se_kappa <- NULL
    } else {
      n_tri <- d_iov * (d_iov + 1L) / 2L
      if (length(result$se_kappa) == d_iov) {
        names(result$se_kappa) <- result$kappa_names
      } else if (length(result$se_kappa) == n_tri) {
        # block kappa: label lower-triangle elements as NAME (diagonal) or COV_i_j
        tri_names <- character(n_tri)
        idx <- 1L
        for (j in seq_len(d_iov)) {
          for (i in j:d_iov) {
            tri_names[idx] <- if (i == j) {
              result$kappa_names[i]
            } else {
              paste0("COV_", result$kappa_names[j], "_", result$kappa_names[i])
            }
            idx <- idx + 1L
          }
        }
        names(result$se_kappa) <- tri_names
      }
    }
    # Per-occasion shrinkage: set column names from kappa_names; NULL when
    # the Rust glue returned an empty/NULL frame (no IOV or unbalanced design).
    if (is.data.frame(result$shrinkage_kappa_by_occ) &&
        nrow(result$shrinkage_kappa_by_occ) > 0L) {
      kn <- result$kappa_names
      # Columns are: occ, <kappa1>, <kappa2>, ... already named by Rust glue
      # but guard against older binaries that may omit kappa column names.
      if (!is.null(kn) && length(kn) == ncol(result$shrinkage_kappa_by_occ) - 1L) {
        colnames(result$shrinkage_kappa_by_occ) <- c("occ", kn)
      }
    } else {
      result$shrinkage_kappa_by_occ <- NULL
    }
    # Sample-size-weighted IOV (ferx-core #1031). The engine leaves both vectors
    # empty unless some kappa carries `weight = <expr>`, so NULL them out on the
    # ordinary IOV path and keep them parallel to kappa_names otherwise.
    result[c("kappa_weights", "kappa_weight_typical")] <-
      .ferx_name_kappa_weights(result$kappa_weights, result$kappa_weight_typical,
                               result$kappa_names, d_iov)
  } else {
    result$omega_iov <- NULL
    result$se_kappa <- NULL
    result$shrinkage_kappa <- NULL
    result$shrinkage_kappa_by_occ <- NULL
    result$kappa_names <- NULL
    result$ebe_kappas <- NULL
    result$kappa_weights <- NULL
    result$kappa_weight_typical <- NULL
  }

  # Reshape omega_param_corr into a square matrix using omega_dim to guard
  # against unexpected vector lengths (NULL when diagonal or absent).
  if (!is.null(result$omega_param_corr) && length(result$omega_param_corr) > 0L) {
    d_pc <- result$omega_dim %||% 0L
    if (d_pc > 0L && length(result$omega_param_corr) == d_pc * d_pc) {
      result$omega_param_corr <- matrix(
        result$omega_param_corr, nrow = d_pc, ncol = d_pc, byrow = TRUE
      )
    } else {
      result$omega_param_corr <- NULL
    }
  } else {
    result$omega_param_corr <- NULL
  }

  # Reshape omega_iov_param_corr into a square matrix (NULL when diagonal or absent)
  if (!is.null(result$omega_iov_param_corr) && length(result$omega_iov_param_corr) > 0L) {
    d_ipc <- result$omega_iov_dim %||% 0L
    if (d_ipc > 0L && length(result$omega_iov_param_corr) == d_ipc * d_ipc) {
      result$omega_iov_param_corr <- matrix(
        result$omega_iov_param_corr, nrow = d_ipc, ncol = d_ipc, byrow = TRUE
      )
    } else {
      result$omega_iov_param_corr <- NULL
    }
  } else {
    result$omega_iov_param_corr <- NULL
  }

  # eta_log_transformed: empty vector -> NULL
  if (is.null(result$eta_log_transformed) || length(result$eta_log_transformed) == 0L) {
    result$eta_log_transformed <- NULL
  }

  # Parameter transform metadata - fall back gracefully for older ferx-core binaries
  # that don't populate these vectors (empty character() from FFI ? treat as absent).
  n_eta_fit   <- if (!is.null(result$omega)) nrow(result$omega) else 0L
  n_theta_fit <- length(result$theta)
  n_sigma_fit <- length(result$sigma)

  if (is.null(result$eta_param_types) || length(result$eta_param_types) != n_eta_fit) {
    result$eta_param_types <- rep("log_normal", n_eta_fit)
  }
  if (is.null(result$eta_linked_theta) || length(result$eta_linked_theta) != n_eta_fit) {
    result$eta_linked_theta <- rep("", n_eta_fit)
  }
  if (is.null(result$theta_transforms) || length(result$theta_transforms) != n_theta_fit) {
    result$theta_transforms <- rep("identity", n_theta_fit)
  }
  if (is.null(result$sigma_types) || length(result$sigma_types) != n_sigma_fit) {
    result$sigma_types <- rep("proportional", n_sigma_fit)
  }
  names(result$theta_transforms) <- names(result$theta)
  names(result$sigma_types)      <- names(result$sigma)

  # Clean up internal fields
  result$theta_names <- NULL
  result$omega_dim <- NULL
  result$omega_iov_dim <- NULL

  # Print mu-referencing detections as informational messages.
  for (eta_names in .ferx_mu_ref_detection_names(result$warnings)) {
    message("Mu-referencing detected for: ", eta_names)
  }

  # Store the dataset name (basename without extension) from the data path
  result$data_name <- tools::file_path_sans_ext(basename(data))

  # Fall back to the model file's basename when the .ferx file declares no name.
  # ferx-core returns "Unnamed" (not NULL/"") when no [model_name] is declared -
  # treat it the same as absent. See fit_result_to_list() in src/rust/src/lib.rs.
  # Long-term fix: engine should return "" so this sentinel check can be removed.
  if (is.null(result$model_name) || !nzchar(result$model_name) || identical(result$model_name, "Unnamed")) {
    result$model_name <- tools::file_path_sans_ext(basename(model))
  }

  # `result$model_structure` is now populated by the Rust engine in
  # `fit_result_to_list()` from the parsed CompiledModel (ferx-core#49), so it
  # reflects exactly what ferx-core ran. The R-side `.ferx_parse_structure()`
  # remains for the pre-fit `ferx_model_inspect(path)` workflow only.

  # Store the requested gradient method, and the method the engine actually
  # resolved to. The engine reports `gradient_method_inner` as one of
  # Map verbose engine labels to the short tokens used
  # by the `gradient` argument so the two values can be compared directly
  # (e.g. requested "auto" -> used "ad"). When the engine omitted the field
  # (older binaries return ""), `gradient_used` is NA.
  # `gradient` defaults to NULL ("defer to the model file / engine default").
  # Assigning that raw NULL would drop the field entirely (and with it the
  # RUN INFO line), so resolve it to the requested token: the explicit arg when
  # given, otherwise the model-file value, otherwise the engine default "auto".
  result$gradient <- if (nzchar(gradient_arg)) {
    gradient_arg
  } else {
    mf_grad <- if ("gradient" %in% names(model_file_opts)) {
      model_file_opts[["gradient"]]
    } else if ("gradient_method" %in% names(model_file_opts)) {
      model_file_opts[["gradient_method"]]
    } else {
      NULL
    }
    if (!is.null(mf_grad)) tolower(as.character(mf_grad)) else "auto"
  }
  result$gradient_used <- .ferx_short_gradient_label(result$gradient_method_inner)

  # Store effective settings (already serialised from settings_parts above)
  result$call_settings <- settings_used

  # Raw [fit_options] from the model file, for inspection and conflict auditing.
  result$model_file_settings <- if (length(model_file_opts) > 0L)
    as.list(model_file_opts) else list()

  # Stash inputs so ferx_save_fit() can embed `model.ferx` and (optionally)
  # `data.csv` into a portable .fitrx bundle. Read failures here must not
  # break the fit, so wrap in tryCatch.
  result$model_source <- tryCatch(
    paste(readLines(model, warn = FALSE), collapse = "\n"),
    error = function(e) ""
  )

  # Source-file provenance. The Rust binding emits the four fields
  # directly (paths as the caller-supplied strings, hashes via
  # `ferx_core::io::hash::sha256_file`). Always run `normalizePath` on
  # the path fields here, regardless of source: the old R-side behavior
  # was to store absolute paths (`normalizePath(data)`), and downstream
  # code (notably ferx_save_fit, where `fit$data_path` is opened from a
  # potentially-different working directory at save time) depends on
  # that. Falling back to NULL ? NA when the binding returns nothing
  # keeps older binaries working.
  empty_to_null <- function(x) {
    if (is.null(x) || length(x) == 0L || !nzchar(x[[1L]])) NULL else as.character(x)
  }
  normalize_or_na <- function(path_str) {
    if (is.null(path_str)) return(NA_character_)
    tryCatch(normalizePath(path_str, mustWork = FALSE),
             error = function(e) NA_character_)
  }
  result$model_path <- normalize_or_na(
    empty_to_null(result$model_path) %||% model
  )
  result$data_path <- normalize_or_na(
    empty_to_null(result$data_path) %||% data
  )
  result$model_hash <- empty_to_null(result$model_hash) %||% NA_character_
  result$data_hash <- empty_to_null(result$data_hash) %||% NA_character_

  # Assemble the structured-warning table. Core supplies severity/category for
  # every warning it emitted (including Durbin-Watson autocorrelation); the R
  # side appends only the diagnostics it computes itself (condition number,
  # ETA normality). No string re-parsing of core messages happens here.
  result$warnings_structured <- .ferx_assemble_structured_warnings(raw, result)

  # Reconstruct the per-iteration IMPMAP parameter trace as a data.frame.
  # The Rust side passes flat vectors + metadata; NULL when not collected.
  # Belt-and-suspenders: only keep it when the caller actually asked for it
  # (settings= arg or [fit_options] in the model file) - guards against the
  # trace leaking onto the fit object from an intermediate stage of a method
  # chain (e.g. c("impmap", "focei")) even when impmap_trace was never
  # requested for that chain.
  impmap_requested <- .ferx_impmap_trace_requested(settings_used, model_file_opts)
  if (!is.null(result$impmap_trace)) {
    result$impmap_trace <- if (impmap_requested) {
      .reconstruct_impmap_trace(result$impmap_trace)
    } else {
      NULL
    }
  }

  # Derived fields (cor_matrix, estimates, eta_cov) - shared with
  # ferx_load_fit() via .ferx_populate_derived_fields() so the two
  # construction paths can't diverge. Uses result$data_path (normalised
  # above), not the raw `data` argument, so a fresh fit and a loaded fit of
  # the same model resolve the dataset identically.
  result <- .ferx_populate_derived_fields(result)

  class(result) <- "ferx_fit"

  result
}
