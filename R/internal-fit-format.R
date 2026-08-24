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

# Print the four structural detail lines (Structural / IIV / IOV / Residual).
# The caller is responsible for any preceding header line.
.ferx_print_structure <- function(ms) {
  cat(sprintf("  Structural:  %s\n", .ferx_format_structural(ms)))
  cat(sprintf("  IIV:         %s\n",
    if (length(ms$iiv) > 0L) paste(ms$iiv, collapse = ", ") else "none"))
  # Annotate a sample-size-weighted kappa (ferx-core #1031) with its weight
  # expression: `KAPPA_EMAX (weight = NARM)`. `iov_weights` is absent (or all
  # NA) for every model that declares no weight. `model_structure` is persisted
  # verbatim under r_extras and read back with `simplifyVector = FALSE`, so
  # after a save/load round-trip both vectors arrive as *lists* whose NA slots
  # are NULL holes - `is.na()` is FALSE for those and `as.character(NULL)` is
  # "NULL", which would print `(weight = NULL)`. Flatten first.
  iov_lbl <- as.character(.fitrx_unwrap_opt_chr_vec(ms$iov) %||% character())
  wts <- as.character(.fitrx_unwrap_opt_chr_vec(ms$iov_weights) %||% character())
  if (length(iov_lbl) > 0L && length(wts) == length(iov_lbl)) {
    has_w <- !is.na(wts) & nzchar(wts)
    iov_lbl[has_w] <- sprintf("%s (weight = %s)", iov_lbl[has_w], wts[has_w])
  }
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
