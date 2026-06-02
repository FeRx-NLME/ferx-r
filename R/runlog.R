#' NONMEM-style run log for a ferx fit
#'
#' Produces a plain-text summary modelled on the information content of a
#' NONMEM \code{.lst} file: model source, data summary, initial vs. final
#' parameter table, and (where available) the final gradient vector with a
#' convergence diagnostic.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} or
#'   \code{\link{ferx_load_fit}}.
#' @param gradient_tol Numeric scalar.  Components of \code{final_gradient}
#'   whose absolute value exceeds this threshold trigger a convergence warning.
#'   Default \code{0.01} (matches typical NONMEM practice).
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
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' ferx_runlog(fit)
#' }
#'
#' @export
ferx_runlog <- function(fit, gradient_tol = 0.01, verbose = TRUE) {
  lines <- character(0)

  .sec <- function(title) {
    c(strrep("=", 60), toupper(title), strrep("-", 60))
  }
  .blank <- function() ""

  # -- Model file -------------------------------------------------------------
  model_text <- fit$model_text
  if (!is.null(model_text) && nzchar(model_text)) {
    lines <- c(lines, .sec("Model file"), model_text, .blank())
  }

  # -- Data summary -----------------------------------------------------------
  lines <- c(lines, .sec("Data summary"))
  n_subj <- fit$n_subjects %||% 0L
  n_obs  <- fit$n_obs      %||% 0L
  tr     <- fit$obs_time_range

  if (!is.null(tr) && length(tr) == 2L && all(is.finite(tr))) {
    lines <- c(lines, sprintf(
      "  %d subjects   %d observations   time range %.4g to %.4g",
      n_subj, n_obs, tr[1], tr[2]
    ))
  } else {
    lines <- c(lines, sprintf("  %d subjects   %d observations", n_subj, n_obs))
  }
  lines <- c(lines, .blank())

  # -- Parameter table: INITIAL vs FINAL -------------------------------------
  has_inits <- !is.null(fit$theta_init) && length(fit$theta_init) > 0L

  lines <- c(lines, .sec("Parameter estimates"))

  if (has_inits) {
    lines <- c(lines, sprintf("  %-28s  %14s  %14s", "PARAMETER", "INITIAL", "FINAL"))
  } else {
    lines <- c(lines, sprintf("  %-28s  %14s", "PARAMETER", "FINAL"))
  }
  lines <- c(lines, sprintf("  %s", strrep("-", if (has_inits) 60 else 46)))

  # Theta rows
  thetas   <- fit$theta       %||% numeric(0)
  th_names <- fit$theta_names %||% paste0("THETA(", seq_along(thetas), ")")
  th_init  <- if (has_inits) fit$theta_init else NULL
  th_fixed <- fit$theta_fixed %||% rep(FALSE, length(thetas))

  for (i in seq_along(thetas)) {
    nm    <- if (i <= length(th_names)) th_names[i] else sprintf("THETA(%d)", i)
    fixed <- isTRUE(th_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", thetas[i])
    if (has_inits) {
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

  for (i in seq_len(n_eta)) {
    nm    <- if (i <= length(eta_nms)) eta_nms[i] else sprintf("OMEGA(%d,%d)", i, i)
    fixed <- isTRUE(om_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", om[i, i])
    if (has_inits && has_om_init) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", om_init[i, i])
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else if (has_inits) {
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, sprintf("%14s", "N/A"), final_str))
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

  for (i in seq_along(sigmas)) {
    nm    <- if (i <= length(sig_nms)) sig_nms[i] else sprintf("SIGMA(%d)", i)
    fixed <- isTRUE(sig_fixed[i])
    final_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", sigmas[i])
    if (has_inits && !is.null(sig_init)) {
      init_str <- if (fixed) sprintf("%14s", "FIXED") else sprintf("%14.6g", sig_init[i])
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, init_str, final_str))
    } else if (has_inits) {
      lines <- c(lines, sprintf("  %-28s  %s  %s", nm, sprintf("%14s", "N/A"), final_str))
    } else {
      lines <- c(lines, sprintf("  %-28s  %s", nm, final_str))
    }
  }
  lines <- c(lines, .blank())

  # -- OFV / AIC / BIC -------------------------------------------------------
  lines <- c(lines, .sec("Objective function"))
  lines <- c(lines,
    sprintf("  OFV   = %.6g", fit$ofv %||% NA_real_),
    sprintf("  AIC   = %.6g", fit$aic %||% NA_real_),
    sprintf("  BIC   = %.6g", fit$bic %||% NA_real_),
    sprintf("  Convergence: %s", if (isTRUE(fit$converged)) "YES" else "NO")
  )
  lines <- c(lines, .blank())

  # -- Final gradient ---------------------------------------------------------
  grad <- fit$final_gradient
  if (!is.null(grad) && length(grad) > 0L) {
    lines <- c(lines, .sec("Final gradient"))
    large <- abs(grad) > gradient_tol
    if (any(large)) {
      lines <- c(lines, sprintf(
        "  WARNING: %d component(s) exceed |gradient| > %.4g -- possible non-convergence.",
        sum(large), gradient_tol
      ))
    } else {
      lines <- c(lines, sprintf(
        "  All %d gradient components satisfy |gradient| <= %.4g.",
        length(grad), gradient_tol
      ))
    }
    # Wrap gradient components at 8 per line
    chunks <- split(grad, ceiling(seq_along(grad) / 8))
    for (chunk in chunks) {
      lines <- c(lines, sprintf("  [ %s ]", paste(sprintf("%.4g", chunk), collapse = "  ")))
    }
    lines <- c(lines, sprintf("  max|gradient| = %.4g", max(abs(grad))))
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

  # -- Assemble and return ----------------------------------------------------
  out <- paste(lines, collapse = "\n")
  if (isTRUE(verbose)) cat(out, "\n")
  invisible(out)
}
