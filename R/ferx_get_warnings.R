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
#' ferx_get_warnings(fit)
#' ferx_get_warnings(fit, as_df = TRUE)
#' @family diagnostics
#' @export
ferx_get_warnings <- function(fit, as_df = FALSE) {
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

# Guidance for the `sir` category. The proposal-conditioning diagnostics
# (ferx-core#1021) are about the *model*, not about SIR tuning: a direction the
# covariance step had to floor is a direction the data do not inform, and no
# amount of extra samples will recover it. Everything else in the category is a
# SIR-availability or tuning issue, which keeps the original advice.
.ferx_sir_guidance <- function(message = "") {
  if (grepl("rank-deficient", message, ignore.case = TRUE)) {
    return(paste0(
      "The covariance matrix has a direction with no uncertainty left after ",
      "FIXed parameters are excluded - the named parameters are not identified ",
      "by the data (they trade off against each other). SIR holds that ",
      "combination at its ML values, so its CIs are not explored. Fix or drop ",
      "one parameter from each named combination, or re-fit a model the data ",
      "can identify."
    ))
  }
  if (grepl("shrunk", message, ignore.case = TRUE)) {
    return(paste0(
      "The proposal was an order of magnitude wider than the room the named ",
      "parameters have between their estimates and their bounds, which means the ",
      "covariance step floored their curvature: they are effectively ",
      "non-identified. SIR shrank those directions so the run could proceed, but ",
      "the CIs along them understate the true uncertainty - treat them as a lower ",
      "bound, and check the covariance-step warning for the same parameters."
    ))
  }
  paste0(
    "SIR uncertainty step issue. Ensure covariance = TRUE and inspect SIR ",
    "tuning (sir_samples / sir_resamples)."
  )
}

# Remediation guidance keyed by the fixed category vocabulary that ferx-core
# (and the R-side additions) emit. Extend this table when core grows a new
# category. Returns NULL for unknown categories so callers can skip printing
# rather than showing a generic placeholder.
#
# The vocabulary is ferx-core's `WarningCode::as_str()` token set; it is pinned
# by the "covers ferx-core's warning vocabulary" test in
# tests/testthat/test-ferx_get_warnings.R, so a core release that adds a code
# fails there rather than silently printing nothing. `general` is deliberately
# unhandled: it is core's fallback bucket for an unrecognised message, where the
# message text is the only guidance there is.
#
# message and uses_sde are used for dw_autocorrelation (positive vs negative
# DW; SDE hint suppressed when already in use) and for the covariance family
# (the specific failure mode embedded in the message text selects targeted
# advice).
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
  # The covariance family. All three of core's covariance codes land here,
  # because the branch that applies is selected by the *message*, not by the
  # code: core routes "Covariance step failed: ..." and "Covariance step:
  # Hessian is not positive definite" to `covariance_failed`, "Covariance step
  # regularized: ..." to `covariance_regularized`, and only the informational
  # N^2-OFV cost note to `covariance_step`.
  #
  # Gating this block on `covariance_step` alone (as it was until this change)
  # made every targeted branch below unreachable: the messages they grep for
  # never arrive under that code. The unit tests did not catch it because they
  # called this function with the category and the message paired by hand, a
  # pairing production never produces.
  if (category %in% c("covariance_failed", "covariance_regularized",
                      "covariance_step")) {
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
    # `covariance_step` is core's Info-level note about the cost of the step
    # (the N^2 OFV evaluations it will run), not a failure. It reaches none of
    # the branches above, so without its own arm it would inherit the
    # "standard errors unavailable" fallback below and tell the user their
    # covariance step had failed when it had not even started.
    if (category == "covariance_step") {
      return(paste0(
        "Informational: the covariance step cost scales with the square of the ",
        "parameter count. No action needed; pass covariance = FALSE to skip it ",
        "during development."
      ))
    }
    # Generic fallback for older or unrecognised covariance failure messages.
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
    bloq_method        = "LOQ censoring note. Set method = \"focei\" explicitly to silence, or review the M3 setup.",
    sir                = .ferx_sir_guidance(message),
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
    eta_shrinkage      = "High ETA shrinkage means the data do not inform that random effect per subject. Individual estimates and any EBE-based diagnostic (ETA plots, covariate screening) are unreliable for it; consider removing the ETA or collecting richer per-subject data.",
    eps_shrinkage      = "High EPS shrinkage means IWRES is pulled toward zero, so residual-based diagnostics understate misfit. Prefer PRED-based or simulation-based checks (VPC/NPDE) for this fit.",
    boundary_estimate  = "A parameter converged onto a bound. The estimate is not a true optimum: widen the bound if the value is plausible, or reconsider the parameterisation if it is not.",
    high_correlation   = "Two parameters are nearly collinear at convergence, so their individual estimates are poorly determined even when the fit looks good. Consider fixing one, removing it, or reparameterising.",
    inflated_rse       = "A relative standard error is implausibly large, which usually signals a parameter the data barely inform. Check identifiability before quoting its precision.",
    flat_parameter     = "A THETA has no effect on the objective at the initial estimate (unmapped, or dropped from the structural/scaling model). It was frozen at its initial value so the rest of the fit could proceed; map it into the model or remove it. Its reported estimate is just its initial value.",
    flip_flop          = "A transit / inverse-Gaussian absorption model entered the flip-flop regime (absorption no faster than elimination). If the model carries an ODE twin the affected subjects were rerouted and the fit is sound; if it does not, that subject's likelihood contribution is degenerate. Check the named subject's MTT/CL (transit) or MAT/CV2/CL (IG) estimates.",
    absorption_twin_declined = "The analytic absorption model kept no ODE fallback, so the fit is a plain closed-form fit but every feature that needs the fallback (time-varying covariates, a TIME-dependent parameter, IOV, steady-state or infusion doses, the flip-flop reroute) is now rejected instead. Read the reason quoted in the message: it is most often an individual parameter named after one of the twin's own compartments (CENTRAL, PERIPH). Renaming it restores the twin.",
    experimental       = "An experimental feature is in use (SDE or neural-network components). Results are usable but should be applied with caution.",
    simulation         = "A simulated subject was handled specially -- a degenerate hazard draw, or an over-large recurrent-event stream that was truncated. The estimated model is unaffected; inspect the named subject's hazard parameters and covariate values.",
    NULL
  )
}
