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

# Guidance for the covariance family. The branch that applies is selected by the
# *message*, not by the code: ferx-core routes "Covariance step failed: ..." and
# "Covariance step: Hessian is not positive definite" to `covariance_failed`,
# "Covariance step regularized: ..." and the off-diagonal FD stencil note to
# `covariance_regularized`, and its informational cost note to `covariance_step`.
# Split out of `.ferx_warning_guidance()` for the same reason as
# `.ferx_sir_guidance()`: a message-dispatching family reads better as its own
# function than as a 100-line block inside a lookup table.
#
# Every branch keys on the message alone, including the informational one -- the
# same text deserves the same answer whichever of the four tokens carries it.
.ferx_covariance_guidance <- function(message = "") {
  # Omega non-PD or near-singular -- checked before the general non-PD branch
  # because omega messages also contain "not positive definite" and
  # "eigenvalue". Matched on "omega matrix is" rather than on a descriptor:
  # ferx-core writes one of two into the same sentence depending on the sign of
  # the smallest eigenvalue ("not positive definite" / "near-singular"), and
  # pinning the grep to the first sent the second to the generic fallback.
  if (grepl("omega matrix is", message, ignore.case = TRUE)) {
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
  # Cancelled part-way. ferx-core returns this from three points inside the FD
  # and score loops (`COV_CANCELLED_MSG`) and it classifies to `general`, so it
  # only reaches this function because `general` is admitted on the message. It
  # matches nothing below, and the generic fallback would tell a user who
  # cancelled to go and question their model's identifiability -- nothing was
  # diagnosed, the step simply did not finish.
  if (grepl("cancelled before completion", message, ignore.case = TRUE)) {
    return(paste0(
      "The covariance step was cancelled before it finished, so no standard ",
      "errors were produced. Nothing was diagnosed about the model: re-run with ",
      "covariance = TRUE and let the step complete if you need them."
    ))
  }
  # Singular / rank-deficient score cross-product. Specific to
  # covariance_method = "s"; the message names the remedy and the generic
  # fallback below does not mention covariance_method at all.
  if (grepl("score cross-product matrix", message, ignore.case = TRUE)) {
    return(paste0(
      "The per-subject score cross-product could not be inverted, which usually ",
      "means fewer subjects than free parameters or collinear per-subject ",
      "scores. Re-run with covariance_method = \"r\" or \"rsr\", or reduce the ",
      "number of free parameters."
    ))
  }
  # An invalid `fd_hessian_step` argument. Core reports it as a covariance-step
  # failure, but the fix is the argument, not the model -- without this arm it
  # inherits the identifiability fallback below.
  if (grepl("fd_hessian_step must be", message, ignore.case = TRUE)) {
    return(paste0(
      "The fd_hessian_step you passed is not a positive, finite number, so the ",
      "covariance step never ran. Pass a positive value, or omit the argument to ",
      "use the default (e.g. ferx_fit(..., fd_hessian_step = 0.05))."
    ))
  }
  # Off-diagonal FD stencil non-finite. Core emits this on the *success* path:
  # the covariance matrix was produced and the SEs exist, they are just missing
  # their cross-partial terms. It reaches none of the branches above (it does not
  # say "regularized"), so without this arm it inherits the "standard errors
  # unavailable" fallback and contradicts the message it is printed under.
  if (grepl("off-diagonal FD stencil", message, ignore.case = TRUE)) {
    return(paste0(
      "Standard errors were produced. The cross-partial terms for the named ",
      "parameter(s) could not be evaluated and were set to zero, so their SEs may ",
      "be over-optimistic -- too narrow rather than missing. Try a different ",
      "fd_hessian_step (e.g. ferx_fit(..., fd_hessian_step = 0.05)), and ",
      "cross-check the affected parameters with ferx_sir()."
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
  # ferx-core's Info-level note about the cost of the step (the n^2 OFV
  # evaluations it will run), not a failure. Without its own arm it would
  # inherit the "standard errors unavailable" fallback below and report a
  # failure that had not happened. Matched on the message, never on the
  # `covariance_step` code: that code's classification ends in a catch-all
  # (`"covariance step:"` + `"parameters"`), so a future message carrying it
  # need not be the cost note -- and answering an unknown one "no action
  # needed" fails in the dangerous direction. Anything else falls through.
  if (grepl("OFV evaluations", message, fixed = TRUE)) {
    return(paste0(
      "Informational: the covariance step cost scales with the square of the ",
      "parameter count. No action needed; pass covariance = FALSE to skip it ",
      "during development."
    ))
  }
  # Generic fallback for older or unrecognised covariance failure messages.
  # Bare rather than return()-wrapped: it is the last expression in the
  # function now, matching `.ferx_sir_guidance()` above (return_linter).
  paste0(
    "Standard errors unavailable. Check identifiability; try a simpler ",
    "omega/sigma structure or covariance = FALSE for development."
  )
}

# Does this warning belong to the covariance family?
#
# Four codes carry covariance-step messages, not one: ferx-core's
# `covariance_failed` / `covariance_regularized` / `covariance_step`, plus
# `"covariance"`, which `ferx_covariance()` assigns to the engine's flat
# covariance warnings when it folds them into the structured table
# (R/ferx_covariance.R). Gating on `covariance_step` alone (as this did until the
# routing fix) made every targeted branch unreachable; omitting `"covariance"`
# left the whole post-hoc `ferx_covariance()` surface unreachable.
#
# `general` is admitted on the message text alone. It is the category
# `ferx_get_warnings()` stamps on every message of a fit saved before structured
# warnings existed, and core's bucket for a covariance message its classifier did
# not recognise. Since every branch is message-keyed, a covariance-step message
# arriving under `general` can be answered exactly as well as one arriving under
# its own code.
.ferx_is_covariance_warning <- function(category, message = "") {
  if (category %in% c("covariance_failed", "covariance_regularized",
                      "covariance_step", "covariance")) {
    return(TRUE)
  }
  category == "general" && grepl("covariance step", message, ignore.case = TRUE)
}



# Remediation guidance keyed by the fixed category vocabulary that ferx-core
# (and the R-side additions) emit. Extend this table when core grows a new
# category. Returns NULL for unknown categories so callers can skip printing
# rather than showing a generic placeholder.
#
# The vocabulary is ferx-core's `WarningCode::as_str()` token set plus the
# categories the R layer adds. Two tests in
# tests/testthat/test-ferx_get_warnings.R cover it, and only one of them is a
# drift guard:
#
#   - "every category in the table returns guidance" walks a hand-written vector
#     of tokens. It cannot notice a code that core added and nobody transcribed,
#     because the loop never visits a token that is not in the vector.
#   - "the table matches ferx-core's WarningCode vocabulary" reads
#     `WarningCode::as_str()` out of a sibling ferx-core checkout and compares
#     the real token set with this table. That one does fail when core grows a
#     code. It skips when the sibling is absent (ferx-r CI builds the pinned
#     crate, not a checkout), so it guards the machine where the drift is
#     introduced rather than the one that consumes it.
#
# `general` is deliberately absent from the table: it is core's bucket for a
# message it did not recognise, where the message text is the only guidance
# there is. It is still routed into the covariance family below when its message
# is a covariance-step message -- see `.ferx_is_covariance_warning()`.
#
# message and uses_sde are used for dw_autocorrelation (positive vs negative
# DW; SDE hint suppressed when already in use). The covariance family branches
# on the message text too, in `.ferx_covariance_guidance()`.
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
  if (.ferx_is_covariance_warning(category, message)) {
    return(.ferx_covariance_guidance(message))
  }
  # ferx-core has no `unused_parameter` code: both of its unused-declaration
  # messages ("theta 'X' is declared in [parameters] but not referenced in any
  # model expression", and the [individual_parameters] "computed but never
  # used" one) fall through `classify_warning` to `general`. So the arm below
  # was unreachable in exactly the way the covariance block was, and is routed
  # the same way -- by the message.
  if (grepl("declared in [parameters] but not referenced", message, fixed = TRUE) ||
      grepl("computed but never used", message, fixed = TRUE)) {
    category <- "unused_parameter"
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
    unused_parameter   = "A declared parameter is never referenced in any model expression, or an [individual_parameters] assignment is computed and then never used. It cannot affect predictions and will not be meaningfully estimated: remove it, or complete the expression that was meant to use it. A commented-out line or a misspelt name is the usual cause.",
    eta_shrinkage      = "High ETA shrinkage means the data do not inform that random effect per subject. Individual estimates and any EBE-based diagnostic (ETA plots, covariate screening) are unreliable for it; consider removing the ETA or collecting richer per-subject data.",
    eps_shrinkage      = "High EPS shrinkage means IWRES is pulled toward zero, so residual-based diagnostics understate misfit. Prefer PRED-based or simulation-based checks (VPC/NPDE) for this fit.",
    boundary_estimate  = "A parameter converged onto a bound. The estimate is not a true optimum: widen the bound if the value is plausible, or reconsider the parameterisation if it is not.",
    high_correlation   = "Two parameters are nearly collinear at convergence, so their individual estimates are poorly determined even when the fit looks good. Consider fixing one, removing it, or reparameterising.",
    inflated_rse       = "A relative standard error is implausibly large, which usually signals a parameter the data barely inform. Check identifiability before quoting its precision.",
    flat_parameter     = "A THETA has no effect on the objective at the initial estimate (unmapped, or dropped from the structural/scaling model). It was frozen at its initial value so the rest of the fit could proceed; map it into the model or remove it. Its reported estimate is just its initial value.",
    flip_flop          = "A transit / inverse-Gaussian absorption model entered the flip-flop regime (absorption no faster than elimination). If the model carries an ODE twin the affected subjects were rerouted and the fit is sound; if it does not, that subject's likelihood contribution is degenerate. Check the named subject's MTT/CL (transit) or MAT/CV2/CL (IG) estimates.",
    absorption_twin_declined = "The analytic absorption model kept no ODE fallback, so the fit is a plain closed-form fit but every feature that needs the fallback (time-varying covariates, a TIME-dependent parameter, IOV, steady-state or infusion doses, the flip-flop reroute) is now rejected instead. The message quotes the reason the fallback could not be built -- the set of reasons is open-ended, so read that text rather than guessing. Resolving it restores the fallback and the features above.",
    experimental       = "An experimental feature is in use (SDE or neural-network components). Results are usable but should be applied with caution.",
    simulation         = "A simulated subject was handled specially -- a degenerate hazard draw, or an over-large recurrent-event stream that was truncated. The estimated model is unaffected; inspect the named subject's hazard parameters and covariate values.",
    NULL
  )
}
