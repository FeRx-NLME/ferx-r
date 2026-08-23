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

# The ill-conditioned-Hessian message carries a per-parameter cause label, and
# ferx-core offers `fd_hessian_step` for only one of the two: a non-finite FD
# stencil is a differencing problem, a zero diagonal is a flat objective and can
# arise on the analytic R-matrix path where nothing is differenced at all
# (`analytic_cov_hessian` defaults to true). Answer the causes present.
.ferx_ill_conditioned_guidance <- function(message = "") {

  overflow <- grepl("FD stencil non-finite", message, fixed = TRUE)
  flat     <- grepl("zero diagonal", message, fixed = TRUE)
  parts <- character(0)
  if (flat) {
    parts <- c(parts, paste0(
      "For the parameter(s) marked \"zero diagonal -- flat objective\": the ",
      "objective does not curve in that direction at all, so the data cannot ",
      "estimate the parameter. Fix it, remove it, or map it into the model if ",
      "it was meant to be used."
    ))
  }
  if (overflow) {
    parts <- c(parts, paste0(
      "For the parameter(s) marked \"FD stencil non-finite\": the objective ",
      "overflowed at the perturbed point. ferx-core's own advice is to tune ",
      "fd_hessian_step (its docs suggest 0.1, or 1e-3 when finite-difference ",
      "noise is the concern); constraining the parameter to a plausible range ",
      "also helps."
    ))
  }
  if (!length(parts)) {
    parts <- paste0(
      "The named parameter(s) have a flat or non-finite Hessian diagonal. The ",
      "message labels the cause for each one; read it before choosing a remedy."
    )
  }
  return(paste(parts, collapse = " "))
}

# Regularisation tiers. Split out of `.ferx_covariance_guidance()` so that
# function stays within a sane branch count; the tier word is embedded in the
# message and `minor` is the only tier ferx-core emits below `moderate`.
.ferx_regularisation_guidance <- function(message = "") {
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
  paste0(
    "Minor Hessian regularisation: standard errors are likely reliable. A ",
    "small eigenvalue floor was applied; this is common on smooth OFV ",
    "surfaces and is usually benign."
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
.ferx_covariance_guidance <- function(message = "", category = "") {
  # Omega non-PD or near-singular -- checked before the general non-PD branch
  # because omega messages also contain "not positive definite" and
  # "eigenvalue". Entered on "omega matrix is" so both descriptors arrive, then
  # split on which one: ferx-core writes "not positive definite" when the
  # smallest eigenvalue is negative and "near-singular" when it is a tiny
  # positive (covariance.rs), a distinction it makes deliberately. Answering
  # both with one text renamed the failure for half of them.
  if (grepl("omega matrix is", message, ignore.case = TRUE)) {
    if (grepl("not positive definite", message, ignore.case = TRUE)) {
      return(paste0(
        "Omega has a negative eigenvalue at convergence, so it is not a valid ",
        "covariance matrix -- some combination of the ETAs carries a negative ",
        "variance. This is usually a block omega whose correlations were driven ",
        "outside the range a positive-definite matrix allows. Simplify the block ",
        "structure (or go diagonal), drop one of the correlated ETAs, or re-fit ",
        "from different starting values."
      ))
    }
    return(paste0(
      "Omega is near-singular at convergence: a variance, or a direction in a ",
      "block, has collapsed towards zero, so the data do not support that much ",
      "random-effect structure. Consider a diagonal omega, fixing the small ",
      "variance to a small positive constant, or removing the corresponding ETA."
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
  # Non-finite or zero Hessian diagonal. ferx-core labels each named parameter
  # with one of two causes inside this single message (covariance.rs), and they
  # have different remedies -- it offers fd_hessian_step for the first only, and
  # the second can arise on the analytic R-matrix path where nothing is
  # differenced at all (`analytic_cov_hessian` defaults to true). So answer the
  # causes actually present rather than merging them.
  if (grepl("ill-conditioned entries", message, ignore.case = TRUE)) {
    return(.ferx_ill_conditioned_guidance(message))
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
      "errors were produced. Nothing was diagnosed about the model -- re-run and ",
      "let the step complete if you need standard errors."
    ))
  }
  # Singular / rank-deficient score cross-product. Specific to
  # covariance_method = "s"; the message names the remedy and the generic
  # fallback below does not mention covariance_method at all.
  if (grepl("score cross-product matrix", message, ignore.case = TRUE)) {
    return(paste0(
      "The per-subject score cross-product could not be inverted, which usually ",
      "means fewer subjects than free parameters or collinear per-subject ",
      "scores. ferx-core's advice in the message is to use covariance_method = ",
      "\"r\" or \"rsr\" instead: pass it as ferx_fit(..., settings = list( ",
      "covariance_method = \"rsr\")), or re-run just the step with ",
      "ferx_covariance(fit, covariance_method = \"rsr\"). Reducing the number of ",
      "free parameters also works."
    ))
  }
  # An invalid `fd_hessian_step` argument. Core reports it as a covariance-step
  # failure, but the fix is the argument, not the model -- without this arm it
  # inherits the identifiability fallback below.
  if (grepl("fd_hessian_step must be", message, ignore.case = TRUE)) {
    return(paste0(
      "The fd_hessian_step you passed is not a positive, finite number, so the ",
      "covariance step never ran. Pass a positive value, or drop the argument ",
      "entirely to use ferx-core's default of 1e-2. (ferx_fit() rejects a ",
      "non-positive value before the engine sees it, so this message normally ",
      "reaches you only from a fit produced outside the R API.)"
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
      "fd_hessian_step -- ferx-core's docs suggest 0.1 when its automatic ",
      "halving was not enough, or 1e-3 when finite-difference noise is the ",
      "concern -- and cross-check the affected parameters with ferx_sir()."
    ))
  }
  # Regularisation path -- severity is embedded in the message.
  if (grepl("covariance step regularized", message, ignore.case = TRUE)) {
    return(.ferx_regularisation_guidance(message))
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
  # Fallback for an unrecognised covariance message. This is the one branch
  # where the code carries something the message does not: `covariance_step` is
  # Info by construction and `covariance_regularized` means the step SUCCEEDED
  # and was merely degraded (types.rs), so telling either that standard errors
  # are unavailable would contradict the row it prints under -- the same
  # invariant the message inventory asserts for the messages we do recognise.
  # Bare rather than return()-wrapped: last expression in the function
  # (return_linter).
  if (category %in% c("covariance_regularized", "covariance_step")) {
    paste0(
      "The covariance step produced a matrix but reported something this ",
      "version does not have specific advice for. Read the message: standard ",
      "errors exist, so treat it as a caveat on their precision rather than as ",
      "a failure, and cross-check with ferx_sir() if it matters."
    )
  } else {
    paste0(
      "Standard errors unavailable. Check identifiability; try a simpler ",
      "omega/sigma structure or covariance = FALSE for development."
    )
  }
}

# --- covariance-family admission -------------------------------------------
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
# Recover the category a `general` row would have carried, from its message.
#
# `general` is core's bucket for a message its classifier did not recognise, and
# -- because `ferx_load_fit()` does not restore the structured table -- it is
# also what every row of a loaded fit arrives as. Two families are recoverable
# from the text alone. Order matters: ferx-core's flat-theta message contains
# "computed but never used", so it must be tested before the unused-declaration
# patterns or it is captured by them.
#
# Returns NULL for a row that is genuinely unrecognised, or for any category
# other than `general` (a live fit's own typed code always wins).
.ferx_resolve_general <- function(category, message = "") {
  if (!identical(category, "general")) return(NULL)
  if (grepl("has no effect on the objective", message, fixed = TRUE)) {
    return("flat_parameter")
  }
  if (grepl("declared in [parameters] but not referenced", message, fixed = TRUE) ||
      grepl("computed but never used", message, fixed = TRUE)) {
    return("unused_parameter")
  }
  NULL
}

.ferx_is_covariance_warning <- function(category, message = "") {
  if (category %in% c("covariance_failed", "covariance_regularized",
                      "covariance_step", "covariance")) {
    return(TRUE)
  }
  # Anchored at the start of the message, optionally after a `[METHOD]` chain
  # prefix, because every covariance warning ferx-core emits begins with
  # "Covariance step". A bare `grepl("covariance step", ...)` also claimed
  # messages that merely mention the step in passing -- notably SIR's
  # "proposal was shrunk ... curvature in the covariance step" and postfit's
  # "SIR requested but the covariance step did not succeed" -- and answered
  # them with "Standard errors unavailable ... try covariance = FALSE", which
  # would remove the matrix SIR needs. Those belong to `sir`.
  #
  # This is not only a legacy-fit path: `ferx_load_fit()` does not restore
  # `warnings_structured`, so every fit read back from disk arrives here with
  # `general` on every row.
  category == "general" &&
    grepl("^(\\[[^]]*\\][[:space:]]*)?Covariance step", message)
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
# there is. Three kinds of `general` row are nonetheless recovered from their
# message text, because `ferx_load_fit()` does not restore the structured table
# and every row of a loaded fit arrives under this category: covariance-step
# messages (see `.ferx_is_covariance_warning()`), and the flat-theta and
# unused-declaration messages (see `.ferx_resolve_general()`).
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
    return(.ferx_covariance_guidance(message, category))
  }
  # ferx-core has no `unused_parameter` code: both of its unused-declaration
  # messages ("theta 'X' is declared in [parameters] but not referenced in any
  # model expression", and the [individual_parameters] "computed but never
  # used" one) fall through `classify_warning` to `general`, so they are routed
  # here by message, the same way the covariance family is.
  #
  # Gating on `general` alone is NOT enough, and the trap is worth naming: the
  # flat-theta message contains the literal phrase "computed but never used",
  # and `ferx_load_fit()` does not restore `warnings_structured`, so on a loaded
  # fit EVERY row arrives under `general` -- flat-theta included. There, the
  # category gate is not an exclusion, it is the admission criterion. The
  # resolver below therefore keys on the distinguishing text and tests
  # flat-theta first, because that message contains the unused phrase but not
  # the other way round.
  resolved <- .ferx_resolve_general(category, message)
  if (!is.null(resolved)) category <- resolved
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
    gradient_fallback  = "SAEM requested HMC sampling but this model is outside HMC's scope, so it fell back to Metropolis-Hastings. The message names which condition applied. The fit is valid; set saem_n_leapfrog = 0 to request MH deliberately and silence the notice.",
    mu_referencing     = if (grepl("not mu-referenced", message, fixed = TRUE)) {
      "The listed individual parameters are not mu-referenced, which ferx-core reports can strongly affect SAEM convergence. Rewrite them in a mu-referenced form (a THETA entering linearly on the transformed scale, e.g. CL = TVCL * exp(ETA_CL)) so the M-step can update them in closed form."
    } else {
      "Mu-referencing was auto-detected for the listed parameters (informational)."
    },
    optimizer_config   = if (grepl("global_search disabled", message, fixed = TRUE)) {
      "Global search was requested but could not be initialised, so the optimiser ran without it and started from your initial estimates alone. The message carries the underlying reason. Treat the result as a single-start fit: check it against a multi-start run before trusting it as the global optimum."
    } else {
      "Optimizer configuration note (informational)."
    },
    multi_start        = "Multi-start information (informational).",
    threads            = "Thread-pool sizing note. Consider matching threads to the subject count.",
    cancelled          = "The fit was cancelled before completion.",
    unused_parameter   = "A declared parameter is never referenced in any model expression, or an [individual_parameters] assignment is computed and then never used. It cannot affect predictions and will not be meaningfully estimated: remove it, or complete the expression that was meant to use it. A commented-out line or a misspelt name is the usual cause.",
    eta_shrinkage      = "High ETA shrinkage means the data do not inform that random effect per subject. Individual estimates and any EBE-based diagnostic (ETA plots, covariate screening) are unreliable for it; consider removing the ETA or collecting richer per-subject data.",
    eps_shrinkage      = "EPS shrinkage is notably negative, meaning mean(IWRES^2) > 1: the residual error model is not absorbing the residuals at the final EBEs, so sigma is under-fit rather than over-fit. Common causes are a SAEM run that stopped at a local optimum (polish with method = c(\"saem\", \"focei\"), or try different starts), model misspecification on a subset of subjects, and sigma sitting at a bound. Inspect the IWRES distribution in the sdtab.",
    boundary_estimate  = "A THETA came to rest on one of its bounds, so the estimate is where the bound put it rather than where the data did. The message names the parameter and the side: widen that bound if the value is plausible, or reconsider the parameterisation if the parameter is drifting because the data do not identify it.",
    high_correlation   = "The named parameter pair(s) are nearly collinear at convergence, so their individual estimates are poorly determined even when the fit as a whole looks good. Consider fixing one of each pair, removing it, or reparameterising.",
    inflated_rse       = "A relative standard error exceeded ferx-core's reporting threshold of 50%, so the named parameter is imprecisely estimated -- the data barely inform it. Check identifiability before quoting its precision, and consider fixing or removing it.",
    flat_parameter     = "A THETA has no effect on the objective at the initial estimate (unmapped, or dropped from the structural/scaling model). It was frozen at its initial value so the rest of the fit could proceed; map it into the model or remove it. Its reported estimate is just its initial value.",
    flip_flop          = "A transit / inverse-Gaussian absorption model entered the flip-flop regime (absorption no faster than elimination). If the model carries an ODE twin the affected subjects were rerouted and the fit is sound; if it does not, that subject's likelihood contribution is degenerate. Check the named subject's MTT/CL (transit) or MAT/CV2/CL (IG) estimates.",
    absorption_twin_declined = "The analytic absorption model kept no ODE fallback, so the fit is a plain closed-form fit but every feature that needs the fallback (time-varying covariates, a TIME-dependent parameter, IOV, steady-state or infusion doses, the flip-flop reroute) is now rejected instead. The message quotes the reason the fallback could not be built -- the set of reasons is open-ended, so read that text rather than guessing. Resolving it restores the fallback and the features above.",
    experimental       = "An experimental feature is in use. For neural-network components ferx-core states that standard errors for the network weights are not reliable; for SDE components it states that standard errors and convergence behaviour are not yet proven. Treat point estimates as provisional and do not quote the precision.",
    simulation         = "A simulated subject was handled specially: a degenerate (non-positive or non-finite) hazard draw, a recurrent-event stream whose expected count exceeded the cap so it was skipped and the subject censored at the horizon, or one truncated at the cap mid-stream. The estimated model is unaffected; inspect the named subject's hazard parameters and covariate values.",
    NULL
  )
}
