#' Structured warnings from a ferx fit
#'
#' Returns or prints the warnings produced by a fit, classified by severity
#' (\code{critical}, \code{warning}, \code{info}) and category. Severity and
#' category are assigned by the ferx-core engine (for engine warnings) or by
#' the R diagnostics layer (for condition number and ETA normality); the R
#' side never re-parses message text to guess severity. A fit read back with
#' \code{\link{ferx_load_fit}} carries only the flat message strings, so each
#' one is run through the engine's own classifier to recover the severity and
#' category the fresh fit had.
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
  df <- .ferx_fit_warnings(fit)
  if (isTRUE(as_df)) {
    return(df)
  }

  use_cli   <- .ferx_use_cli()
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
    guide <- .ferx_warning_guidance(row$category, message = row$message)
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

# Severity and category for flat warning strings, as ferx-core assigns them.
#
# Delegates to the engine's `classify_warning` (via ferx_rust_classify_warnings)
# so the R side keeps no second copy of its message patterns. The one R-side
# message the engine does not know - the dropped-`[output]` finding, which goes
# into the flat vector so it survives a `.fitrx` round trip - is recovered here.
# Returns the same four-column data frame as `fit$warnings_structured`.
.ferx_classify_flat_warnings <- function(msgs) {
  msgs <- as.character(msgs)
  msgs[is.na(msgs)] <- ""
  if (length(msgs) == 0L) {
    return(data.frame(severity = character(0), category = character(0),
                      message = character(0), source_method = character(0),
                      stringsAsFactors = FALSE))
  }
  cl <- ferx_rust_classify_warnings(msgs)
  df <- data.frame(
    severity      = as.character(cl$severity),
    category      = as.character(cl$category),
    message       = as.character(cl$message),
    source_method = as.character(cl$source_method),
    stringsAsFactors = FALSE
  )
  is_output <- df$category == "general" &
    grepl("named in [output]", df$message, fixed = TRUE)
  df$category[is_output] <- "output"
  df
}

# The fit's warnings as one structured table: `fit$warnings_structured` plus
# every flat message it does not already carry, classified by the engine.
#
# `ferx_load_fit()` does not restore `warnings_structured`, so a loaded fit has
# only the flat strings - and `ferx_sir()` / `ferx_covariance()` run on it then
# install a table holding just their own rows. Reading the table alone would
# drop every older warning on that fit; reading it only when absent (as this
# once did) dropped them the moment either post-hoc step ran (#308 review). On
# a fresh fit every flat message is already in the table, so it is returned
# unchanged. The fit itself is never modified, so check_strictness(), which
# reads `warnings_structured` directly, is unaffected.
.ferx_fit_warnings <- function(fit) {
  cols <- c("severity", "category", "message", "source_method")
  ws <- fit$warnings_structured
  if (is.data.frame(ws) && all(cols[1:3] %in% names(ws))) {
    if (!"source_method" %in% names(ws)) ws$source_method <- rep("", nrow(ws))
    ws <- ws[, cols, drop = FALSE]
  } else {
    ws <- NULL
  }
  msgs <- unique(as.character(unlist(fit$warnings %||% character(0),
                                     use.names = FALSE)))
  msgs <- msgs[!is.na(msgs) & nzchar(msgs)]
  flat <- .ferx_classify_flat_warnings(msgs)
  if (is.null(ws)) return(flat)
  missing <- !(flat$message %in% ws$message)
  if (!any(missing)) return(ws)
  out <- rbind(ws, flat[missing, , drop = FALSE])
  rownames(out) <- NULL
  out
}

# Guidance for the `mu_referencing` category. ferx-core's classifier gives this
# code two severities: Warning for "not mu-referenced" (a parameter the
# closed-form updates cannot use) and Info for every other mu-referencing note.
.ferx_mu_referencing_guidance <- function(message = "") {
  if (grepl("not mu-referenced", message, ignore.case = TRUE)) {
    return(paste0(
      "The listed individual parameters are not mu-referenced, so the ",
      "closed-form updates SAEM and IMPMAP use for mu-referenced parameters do ",
      "not apply to them and convergence can be slower or less stable. Where ",
      "possible, write them as typical value times exp(ETA), e.g. ",
      "CL = TVCL * exp(ETA_CL)."
    ))
  }
  "Mu-referencing note (informational): it reports how mu-referencing was detected or applied for the listed parameters. No action needed."
}

# Guidance for the `optimizer_config` category. Warning for "global_search
# disabled" - the CRS2-LM pre-search failed at runtime and the fit ran without
# it - and Info for the other global_search notes.
.ferx_optimizer_config_guidance <- function(message = "") {
  if (grepl("global_search disabled", message, ignore.case = TRUE)) {
    return(paste0(
      "The global pre-search (CRS2-LM) could not start, so the fit ran from the ",
      "declared initial values without it and may have settled in a local ",
      "optimum. Check the reason quoted in the message, or cover the parameter ",
      "space another way, e.g. settings = list(n_starts = 4L)."
    ))
  }
  "Optimizer configuration note (informational)."
}

# Guidance for the `flip_flop` category. Two emitters share it: the auto-reroute
# note (the ODE twin evaluates flip-flop parameters, so the profile is correct)
# and a twin-less model whose subjects silently degenerate at their EBEs.
.ferx_flip_flop_guidance <- function(message = "") {
  if (grepl("automatically evaluates", message, fixed = TRUE)) {
    return(paste0(
      "Informational: the flip-flop parameters are evaluated with the equivalent ",
      "ODE model, so predictions are correct but slower than the closed form. ",
      "Check the absorption starting estimates if flip-flop kinetics are not ",
      "expected."
    ))
  }
  paste0(
    "The listed subjects fall where the analytic absorption closed form returns ",
    "an all-zero concentration profile, and this model has no ODE equivalent to ",
    "fall back on, so their likelihood contributions are silently degenerate. ",
    "Rewrite the absorption as an explicit [odes] model, or revisit the ",
    "absorption and clearance starting estimates."
  )
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

# Message phrases that say which clauses an `ode_solver` statistics warning
# carries. ferx-core builds that warning as one message with a clause per
# non-zero solver counter (ode_solver_diagnostics_warning() in
# src/api/postfit.rs), and the guidance below depends on which clauses are
# there.
#
# The warning's `details` payload names every counter, but it never reaches R:
# fit_result_to_list() in src/rust/src/lib.rs passes each warning's severity,
# category, message and source_method, and drops `details`. So the clauses have
# to be told apart by their text. Each phrase below is the one that says what
# its counter measures, and every one except the rejected-escalation phrase is
# also asserted on by ferx-core's own tests
# (src/api/tests/ode_solver_diagnostics_tests.rs), so rewording it breaks a
# core test too. test-ferx_get_warnings.R checks every phrase against the
# engine source at the pinned revision.
.ferx_ode_solver_anchors <- function() {
  list(
    # Counters no solver setting can fix.
    #   abandoned_non_finite_timeline (ferx-core #1234)
    abandoned_timeline = "could not be ordered",
    #   auto_stiff_rejected_jets (ferx-core #1204)
    sensitivity_overflow = "analytic-sensitivity solve",
    # Counters for an integration that ran but not cleanly, where a solver
    # setting is the remedy: kept_clamped_steps, auto_stiff_rejected,
    # auto_fallback_failed, kept_unfinished_segments, stiff_aborted_segments.
    unclean_integration = c(
      "clamped at the minimum step size",
      "stiff escalation(s) chosen by",
      "had both attempts fail",
      "returned segment(s) stopped",
      "were abandoned early"
    )
  )
}

# Guidance for the `ode_solver` category.
#
# Two emitters share this code and disagree about severity: the `auto` guard's
# escalation note is informational - the stiff method coped - while the
# post-fit statistics pass is not. ferx-core prefixes its own token to each
# message and its classifier matches on that, so the token picks between them.
#
# The token cannot tell the statistics warning's clauses apart, and two of them
# are not solver-setting problems. A walk abandoned because the subject's
# timeline could not be ordered never integrated at all (ferx-core #1234), and a
# segment whose analytic sensitivities overflowed was integrated by the stiff
# method without trouble (ferx-core #1204). The engine's clause says so for
# each - no solver setting changes an abandoned walk, and a different
# ode_method will not help an overflow - so advice to change them, printed
# directly beneath, contradicted it. Those two now get their own advice. The
# solver-setting advice is kept for the clauses it does apply to. A message
# carrying neither of the two - including one whose wording is not recognised -
# gets the solver-setting advice as before.
.ferx_ode_solver_guidance <- function(message = "") {
  if (grepl("W_ODE_SOLVER_ESCALATION_NOTE", message, fixed = TRUE)) {
    return(paste0(
      "ode_method = \"auto\" escalated to a stiff solver and the stiff ",
      "method coped (informational). Set ode_method explicitly to skip the ",
      "non-stiff attempt."
    ))
  }
  anchors <- .ferx_ode_solver_anchors()
  carries <- function(phrases) {
    any(vapply(phrases, grepl, logical(1), x = message, fixed = TRUE))
  }
  abandoned <- carries(anchors$abandoned_timeline)
  overflow  <- carries(anchors$sensitivity_overflow)
  if (!abandoned && !overflow) {
    return(paste0(
      "Integration under the final estimates was not clean: steps clamped at ",
      "the minimum step size, a stiff escalation the auto guard discarded, or ",
      "a segment cut short by ode_stiff_abort_after. The optimizer may still ",
      "have converged. Set ode_method explicitly, or adjust ode_abstol / ",
      "ode_reltol / ode_max_steps."
    ))
  }
  parts <- character(0)
  if (abandoned) {
    parts <- c(parts, paste0(
      "Nothing was integrated for some subjects: their timeline could not be ",
      "ordered (a NaN or infinite dose time, lagtime, route lag or infusion ",
      "duration), so their predictions are NaN. No ode_method or tolerance ",
      "setting changes that. Check those dose records, and any covariate model ",
      "on ALAG / F / D / R that can overflow."
    ))
  }
  if (overflow) {
    parts <- c(parts, paste0(
      "The analytic sensitivities (the gradients FOCE/FOCEI use) overflowed on ",
      "some segments although the predictions stayed finite. A different ",
      "ode_method will not help. Check the model's units and scaling (a state ",
      "in ng rather than mg, an unbounded growth term) before trusting the ",
      "estimates."
    ))
  }
  if (carries(anchors$unclean_integration)) {
    parts <- c(parts, paste0(
      "For the clamped, discarded, unfinished or aborted segments the message ",
      "also reports, set ode_method explicitly, or adjust ode_abstol / ",
      "ode_reltol / ode_max_steps."
    ))
  }
  paste(parts, collapse = " ")
}

# Does this warning belong to the covariance family?
#
# Four categories carry covariance-step messages, not one: ferx-core's
# `covariance_failed` / `covariance_regularized` / `covariance_step`, plus
# `"covariance"`, which `ferx_covariance()` assigns to the engine's flat
# covariance warnings (R/ferx_covariance.R). Gating on `covariance_step` alone
# made every targeted branch below unreachable, and omitting `"covariance"` left
# the whole post-hoc covariance surface without guidance.
#
# `general` is admitted on the message text alone, anchored at the start of the
# message (optionally after a `[METHOD]` chain prefix) because every covariance
# warning ferx-core emits begins with "Covariance step". That matters beyond
# legacy fits: `ferx_load_fit()` does not restore `warnings_structured`, so
# every row of a fit read back from disk arrives under `general`. Anchoring
# keeps SIR's diagnostics -- which mention the covariance step in passing -- out
# of this family, where they would be told to re-run with covariance = FALSE,
# the one setting that removes the matrix SIR needs.
.ferx_is_covariance_warning <- function(category, message = "") {
  if (category %in% c("covariance_failed", "covariance_regularized",
                      "covariance_step", "covariance")) {
    return(TRUE)
  }
  category == "general" &&
    grepl("^(\\[[^]]*\\][[:space:]]*)?Covariance step", message)
}

# Remediation guidance keyed by the fixed category vocabulary that ferx-core
# (and the R-side additions) emit. Extend this table when core grows a new
# category. Returns NULL for unknown categories so callers can skip printing
# rather than showing a generic placeholder.
#
# message picks the dw_autocorrelation arm (positive vs negative DW) and, for
# covariance_step, selects targeted advice from the specific failure mode
# embedded in the message text.
#
# No [diffusion] arm here, on purpose. The positive-autocorrelation guidance
# used to append "For ODE models, also consider SDE process noise ([diffusion]
# block)." on every ODE fit. ferx-core dropped the equivalent sentence from its
# own Durbin-Watson warning in ferx-core #1426 (refs ferx-core #1285): ferx
# implements the covariance half of the EKF only, so the state mean stays the
# deterministic ODE solution and is never corrected by the observed values. A
# [diffusion] term re-weights the fit rather than following a subject's drift,
# and is not a remedy for IWRES autocorrelation.
.ferx_warning_guidance <- function(category, message = "") {
  # A `general` row is recovered to the category the engine gives its message
  # first. `ferx_load_fit()` does not restore `warnings_structured`, so a loaded
  # fit's rows would otherwise all arrive here as `general`, and every arm
  # below keyed on a category would be dead for them (#308). A message the
  # engine itself files under `general` is unchanged by this.
  if (identical(category, "general") && is.character(message) &&
      length(message) == 1L && !is.na(message) && nzchar(message)) {
    category <- .ferx_classify_flat_warnings(message)$category
  }
  if (category == "dw_autocorrelation") {
    if (grepl("egative", message, ignore.case = TRUE)) {
      return("Negative IWRES autocorrelation suggests over-parameterisation or a misspecified error model. Consider removing a parameter or simplifying the residual model.")
    }
    return(paste0(
      "Positive IWRES autocorrelation suggests missing structural dynamics.",
      " Consider transit absorption, an extra compartment, or IOV on ka/F."
    ))
  }
  if (category == "ode_solver") {
    return(.ferx_ode_solver_guidance(message))
  }
  if (.ferx_is_covariance_warning(category, message)) {
    # Omega non-PD -- checked before general non-PD because omega messages also
    # contain "not positive definite" and "eigenvalue".
    #
    # Matched on "omega matrix is" rather than on a descriptor: ferx-core
    # interpolates one of two into this sentence depending on the sign of the
    # smallest eigenvalue ("not positive definite" / "near-singular",
    # covariance.rs), and pinning the grep to the first sent the second to the
    # generic "standard errors unavailable" fallback instead of here. The text
    # below is unchanged and still speaks of a near-singular omega; wording that
    # distinguishes the two cases is a separate change.
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
    # Off-diagonal FD stencil non-finite. ferx-core emits this on its *success*
    # path: the covariance matrix was produced and the standard errors exist,
    # they are just missing their cross-partial terms. It matches none of the
    # branches above, so without this arm it inherits the "standard errors
    # unavailable" fallback and contradicts the message it prints under.
    if (grepl("off-diagonal FD stencil", message, ignore.case = TRUE)) {
      return(paste0(
        "Standard errors were produced. The cross-partial terms for the named ",
        "parameter(s) could not be evaluated and were set to zero, so their SEs ",
        "may be over-optimistic. The message suggests tuning fd_hessian_step; ",
        "cross-check the affected parameters with ferx_sir()."
      ))
    }
    # The S matrix is the sum of per-subject score outer products, so one
    # subject whose score cannot be evaluated takes the whole matrix with it -
    # and the message names that subject, which is the useful half of it.
    if (grepl("quadrature score for subject", message, ignore.case = TRUE)) {
      return(paste0(
        "The score for the named subject could not be evaluated, so the S ",
        "matrix (and any covariance built on it) is unavailable. Inspect that ",
        "subject: very few observations, all-BLOQ records or an extreme ",
        "covariate are the usual causes. The R-matrix covariance does not use ",
        "per-subject scores - ferx_covariance(fit, covariance_method = \"r\") ",
        "- and ferx_sir() estimates uncertainty without a covariance step at ",
        "all."
      ))
    }
    # The sum itself came out non-finite: no single subject is named, so the
    # advice is about the matrix rather than about a record.
    if (grepl("non-finite score cross-product", message, ignore.case = TRUE)) {
      return(paste0(
        "The score cross-product overflowed or contained a non-finite entry, ",
        "so the S matrix is unavailable and no standard errors were produced ",
        "from it. This usually follows extreme parameter values at ",
        "convergence. Check the estimates against their bounds, then try the ",
        "R-matrix covariance (ferx_covariance(fit, covariance_method = \"r\")) ",
        "or ferx_sir()."
      ))
    }
    # An invalid step never reaches the Hessian: nothing about the model was
    # diagnosed, so the identifiability advice in the fallback does not apply.
    if (grepl("fd_hessian_step must be positive and finite", message, fixed = TRUE)) {
      return(paste0(
        "The covariance step did not run because fd_hessian_step is not a ",
        "positive finite number. Set a positive value (the default is 0.01) in ",
        "ferx_fit(..., fd_hessian_step = ) or in [fit_options], and re-run."
      ))
    }
    # covariance_method = "s" alone: the S matrix itself cannot be inverted.
    # The R-matrix estimators do not need it inverted, so they are the remedy.
    if (grepl("score cross-product matrix S is singular", message, fixed = TRUE)) {
      return(paste0(
        "The score cross-product S is singular, which usually means fewer ",
        "subjects than free parameters or collinear per-subject scores, so the ",
        "S-only covariance cannot be formed. Use an estimator that does not ",
        "invert S: ferx_covariance(fit, covariance_method = \"r\") or ",
        "covariance_method = \"rsr\"."
      ))
    }
    # Cancelled part-way (`COV_CANCELLED_MSG`). No standard errors, but nothing
    # was diagnosed about the model either, so the identifiability advice in the
    # fallback would report a finding the engine never made.
    if (grepl("cancelled before completion", message, ignore.case = TRUE)) {
      return(paste0(
        "The covariance step was cancelled before it finished, so no standard ",
        "errors were produced and nothing was diagnosed about the model. Re-run ",
        "and let the step complete if you need them."
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
    # ferx-core's Info-level note about the cost of the step, emitted BEFORE it
    # runs. Matched on the message, never on the `covariance_step` code alone:
    # that code's classification ends in a catch-all, so a future message
    # carrying it need not be the cost note. Without this arm the note inherits
    # the failure fallback and reports a failure that has not happened.
    if (grepl("OFV evaluations", message, fixed = TRUE)) {
      return(paste0(
        "Informational: the covariance step cost scales with the square of the ",
        "parameter count. No action needed; pass covariance = FALSE to skip it ",
        "during development."
      ))
    }
    # Generic fallback for older or unrecognised covariance messages.
    return(paste0(
      "Standard errors unavailable. Check identifiability; try a simpler ",
      "omega/sigma structure or covariance = FALSE for development."
    ))
  }
  # ferx-core has no `unused_parameter` code -- its unused-declaration messages
  # fall through `classify_warning` to `general`, so the arm below was
  # unreachable in production. Recover it from the message.
  #
  # The flat-theta message contains the literal phrase "computed but never
  # used" too, so it must be excluded explicitly: gating on `general` is NOT
  # enough, because `ferx_load_fit()` does not restore `warnings_structured`
  # and every row of a loaded fit arrives under that category.
  if (identical(category, "general") &&
      !grepl("has no effect on the objective", message, fixed = TRUE) &&
      (grepl("declared in [parameters] but not referenced", message, fixed = TRUE) ||
       grepl("computed but never used", message, fixed = TRUE))) {
    category <- "unused_parameter"
  }
  switch(category,
    convergence        = "Optimizer did not reach convergence. Try different initial values, method = c(\"saem\", \"focei\"), or settings = list(n_starts = 4L).",
    # ferx-core #1380. Emitted even when `converged` is TRUE: a fit that never
    # moved has a flat objective trace to plateau on, so the flag cannot be
    # trusted by itself. The engine suppresses it for an evaluation-only run.
    stalled_at_init    = "The fit never left its initial estimates, so the reported OFV is the objective of the starting values and says nothing about the model, even if converged is TRUE. Check fit$final_gradient, try different initial estimates, or a gradient-based optimizer if this fit used a derivative-free one.",
    # ferx-core #1386. The reported fit already uses the better of the two EBE
    # sets; the warning is about how far the EBE-derived output can be trusted.
    ebe_start_dependent = "The empirical Bayes estimates depend on where the inner loop starts, so IPRED, IWRES, CWRES, shrinkage and the covariance step do too. Raise settings = list(inner_maxiter = 500L) to rule out an inner budget a cold start cannot finish within, or inner_restarts for a suspected second mode.",
    condition_number   = "Parameters are correlated/ill-scaled. Consider fixing or removing a parameter, or reparameterising.",
    optimizer_health   = "Optimizer struggled (trust region / Hessian). Inspect the trace and consider better starting values.",
    vi_bad_basin       = "VI's final ELBO check found that the flat objective is a bad basin, not a usable variational approximation. Refit from different initial values, raise settings = list(n_starts = 4L), or use method = \"focei\".",
    boundary_estimate  = "A THETA estimate finished on one of its declared bounds, so it is not an interior optimum and its standard error is not meaningful. Widen the bound if it is too tight, FIX the parameter if the boundary value is intended, or simplify the model if the data do not inform it.",
    parameter_at_runaway_guard = "A coordinate is pinned to an internal safety limit (an implicit theta cap, or an omega/sigma guard), so this is not an interior optimum. Give the parameter explicit bounds, fix it, or remove the term it belongs to.",
    # The start-side twin of parameter_at_runaway_guard above: same internal
    # rails, but this one fires before the first objective evaluation rather
    # than at convergence. It is deliberately NOT boundary_estimate, which is
    # about where a fit ENDED and which drives three default-on rejection
    # filters (bootstrap's skip_estimate_near_boundary, reject_on_boundary,
    # and .ferx_boundary_detail() in check_strictness.R) - a clamped start
    # wearing that category would silently drop bootstrap replicates.
    #
    # The advice names omega and sigma specifically because that is the whole
    # population of this arm, measured at pin 944cbf1e: a theta start past the
    # hidden 1e9 cap does NOT arrive here, it is the error
    # E_THETA_INIT_OUTSIDE_BOUNDS ("above its own declared upper bound of 1e9")
    # and stops the fit. Only the omega variance rail and the sigma SD rail
    # produce a W_INIT_OUTSIDE_BOUNDS warning row. Neither declaration form
    # (`omega X ~ v`, `sigma X ~ v`) accepts bounds, so advising the reader to
    # declare some -- as this arm first did -- is an instruction they cannot
    # follow, printed directly beneath an engine message that already gives the
    # right one.
    init_outside_bounds = "A start value was clamped onto one of the optimizer's internal rails before the first objective evaluation, so the fit did not begin from what the model file declares. This is an omega variance or a sigma SD, neither of which takes explicit bounds: move the start inside the rail quoted above, or FIX the parameter to hold the declared value.",
    inflated_rse       = "The listed THETAs are imprecisely estimated: the data barely inform them. Check that the design covers them (e.g. absorption-phase samples for KA), consider fixing or removing them, and cross-check their intervals with ferx_sir().",
    high_correlation   = "Highly correlated estimates are not separated by the data. Inspect fit$cor_matrix, then fix or remove one parameter of each pair, or reparameterise.",
    eta_shrinkage      = "EBE-based diagnostics (ETA-versus-covariate plots, IPRED, IWRES) are unreliable for the listed ETAs - do not screen covariates on them. Consider removing that IIV term, or a design that informs it.",
    eps_shrinkage      = "Negative EPS shrinkage means the residuals at the EBEs are larger than the residual error model allows. Inspect IWRES in fit$sdtab; if a SAEM fit, polish with method = c(\"saem\", \"focei\"), otherwise revisit the error model or the subjects it fits poorly.",
    flat_parameter     = "The THETA had no effect on the objective at its initial value, so it was FIXed there and not estimated. Map it into the model (e.g. in [structural_model]) or remove it from [parameters].",
    experimental       = "This model uses a feature marked experimental (SDE [diffusion] or neural-network components), validated on few datasets. Cross-check the result against a model without the feature, and treat its standard errors with caution.",
    absorption_twin_declined = "The absorption model stays closed-form with no ODE fallback, so subjects that need one (time-varying covariates, IOV, steady-state or infusion doses, flip-flop kinetics) will be rejected with an error. Fix the reason quoted after 'Reason:', or write the absorption as an explicit [odes] model.",
    flip_flop          = .ferx_flip_flop_guidance(message),
    simulation         = "A simulated subject drew a degenerate hazard; it was censored at the end of its observation window instead of producing events. Check the hazard parameters and that subject's covariate values.",
    eta_normality      = "ETA distribution may be non-normal. High shrinkage or sparse data can cause this; prefer QQ-plots for diagnosis.",
    bloq_method        = "LOQ censoring note. Set method = \"focei\" explicitly to silence, or review the M3 setup.",
    sir                = .ferx_sir_guidance(message),
    importance_sampling = "Importance-sampling ESS collapsed for some subjects. Raise imp_samples / imp_proposal_df or check EBE quality.",
    data_quality       = "Data issue detected. Review the flagged observations in the dataset.",
    omega_structure    = "Mixed parameterisation in a block omega. Check the [individual_parameters] forms for the correlated etas.",
    ebe_convergence    = "Some subjects' inner EBE search did not converge. Inspect those subjects or relax inner_tol / max_unconverged_frac.",
    gradient_fallback  = "Gradient method fell back (e.g. AD -> FD or HMC -> MH). The fit is valid; expect a longer runtime.",
    mu_referencing     = .ferx_mu_referencing_guidance(message),
    optimizer_config   = .ferx_optimizer_config_guidance(message),
    multi_start        = "Multi-start information (informational).",
    threads            = "Thread-pool sizing note. Consider matching threads to the subject count.",
    cancelled          = "The fit was cancelled before completion.",
    unused_parameter   = "A declared parameter is never referenced in [individual_parameters] or [error_model]. Remove it from [parameters] or complete the expression that uses it.",
    output             = "A name in [output] produced no sdtab column. Drop it from the block, or read the quantity where it does live (etas: fit$ebe_etas).",
    NULL
  )
}
