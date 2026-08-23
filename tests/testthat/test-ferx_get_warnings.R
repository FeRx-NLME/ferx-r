test_that("ferx_get_warnings(as_df = TRUE) returns the underlying data frame", {
  fit <- warfarin_fit_cov()
  df <- ferx_get_warnings(fit, as_df = TRUE)
  expect_identical(df, fit$warnings_structured)
})
test_that("ferx_get_warnings() prints a grouped summary", {
  fit <- warfarin_fit_cov()
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("ferx fit warnings", out)))
  # Footer tallies are always present
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})
test_that("ferx_get_warnings() rejects non-fit input", {
  expect_error(ferx_get_warnings(list()), "ferx_fit")
})
test_that("ferx_get_warnings() falls back to flat warnings when structured is absent", {
  fake <- structure(
    list(
      model_name = "legacy",
      warnings = c("Outer optimization did not converge", "something else"),
      warnings_structured = NULL
    ),
    class = "ferx_fit"
  )
  df <- ferx_get_warnings(fake, as_df = TRUE)
  expect_equal(nrow(df), 2L)
  expect_true(all(df$severity == "warning"))
  expect_true(all(df$category == "general"))
})
test_that(".ferx_warning_guidance dispatches the covariance family by message content", {
  # Each message is paired with the category ferx-core's `classify_warning`
  # actually assigns it, not with a hand-picked one. Every message below carries
  # "covariance step failed" or "covariance step" + "not positive definite", so
  # core codes it `covariance_failed`; the regularisation messages carry
  # "covariance step regularized", so core codes those `covariance_regularized`.
  # Neither reaches the code `covariance_step`, which core reserves for its
  # Info-level cost note -- so pinning these to "covariance_step" (as this test
  # did until the routing fix) asserted a pairing production never emits, and
  # every branch under test was dead in the real call path.
  g <- function(msg, category = "covariance_failed") {
    ferx:::.ferx_warning_guidance(category, message = msg)
  }

  # NonPdHessian path: eigenvalue list in message.
  msg_npd <- paste0(
    "Covariance step: Hessian is not positive definite. ",
    "Eigenvalues: [8.4000, 2.1000, -0.0100]. SE estimates not available."
  )
  expect_match(g(msg_npd), "eigenvalue", ignore.case = TRUE)
  expect_match(g(msg_npd), "near-zero|negative", ignore.case = TRUE, perl = TRUE)

  # Ill-conditioned Hessian entries: names a parameter.
  msg_ic <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) -- theta[CL] (non-finite diagonal). ",
    "SE estimates not available."
  )
  expect_match(g(msg_ic), "fd_hessian_step", ignore.case = TRUE)

  # Omega non-PD.
  msg_omega <- paste0(
    "Covariance step failed: Omega matrix is not positive definite at ",
    "convergence (min eigenvalue = 1.2e-10; eigenvalues: [0.5000, 1.2e-10]). ",
    "SE estimates not available."
  )
  expect_match(g(msg_omega), "near-singular", ignore.case = TRUE)
  expect_false(grepl("eigenvalue list", g(msg_omega), ignore.case = TRUE))

  # Non-finite OFV.
  msg_ofv <- paste0(
    "Covariance step failed: base OFV is non-finite at convergence ",
    "(likely numerical overflow or underflow in model evaluation). ",
    "SE estimates not available."
  )
  expect_match(g(msg_ofv), "overflow|underflow", ignore.case = TRUE, perl = TRUE)

  # Regularisation: minor, moderate, severe.
  base_reg <- function(sev) paste0(
    "Covariance step regularized: eigenvalue floor applied to FD Hessian ",
    "(1 of 3 free-block eigenvalues clipped; min eig = 1.2e-6, floor = 8.4e-14; ",
    "severity: ", sev, "). Standard errors are likely reliable."
  )
  reg <- function(sev) g(base_reg(sev), category = "covariance_regularized")
  expect_match(reg("minor"),    "benign|reliable",     ignore.case = TRUE, perl = TRUE)
  expect_match(reg("moderate"), "ferx_sir|caution",    ignore.case = TRUE, perl = TRUE)
  expect_match(reg("severe"),   "ferx_sir|unreliable", ignore.case = TRUE, perl = TRUE)
  # Minor guidance should not suggest SIR.
  expect_false(grepl("ferx_sir", reg("minor"), ignore.case = TRUE))

  # Generic fallback for unrecognised message.
  expect_match(g("Covariance step failed"), "identifiability", ignore.case = TRUE)
})

test_that("the informational covariance_step note does not read as a failure", {
  # Core codes only its Info-level cost note `covariance_step`. It matches none
  # of the targeted branches, so before this arm existed it fell through to the
  # "standard errors unavailable" fallback and told the user the step had failed
  # when it had not yet run.
  msg <- "Covariance step: 132 N^2 OFV evaluations for 12 parameters."
  g <- ferx:::.ferx_warning_guidance("covariance_step", message = msg)
  expect_match(g, "Informational", ignore.case = TRUE)
  expect_false(grepl("unavailable|identifiability", g, ignore.case = TRUE))
})

test_that("an unrecognised covariance_step message is not called benign", {
  # The informational arm is matched on the message, not on the code alone.
  # Core's `CovarianceStep` classification ends in a catch-all -- any otherwise
  # unmatched message carrying "covariance step:" and "parameters" lands there
  # -- so the code is not proof the message is the cost note. Answering an
  # unknown one "no action needed" fails in the dangerous direction, so it must
  # fall through to the failure fallback instead.
  g <- ferx:::.ferx_warning_guidance(
    "covariance_step",
    message = "Covariance step: something new about 12 parameters."
  )
  expect_match(g, "identifiability", ignore.case = TRUE)
  expect_false(grepl("Informational|No action needed", g, ignore.case = TRUE))
})

test_that("a covariance step that succeeded is not reported as unavailable", {
  # ferx-core emits this on the CovarianceStepResult::Success path: the matrix
  # was produced and the SEs exist, they are merely missing their cross-partial
  # terms. `classify_warning` codes it `covariance_regularized` (the
  # "off-diagonal fd stencil" arm), so the routing fix admits it to the
  # covariance block for the first time -- where, before this arm, it matched
  # nothing and inherited "Standard errors unavailable", contradicting the
  # message it prints under.
  msg <- paste0(
    "Covariance step: off-diagonal FD stencil(s) non-finite for theta[CL], ",
    "sigma[1]. Cross-partial correlation set to 0; SE for these parameter(s) ",
    "may be over-optimistic. Try tuning fd_hessian_step."
  )
  g <- ferx:::.ferx_warning_guidance("covariance_regularized", message = msg)
  expect_match(g, "over-optimistic", ignore.case = TRUE)
  expect_match(g, "fd_hessian_step", fixed = TRUE)
  expect_false(grepl("unavailable", g, ignore.case = TRUE))
})

test_that("an invalid fd_hessian_step is answered as a bad argument", {
  # Core reports it as a covariance-step failure, but the fix is the argument
  # the user passed, not the model -- the identifiability fallback would send
  # them to restructure omega over a typo.
  msg <- paste0(
    "Covariance step failed: fd_hessian_step must be positive and finite, ",
    "got 0. SE estimates not available."
  )
  g <- ferx:::.ferx_warning_guidance("covariance_failed", message = msg)
  expect_match(g, "fd_hessian_step", fixed = TRUE)
  expect_false(grepl("identifiability|omega", g, ignore.case = TRUE))
})

test_that("ferx_covariance()'s `covariance` category reaches the guidance", {
  # `ferx_covariance()` folds the engine's flat covariance warnings into the
  # structured table under `category = "covariance"` (R/ferx_covariance.R) --
  # a fourth token, not one of core's three. Omitting it left the entire
  # post-hoc covariance surface without guidance even after the routing fix.
  msg <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) -- theta[CL] (non-finite diagonal)."
  )
  g <- ferx:::.ferx_warning_guidance("covariance", message = msg)
  expect_match(g, "Hessian diagonal", ignore.case = TRUE)
  expect_match(g, "fd_hessian_step", fixed = TRUE)
})

test_that("a legacy fit's `general` covariance message still gets guidance", {
  # `ferx_get_warnings()` stamps `general` on every message of a fit saved
  # before structured warnings existed. Every branch is keyed on the message,
  # so such a message can be answered exactly as well as one arriving under its
  # own code -- routing it around the block re-created, for those fits, the
  # silence this whole change removes.
  msg <- paste0(
    "Covariance step failed: Omega matrix is not positive definite at ",
    "convergence."
  )
  g <- ferx:::.ferx_warning_guidance("general", message = msg)
  expect_match(g, "Omega is near-singular", ignore.case = TRUE)
  # A `general` message that is not about the covariance step stays unhandled.
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Something else"))
})
test_that("ferx_get_warnings() shows guidance for an unused declaration", {
  # End to end through the printer, with the category ferx-core actually
  # assigns: its unused-declaration messages classify to `general`, not to
  # `unused_parameter`. Pinning this to the hand-made pairing (as it did until
  # the routing fix) asserted a row production never produces.
  msg <- paste0(
    "theta 'TVCL' is declared in [parameters] but not referenced in any model ",
    "expression -- it will not affect predictions or be meaningfully estimated"
  )
  fake <- structure(
    list(
      model_name = "m",
      warnings_structured = data.frame(
        severity      = "warning",
        category      = "general",
        message       = msg,
        source_method = "",
        stringsAsFactors = FALSE
      ),
      warnings = msg,
      condition_number = NULL,
      eta_normality = NULL,
      uses_sde = FALSE
    ),
    class = "ferx_fit"
  )
  out <- capture.output(ferx_get_warnings(fake))
  expect_true(
    any(grepl("cannot affect predictions", out, fixed = TRUE)),
    info = paste("Expected guidance text not found in output:\n",
                 paste(out, collapse = "\n"))
  )
})

# ---- header from test-diagnostics-helpers.R ----
# Unit tests for the pure-R diagnostic/formatting helpers in diagnostics.R:
# the warning-guidance lookup and the cli styling shims. No model fit required.

# ---------------------------------------------------------------------------
# .ferx_warning_guidance
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# .ferx_use_cli / .ferx_style
# ---------------------------------------------------------------------------




test_that(".ferx_warning_guidance toggles the SDE hint for positive autocorrelation", {
  pos_ode <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = FALSE)
  expect_match(pos_ode, "Positive IWRES autocorrelation")
  expect_match(pos_ode, "SDE process noise")            # hint shown when not already using SDE
  pos_sde <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive", uses_sde = TRUE)
  expect_match(pos_sde, "Positive IWRES autocorrelation")
  expect_false(grepl("SDE process noise", pos_sde))     # hint suppressed
})
test_that(".ferx_warning_guidance returns negative-autocorrelation guidance", {
  neg <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "Negative autocorrelation")
  expect_match(neg, "Negative IWRES autocorrelation")
})
# Every token ferx-core's `WarningCode::as_str()` can emit, minus `general`.
# Hand-transcribed, and deliberately so: `.ferx_core_warning_codes()` below
# checks it against the real thing whenever a sibling ferx-core checkout is
# present, which is the only mechanism here that can notice core growing a code.
.core_warning_cats <- function() {
  c(
    "absorption_twin_declined", "bloq_method", "boundary_estimate", "cancelled",
    "condition_number", "convergence", "covariance_failed",
    "covariance_regularized", "covariance_step", "data_quality",
    "dw_autocorrelation", "eps_shrinkage", "eta_normality", "eta_shrinkage",
    "experimental", "flat_parameter", "flip_flop", "gradient_fallback",
    "high_correlation", "importance_sampling", "inflated_rse", "mu_referencing",
    "multi_start", "omega_structure", "optimizer_config", "optimizer_health",
    "simulation", "sir", "threads"
  )
}

# Categories that reach `.ferx_warning_guidance()` without being ferx-core
# `WarningCode` tokens.
#
#   covariance        - assigned by ferx_covariance() (R/ferx_covariance.R) to
#                       the engine's flat covariance warnings.
#   ebe_convergence   - no current emitter. ferx-r surfaces EBE convergence as
#   unused_parameter    the count field `ebe_convergence_warnings`, not as a
#                       warning row; core has no `WarningCode` for either, and
#                       its "computed but never used" parser message classifies
#                       to `general`. The arms are kept because they cost
#                       nothing and read as the intended taxonomy, but they are
#                       unreachable today -- see the note on the roxygen for
#                       ferx_fit(), which still describes `unused_parameter` as
#                       parser-emitted.
.extra_warning_cats <- function() {
  c("covariance", "ebe_convergence", "unused_parameter")
}

test_that("every category in the guidance table returns guidance", {
  # A completeness walk over the table, NOT a drift guard: the loop only visits
  # tokens listed above, so a code core added and nobody transcribed is invisible
  # to it. `.ferx_warning_guidance() matches ferx-core's WarningCode vocabulary`
  # is the test that catches that.
  for (cat in c(.core_warning_cats(), .extra_warning_cats())) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
})

test_that(".ferx_warning_guidance matches ferx-core's WarningCode vocabulary", {
  # The real drift guard. Reads `WarningCode::as_str()` out of a sibling
  # ferx-core checkout instead of trusting a second hand-maintained copy, so a
  # code added in core fails here even though nobody edited this file.
  #
  # Skipped when the sibling is absent: ferx-r CI builds the pinned crate, not a
  # checkout. That is the right trade -- drift is introduced on a developer
  # machine, which per CLAUDE.md always has ../ferx-core, and this fails there.
  candidates <- c(
    file.path("..", "..", "..", "ferx-core", "src", "types.rs"),
    file.path("~", "ferx-core", "src", "types.rs")
  )
  types_rs <- Filter(function(f) file.exists(path.expand(f)), candidates)
  skip_if(length(types_rs) == 0L, "no sibling ferx-core checkout to read")

  src <- readLines(path.expand(types_rs[[1]]), warn = FALSE)
  hits <- regmatches(src, regexpr('WarningCode::[A-Za-z]+ => "[a-z_]+"', src))
  tokens <- gsub('.*=> "|"$', "", hits)
  expect_gt(length(tokens), 20L)          # the file was found and understood

  # `general` is core's unrecognised-message bucket and is deliberately absent
  # from the table -- the message text is the only guidance there is.
  expect_setequal(setdiff(tokens, "general"), .core_warning_cats())
  for (cat in setdiff(tokens, "general")) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
})

test_that(".ferx_warning_guidance gives `general` no category-level guidance", {
  # `general` is core's bucket for a message its classifier did not recognise,
  # so there is no category-level remediation to give -- the message text is the
  # guidance. Asserted so it reads as a decision rather than an omission. The
  # one exception is a covariance-step message arriving under `general`, which
  # the message-keyed covariance family still answers (tested above).
  expect_null(ferx:::.ferx_warning_guidance("general"))
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Anything at all."))
})
test_that("core's unused-declaration messages reach the unused_parameter guidance", {
  # ferx-core has no `unused_parameter` WarningCode -- both messages fall
  # through `classify_warning` to `general`, so the arm was unreachable in
  # production for exactly the reason the covariance block was. The pre-existing
  # test passed only because it fed the category in by hand.
  theta_msg <- paste0(
    "theta 'TVCL' is declared in [parameters] but not referenced in any model ",
    "expression -- it will not affect predictions or be meaningfully estimated"
  )
  ip_msg <- paste0(
    "[individual_parameters] KA is computed but never used -- not mapped into ",
    "the `pk(...)` model and not referenced in any other block, so it has no ",
    "effect. Map KA in [structural_model] (e.g. `f=F`) or remove KA."
  )
  for (msg in c(theta_msg, ip_msg)) {
    g <- ferx:::.ferx_warning_guidance("general", message = msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = msg)
    expect_match(g, "never referenced in any model expression|never used",
                 ignore.case = TRUE, perl = TRUE, info = msg)
  }
  # An unrelated `general` message is still left to its own text.
  expect_null(ferx:::.ferx_warning_guidance("general", message = "Something else."))
})

test_that("every covariance message ferx-core emits gets non-contradictory guidance", {
  # The inventory test. Every message below is copied verbatim from a ferx-core
  # emit site and paired with the category `classify_warning` assigns it -- not
  # with a category chosen by hand. The previous tests enumerated only messages
  # the block was already known to handle, which is why a whole dead code path
  # had green tests, and why the routing fix could newly admit messages the
  # block answers wrongly without anything failing.
  #
  # `ok` is what the guidance must say; `never` is what it must not. The
  # standing invariant is the `never` column on the two success-path messages:
  # guidance must not tell a user their standard errors are unavailable when the
  # engine produced them.
  cases <- list(
    list(cat = "covariance_failed",
         # covariance.rs format_non_pd_warning
         msg = paste0("Covariance step: Hessian is not positive definite. ",
                      "Eigenvalues: [-1.2e-03, 4.5e-01]. SE estimates not available."),
         ok = "eigenvalue list", never = NULL),
    list(cat = "covariance_failed",
         # covariance.rs ill-conditioned entries
         msg = paste0("Covariance step failed: Hessian has ill-conditioned entries ",
                      "for the following parameter(s) -- theta[CL] (non-finite ",
                      "diagonal). SE estimates not available."),
         ok = "Hessian diagonal", never = NULL),
    list(cat = "covariance_failed",
         # covariance.rs omega branch, descriptor = "not positive definite"
         msg = paste0("Covariance step failed: Omega matrix is not positive definite ",
                      "at convergence (min eigenvalue = -1.0e-09; eigenvalues: ",
                      "[-1.0e-09, 9.0e-02]). SE estimates not available."),
         ok = "Omega is near-singular", never = NULL),
    list(cat = "covariance_failed",
         # SAME emit site, descriptor = "near-singular". Pinning the grep to the
         # other descriptor sent this one to the generic fallback.
         msg = paste0("Covariance step failed: Omega matrix is near-singular at ",
                      "convergence (min eigenvalue = 1.0e-12; eigenvalues: ",
                      "[1.0e-12, 9.0e-02]). SE estimates not available."),
         ok = "Omega is near-singular", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: base OFV is non-finite at convergence ",
                      "(likely numerical overflow or underflow in model evaluation). ",
                      "SE estimates not available."),
         ok = "overflowed", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: fd_hessian_step must be positive and ",
                      "finite, got 0. SE estimates not available."),
         ok = "fd_hessian_step", never = "identifiability"),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: the score cross-product matrix S is ",
                      "singular or rank-deficient (covariance_method = s); typically ",
                      "fewer subjects than free parameters, or collinear per-subject ",
                      "scores. Use covariance_method = r or rsr. SE estimates not ",
                      "available."),
         ok = "covariance_method", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: could not compute eigenvalues of the ",
                      "FD Hessian (Hessian may contain NaN or Inf). SE estimates ",
                      "not available."),
         ok = "identifiability", never = NULL),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (1 of 3 free-block eigenvalues clipped; min eig = ",
                      "1.2e-06, floor = 8.4e-14; severity: minor). Standard errors ",
                      "are likely reliable."),
         ok = "benign", never = "ferx_sir"),
    list(cat = "covariance_regularized",
         # SUCCESS path -- the covariance matrix exists.
         msg = paste0("Covariance step: off-diagonal FD stencil(s) non-finite for ",
                      "theta[CL], sigma[1]. Cross-partial correlation set to 0; SE ",
                      "for these parameter(s) may be over-optimistic. Try tuning ",
                      "fd_hessian_step."),
         ok = "over-optimistic", never = "unavailable"),
    list(cat = "covariance_step",
         # Info-level cost note -- the step has not even run.
         msg = "Covariance step: 35 parameters -> n^2 OFV evaluations",
         ok = "Informational", never = "unavailable"),
    list(cat = "general",
         # COV_CANCELLED_MSG. Classifies to `general`, so it reaches the block
         # only because `general` is admitted on the message -- which means the
         # routing that admits it also owes it an answer. "Check identifiability"
         # is wrong here: the user cancelled, nothing was diagnosed.
         msg = paste0("Covariance step cancelled before completion; standard ",
                      "errors not available."),
         ok = "cancelled", never = "identifiability")
  )

  for (case in cases) {
    g <- ferx:::.ferx_warning_guidance(case$cat, message = case$msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = case$msg)
    expect_match(g, case$ok, ignore.case = TRUE, info = case$msg)
    if (!is.null(case$never)) {
      expect_false(grepl(case$never, g, ignore.case = TRUE), info = case$msg)
    }
    # Same message under ferx_covariance()'s token must give the same answer.
    expect_identical(ferx:::.ferx_warning_guidance("covariance", message = case$msg),
                     ferx:::.ferx_warning_guidance(case$cat, message = case$msg),
                     info = case$msg)
  }
})

test_that(".ferx_warning_guidance returns NULL for an unknown category", {
  expect_null(ferx:::.ferx_warning_guidance("not_a_real_category"))
})

# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   ferx_estimates(), .ferx_est_row(), ferx_eta_cov(), ferx_cor_matrix(),
#   ferx_get_warnings(), and .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row       <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm  <- getFromNamespace(".ferx_compute_eta_normality", "ferx")

# ---------------------------------------------------------------------------
# ferx_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_eta_cov — message branches plus the correlation path
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# ferx_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that("ferx_get_warnings prints 'No warnings.' for an empty structured table", {
  empty <- data.frame(severity = character(0), category = character(0),
                      message = character(0), source_method = character(0),
                      stringsAsFactors = FALSE)
  fit <- make_fake_fit(warnings_structured = empty)
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("No warnings", out)))
})
test_that("ferx_get_warnings labels each severity level", {
  df <- data.frame(
    severity      = c("critical", "warning", "info"),
    category      = c("convergence", "condition_number", "mu_referencing"),
    message       = c("did not converge", "ill-conditioned", "mu detected"),
    source_method = c("focei", "focei", "focei"),
    stringsAsFactors = FALSE
  )
  fit <- make_fake_fit(warnings_structured = df)
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})

# ferx-core#1021: the `sir` category now dispatches on message content. The
# proposal-conditioning diagnostics are a statement about the model (parameters
# the data do not identify), not about SIR tuning, so they must not be answered
# with "raise sir_samples".
test_that(".ferx_warning_guidance dispatches sir by message content", {
  g <- function(msg) ferx:::.ferx_warning_guidance("sir", message = msg)

  rank_def <- g(paste0(
    "SIR: proposal covariance is rank-deficient beyond the FIX-ed parameters: ",
    "1 direction(s) carry no uncertainty [TVCL +0.71, TVV -0.70]."
  ))
  expect_match(rank_def, "not identified", fixed = TRUE)
  expect_false(grepl("sir_samples", rank_def, fixed = TRUE))

  shrunk <- g(paste0(
    "SIR: proposal was shrunk in 3 direction(s) so draws stay inside the ",
    "parameter bounds [TVKA -0.74 (sd 6.15e3 -> 2.89e0)]."
  ))
  expect_match(shrunk, "understate", fixed = TRUE)
  expect_false(grepl("sir_samples", shrunk, fixed = TRUE))

  # Availability / tuning messages keep the original advice.
  tuning <- g("SIR requested but covariance matrix is not available")
  expect_match(tuning, "sir_samples", fixed = TRUE)
})
