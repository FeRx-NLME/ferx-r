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
  # The flat strings are re-classified by the engine's own classifier (#308),
  # so a recognised message gets the severity and category a fresh fit would
  # have carried, and an unrecognised one stays `general`.
  expect_identical(df$severity, c("critical", "warning"))
  expect_identical(df$category, c("convergence", "general"))
  expect_identical(df$source_method, c("", ""))
})

test_that("the flat-warning fallback splits a [METHOD] prefix like a fresh fit", {
  fake <- structure(
    list(model_name = "legacy",
         warnings = "[FOCEI] Outer optimization did not converge",
         warnings_structured = NULL),
    class = "ferx_fit"
  )
  df <- ferx_get_warnings(fake, as_df = TRUE)
  expect_identical(df$source_method, "FOCEI")
  expect_identical(df$message, "Outer optimization did not converge")
  expect_identical(df$category, "convergence")
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
  # Every parameter in this message carries one of ferx-core's two cause
  # labels; "non-finite diagonal" was never one of them.
  msg_ic <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) \u2014 theta[CL] (FD stencil non-finite; model may ",
    "overflow at perturbation \u2014 try tuning fd_hessian_step). ",
    "SE estimates not available."
  )
  expect_match(g(msg_ic), "fd_hessian_step", ignore.case = TRUE)

  # Omega non-PD.
  # A tiny POSITIVE minimum eigenvalue: ferx-core writes "near-singular" for
  # this case, not "not positive definite" (covariance.rs picks the descriptor
  # from the sign), so the fixture must too.
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
  # Anchored on the tier word itself. The previous alternations did not
  # discriminate: "unreliable" contains "reliable", so the severe text matched
  # the minor pattern, and the moderate text matched the severe pattern --
  # swapping the moderate and severe return values kept the suite green.
  expect_match(reg("minor"),    "Minor Hessian regularisation",    fixed = TRUE)
  expect_match(reg("moderate"), "Moderate Hessian regularisation", fixed = TRUE)
  expect_match(reg("severe"),   "Severe Hessian regularisation",   fixed = TRUE)
  # Each tier's own word, and nobody else's.
  expect_false(grepl("Moderate|Severe", reg("minor")))
  expect_false(grepl("Minor|Severe",    reg("moderate")))
  expect_false(grepl("Minor|Moderate",  reg("severe")))
  # Only minor is benign enough to omit the SIR cross-check.
  expect_false(grepl("ferx_sir", reg("minor"), ignore.case = TRUE))
  expect_match(reg("moderate"), "ferx_sir", fixed = TRUE)
  expect_match(reg("severe"),   "ferx_sir", fixed = TRUE)
  # These are count-graded messages (pre ferx-core #1508): the tier said what
  # FRACTION of the eigenvalues were clipped, not how far the SEs moved, so the
  # guidance must not volunteer how common or how benign "minor" is.
  expect_false(grepl("benign|no action", reg("minor"), ignore.case = TRUE))
  expect_match(reg("minor"), "%RSE", fixed = TRUE)

  # Generic fallback for unrecognised message.
  expect_match(g("Covariance step failed"), "identifiability", ignore.case = TRUE)

})

test_that("magnitude-graded regularisation guidance gives the action, not the mechanism (#395)", {
  # Message shape from ferx-core #1508 (cov_diagnostics.rs
  # format_regularized_warning): the tier is an OR over the worst variance
  # inflation and |min eig| / max eig, and core's own sentence already says
  # which leg fired. `interp` is that sentence.
  reg_msg <- function(sev, interp, tail = "", source = "the FD Hessian") {
    paste0(
      "Covariance step regularized: eigenvalue floor applied to ", source,
      " (1 of 13 free-block eigenvalues clipped; min eig = -3.100e-2, ",
      "max eig = 4.200e+3, |min eig|/max eig = 7.38e-6, floor = 4.200e-11; ",
      "worst inflation of a reported variance = 4.400e3x; severity: ", sev,
      "). ", interp, tail
    )
  }
  g <- function(msg) ferx:::.ferx_warning_guidance("covariance_regularized",
                                                   message = msg)
  sev_infl <- paste0("Standard errors for the affected parameters come mostly ",
                     "from the floor rather than from the data and are not ",
                     "reliable; SIR-based confidence intervals are recommended.")
  sev_indef <- paste0("The Hessian was materially altered by the floor and ",
                      "these standard errors are not reliable; SIR-based ",
                      "confidence intervals are recommended.")
  mod <- paste0("Part of the reported standard errors for the affected ",
                "parameters comes from the floor rather than from the data.")

  # Severe: same action on both legs, and no claim that the SEs "come from the
  # floor" -- false on the indefiniteness-only cell.
  for (interp in c(sev_infl, sev_indef)) {
    s <- g(reg_msg("severe", interp))
    expect_match(s, "Severe Hessian regularisation", fixed = TRUE)
    expect_match(s, "ferx_sir()", fixed = TRUE)
    expect_false(grepl("from the floor|mostly", s, ignore.case = TRUE))
  }
  expect_identical(g(reg_msg("severe", sev_infl)), g(reg_msg("severe", sev_indef)))

  # Moderate now fires on fits that used to print minor: "worth a look", not a
  # failure, and still pointing at SIR.
  m <- g(reg_msg("moderate", mod))
  expect_match(m, "Moderate Hessian regularisation", fixed = TRUE)
  expect_match(m, "worth a look", fixed = TRUE)
  expect_match(m, "ferx_sir()", fixed = TRUE)

  # Minor under magnitude grading states the measured bound and asks for
  # nothing; the old hedge is gone.
  n <- g(reg_msg("minor", "Standard errors are likely reliable.",
                 source = "the analytic R-matrix"))
  expect_match(n, "Minor Hessian regularisation", fixed = TRUE)
  expect_match(n, "no action", fixed = TRUE)
  expect_match(n, "less than 1%", fixed = TRUE)
  expect_false(grepl("benign|ferx_sir|%RSE", n))

  # [scaling] obs_scale named as the declining clause: the one-line rewrite,
  # and whether it alone moves the route follows what core says.
  decl <- paste0(" The exact analytic covariance R-matrix was declined because ",
                 "[scaling] obs_scale = ... is in use.")
  alone <- g(reg_msg("severe", sev_infl, paste0(
    decl, " Writing the readout as an explicit expression ([scaling] y = ",
    "central / V) instead of obs_scale moves the fit onto the analytic route.")))
  expect_match(alone, "[scaling] y = central / V", fixed = TRUE)
  expect_match(alone, "moves the fit onto the analytic route", fixed = TRUE)
  expect_match(alone, "Severe Hessian regularisation", fixed = TRUE)

  blocked <- g(reg_msg("severe", sev_infl, paste0(
    " The exact analytic covariance R-matrix was declined because [scaling] ",
    "obs_scale = ... is in use and the model is a mixture model. Writing the ",
    "readout as an explicit expression ([scaling] y = central / V) instead of ",
    "obs_scale clears that clause, but the remaining clause has no one-line ",
    "remedy, so the fit stays on the finite-difference route until all of them ",
    "are cleared.")))
  expect_match(blocked, "[scaling] y = central / V", fixed = TRUE)
  expect_match(blocked, "keep the fit on the finite-difference route", fixed = TRUE)
  expect_false(grepl("moves the fit onto", blocked, fixed = TRUE))

  together <- g(reg_msg("moderate", mod, paste0(
    " The exact analytic covariance R-matrix was declined because gradient = ",
    "fd and [scaling] obs_scale = ... is in use. Dropping gradient = fd and ",
    "writing the readout as an explicit expression ([scaling] y = central / V) ",
    "instead of obs_scale together move the fit onto the analytic route.")))
  expect_match(together, "together with the other changes", fixed = TRUE)

  # ODE tolerances looser than the FD stencil's plateau.
  ode <- g(reg_msg("moderate", mod, paste0(
    " Note that this model integrates ODEs at ode_reltol = 1e-4 / ode_abstol ",
    "= 1e-6, and the FD covariance stencil amplifies integration noise by ",
    "1/h\u00b2; ode_reltol = 1e-6 / ode_abstol = 1e-8 is on the measured ",
    "accuracy plateau for this stencil (#520).")))
  expect_match(ode, "settings = list(ode_reltol = 1e-6, ode_abstol = 1e-8)",
               fixed = TRUE)
  # Neither extra appears unless core named it.
  expect_false(grepl("obs_scale|ode_reltol", g(reg_msg("severe", sev_infl))))
  # A closed-form twin gets the same pointer.
  twin <- g(reg_msg("minor", "Standard errors are likely reliable.", paste0(
    " Note that this model's closed-form absorption ODE twin integrates at ",
    "ode_reltol = 1e-4 / ode_abstol = 1e-6, and the FD covariance stencil ",
    "amplifies integration noise by 1/h2; ode_reltol = 1e-6 / ode_abstol = ",
    "1e-8 is on the measured accuracy plateau for this stencil (#520).")))
  expect_match(twin, "ode_reltol = 1e-6", fixed = TRUE)
  expect_match(twin, "Minor Hessian regularisation", fixed = TRUE)
  # Pure FD route: the rewrite is the whole fit's route change.
  expect_false(grepl("finite-differenced subjects", alone, fixed = TRUE))
})

test_that("hybrid analytic/FD regularisation keeps the route change on those subjects", {
  # ferx-core #1514 / #1516: on the hybrid route the declining clauses belong
  # to the finite-differenced subjects only; everyone else is already
  # analytic, so "moves the fit" / "keeps the fit on FD" would be false.
  # Exact shapes from cov_diagnostics.rs decline_sentences().
  hyb <- function(tail) paste0(
    "Covariance step regularized: eigenvalue floor applied to the hybrid ",
    "analytic/FD R-matrix (1 of 13 free-block eigenvalues clipped; min eig = ",
    "-3.100e-2, max eig = 4.200e+3, |min eig|/max eig = 7.38e-6, floor = ",
    "4.200e-11; worst inflation of a reported variance = 4.400e3x; severity: ",
    "severe). Standard errors for the affected parameters come mostly from the ",
    "floor rather than from the data and are not reliable; SIR-based confidence ",
    "intervals are recommended.", tail)
  g <- function(msg) ferx:::.ferx_warning_guidance("covariance_regularized",
                                                   message = msg)
  alone <- g(hyb(paste0(
    " The finite-differenced subjects declined the exact analytic covariance ",
    "R-matrix because [scaling] obs_scale = ... is in use. Writing the readout ",
    "as an explicit expression ([scaling] y = central / V) instead of obs_scale ",
    "moves those subjects onto the analytic route.")))
  expect_match(alone, "moves the finite-differenced subjects onto the analytic route",
               fixed = TRUE)
  expect_false(grepl("the fit", alone, fixed = TRUE))

  together <- g(hyb(paste0(
    " The finite-differenced subjects declined the exact analytic covariance ",
    "R-matrix because gradient = fd and [scaling] obs_scale = ... is in use. ",
    "Dropping gradient = fd and writing the readout as an explicit expression ",
    "([scaling] y = central / V) instead of obs_scale together move those ",
    "subjects onto the analytic route.")))
  expect_match(together, "moves the finite-differenced subjects onto the analytic route",
               fixed = TRUE)
  expect_false(grepl("the fit", together, fixed = TRUE))

  blocked <- g(hyb(paste0(
    " The finite-differenced subjects declined the exact analytic covariance ",
    "R-matrix because [scaling] obs_scale = ... is in use and the model is a ",
    "mixture model. Writing the readout as an explicit expression ([scaling] ",
    "y = central / V) instead of obs_scale clears that clause, but the ",
    "remaining clause has no one-line remedy, so those subjects stay on the ",
    "finite-difference route until all of them are cleared.")))
  expect_match(blocked,
               "keep the finite-differenced subjects on the finite-difference route",
               fixed = TRUE)
  expect_false(grepl("the fit", blocked, fixed = TRUE))
})

test_that("off-diagonal FD stencil guidance distinguishes pure-FD from hybrid", {
  g <- function(msg) ferx:::.ferx_warning_guidance("covariance_regularized",
                                                   message = msg)
  fd <- g(paste0(
    "Covariance step: off-diagonal FD stencil(s) non-finite for theta[CL]. ",
    "Cross-partial correlation set to 0; SE for these parameter(s) may be ",
    "over-optimistic. Try tuning fd_hessian_step."))
  expect_match(fd, "were set to zero", fixed = TRUE)
  expect_false(grepl("finite-differenced subjects", fd, fixed = TRUE))

  hy <- g(paste0(
    "Covariance step: off-diagonal FD stencil(s) non-finite for theta[CL]. ",
    "Those cross-partials keep only the analytically assembled subjects' ",
    "contribution \u2014 the finite-differenced subjects' share of them is ",
    "missing \u2014 so SE for these parameter(s) may be over-optimistic. Try ",
    "tuning fd_hessian_step."))
  expect_match(hy, "finite-differenced subjects' share", fixed = TRUE)
  expect_match(hy, "analytically assembled subjects' share is kept", fixed = TRUE)
  expect_false(grepl("set to zero", hy, fixed = TRUE))
  for (x in c(fd, hy)) {
    expect_match(x, "Standard errors were produced", fixed = TRUE)
    expect_match(x, "fd_hessian_step", fixed = TRUE)
  }
})

test_that("W_COV_ANALYTIC_SALVAGE is informational, not a failure (ferx-core #1514)", {
  # ferx-core classifies it (Info, CovarianceStep) on its token.
  msg <- paste0(
    "W_COV_ANALYTIC_SALVAGE: 2 of 30 subjects (IDs 7 and 12) are outside the ",
    "exact analytic covariance R-matrix scope; their information terms were ",
    "finite-differenced from their own marginals, and the remaining 28 subjects ",
    "were assembled analytically. Each subject contributes its own term to the ",
    "information matrix, so only the named subjects' terms use a different ",
    "estimator.")
  g <- ferx:::.ferx_warning_guidance("covariance_step", message = msg)
  expect_match(g, "Informational", fixed = TRUE)
  expect_match(g, "complete", fixed = TRUE)
  expect_match(g, "only their terms were finite-differenced", fixed = TRUE)
  expect_false(grepl("unavailable|identifiability", g, ignore.case = TRUE))
  expect_identical(ferx:::.ferx_warning_guidance("covariance", message = msg), g)
})

test_that("the informational covariance_step note does not read as a failure", {
  # Core codes only its Info-level cost note `covariance_step`. It matches none
  # of the targeted branches, so before this arm existed it fell through to the
  # "standard errors unavailable" fallback and told the user the step had failed
  # when it had not yet run.
  # ferx-core's primary cost note, verbatim (api/fit.rs).
  msg <- paste0("Covariance step: 35 parameters \u2192 1225 OFV evaluations ",
                "(finite-difference Hessian). This may take several minutes on ",
                "complex models.")
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
  # Not the benign cost note (the message is not the cost note) and not the
  # failure fallback either: `covariance_step` is Info by construction, so
  # "standard errors unavailable" would contradict the row it prints under.
  # Not the benign cost note ...
  # Not the benign cost note: that arm matches on the message, never on the
  # code alone, because `covariance_step`'s classification ends in a catch-all.
  # An unrecognised message under it falls through to the generic fallback.
  expect_false(grepl("No action needed", g, fixed = TRUE))
  expect_match(g, "identifiability", ignore.case = TRUE)
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

test_that("ferx_covariance()'s `covariance` category reaches the guidance", {
  # `ferx_covariance()` folds the engine's flat covariance warnings into the
  # structured table under `category = "covariance"` (R/ferx_covariance.R) --
  # a fourth token, not one of core's three. Omitting it left the entire
  # post-hoc covariance surface without guidance even after the routing fix.
  msg <- paste0(
    "Covariance step failed: Hessian has ill-conditioned entries for the ",
    "following parameter(s) \u2014 theta[CL] (zero diagonal \u2014 flat objective)."
  )
  g <- ferx:::.ferx_warning_guidance("covariance", message = msg)
  expect_match(g, "Hessian diagonal", fixed = TRUE)
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
  expect_match(g, "Omega is near-singular", fixed = TRUE)
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
  # Collapsed and whitespace-normalised before matching: the printer wraps at
  # 70 columns, so asserting on a phrase in a single output line passes or fails
  # on where the wrap happens to fall rather than on the behaviour under test.
  out <- capture.output(ferx_get_warnings(fake))
  flat <- gsub("[[:space:]]+", " ", paste(out, collapse = " "))
  expect_true(
    grepl("Remove it from", flat, fixed = TRUE),
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




test_that(".ferx_warning_guidance never recommends SDE process noise for positive autocorrelation", {
  pos <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "positive")
  # The remedies the guidance does still name, so this is not a pure negative.
  expect_match(pos, "Positive IWRES autocorrelation")
  expect_match(pos, "transit absorption", fixed = TRUE)
  expect_match(pos, "an extra compartment", fixed = TRUE)
  expect_match(pos, "IOV on ka/F", fixed = TRUE)
  # ferx-core #1426 dropped the equivalent sentence from the engine's own
  # Durbin-Watson warning; R must not contradict it.
  expect_false(grepl("SDE", pos, fixed = TRUE))
  expect_false(grepl("diffusion", pos, fixed = TRUE))
  # The hint used to be toggled by a `uses_sde` argument, so it reached every
  # ODE model without a [diffusion] block - exactly the population it should no
  # longer reach. There is no longer a side of that toggle that shows it.
  expect_false("uses_sde" %in% names(formals(ferx:::.ferx_warning_guidance)))
})
test_that(".ferx_warning_guidance returns negative-autocorrelation guidance", {
  neg <- ferx:::.ferx_warning_guidance("dw_autocorrelation", "Negative autocorrelation")
  expect_match(neg, "Negative IWRES autocorrelation")
})
# Every token ferx-core's `WarningCode::as_str()` can emit, minus `general`.
# Hand-transcribed, and deliberately so: `.core_warning_cats()` below
# checks it against the real thing whenever a sibling ferx-core checkout is
# present, which is the only mechanism here that can notice core growing a code.
.core_warning_cats <- function() {
  c(
    "absorption_twin_declined", "bloq_method", "boundary_estimate", "cancelled",
    "condition_number", "convergence", "covariance_failed",
    "covariance_regularized", "covariance_step", "data_quality",
    "dw_autocorrelation", "ebe_start_dependent", "eps_shrinkage",
    "eta_normality", "eta_shrinkage",
    "experimental", "flat_parameter", "flip_flop", "gradient_fallback",
    "high_correlation", "importance_sampling", "inflated_rse",
    "init_outside_bounds", "mu_referencing",
    "multi_start", "ode_solver", "omega_structure", "optimizer_config",
    "optimizer_health", "parameter_at_runaway_guard", "simulation", "sir",
    "stalled_at_init", "threads", "vi_bad_basin"
  )
}

# Categories that reach `.ferx_warning_guidance()` without being ferx-core
# `WarningCode` tokens.
#
#   covariance        - assigned by ferx_covariance() (R/ferx_covariance.R) to
#                       the engine's flat covariance warnings.
#   unused_parameter  - no core `WarningCode`; both of core's unused-declaration
#                       messages classify to `general` and are rerouted here by
#                       message, so the arm is reachable (see the test below
#                       that drives it with core's real message texts).
#
# `ebe_convergence` used to sit here too and has been removed: nothing emits it.
# Core has no such `WarningCode`, and ferx-r reports EBE convergence as the
# integer field `ebe_convergence_warnings`, never as a warning row. Keeping a
# guidance arm for it in a change whose premise is that unreachable branches are
# a defect would have been the same mistake in miniature.
.extra_warning_cats <- function() {
  c("covariance", "unused_parameter")
}

# ferx-core tokens this package deliberately answers with nothing. Empty since
# #308 answered the last ten. Kept, rather than deleted, so a code that has to
# go unanswered for a while has one explicit place to be listed: a code
# ferx-core adds later lands in neither list and fails the drift test rather
# than silently printing nothing.
.unanswered_warning_cats <- function() {
  character(0)
}

test_that("every category in the guidance table returns guidance", {
  # A completeness walk, and nothing more. Two things it deliberately does NOT
  # establish, both covered elsewhere:
  #   - it is not a drift guard. The loop visits only the tokens listed above,
  #     so a code core added and nobody transcribed is invisible to it; that is
  #     what `matches ferx-core's WarningCode vocabulary` is for.
  #   - it says nothing about the covariance family. Called with an empty
  #     message, all four covariance tokens land on the generic fallback, so
  #     this would still pass if `.ferx_covariance_guidance()` were reduced to
  #     its last line. The message inventory below is what pins those branches.
  answered <- setdiff(c(.core_warning_cats(), .extra_warning_cats()),
                      .unanswered_warning_cats())
  for (cat in answered) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
  # ... and the deliberately unanswered ones really are unanswered, so the list
  # cannot rot into a place where handled codes hide.
  for (cat in .unanswered_warning_cats()) {
    expect_null(ferx:::.ferx_warning_guidance(cat), info = cat)
  }
})

# Locate a sibling ferx-core checkout. Walks upward rather than using a fixed
# relative depth: testthat's wd is tests/testthat, so `../../../ferx-core` only
# resolves for a plain checkout. ferx-r's CLAUDE.md mandates working in a git
# worktree under <repo>/.claude/worktrees/<name>, where that path points at
# <repo>/.claude/worktrees/ferx-core and never exists -- so the guards that
# depend on it silently skipped in exactly the workflow they were written for.
.find_core_src <- function() {
  d <- normalizePath(".", mustWork = FALSE)
  for (i in seq_len(10)) {
    cand <- file.path(d, "ferx-core", "src")
    if (dir.exists(cand)) return(cand)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NULL
}

# Materialise ferx-core's `src/` at the revision this package is PINNED to, and
# return that directory (or NULL).
#
# The sibling checkout's working tree is NOT the right source to check against:
# `Cargo.toml` tracks `branch = "main"`, but CI builds the revision recorded in
# `src/rust/Cargo.lock`, and a local checkout can sit either side of it. Reading
# the pin means these guards assert against the vocabulary and the message texts
# that actually ship. Falls back to NULL (skip) rather than silently reading a
# different revision.
.core_src_at_pin <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(if (identical(cached, NA_character_)) NULL else cached)
    fail <- function() { cached <<- NA_character_; NULL }
    core_src <- .find_core_src(); if (is.null(core_src)) return(fail())
    repo <- dirname(core_src)
    lock <- file.path("..", "..", "src", "rust", "Cargo.lock")
    if (!file.exists(lock)) return(fail())
    ln <- readLines(lock, warn = FALSE)
    i <- grep('^name = "ferx-core"$', ln)
    if (!length(i)) return(fail())
    src_line <- grep("ferx-core.*#", ln[i[1]:min(length(ln), i[1] + 5L)], value = TRUE)
    if (!length(src_line)) return(fail())
    rev <- sub('".*$', "", sub("^.*#", "", src_line[1]))
    if (!grepl("^[0-9a-f]{7,40}$", rev)) return(fail())
    have <- suppressWarnings(system2("git", c("-C", repo, "cat-file", "-e",
                                              paste0(rev, "^{commit}")),
                                     stdout = FALSE, stderr = FALSE))
    if (!identical(have, 0L)) return(fail())
    dest <- file.path(tempdir(), paste0("ferx-core-", substr(rev, 1, 12)))
    if (!dir.exists(dest)) {
      dir.create(dest, recursive = TRUE)
      tar <- file.path(tempdir(), paste0(substr(rev, 1, 12), ".tar"))
      ok <- suppressWarnings(system2("git", c("-C", repo, "archive", "--format=tar",
                                              "-o", tar, rev, "src"),
                                     stdout = FALSE, stderr = FALSE))
      if (!identical(ok, 0L)) return(fail())
      untar(tar, exdir = dest)
    }
    out <- file.path(dest, "src")
    if (!dir.exists(out)) return(fail())
    cached <<- out
    out
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
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")
  types_rs <- file.path(core_src, "types.rs")
  skip_if(!file.exists(types_rs), "sibling ferx-core has no src/types.rs")

  src <- readLines(types_rs, warn = FALSE)
  hits <- regmatches(src, regexpr('WarningCode::[A-Za-z]+ => "[a-z_]+"', src))
  tokens <- gsub('.*=> "|"$', "", hits)
  expect_gt(length(tokens), 20L)          # the file was found and understood

  # `general` is core's unrecognised-message bucket and is deliberately absent
  # from the table -- the message text is the only guidance there is.
  expect_setequal(setdiff(tokens, "general"), .core_warning_cats())
  for (cat in setdiff(tokens, c("general", .unanswered_warning_cats()))) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
  }
})

test_that("init_outside_bounds guidance is start-side and distinct from boundary_estimate", {
  # ferx-core #1251 gave a clamped START its own WarningCode rather than reusing
  # BoundaryEstimate, because that category drives three default-on rejection
  # filters (bootstrap's skip_estimate_near_boundary, reject_on_boundary, and
  # .ferx_boundary_detail() in check_strictness.R) and a start wearing it would
  # silently drop bootstrap replicates. The guidance must keep the two apart:
  # this one is about where the fit BEGAN.
  # Only the phrase anchor lives here. That the arm returns a non-empty string
  # at all is already covered by the completeness walk above (the token is in
  # .core_warning_cats()), and that `boundary_estimate` stays unanswered is
  # already covered by the .unanswered_warning_cats() loop - asserting either
  # again would just add another place to edit when the vocabulary moves.
  expect_match(ferx:::.ferx_warning_guidance("init_outside_bounds"),
               "before the first objective evaluation", fixed = TRUE)
})

test_that("a real fit produces an init_outside_bounds row that carries guidance", {
  # The end-to-end path, which the table test above cannot reach: it hands the
  # category in by hand, so it asserts the entry exists, never that any fit
  # produces it. This drives the whole chain - engine emits
  # W_INIT_OUTSIDE_BOUNDS, the glue maps the code to `init_outside_bounds`, the
  # row survives into fit$warnings_structured, and ferx_get_warnings() prints
  # the guidance under it.
  #
  # An omega over the variance rail is the cheapest trigger, and one of only two
  # things that can reach this arm at all (the other is a sigma over the SD
  # rail). A theta past the hidden cap cannot: that is the *error*
  # E_THETA_INIT_OUTSIDE_BOUNDS, which stops the fit before any warning exists.
  ex   <- ferx_example("warfarin")
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(0.134, 0.001, 10.0)",
    "  omega ETA_CL ~ 1e6",          # far above the optimizer's variance rail
    "  sigma PROP_ERR ~ 0.01",
    "",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=10.0, ka=1.0)",
    "",
    "[error_model]",
    "  DV ~ proportional(PROP_ERR)"
  ), path)
  on.exit(unlink(path))

  fit <- ferx_fit(path, ex$data, settings = list(maxiter = 2L),
                  covariance = FALSE, verbose = FALSE)
  ws <- ferx_get_warnings(fit, as_df = TRUE)
  expect_true("init_outside_bounds" %in% ws$category)

  # The guidance must actually print, not merely exist in the table.
  out <- capture.output(ferx_get_warnings(fit))
  expect_true(any(grepl("clamped onto one of the optimizer's internal rails",
                        out, fixed = TRUE)))
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
    "[individual_parameters] `KA` is computed but never used \u2014 not mapped ",
    "into the `pk(...)` model and not referenced in any other block, so it has ",
    "no effect. Map `KA` in [structural_model] (e.g. `f=F`) or remove `KA`."
  )
  for (msg in c(theta_msg, ip_msg)) {
    g <- ferx:::.ferx_warning_guidance("general", message = msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = msg)
    expect_match(g, "never referenced", fixed = TRUE, info = msg)
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
  # Each `msg` is ferx-core's own text with its `{}` placeholders filled in --
  # em dashes, arrows and the superscript are the real characters, not ASCII
  # stand-ins. `ok` is what the guidance must say; `never` is what it must not.
  # The standing invariant is the `never` column on the success-path messages:
  # guidance must not report standard errors unavailable when the engine
  # produced them.
  # Each `msg` is ferx-core's own text with its `{}` placeholders filled in.
  # `ok` is what the guidance must say; `never` what it must not. The standing
  # invariant is the `never` column on the success-path messages: guidance must
  # not report standard errors unavailable when the engine produced them.
  cases <- list(
    list(cat = "covariance_failed",
         msg = paste0("Covariance step: Hessian is not positive definite. ",
                      "Eigenvalues: [-0.0012, 0.4500]. SE estimates not available."),
         ok = "eigenvalue list", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: Hessian has ill-conditioned entries ",
                      "for the following parameter(s) \u2014 theta[CL] (FD stencil ",
                      "non-finite; model may overflow at perturbation \u2014 try ",
                      "tuning fd_hessian_step). SE estimates not available."),
         ok = "Hessian diagonal", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: Omega matrix is not positive definite ",
                      "at convergence (min eigenvalue = -1.000e-9; eigenvalues: ",
                      "[-1.000e-9, 0.0900]). SE estimates not available."),
         ok = "Omega is near-singular", never = NULL),
    list(cat = "covariance_failed",
         # Same emit site, the other descriptor. ferx-core picks it from the
         # sign of the smallest eigenvalue, so both must reach the omega arm
         # rather than one of them falling to the generic fallback.
         msg = paste0("Covariance step failed: Omega matrix is near-singular at ",
                      "convergence (min eigenvalue = 1.000e-12; eigenvalues: ",
                      "[1.000e-12, 0.0900]). SE estimates not available."),
         ok = "Omega is near-singular", never = "identifiability"),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: base OFV is non-finite at convergence ",
                      "(likely numerical overflow or underflow in model evaluation). ",
                      "SE estimates not available."),
         ok = "overflow", never = NULL),
    list(cat = "covariance_failed",
         msg = paste0("Covariance step failed: could not compute eigenvalues of the ",
                      "FD Hessian (Hessian may contain NaN or Inf). SE estimates ",
                      "not available."),
         ok = "identifiability", never = NULL),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (1 of 3 free-block eigenvalues clipped; min eig = ",
                      "1.200e-6, floor = 8.400e-14; severity: minor). Standard ",
                      "errors are likely reliable."),
         ok = "Minor Hessian regularisation", never = "ferx_sir"),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (2 of 3 free-block eigenvalues clipped; min eig = ",
                      "1.200e-8, floor = 8.400e-14; severity: moderate). Standard ",
                      "errors should be interpreted with caution."),
         ok = "Moderate Hessian regularisation", never = "Severe"),
    list(cat = "covariance_regularized",
         msg = paste0("Covariance step regularized: eigenvalue floor applied to FD ",
                      "Hessian (3 of 3 free-block eigenvalues clipped; min eig = ",
                      "-1.200e-8, floor = 8.400e-14; severity: severe). SIR-based ",
                      "confidence intervals are recommended."),
         ok = "Severe Hessian regularisation", never = "Moderate"),
    list(cat = "covariance_regularized",
         # SUCCESS path -- the covariance matrix exists.
         msg = paste0("Covariance step: off-diagonal FD stencil(s) non-finite for ",
                      "theta[CL], sigma[1]. Cross-partial correlation set to 0; SE ",
                      "for these parameter(s) may be over-optimistic. Try tuning ",
                      "fd_hessian_step."),
         ok = "over-optimistic", never = "unavailable"),
    list(cat = "covariance_regularized",
         # Hybrid analytic/FD route (ferx-core #1514): only the
         # finite-differenced subjects' share is missing.
         msg = paste0("Covariance step: off-diagonal FD stencil(s) non-finite for ",
                      "theta[CL]. Those cross-partials keep only the analytically ",
                      "assembled subjects' contribution \u2014 the finite-differenced ",
                      "subjects' share of them is missing \u2014 so SE for these ",
                      "parameter(s) may be over-optimistic. Try tuning ",
                      "fd_hessian_step."),
         ok = "over-optimistic", never = "set to zero"),
    list(cat = "covariance_step",
         # Informational salvage note after a successful hybrid covariance.
         msg = paste0("W_COV_ANALYTIC_SALVAGE: 1 of 30 subjects (ID 7) is outside ",
                      "the exact analytic covariance R-matrix scope; its information ",
                      "term was finite-differenced from its own marginal, and the ",
                      "remaining 29 subjects were assembled analytically. Each ",
                      "subject contributes its own term to the information matrix, ",
                      "so only the named subjects' terms use a different estimator."),
         ok = "Informational", never = "unavailable"),
    list(cat = "covariance_step",
         # Info-level cost note, primary form -- the step has not even run.
         msg = paste0("Covariance step: 35 parameters \u2192 1225 OFV evaluations ",
                      "(finite-difference Hessian). This may take several minutes ",
                      "on complex models."),
         ok = "Informational", never = "unavailable"),
    list(cat = "covariance_step",
         # Second cost-note form: the evaluation count overflows usize.
         msg = paste0("Covariance step: 4294967296 parameters \u2192 n\u00b2 OFV ",
                      "evaluations (finite-difference Hessian). Estimate exceeds ",
                      "usize range; expect this to be very slow."),
         ok = "Informational", never = "unavailable"),
    list(cat = "general",
         # COV_CANCELLED_MSG, which classifies to `general` and so reaches the
         # block only because `general` is admitted on the message.
         msg = paste0("Covariance step cancelled before completion; standard ",
                      "errors not available."),
         ok = "nothing was diagnosed", never = "identifiability"),
    list(cat = "covariance_failed",
         # S-matrix path: one subject's quadrature score could not be
         # evaluated, and the message names that subject. The guidance has to
         # point at the subject rather than at the model, because the engine
         # diagnosed a record, not an identifiability problem.
         msg = paste0("Covariance step failed: could not obtain a converged, ",
                      "finite quadrature score for subject 42. SE estimates ",
                      "not available."),
         ok = "named subject", never = NULL),
    list(cat = "covariance_failed",
         # The same path, one level up: the summed cross-product is non-finite,
         # so no subject is named and the advice is about the estimates.
         msg = paste0("Covariance step failed: non-finite score cross-product. ",
                      "SE estimates not available."),
         ok = "score cross-product", never = NULL),
    list(cat = "covariance_failed",
         # An invalid step: the Hessian was never attempted, so nothing about
         # the model was diagnosed (#308).
         msg = paste0("Covariance step failed: fd_hessian_step must be positive ",
                      "and finite, got -1. SE estimates not available."),
         ok = "positive finite", never = "identifiability"),
    list(cat = "covariance_failed",
         # covariance_method = "s": S itself is singular, and the remedy is an
         # estimator that does not invert it (#308).
         msg = paste0("Covariance step failed: the score cross-product matrix S ",
                      "is singular or rank-deficient (covariance_method = s); ",
                      "typically fewer subjects than free parameters, or ",
                      "collinear per-subject scores. Use covariance_method = r ",
                      "or rsr. SE estimates not available."),
         ok = "covariance_method = \"rsr\"", never = "identifiability")
  )
  for (case in cases) {
    g <- ferx:::.ferx_warning_guidance(case$cat, message = case$msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = case$msg)
    expect_match(g, case$ok, ignore.case = TRUE, info = case$msg)
    if (!is.null(case$never)) {
      expect_false(grepl(case$never, g, ignore.case = TRUE), info = case$msg)
    }
    # Same message under ferx_covariance()'s token must give the same answer --
    # true for every case here because all of them hit a targeted branch. It is
    # deliberately NOT true in the fallback, where the code carries information
    # the message does not; see "the fallback does not report failure under a
    # code that means success".
    expect_identical(ferx:::.ferx_warning_guidance("covariance", message = case$msg),
                     ferx:::.ferx_warning_guidance(case$cat, message = case$msg),
                     info = case$msg)
  }
})


test_that("a flat THETA is not answered as an unused parameter, by either door", {
  # ferx-core's flat-theta message contains the literal phrase "computed but
  # never used", so the unused-parameter reroute captured it. Gating the reroute
  # on `general` did NOT fix that: `ferx_load_fit()` does not restore
  # `warnings_structured`, so on a loaded fit every row arrives under `general`
  # and the gate becomes the admission criterion rather than an exclusion. The
  # message must be answered the same way through both doors.
  msg <- paste0("[parameters] `TVCL` has no effect on the objective (gradient ",
                "\u2248 0 at the initial estimate) \u2014 it is likely computed but ",
                "never used (unmapped, or dropped from the structural / scaling ",
                "model). Freezing it at its initial value (1.5) so the remaining ",
                "parameters can be estimated; map or remove `TVCL` to silence this.")
  # Through the `general` door the engine's classifier recovers
  # `flat_parameter`, so both doors now give the flat-THETA answer - and
  # neither gives the unused-parameter one.
  loaded <- ferx:::.ferx_warning_guidance("general", message = msg)
  fresh  <- ferx:::.ferx_warning_guidance("flat_parameter", message = msg)
  expect_identical(loaded, fresh)
  expect_match(loaded, "not estimated", fixed = TRUE)
  expect_false(identical(loaded, ferx:::.ferx_warning_guidance("unused_parameter")))

  # The genuine unused-declaration messages still resolve, and a live fit's own
  # typed code still wins over anything the text might suggest.
  unused <- paste0("theta 'TVX' is declared in [parameters] but not referenced in ",
                   "any model expression \u2014 it will not affect predictions")
  expect_match(ferx:::.ferx_warning_guidance("general", message = unused),
               "never referenced", fixed = TRUE)
  expect_false(grepl("never referenced",
                     ferx:::.ferx_warning_guidance("flat_parameter", message = unused),
                     fixed = TRUE))
})

test_that("`general` does not claim messages that only mention the covariance step", {
  # `ferx_load_fit()` does not restore `warnings_structured`, so EVERY row of a
  # fit read back from disk arrives under `general`. Admitting any message that
  # merely contains "covariance step" therefore handed SIR's own diagnostics the
  # covariance fallback -- which tells the user to pass covariance = FALSE, the
  # one setting that removes the matrix SIR needs.
  # Verbatim from sir.rs and postfit.rs -- both clauses these previously dropped
  # ("so draws mostly stay inside the parameter bounds", "see the covariance
  # warning above for the cause") are the parts that mention the covariance
  # step, so trimming them weakened the very thing under test.
  shrunk <- paste0("SIR: proposal was shrunk in 1 direction(s) so draws mostly stay ",
                   "inside the parameter bounds [KA +1.00]. Those directions come ",
                   "from eigenvalue-floored (non-identified) curvature in the ",
                   "covariance step; the SIR CIs along them understate the true ",
                   "uncertainty.")
  requested <- paste0("SIR requested but the covariance step did not succeed and no ",
                      "usable SIR proposal could be built from it, so SIR could not ",
                      "run \u2014 see the covariance warning above for the cause.")
  # Since #308 a `general` row is first recovered to the engine's own category,
  # so these now get SIR's guidance - which is the point: never the covariance
  # fallback and its covariance = FALSE advice.
  for (msg in c(shrunk, requested)) {
    g <- ferx:::.ferx_warning_guidance("general", message = msg)
    expect_identical(g, ferx:::.ferx_warning_guidance("sir", message = msg), info = msg)
    expect_false(grepl("covariance = FALSE", g, fixed = TRUE), info = msg)
  }
  # ... while a real covariance message under `general` is still answered, with
  # or without a [METHOD] chain prefix.
  cancelled <- "Covariance step cancelled before completion; standard errors not available."
  expect_match(ferx:::.ferx_warning_guidance("general", message = cancelled),
               "cancelled", ignore.case = TRUE)
  expect_match(ferx:::.ferx_warning_guidance("general", message = paste0("[FOCEI] ", cancelled)),
               "cancelled", ignore.case = TRUE)
})

test_that("no covariance message ferx-core emits is missing from the inventory", {
  # The completeness guard for the inventory test above.
  #
  # That test is only as good as the list of messages someone hand-collected for
  # it, and a hand-collected list cannot notice a message added later. This one
  # derives the list from ferx-core's source instead: every string literal that
  # *begins* with "Covariance step" -- which is the shape of every covariance
  # warning the engine emits -- must be accounted for, either by a case in the
  # inventory or by an explicit exemption naming why it is not a warning.
  #
  # Skipped without a sibling ferx-core checkout, for the same reason as the
  # WarningCode vocabulary test: ferx-r CI builds the pinned crate rather than a
  # checkout, and drift is introduced on a developer machine.
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")

  files <- list.files(core_src, pattern = "[.]rs$",
                      recursive = TRUE, full.names = TRUE)
  # Test-only sources: sibling `*_tests.rs` files and the `src/api/tests/` dir.
  files <- files[!grepl("_tests[.]rs$|/tests/", files)]
  lits <- character(0)
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    # Drop whole-line `//` comments before joining: ferx-core quotes its own
    # message text inside comments, and those quotes are not emit sites.
    ln <- ln[!grepl("^[[:space:]]*//", ln)]
    # Rust continues a string literal with a trailing backslash; rejoin so a
    # wrapped message is recovered as one literal.
    joined <- gsub("\\\\\n[[:space:]]*", "", paste(ln, collapse = "\n"))
    m <- gregexpr('"(?:[^"\\\\]|\\\\.)*"', joined, perl = TRUE)
    for (s in regmatches(joined, m)[[1]]) {
      body <- substr(s, 2, nchar(s) - 1)
      if (grepl("^Covariance step", body)) lits <- c(lits, body)
    }
  }
  lits <- sort(unique(lits))
  expect_gt(length(lits), 8L)          # the source was found and understood

  # Not warnings, so no guidance is owed. Keyed on a distinctive fragment.
  exempt <- c(
    # Inline `#[cfg(test)]` fixture in run_sir.rs (the file itself is not a
    # test-only source, so the path filter above cannot drop it).
    "matrix was not positive definite"
  )
  # Fragments the inventory test covers. Not one per case: the two Omega
  # descriptors share "Omega matrix is" and the two cost-note forms share
  # "OFV evaluations", so 13 fragments cover 17 cases.
  covered <- c(
    "cancelled before completion",
    "ill-conditioned entries",
    "Omega matrix is",
    "base OFV is non-finite",
    "could not compute eigenvalues",
    "regularized: eigenvalue floor",
    "Hessian is not positive definite",
    "off-diagonal FD stencil",
    "OFV evaluations",
    # The two S-matrix failures (ferx-core covariance.rs): one names the
    # subject whose quadrature score could not be evaluated, the other reports
    # the summed cross-product coming out non-finite.
    "quadrature score for subject",
    "non-finite score cross-product",
    # Given targeted arms in #308; answered by the generic fallback before.
    "fd_hessian_step must be",
    "score cross-product matrix"
  )
  accounted <- function(lit) any(vapply(c(exempt, covered),
                                        function(k) grepl(k, lit, fixed = TRUE),
                                        logical(1)))
  unaccounted <- lits[!vapply(lits, accounted, logical(1))]
  expect_identical(
    unaccounted, character(0),
    info = paste0(
      "ferx-core emits covariance message(s) with no case in the inventory test ",
      "and no exemption. Add a case to \"every covariance message ferx-core ",
      "emits gets non-contradictory guidance\", or exempt it here with the ",
      "reason it is not a warning:\n  ",
      paste(unaccounted, collapse = "\n  ")
    )
  )
})

# ---------------------------------------------------------------------------
# End-to-end anchor: the engine's own warnings, through the real classifier
# ---------------------------------------------------------------------------
# Every other test here feeds `.ferx_warning_guidance()` a message typed by
# hand. That exercises the dispatch but not the seam that actually matters --
# ferx-core emitting a message, `classify_warning` assigning it a category, and
# this package answering the pair. Nothing pinned that seam, and it is exactly
# where the routing defect this file exists to fix lived: the categories were
# hand-paired in the tests, so a pairing production never emits looked correct.
#
# The model is deliberately pathological: one theta declared and never
# referenced makes ferx-core emit both a real unused-declaration warning and a
# real covariance-step failure, on a fit that settles in a couple of seconds.
warning_anchor_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      mf <- tempfile(fileext = ".ferx")
      writeLines(c(
        "[parameters]",
        "  theta TVCL(0.134, 0.001, 10.0)",
        "  theta TVV(8.1, 0.1, 500.0)",
        "  theta TVKA(1.0, 0.01, 50.0)",
        "  theta TVUNUSED(1.0, 0.1, 10.0)",
        "  omega ETA_CL ~ 0.07",
        "  omega ETA_V  ~ 0.02",
        "  omega ETA_KA ~ 0.40",
        "  sigma PROP_ERR ~ 0.01 (sd)",
        "",
        "[individual_parameters]",
        "  CL = TVCL * exp(ETA_CL)",
        "  V  = TVV  * exp(ETA_V)",
        "  KA = TVKA * exp(ETA_KA)",
        "",
        "[structural_model]",
        "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
        "",
        "[error_model]",
        "  DV ~ proportional(PROP_ERR)"
      ), mf)
      ex <- ferx_example("warfarin")
      fit <<- ferx_fit(mf, ex$data, method = "gn", verbose = FALSE,
                       covariance = TRUE, settings = list(maxiter = 30L))
    }
    fit
  }
})

test_that("the engine's own warnings reach the guidance they were routed for", {
  ws <- warning_anchor_fit()$warnings_structured
  expect_true(is.data.frame(ws) && nrow(ws) > 0L)
  guide <- function(i) ferx:::.ferx_warning_guidance(ws$category[i], message = ws$message[i])

  # 1. The unused theta. ferx-core has no `unused_parameter` code, so this
  #    arrives under `general` and is recovered from the message. If core gains
  #    a code for it, or rewords it, this fails instead of silently printing
  #    nothing.
  i_unused <- which(grepl("TVUNUSED", ws$message, fixed = TRUE) &
                    !grepl("^Covariance step", ws$message))
  expect_length(i_unused, 1L)
  expect_identical(ws$category[i_unused], "general")
  expect_match(guide(i_unused), "never referenced", fixed = TRUE)

  # 2. The covariance failure: a CRITICAL row that printed no guidance at all
  #    before the routing fix, because the block was gated on `covariance_step`
  #    while `classify_warning` codes this `covariance_failed`.
  i_cov <- which(ws$category == "covariance_failed")
  expect_length(i_cov, 1L)
  expect_match(ws$message[i_cov], "ill-conditioned entries", fixed = TRUE)
  # The engine labels TVUNUSED "(zero diagonal -- flat objective)", so the
  # guidance must answer THAT cause -- not the finite-difference one, which
  # would be wrong here twice over: the label says the stencil succeeded, and
  # this fit ran on the analytic R-matrix path where nothing is differenced.
  g_cov <- guide(i_cov)
  expect_match(ws$message[i_cov], "zero diagonal", fixed = TRUE)
  expect_match(g_cov, "Hessian diagonal", fixed = TRUE)
  # ... and specifically not the generic fallback, which is where a message the
  # targeted branches failed to recognise would land.
  expect_false(grepl("Check identifiability", g_cov, fixed = TRUE))

  # 3. No row in a real table may produce an empty guidance block.
  for (i in seq_len(nrow(ws))) {
    g <- guide(i)
    expect_true(is.null(g) || (is.character(g) && nzchar(g)), info = ws$message[i])
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

# ---------------------------------------------------------------------------
# ode_solver: which clauses the statistics warning carries
# ---------------------------------------------------------------------------
# ferx-core's post-fit `ode_solver` warning, built the way
# ode_solver_diagnostics_warning() in src/api/postfit.rs builds it: one clause
# per non-zero counter, joined by "; ", then a lead-in and trailing advice
# picked by which counters are non-zero. Only the four counters these tests need
# are modelled. Clause texts are verbatim as of ferx-core 944cbf1e; the drift
# guard at the end of this block checks the phrases the guidance keys on against
# the engine source at whatever revision is pinned.
.ode_solver_warning <- function(abandoned = 0L, clamped = 0L, jets = 0L,
                                aborted = 0L) {
  parts <- character(0)
  if (abandoned > 0L) parts <- c(parts, paste0(
    abandoned, " solver walk(s) were abandoned before integrating because the ",
    "subject's timeline could not be ordered — a NaN or infinite dose time, ",
    "lagtime, route lag, or infusion duration at the final estimates — so ",
    "those subjects' predictions are NaN by construction, and they contributed ",
    "nothing to any other counter in this payload because nothing was integrated ",
    "for them. This counts walks, not subjects: one subject reaches more than one ",
    "engine in this pass (its predictions and its [odes] state readout are ",
    "separate walks), so it contributes more than one. Check the dose records and ",
    "any exponential covariate model on ALAG / F / D / R for a value that ",
    "overflows at typical covariates"
  ))
  if (clamped > 0L) parts <- c(parts, paste0(
    clamped, " step(s) clamped at the minimum step size — the local-error ",
    "test failed and the step was accepted anyway because dt could not shrink ",
    "further, so those segments are stability-limited rather than ",
    "accuracy-limited, and any output times left in a segment the solver could ",
    "not finish are freeze-padded with the last state (finite, but not ",
    "integrated)"
  ))
  if (jets > 0L) parts <- c(parts, paste0(
    jets, " segment(s) of the analytic-sensitivity solve were discarded and ",
    "re-solved with rk45: their predicted values were all finite while their ",
    "analytic derivatives — the gradients FOCE/FOCEI differentiate — ",
    "had overflowed to inf/NaN. The stiff method integrated those segments; the ",
    "trajectory simply reached a magnitude the sensitivities cannot represent, so ",
    "naming a different ode_method will not help. Check the model's units and ",
    "scaling (a state in ng rather than mg, an unbounded growth term, a rate ",
    "constant on the wrong clock) before trusting the estimates. This count comes ",
    "from the sensitivity sweep and is separate from the escalation counts ",
    "reported above"
  ))
  if (aborted > 0L) parts <- c(parts, paste0(
    aborted, " segment(s) were abandoned early by ode_stiff_abort_after = 2, ",
    "which bounds their cost and freeze-pads their tails"
  ))
  ran <- clamped > 0L || jets > 0L || aborted > 0L
  paste0(
    "W_ODE_SOLVER_DIAGNOSTICS: the ODE solver ",
    if (ran) "did not integrate cleanly" else "did not produce a usable integration",
    " at the final estimates (ode_method = auto): ",
    paste(parts, collapse = "; "),
    ". Counters are from the post-fit prediction pass over all subjects.",
    if (ran) paste0(
      " For the segments that did integrate, consider a different ode_method, a ",
      "looser ode_reltol / ode_abstol, or checking the parameter estimates that ",
      "produce these dynamics."
    ),
    if (abandoned > 0L) paste0(
      " The abandoned walk(s) are not an ode_method or tolerance problem — ",
      "nothing was integrated for them, so no solver setting changes the outcome; ",
      "fix the record or the parameter that produces the non-finite time."
    )
  )
}

# Every piece of solver-setting advice the guidance gives names ode_abstol, and
# nothing else it says does.
.recommends_solver_settings <- function(g) grepl("ode_abstol", g, fixed = TRUE)

test_that("an abandoned ODE timeline is not answered with solver settings", {
  # ferx-core #1234. The walk never integrated - the subject's timeline could
  # not be ordered - so the engine's own clause says no solver setting changes
  # the outcome. The guidance used to tell the user to change them anyway,
  # directly beneath that clause.
  g <- ferx:::.ferx_warning_guidance("ode_solver",
                                     message = .ode_solver_warning(abandoned = 3L))
  expect_match(g, "Nothing was integrated for some subjects", fixed = TRUE)
  expect_match(g, "dose records", fixed = TRUE)
  expect_false(.recommends_solver_settings(g))
})

test_that("an analytic-sensitivity overflow is not answered with solver settings", {
  # ferx-core #1204. The stiff method integrated the segment; the derivatives
  # overflowed. The clause says a different ode_method will not help. The
  # engine's lead-in still counts this counter as an integration that ran
  # badly, so its message ends with its general solver-setting sentence
  # anyway. The guidance follows the clause, not that sentence.
  g <- ferx:::.ferx_warning_guidance("ode_solver",
                                     message = .ode_solver_warning(jets = 2L))
  expect_match(g, "will not help", fixed = TRUE)
  expect_match(g, "units and scaling", fixed = TRUE)
  expect_false(.recommends_solver_settings(g))
})

test_that("solver-setting advice stays with the ode_solver clauses it fits", {
  g <- function(...) {
    ferx:::.ferx_warning_guidance("ode_solver", message = .ode_solver_warning(...))
  }
  # A plain integration problem keeps the advice it always had.
  clamps <- g(clamped = 2L)
  expect_match(clamps, "Integration under the final estimates was not clean",
               fixed = TRUE)
  expect_true(.recommends_solver_settings(clamps))

  # ode_stiff_abort_after's clause says "abandoned" too, but that segment did
  # integrate, and the budget is a solver setting.
  aborts <- g(aborted = 3L)
  expect_true(.recommends_solver_settings(aborts))
  expect_false(grepl("Nothing was integrated", aborts, fixed = TRUE))

  # Several clauses in one warning: each gets its own advice, and the
  # solver-setting advice is scoped to the segments it can help.
  walks_and_clamps <- g(abandoned = 3L, clamped = 2L)
  expect_match(walks_and_clamps, "Nothing was integrated", fixed = TRUE)
  expect_match(walks_and_clamps, "For the clamped, discarded", fixed = TRUE)

  jets_and_clamps <- g(jets = 2L, clamped = 2L)
  expect_match(jets_and_clamps, "units and scaling", fixed = TRUE)
  expect_match(jets_and_clamps, "For the clamped, discarded", fixed = TRUE)
})

test_that("the ode_solver escalation note keeps its informational guidance", {
  # Matched on its token before any clause phrase: this note contains "clamped
  # at the minimum step size" as well ("no step clamped ...").
  note <- paste0(
    "W_ODE_SOLVER_ESCALATION_NOTE: ode_method = auto escalated 240 integration ",
    "segment(s) to a stiff stepper at the final estimates; every other segment ",
    "used rk45, no escalation was rejected, and no step clamped at the minimum ",
    "step size. Informational — set ode_method = rk45 to pin the explicit ",
    "stepper, or name a stiff method to pin the other half."
  )
  g <- ferx:::.ferx_warning_guidance("ode_solver", message = note)
  expect_match(g, "(informational)", fixed = TRUE)
  expect_false(.recommends_solver_settings(g))
})

test_that("a real fit over an unorderable ODE timeline prints the timeline guidance", {
  # End to end, like the init_outside_bounds test above: the engine emits the
  # abandoned-walk clause (ferx-core #1234), the glue carries it into
  # fit$warnings_structured under `ode_solver`, and ferx_get_warnings() prints
  # the guidance under it. Unlike the transcribed fixtures, this message comes
  # from the pinned engine itself, so it runs in CI as well.
  #
  # The fixture is ferx-core's own (a_fit_over_an_unorderable_timeline_says_so):
  # a two-state ODE model and one subject whose dose TIME is NaN. Nothing
  # upstream rejects that - the engine's data checks cover dose attributes
  # (ALAG / F / D / R), not record times. The CSV is written with na = "NaN"
  # because write.csv() writes NaN as NA by default, which the reader treats as
  # a missing value, so the timeline never becomes NaN.
  model <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(1.0, 0.1, 50.0)",
    "  theta TVV(10.0, 1.0, 500.0)",
    "  theta KFAST(1.0, 1e-6, 1e6)",
    "  omega ETA_CL ~ 0.04",
    "  sigma PROP ~ 0.04",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "  KF = KFAST",
    "[structural_model]",
    "  ode(obs_cmt=central, states=[central, periph])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central - KF * central + KF * periph",
    "  d/dt(periph)  = KF * central - KF * periph",
    "[error_model]",
    "  DV ~ proportional(PROP)"
  ), model)
  obs_t <- c(0.5, 2, 8, 24)
  subject <- function(id, dose_time, scale) {
    data.frame(ID = id, TIME = c(dose_time, obs_t),
               DV = c(0, scale * 50 / (1 + obs_t)),
               EVID = c(1, 0, 0, 0, 0), AMT = c(100, 0, 0, 0, 0),
               CMT = 1, MDV = c(1, 0, 0, 0, 0))
  }
  data <- tempfile(fileext = ".csv")
  utils::write.csv(rbind(subject(1, NaN, 1.0), subject(2, 0, 1.2)), data,
                   row.names = FALSE, na = "NaN")
  on.exit(unlink(c(model, data)))

  fit <- ferx_fit(model, data, method = "focei", settings = list(maxiter = 1L),
                  covariance = FALSE, verbose = FALSE)
  ws <- ferx_get_warnings(fit, as_df = TRUE)
  i <- which(ws$category == "ode_solver")
  expect_length(i, 1L)
  expect_match(ws$message[i], "could not be ordered", fixed = TRUE)

  g <- ferx:::.ferx_warning_guidance(ws$category[i], message = ws$message[i])
  expect_match(g, "Nothing was integrated for some subjects", fixed = TRUE)
  expect_false(.recommends_solver_settings(g))

  # ... and it actually prints. Whitespace-normalised: the printer wraps at 70.
  out <- capture.output(ferx_get_warnings(fit))
  flat <- gsub("[[:space:]]+", " ", paste(out, collapse = " "))
  expect_true(grepl("Nothing was integrated for some subjects", flat, fixed = TRUE))
})

test_that("the ode_solver phrases the guidance keys on are still in the engine", {
  # The fixtures above are transcribed, so on their own they cannot notice
  # ferx-core rewording a clause. This reads the string literals of the function
  # that builds the warning, at the pinned revision, and requires every phrase
  # `.ferx_ode_solver_anchors()` matches on to still be in one of them.
  #
  # Skipped without a sibling ferx-core checkout, like the other source guards
  # in this file: CI builds the pinned crate, and a pin bump happens on a
  # developer machine, where this runs.
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")
  postfit <- file.path(core_src, "api", "postfit.rs")
  skip_if(!file.exists(postfit), "ferx-core at the pin has no src/api/postfit.rs")

  ln <- readLines(postfit, warn = FALSE)
  start <- grep("fn ode_solver_diagnostics_warning(", ln, fixed = TRUE)
  expect_length(start, 1L)
  if (length(start) == 1L) {
    # A top-level fn closes on a bare `}` in column one.
    end <- start + which(ln[(start + 1L):length(ln)] == "}")[1]
    body <- ln[start:end]
    # Drop `//` comments (they quote the message text too), then rejoin Rust's
    # backslash-continued literals, as the covariance inventory guard does.
    body <- body[!grepl("^[[:space:]]*//", body)]
    joined <- gsub("\\\\\n[[:space:]]*", "", paste(body, collapse = "\n"))
    lits <- regmatches(joined,
                       gregexpr('"(?:[^"\\\\]|\\\\.)*"', joined, perl = TRUE))[[1]]
    expect_gt(length(lits), 5L)            # the function was found and understood

    for (phrase in unlist(ferx:::.ferx_ode_solver_anchors(), use.names = FALSE)) {
      expect_true(any(grepl(phrase, lits, fixed = TRUE)),
                  info = paste0("no longer in ode_solver_diagnostics_warning(): ",
                                phrase))
    }
  }
  # The escalation note's token is a const outside the function.
  expect_true(any(grepl('"W_ODE_SOLVER_ESCALATION_NOTE"', ln, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# #308: the codes that used to print no guidance, and the two split codes
# ---------------------------------------------------------------------------
# Each `msg` is ferx-core's own text at the pinned revision with its `{}`
# placeholders filled in. Every case is routed through the engine's real
# classifier first, so a pairing production never emits cannot pass here: the
# case asserts the category ferx-core gives the message, then the guidance that
# category earns, then that the `general` door (a loaded fit) gives the same.
.engine_routed_cases <- function() {
  list(
    list(cat = "boundary_estimate", sev = "warning",
         msg = paste0("Parameter estimate(s) pinned to an optimizer bound: TVKA ",
                      "(50.0000 at upper bound). This often indicates ",
                      "non-identifiability or a too-tight bound; inspect the ",
                      "affected parameter(s) and consider relaxing the bound or ",
                      "simplifying the model."),
         ok = "declared bounds"),
    list(cat = "inflated_rse", sev = "warning",
         msg = paste0("High relative standard error (RSE > 50%): TVKA (84.2%). ",
                      "These parameter(s) are imprecisely estimated — often a ",
                      "sign of over-parameterization or data that do not inform ",
                      "them; consider simplifying the model."),
         ok = "imprecisely estimated"),
    list(cat = "high_correlation", sev = "warning",
         msg = paste0("Highly correlated parameter pair(s) (|r| >= 0.95): TVCL ~ ",
                      "TVV (0.97). Highly correlated estimates indicate ",
                      "over-parameterization or non-identifiability; consider ",
                      "fixing or removing one of each pair."),
         ok = "fit$cor_matrix"),
    list(cat = "eta_shrinkage", sev = "warning",
         msg = paste0("High ETA shrinkage (≥ 30%): ETA_KA (54%). EBE-based ",
                      "diagnostics for these random effects are unreliable and ",
                      "the data poorly inform their individual estimates; ",
                      "consider removing the IIV on the affected parameter(s) or ",
                      "collecting more informative data."),
         ok = "do not screen covariates"),
    list(cat = "eps_shrinkage", sev = "warning",
         msg = paste0("EPS shrinkage is notably negative (-35.0%): mean(IWRES^2) ",
                      "> 1, which means the residual error model does not absorb ",
                      "the residuals at the final EBE etas. Common causes: SAEM ",
                      "converged to a local optimum with under-fit sigma (try ",
                      "`method = [saem, focei]` to polish with FOCEI, or ",
                      "different starts); model misspecification on a subset of ",
                      "subjects; sigma at a bound. Inspect the IWRES distribution ",
                      "in the sdtab."),
         ok = "IWRES"),
    list(cat = "flat_parameter", sev = "warning",
         msg = paste0("[parameters] `TVX` has no effect on the objective ",
                      "(gradient ≈ 0 at the initial estimate) — it is ",
                      "likely computed but never used (unmapped, or dropped from ",
                      "the structural / scaling model). Freezing it at its ",
                      "initial value (1) so the remaining parameters can be ",
                      "estimated; map or remove `TVX` to silence this."),
         ok = "not estimated"),
    list(cat = "experimental", sev = "warning",
         msg = paste0("Stochastic differential equations ([diffusion] / Extended ",
                      "Kalman Filter) are an EXPERIMENTAL feature: validated only ",
                      "on a small set of toy examples, with estimator support ",
                      "limited to FOCE/FOCEI. Standard errors and convergence ",
                      "behaviour are not yet proven across diverse datasets ",
                      "— validate results carefully before relying on them. ",
                      "See the Feature Maturity page in the documentation."),
         ok = "experimental"),
    list(cat = "absorption_twin_declined", sev = "warning",
         msg = paste0("This absorption model's ODE equivalent could not be built, ",
                      "so the model stays closed-form. Subjects needing the ODE ",
                      "fallback (time-varying covariates, a `TIME`-dependent ",
                      "parameter, IOV, steady-state or infusion doses, or the ",
                      "flip-flop regime) will be rejected with an explicit error ",
                      "instead of silently rerouting (W_ABSORPTION_TWIN_DECLINED, ",
                      "#1008). Reason: unknown symbol `Q`."),
         ok = "no ODE fallback"),
    # flip_flop, both emitters.
    list(cat = "flip_flop", sev = "warning",
         msg = paste0("one_cpt_transit disposition rate (2.5000) ≥ transit ",
                      "rate KTR = (n+1)/mtt (2.0000) at typical values (subject ",
                      "1): the flip-flop regime, outside the analytic absorption ",
                      "closed form's convergence domain. ferx automatically ",
                      "evaluates the equivalent ODE transit() model for such ",
                      "parameters (correct, but slower than the closed form) ",
                      "— check the MTT / CL starting estimates if the ",
                      "flip-flop is unexpected."),
         ok = "Informational", never = "degenerate"),
    list(cat = "flip_flop", sev = "warning",
         msg = paste0("2 subject(s) [3, 7] enter the analytic absorption closed ",
                      "form's clamp region — the flip-flop regime ",
                      "(disposition rate ≥ transit rate KTR = (n+1)/mtt), or ",
                      "(for a 2-cpt model) coincident disposition eigenvalues ",
                      "— at their fitted empirical-Bayes estimates, where the ",
                      "closed form returns an identically-zero concentration ",
                      "profile, silently degenerating those subjects' likelihood ",
                      "contributions."),
         ok = "degenerate", never = "Informational"),
    list(cat = "simulation", sev = "warning",
         msg = paste0("W_TTE_DEGENERATE_HAZARD: subject '4' (CMT=2) drew a ",
                      "non-positive / non-finite effective hazard rate; no event ",
                      "was generated and the subject is censored at the ",
                      "observation window (#763). Check the hazard parameters / ",
                      "covariate values."),
         ok = "degenerate hazard"),
    list(cat = "simulation", sev = "warning",
         msg = paste0("W_RTTE_DEGENERATE: subject '4' (CMT=2) has a degenerate ",
                      "recurrent hazard (~1.000e7 events expected before horizon ",
                      "100); its event stream was skipped and the subject ",
                      "censored at the horizon (#762). Check the hazard ",
                      "parameters / covariate values."),
         ok = "degenerate hazard"),
    # mu_referencing: one code, two severities, two answers.
    list(cat = "mu_referencing", sev = "warning",
         msg = paste0("individual parameter(s) not mu-referenced: KA. This can ",
                      "strongly affect convergence; prefer forms such as ",
                      "`CL = TVCL * exp(ETA_CL)` when possible."),
         ok = "not mu-referenced", never = "informational"),
    list(cat = "mu_referencing", sev = "info",
         msg = paste0("covariate mu-reference on CL (reads WT) is not used: its ",
                      "random effect has no variance."),
         ok = "informational", never = "not mu-referenced"),
    # optimizer_config: likewise.
    list(cat = "optimizer_config", sev = "warning",
         msg = "global_search disabled: CRS2-LM initialisation failed",
         ok = "could not start", never = "informational"),
    list(cat = "optimizer_config", sev = "info",
         msg = paste0("global_search = true with n_starts = 4: CRS2-LM only runs ",
                      "on start 0 (it ignores the starting point and would ",
                      "override the theta perturbation on starts 1..4)"),
         ok = "informational", never = "could not start")
  )
}

test_that("the engine routes each #308 message to the category its guidance keys on", {
  cases <- .engine_routed_cases()
  cl <- ferx:::.ferx_classify_flat_warnings(vapply(cases, `[[`, "", "msg"))
  for (i in seq_along(cases)) {
    case <- cases[[i]]
    expect_identical(cl$category[i], case$cat, info = case$msg)
    expect_identical(cl$severity[i], case$sev, info = case$msg)
  }
})

test_that("each #308 message gets its own guidance, through either door", {
  for (case in .engine_routed_cases()) {
    g <- ferx:::.ferx_warning_guidance(case$cat, message = case$msg)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = case$msg)
    expect_match(g, case$ok, fixed = TRUE, info = case$msg)
    if (!is.null(case$never)) {
      expect_false(grepl(case$never, g, fixed = TRUE), info = case$msg)
    }
    # A loaded fit hands every row in as `general`; the answer must not change.
    expect_identical(ferx:::.ferx_warning_guidance("general", message = case$msg),
                     g, info = case$msg)
  }
})

test_that("the #308 message texts are still what the engine emits", {
  # The cases above are hand-copied; this reads a fixed fragment of each back
  # out of ferx-core at the pinned revision, so a reworded emitter fails here
  # rather than leaving the cases testing a message nobody sends.
  core_src <- .core_src_at_pin()
  skip_if(is.null(core_src), "cannot materialise ferx-core at the pinned revision")
  rs <- list.files(core_src, pattern = "\\.rs$", recursive = TRUE, full.names = TRUE)
  rs <- rs[!grepl("tests?(/|\\.rs$)|_tests\\.rs$", rs)]
  # Join backslash-continued string literals so a fragment spanning a line
  # break in the source still matches.
  src <- paste(vapply(rs, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
                      ""), collapse = "\n")
  src <- gsub("\\\\\n[[:space:]]*", "", src)
  fragments <- c(
    "pinned to an optimizer bound: {list}",
    "High relative standard error (RSE > {:.0}%)",
    "Highly correlated parameter pair(s)",
    "EBE-based diagnostics for these random effects are unreliable",
    "EPS shrinkage is notably negative",
    "has no effect on the objective",
    "are an EXPERIMENTAL feature",
    "(W_ABSORPTION_TWIN_DECLINED, #1008)",
    "ferx automatically evaluates the equivalent ODE",
    "identically-zero concentration profile, silently degenerating",
    "W_TTE_DEGENERATE_HAZARD: subject",
    "W_RTTE_DEGENERATE: subject",
    "individual parameter(s) not mu-referenced",
    "global_search disabled: {}",
    "CRS2-LM only runs on start 0",
    "fd_hessian_step must be positive and finite",
    "the score cross-product matrix S is singular"
  )
  for (f in fragments) {
    expect_true(grepl(f, src, fixed = TRUE), info = f)
  }
})

test_that("a loaded fit's warnings get the same categories and guidance as the fresh fit", {
  # ferx_load_fit() does not restore warnings_structured, so before #308 every
  # row of a reloaded fit arrived as `general` and all but three guidance arms
  # were dead for it. This drives the real round trip.
  fit <- warning_anchor_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_null(loaded$warnings_structured)

  fresh_ws  <- ferx_get_warnings(fit, as_df = TRUE)
  loaded_ws <- ferx_get_warnings(loaded, as_df = TRUE)
  expect_gt(nrow(loaded_ws), 0L)
  # Something other than `general` came back - the recovery actually ran.
  expect_true(any(loaded_ws$category != "general"))
  for (i in seq_len(nrow(loaded_ws))) {
    j <- match(loaded_ws$message[i], fresh_ws$message)
    if (is.na(j)) next   # an R-side row the flat vector does not carry
    expect_identical(loaded_ws$category[i], fresh_ws$category[j],
                     info = loaded_ws$message[i])
    expect_identical(loaded_ws$severity[i], fresh_ws$severity[j],
                     info = loaded_ws$message[i])
    expect_identical(
      ferx:::.ferx_warning_guidance(loaded_ws$category[i], loaded_ws$message[i]),
      ferx:::.ferx_warning_guidance(fresh_ws$category[j], fresh_ws$message[j]),
      info = loaded_ws$message[i]
    )
  }
})

test_that("a dropped-[output] row keeps its guidance on a loaded fit", {
  # The one R-side message that rides in the flat vector; the engine files it
  # under `general`, so it is recovered by its text.
  msg <- "`FOO` named in [output] was not written to sdtab."
  expect_identical(ferx:::.ferx_classify_flat_warnings(msg)$category, "output")
  expect_identical(ferx:::.ferx_warning_guidance("general", message = msg),
                   ferx:::.ferx_warning_guidance("output"))
})

# ---------------------------------------------------------------------------
# #308 review: a partial structured table must not hide the flat warnings
# ---------------------------------------------------------------------------
test_that("flat warnings missing from a partial structured table are merged in", {
  # The shape a loaded fit has after ferx_sir() / ferx_covariance(): two older
  # flat warnings, and a structured table holding only the post-hoc row.
  fake <- structure(
    list(
      model_name = "loaded",
      warnings = c("Outer optimization did not converge", "something else",
                   "SIR failed: no usable proposal"),
      warnings_structured = data.frame(
        severity = "warning", category = "sir",
        message = "SIR failed: no usable proposal", source_method = "",
        stringsAsFactors = FALSE
      )
    ),
    class = "ferx_fit"
  )
  df <- ferx_get_warnings(fake, as_df = TRUE)
  expect_equal(nrow(df), 3L)
  expect_setequal(df$message, fake$warnings)
  # The table's own row keeps the category it was given ...
  expect_identical(df$category[df$message == "SIR failed: no usable proposal"], "sir")
  # ... and the recovered ones get the engine's.
  expect_identical(df$category[df$message == "Outer optimization did not converge"],
                   "convergence")
  # The fit is read, never written: strictness still sees the stored table.
  expect_equal(nrow(fake$warnings_structured), 1L)
})

test_that("a fresh fit's warnings come back exactly as its structured table", {
  # On a fresh fit every flat message is already in the table, so the merge
  # must add nothing - no duplicate rows under a second category.
  fit <- warning_anchor_fit()
  expect_identical(ferx_get_warnings(fit, as_df = TRUE),
                   fit$warnings_structured[, c("severity", "category", "message",
                                               "source_method")])
})

test_that("a loaded fit keeps its older warnings after ferx_covariance()", {
  fit <- warning_anchor_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  before <- ferx_get_warnings(loaded, as_df = TRUE)
  expect_gt(nrow(before), 0L)

  recov <- suppressWarnings(ferx_covariance(loaded))
  expect_true(is.data.frame(recov$warnings_structured))
  after <- ferx_get_warnings(recov, as_df = TRUE)
  # Every warning the loaded fit showed is still shown, under its category.
  for (i in seq_len(nrow(before))) {
    j <- match(before$message[i], after$message)
    expect_false(is.na(j), info = before$message[i])
    if (!is.na(j) && !identical(after$category[j], "covariance")) {
      expect_identical(after$category[j], before$category[i], info = before$message[i])
    }
  }
})

test_that("a loaded fit keeps its older warnings after ferx_sir()", {
  ex <- ferx_example("warfarin")
  fit <- ferx_fit(ex$model, ex$data, method = "focei", covariance = TRUE,
                  verbose = FALSE, settings = list(maxiter = 30L))
  skip_if(is.null(fit$cov_matrix), "covariance step did not produce a matrix")
  # Give the loaded fit an older warning of a known category to track.
  old <- "Outer optimization did not converge"
  fit$warnings <- unique(c(fit$warnings, old))
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_null(loaded$warnings_structured)

  sired <- suppressWarnings(ferx_sir(loaded, sir_samples = 200L,
                                     sir_resamples = 50L, sir_seed = 1L))
  ws <- ferx_get_warnings(sired, as_df = TRUE)
  expect_true(old %in% ws$message)
  expect_identical(ws$category[ws$message == old], "convergence")
  expect_identical(ws$severity[ws$message == old], "critical")
})

test_that("print() of a reloaded fit tallies the recovered severities", {
  fit <- warning_anchor_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_null(loaded$warnings_structured)
  out <- paste(capture.output(print(loaded)), collapse = "\n")
  # The anchor fit's covariance failure is critical; before the fix a loaded
  # fit printed only "N warning(s) -- inspect fit$warnings".
  expect_match(out, "critical", fixed = TRUE)
  expect_match(out, "covariance_failed", fixed = TRUE)
  expect_false(grepl("inspect fit$warnings", out, fixed = TRUE))
})
