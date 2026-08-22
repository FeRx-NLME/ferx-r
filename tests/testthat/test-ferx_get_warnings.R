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
test_that(".ferx_warning_guidance dispatches covariance_step by message content", {
  g <- function(msg) ferx:::.ferx_warning_guidance("covariance_step", message = msg)

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
  expect_match(g(base_reg("minor")),    "benign|reliable",   ignore.case = TRUE, perl = TRUE)
  expect_match(g(base_reg("moderate")), "ferx_sir|caution",  ignore.case = TRUE, perl = TRUE)
  expect_match(g(base_reg("severe")),   "ferx_sir|unreliable", ignore.case = TRUE, perl = TRUE)
  # Minor guidance should not suggest SIR.
  expect_false(grepl("ferx_sir", g(base_reg("minor")), ignore.case = TRUE))

  # Generic fallback for unrecognised message.
  expect_match(g("Covariance step failed"), "identifiability", ignore.case = TRUE)
})
test_that("ferx_get_warnings() shows guidance for unused_parameter category", {
  fake <- structure(
    list(
      model_name = "m",
      warnings_structured = data.frame(
        severity      = "warning",
        category      = "unused_parameter",
        message       = "TVCL is declared but never referenced",
        source_method = "",
        stringsAsFactors = FALSE
      ),
      warnings = "TVCL is declared but never referenced",
      condition_number = NULL,
      eta_normality = NULL,
      uses_sde = FALSE
    ),
    class = "ferx_fit"
  )
  out <- capture.output(ferx_get_warnings(fake))
  # Guidance for unused_parameter must appear in the output
  expect_true(
    any(grepl("Remove it from", out, fixed = TRUE)),
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
test_that(".ferx_warning_guidance maps every known category to non-empty text", {
  cats <- c("convergence", "covariance_step", "condition_number", "optimizer_health",
            "eta_normality", "bloq_method", "sir", "importance_sampling", "data_quality",
            "omega_structure", "ebe_convergence", "gradient_fallback", "mu_referencing",
            "optimizer_config", "multi_start", "threads", "cancelled", "unused_parameter")
  for (cat in cats) {
    g <- ferx:::.ferx_warning_guidance(cat)
    expect_true(is.character(g) && length(g) == 1L && nzchar(g), info = cat)
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
