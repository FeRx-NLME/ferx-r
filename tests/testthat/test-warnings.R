test_that("fit carries a structured warnings data frame", {
  fit <- warfarin_fit_cov()
  expect_true(is.data.frame(fit$warnings_structured))
  expect_named(fit$warnings_structured, c("severity", "category", "message", "source_method"))
  # Severity values are always from the fixed vocabulary
  if (nrow(fit$warnings_structured) > 0L) {
    expect_true(all(fit$warnings_structured$severity %in%
      c("critical", "warning", "info")))
  }
})

test_that("ferx_warnings(as_df = TRUE) returns the underlying data frame", {
  fit <- warfarin_fit_cov()
  df <- ferx_warnings(fit, as_df = TRUE)
  expect_identical(df, fit$warnings_structured)
})

test_that("ferx_warnings() prints a grouped summary", {
  fit <- warfarin_fit_cov()
  out <- capture.output(ferx_warnings(fit))
  expect_true(any(grepl("ferx fit warnings", out)))
  # Footer tallies are always present
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})

test_that("ferx_warnings() rejects non-fit input", {
  expect_error(ferx_warnings(list()), "ferx_fit")
})

test_that("ferx_warnings() falls back to flat warnings when structured is absent", {
  fake <- structure(
    list(
      model_name = "legacy",
      warnings = c("Outer optimization did not converge", "something else"),
      warnings_structured = NULL
    ),
    class = "ferx_fit"
  )
  df <- ferx_warnings(fake, as_df = TRUE)
  expect_equal(nrow(df), 2L)
  expect_true(all(df$severity == "warning"))
  expect_true(all(df$category == "general"))
})

test_that("structured warnings include R-side condition number flag", {
  fake <- structure(
    list(
      model_name = "m",
      condition_number = 5000,
      eta_normality = NULL
    ),
    class = "ferx_fit"
  )
  raw <- list(
    warnings_severity = character(0),
    warnings_category = character(0),
    warnings_message  = character(0)
  )
  df <- ferx:::.ferx_assemble_structured_warnings(raw, fake)
  expect_true(any(df$category == "condition_number"))
  expect_equal(df$severity[df$category == "condition_number"], "critical")
})

test_that("multiple non-normal ETAs fold into one structured warning (#163)", {
  eta_norm <- data.frame(
    eta   = c("ETA_CL", "ETA_V", "ETA_KA"),
    W     = c(0.90, 0.80, 0.99),
    p_val = c(0.0001, 0.0000, 0.6000),
    flag  = c("[!]", "[!]", ""),
    stringsAsFactors = FALSE
  )
  fake <- structure(
    list(model_name = "m", condition_number = NULL, eta_normality = eta_norm),
    class = "ferx_fit"
  )
  raw <- list(
    warnings_severity = character(0),
    warnings_category = character(0),
    warnings_message  = character(0)
  )
  df <- ferx:::.ferx_assemble_structured_warnings(raw, fake)
  norm_rows <- df[df$category == "eta_normality", , drop = FALSE]
  # Exactly one eta_normality warning, regardless of how many ETAs were flagged
  expect_equal(nrow(norm_rows), 1L)
  # The single message lists every flagged ETA but not the normal one
  expect_match(norm_rows$message, "ETA_CL")
  expect_match(norm_rows$message, "ETA_V")
  expect_false(grepl("ETA_KA", norm_rows$message))
  # Plural noun for >1 flagged ETA
  expect_match(norm_rows$message, "2 ETAs:")
})

test_that(".ferx_eta_normality_warning uses singular noun for 1 ETA", {
  one <- data.frame(
    eta = "ETA_CL", W = 0.9, p_val = 0.0001, flag = "[!]",
    stringsAsFactors = FALSE
  )
  msg <- ferx:::.ferx_eta_normality_warning(one)
  expect_match(msg, "1 ETA:")
  expect_false(grepl("ETA(s)", msg, fixed = TRUE))
})

test_that(".ferx_eta_normality_warning returns NULL when nothing is flagged", {
  expect_null(ferx:::.ferx_eta_normality_warning(NULL))
  all_ok <- data.frame(
    eta = c("ETA_CL", "ETA_V"), W = c(0.99, 0.98),
    p_val = c(0.7, 0.6), flag = c("", ""), stringsAsFactors = FALSE
  )
  expect_null(ferx:::.ferx_eta_normality_warning(all_ok))
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

test_that("ferx_warnings() shows guidance for unused_parameter category", {
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
  out <- capture.output(ferx_warnings(fake))
  # Guidance for unused_parameter must appear in the output
  expect_true(
    any(grepl("Remove it from", out, fixed = TRUE)),
    info = paste("Expected guidance text not found in output:\n",
                 paste(out, collapse = "\n"))
  )
})
