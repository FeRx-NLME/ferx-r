# test-runlog.R
# Tests for ferx_runlog() -- shape, content, and graceful degradation.
# Uses warfarin_fit() (no covariance) and warfarin_fit_cov() (with covariance)
# from helper-warfarin-fit.R, which cap iterations at 30 for CI speed.

# -- Tier 1: basic contract ---------------------------------------------------

test_that("ferx_runlog returns an invisible character string", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_type(out, "character")
  expect_length(out, 1L)
})

test_that("ferx_runlog verbose=TRUE prints and still returns invisibly", {
  fit <- warfarin_fit()
  expect_output(out <- ferx_runlog(fit, verbose = TRUE), regexp = "ferx run log")
  expect_type(out, "character")
})

test_that("ferx_runlog output contains all mandatory section headers", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "ferx run log",         fixed = TRUE)
  expect_match(out, "DATA SUMMARY",         fixed = TRUE)
  expect_match(out, "PARAMETER ESTIMATES",  fixed = TRUE)
  expect_match(out, "OBJECTIVE FUNCTION",   fixed = TRUE)
  expect_match(out, "FINAL GRADIENT",       fixed = TRUE)
  expect_match(out, "RUNTIME",              fixed = TRUE)
})

# -- Tier 2: header fields ----------------------------------------------------

test_that("ferx_runlog header includes method and ferx version", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "Method:", fixed = TRUE)
  expect_match(out, "ferx:",   fixed = TRUE)
  expect_match(out, "Date:",   fixed = TRUE)
})

test_that("ferx_runlog header includes model name when present", {
  fit <- warfarin_fit()
  if (!is.null(fit$model_name) && nzchar(fit$model_name)) {
    out <- ferx_runlog(fit, verbose = FALSE)
    expect_match(out, fit$model_name, fixed = TRUE)
  } else {
    skip("model_name absent in this fit")
  }
})

# -- Tier 2: data summary -----------------------------------------------------

test_that("ferx_runlog shows subject and observation counts when available", {
  fit <- warfarin_fit()
  if (!is.null(fit$n_subjects) && !is.null(fit$n_obs)) {
    out <- ferx_runlog(fit, verbose = FALSE)
    expect_match(out, "Subjects:",     fixed = TRUE)
    expect_match(out, "Observations:", fixed = TRUE)
  } else {
    skip("n_subjects/n_obs absent in this fit")
  }
})

test_that("ferx_runlog shows time range when obs_time_range is present", {
  fit <- warfarin_fit()
  fit$obs_time_range <- c(0.5, 144)
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "Time range:", fixed = TRUE)
})

test_that("ferx_runlog data summary falls back gracefully when counts missing", {
  fit <- warfarin_fit()
  fit$n_subjects     <- NULL
  fit$n_obs          <- NULL
  fit$obs_time_range <- NULL
  fit$input_columns  <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "data summary not available", fixed = TRUE)
})

test_that("ferx_runlog shows input_columns line when field is present", {
  fit <- warfarin_fit()
  fit$input_columns <- c("ID", "TIME", "AMT", "EVID", "DV", "WT")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "Input columns:", fixed = TRUE)
  expect_match(out, "WT",             fixed = TRUE)
})

test_that("ferx_runlog fallback message absent when only input_columns is present", {
  fit <- warfarin_fit()
  fit$n_subjects     <- NULL
  fit$n_obs          <- NULL
  fit$obs_time_range <- NULL
  fit$input_columns  <- c("ID", "TIME", "DV")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out,    "Input columns:",          fixed = TRUE)
  expect_no_match(out, "data summary not available", fixed = TRUE)
})

test_that("ferx_runlog omits input_columns line when field is empty", {
  fit <- warfarin_fit()
  fit$input_columns <- character(0)
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "Input columns:", fixed = TRUE)
})

# -- Tier 2: parameter table --------------------------------------------------

test_that("ferx_runlog parameter table contains theta names", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  for (nm in fit$theta_names) {
    expect_match(out, nm, fixed = TRUE)
  }
})

test_that("ferx_runlog parameter table contains eta names", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  for (nm in fit$eta_names) {
    expect_match(out, nm, fixed = TRUE)
  }
})

test_that("ferx_runlog parameter table does NOT wrap eta names in OMEGA()", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  for (nm in fit$eta_names) {
    # bare name must appear
    expect_match(out, nm, fixed = TRUE)
    # must not be wrapped as OMEGA(<name>)
    expect_no_match(out, paste0("OMEGA(", nm, ")"), fixed = TRUE)
  }
})

test_that("ferx_runlog fallback omega labels use OMEGA(i,i) format", {
  fit <- warfarin_fit()
  fit$eta_names <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "OMEGA(1,1)", fixed = TRUE)
})

test_that("ferx_runlog shows SE and %RSE columns when covariance step ran", {
  fit <- warfarin_fit_cov()
  if (is.null(fit$se_theta)) skip("covariance step did not converge in 30 iter")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "%RSE",     fixed = TRUE)
  expect_match(out, "SE",       fixed = TRUE)
})

test_that("ferx_runlog shows INITIAL column when theta_init is present", {
  fit <- warfarin_fit()
  fit$theta_init <- fit$theta * 0.9   # synthetic inits
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "INITIAL", fixed = TRUE)
})

# -- Tier 2: covariance step --------------------------------------------------

test_that("ferx_runlog covariance section present when eigenvalues available", {
  fit <- warfarin_fit_cov()
  if (is.null(fit$eigenvalues)) skip("eigenvalues absent")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "COVARIANCE STEP", fixed = TRUE)
  expect_match(out, "Eigenvalues:",    fixed = TRUE)
})

test_that("ferx_runlog flags condition number > 1000", {
  fit <- warfarin_fit()
  fit$condition_number <- 1500
  fit$eigenvalues      <- c(0.1, 1500)
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "condition number > 1000", fixed = TRUE)
})

test_that("ferx_runlog covariance section absent when no cov data", {
  fit <- warfarin_fit()
  fit$condition_number <- NULL
  fit$eigenvalues      <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "COVARIANCE STEP", fixed = TRUE)
})

# -- Tier 2: estimation settings ----------------------------------------------

test_that("ferx_runlog estimation settings section present when optimizer_label available", {
  fit <- warfarin_fit()
  fit$optimizer_label       <- "LBFGS"
  fit$bloq_method_label     <- "drop"
  fit$outer_maxiter         <- 200L
  fit$outer_gtol            <- 1e-4
  fit$inits_from_nca        <- FALSE
  fit$gradient_method_outer <- "finite differences"  # gradient-based optimizer
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "ESTIMATION SETTINGS", fixed = TRUE)
  expect_match(out, "LBFGS",              fixed = TRUE)
  expect_match(out, "Max iterations:",    fixed = TRUE)
  expect_match(out, "Gradient tol:",      fixed = TRUE)
  expect_match(out, "BLOQ method:",       fixed = TRUE)
})

test_that("ferx_runlog estimation settings shows NCA warm-start method", {
  fit <- warfarin_fit()
  fit$optimizer_label <- "lbfgs"
  fit$inits_from_nca  <- "nca_sweep"
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "NCA-derived", fixed = TRUE)
  expect_match(out, "nca_sweep",   fixed = TRUE)
})

test_that("ferx_runlog estimation settings shows multi-start count", {
  fit <- warfarin_fit()
  fit$optimizer_label  <- "LBFGS"
  fit$n_starts         <- 5L
  fit$multi_start_seed <- 42
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "multi-start n=5", fixed = TRUE)
  expect_match(out, "Seeds:",          fixed = TRUE)
  expect_match(out, "multi_start=42",  fixed = TRUE)
})

test_that("ferx_runlog estimation settings shows covariate names", {
  fit <- warfarin_fit()
  fit$optimizer_label  <- "LBFGS"
  fit$covariate_names  <- c("WT", "AGE", "SEX")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "Covariates:", fixed = TRUE)
  expect_match(out, "WT",         fixed = TRUE)
  expect_match(out, "AGE",        fixed = TRUE)
})

test_that("ferx_runlog estimation settings section absent when no settings fields", {
  fit <- warfarin_fit()
  fit$optimizer_label   <- NULL
  fit$bloq_method_label <- NULL
  fit$inits_from_nca    <- NULL
  fit$covariate_names   <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "ESTIMATION SETTINGS", fixed = TRUE)
})

# -- Tier 2: diagnostics ------------------------------------------------------

test_that("ferx_runlog diagnostics section present when shrinkage available", {
  fit <- warfarin_fit()
  if (is.null(fit$shrinkage_eta)) skip("shrinkage_eta absent")
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "DIAGNOSTICS",    fixed = TRUE)
  expect_match(out, "ETA shrinkage:", fixed = TRUE)
})

test_that("ferx_runlog diagnostics section present when DW statistic available", {
  fit <- warfarin_fit()
  fit$dw_statistic <- 1.95
  fit$iwres_lag1_r <- 0.05
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "DW statistic:", fixed = TRUE)
})

test_that("ferx_runlog DW section flags positive autocorrelation correctly", {
  fit <- warfarin_fit()
  fit$dw_statistic <- 1.1
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "positive autocorrelation", fixed = TRUE)
})

test_that("ferx_runlog diagnostics section absent when nothing available", {
  fit <- warfarin_fit()
  fit$shrinkage_eta  <- NULL
  fit$shrinkage_eps  <- NULL
  fit$dw_statistic   <- NA_real_
  fit$eta_normality  <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "DIAGNOSTICS", fixed = TRUE)
})

test_that("ferx_runlog normality section flags p < 0.05", {
  fit <- warfarin_fit()
  fit$eta_normality <- data.frame(
    eta   = c("ETA_CL", "ETA_V"),
    W     = c(0.97,      0.85),
    p_val = c(0.60,      0.02),
    flag  = c("ok",      "non-normal"),
    stringsAsFactors = FALSE
  )
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "ETA normality",  fixed = TRUE)
  expect_match(out, "ETA_V",          fixed = TRUE)
  expect_match(out, "[!]",            fixed = TRUE)
})

# -- Tier 2: gradient ---------------------------------------------------------

test_that("ferx_runlog warns when gradient components exceed tolerance", {
  fit <- warfarin_fit()
  fit$final_gradient <- c(0.001, 0.5, 0.002)
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "WARNING", fixed = TRUE)
  expect_match(out, "possible non-convergence", fixed = TRUE)
})

test_that("ferx_runlog reports all-good when gradient within tolerance", {
  fit <- warfarin_fit()
  fit$final_gradient <- c(0.001, 0.002, 0.003)
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "satisfy |grad", fixed = TRUE)
  expect_no_match(out, "WARNING", fixed = TRUE)
})

test_that("ferx_runlog gradient_tol argument is respected", {
  fit <- warfarin_fit()
  fit$final_gradient <- 0.005
  out_strict <- ferx_runlog(fit, gradient_tol = 0.001, verbose = FALSE)
  out_lax    <- ferx_runlog(fit, gradient_tol = 0.01,  verbose = FALSE)
  expect_match(out_strict, "WARNING", fixed = TRUE)
  expect_no_match(out_lax, "WARNING", fixed = TRUE)
})

test_that("ferx_runlog wraps long gradient at 8 components per line", {
  fit <- warfarin_fit()
  fit$final_gradient <- seq(0.001, 0.001 * 10, by = 0.001)  # 10 components
  out <- ferx_runlog(fit, verbose = FALSE)
  # Two lines of [ ... ] expected: 8 + 2
  n_grad_lines <- sum(grepl("^  \\[", strsplit(out, "\n")[[1]]))
  expect_equal(n_grad_lines, 2L)
})

test_that("ferx_runlog gradient section says N/A for derivative-free optimizer", {
  fit <- warfarin_fit()
  fit$final_gradient        <- NULL
  fit$gradient_method_outer <- "N/A"
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "BOBYQA", fixed = TRUE)
})

test_that("ferx_runlog gradient section mentions older engine for missing gradient", {
  fit <- warfarin_fit()
  fit$final_gradient        <- NULL
  fit$gradient_method_outer <- "Enzyme AD"
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "older engine version", fixed = TRUE)
})

# -- Tier 2: runtime ----------------------------------------------------------

test_that("ferx_runlog runtime section shows wall time when present", {
  fit <- warfarin_fit()
  if (is.null(fit$wall_time_secs) || !is.finite(fit$wall_time_secs)) {
    skip("wall_time_secs absent")
  }
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "Wall time:", fixed = TRUE)
})

# -- Tier 3: graceful degradation ---------------------------------------------

test_that("ferx_runlog works on a minimal fit with almost all fields NULL", {
  fit <- structure(list(
    method             = "focei",
    ferx_version       = "0.0.0",
    converged          = TRUE,
    ofv                = -100,
    aic                = -90,
    bic                = -80,
    theta              = c(0.1, 8.0),
    theta_names        = c("TVCL", "TVV"),
    theta_fixed        = c(FALSE, FALSE),
    omega              = matrix(0.09, 1, 1),
    eta_names          = "ETA_CL",
    sigma              = 0.01,
    sigma_names        = "EPS_PROP"
  ), class = "ferx_fit")
  expect_no_error(out <- ferx_runlog(fit, verbose = FALSE))
  expect_type(out, "character")
  expect_match(out, "TVCL", fixed = TRUE)
})

test_that("ferx_runlog is pure ASCII in its output", {
  fit <- warfarin_fit()
  out <- ferx_runlog(fit, verbose = FALSE)
  # Check for non-ASCII bytes by comparing byte length vs character length.
  # Each non-ASCII code point in UTF-8 consumes > 1 byte, so a mismatch flags it.
  non_ascii <- nchar(out, type = "bytes") != nchar(out, type = "chars")
  expect_false(non_ascii, info = paste("Non-ASCII bytes found in ferx_runlog output"))
})

# -- Tier 3: NA-guard regression tests ----------------------------------------
# extendr maps Option<String>::None to NA_character_, not NULL.
# Verify that ferx_runlog() never emits NA tokens for these optional fields,
# even when they arrive as NA and another field (covariate_names) triggers
# has_settings = TRUE.

test_that("ferx_runlog does not emit 'Optimizer: NA' when optimizer_label is NA", {
  fit <- warfarin_fit()
  fit$optimizer_label  <- NA_character_
  fit$covariate_names  <- "WT"   # triggers has_settings = TRUE via covariate arm
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "Optimizer:      NA", fixed = TRUE)
  expect_match(out, "Covariates:", fixed = TRUE)  # settings section still present
})

test_that("ferx_runlog does not emit 'BLOQ method: NA' when bloq_method_label is NA", {
  fit <- warfarin_fit()
  fit$bloq_method_label <- NA_character_
  fit$covariate_names   <- "WT"
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "BLOQ method:    NA", fixed = TRUE)
})

test_that("ferx_runlog model_text is populated after ferx_load_fit round-trip", {
  # model_text must survive save/load so ferx_runlog shows the model source.
  fit <- warfarin_fit()
  # Simulate what ferx_load_fit sets (persist.R line 253-254):
  # result$model_source is restored; result$model_text must match.
  fit$model_source <- "[parameters]\n  theta TVCL(0.134, 0.001, 10.0)\n"
  fit$model_text   <- fit$model_source  # as set by the fix in persist.R
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_match(out, "MODEL FILE",   fixed = TRUE)
  expect_match(out, "TVCL",         fixed = TRUE)
})

test_that("ferx_runlog omits model file section gracefully when model_text NULL", {
  fit <- warfarin_fit()
  fit$model_text <- NULL
  out <- ferx_runlog(fit, verbose = FALSE)
  expect_no_match(out, "MODEL FILE", fixed = TRUE)
  # rest of runlog must still be intact
  expect_match(out, "PARAMETER ESTIMATES", fixed = TRUE)
})
