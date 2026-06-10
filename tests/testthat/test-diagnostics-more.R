# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   ferx_estimates(), .ferx_est_row(), ferx_eta_cov(), ferx_cor_matrix(),
#   ferx_warnings(), and .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row       <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm  <- getFromNamespace(".ferx_compute_eta_normality", "ferx")

# ---------------------------------------------------------------------------
# ferx_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------
test_that("ferx_estimates builds a row per theta / omega / sigma / kappa", {
  fit <- make_fake_fit(
    theta      = c(1.0, 10.0),          # unnamed -> THETA1/THETA2 fallback
    se_theta   = c(0.1, 0.5),
    omega      = 0.09,                  # scalar -> coerced to 1x1
    se_omega   = 0.01,
    sigma      = c(0.05),
    sigma_types = "proportional",
    se_sigma   = 0.005,
    omega_iov  = 0.04,                  # scalar IOV -> KAPPA1
    se_kappa   = 0.004
  )
  tab <- ferx_estimates(fit)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("THETA1", "THETA2") %in% tab$param))
  expect_true(any(grepl("OMEGA", tab$param)))
  expect_true(any(grepl("SIGMA", tab$param)))
  expect_true(any(grepl("KAPPA", tab$param)))
})

test_that(".ferx_est_row back-transforms logit parameters", {
  row <- .est_row("LOGIT_P", estimate = 0.0, se = 0.2, transform = "logit",
                  init_as_sd = FALSE)
  expect_s3_class(row, "data.frame")
  expect_equal(row$estimate_natural, 0.5)            # inv_logit(0) == 0.5
  expect_true(row$lower_95_natural < row$upper_95_natural)
})

# ---------------------------------------------------------------------------
# ferx_eta_cov — message branches plus the correlation path
# ---------------------------------------------------------------------------
test_that("ferx_eta_cov messages when there are no ETA columns", {
  expect_message(
    res <- ferx_eta_cov(list(ebe_etas = data.frame(ID = 1:3)),
                        data.frame(ID = 1:3, WT = c(70, 80, 90))),
    "No ETA columns"
  )
  expect_null(res)
})

test_that("ferx_eta_cov messages when data has no numeric covariates", {
  expect_message(
    ferx_eta_cov(list(ebe_etas = data.frame(ID = 1:3, ETA_CL = c(.1, .2, .3))),
                 data.frame(ID = 1:3, SEX = c("M", "F", "M"))),
    "No numeric covariate"
  )
})

test_that("ferx_eta_cov messages when no covariate is constant per subject", {
  expect_message(
    ferx_eta_cov(
      list(ebe_etas = data.frame(ID = c(1, 2), ETA_CL = c(.1, .2))),
      data.frame(ID = c(1, 1, 2, 2), FOO = c(1, 2, 3, 4))  # varies within subject
    ),
    "constant-per-subject"
  )
})

test_that("ferx_eta_cov returns a correlation table for constant covariates", {
  fit  <- list(ebe_etas = data.frame(ID = 1:5, ETA_CL = c(.1, .2, .3, .4, .5)))
  data <- data.frame(ID = 1:5, WT = c(70, 80, 90, 100, 110))
  out <- capture.output(res <- ferx_eta_cov(fit, data))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("eta", "covariate", "r", "p_val", "flag") %in% names(res)))
})

test_that("ferx_eta_cov yields NA stats when fewer than 3 pairs are usable", {
  fit  <- list(ebe_etas = data.frame(ID = c(1, 2), ETA_CL = c(.1, .2)))
  data <- data.frame(ID = c(1, 2), WT = c(70, 80))
  out <- capture.output(res <- ferx_eta_cov(fit, data))
  expect_true(is.na(res$r))
})

test_that("ferx_eta_cov validates its inputs", {
  expect_error(ferx_eta_cov(list(ebe_etas = NULL), data.frame(x = 1)), "ebe_etas")
  expect_error(ferx_eta_cov(list(ebe_etas = data.frame(ID = 1)), NULL), "data.frame")
})

# ---------------------------------------------------------------------------
# ferx_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------
test_that("ferx_cor_matrix warns on a non-positive variance diagonal", {
  # A zero on the diagonal gives sqrt() == 0, hitting the `se <= 0` branch
  # (a negative diagonal would instead yield NaN and slip past the guard).
  fit <- make_fake_fit(cov_matrix = matrix(c(0, 0, 0, 4), 2, 2))
  expect_warning(out <- capture.output(cm <- ferx_cor_matrix(fit)), "non-positive")
})

test_that("ferx_cor_matrix errors when no covariance matrix is present", {
  expect_error(ferx_cor_matrix(make_fake_fit(cov_matrix = NULL)),
               "Covariance matrix not available")
})

# ---------------------------------------------------------------------------
# ferx_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------
test_that("ferx_warnings prints 'No warnings.' for an empty structured table", {
  empty <- data.frame(severity = character(0), category = character(0),
                      message = character(0), source_method = character(0),
                      stringsAsFactors = FALSE)
  fit <- make_fake_fit(warnings_structured = empty)
  out <- capture.output(ferx_warnings(fit))
  expect_true(any(grepl("No warnings", out)))
})

test_that("ferx_warnings labels each severity level", {
  df <- data.frame(
    severity      = c("critical", "warning", "info"),
    category      = c("convergence", "condition_number", "mu_referencing"),
    message       = c("did not converge", "ill-conditioned", "mu detected"),
    source_method = c("focei", "focei", "focei"),
    stringsAsFactors = FALSE
  )
  fit <- make_fake_fit(warnings_structured = df)
  out <- capture.output(ferx_warnings(fit))
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})

# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------
test_that(".ferx_compute_eta_normality returns NULL for missing or eta-less input", {
  expect_null(.compute_norm(NULL))
  expect_null(.compute_norm("not a df"))
  expect_null(.compute_norm(data.frame(ID = 1:3)))   # only ID, no eta columns
})

test_that(".ferx_compute_eta_normality subsamples and warns above the 5000 cap", {
  set.seed(1)
  big <- data.frame(ID = seq_len(5001L), ETA_CL = rnorm(5001L))
  expect_warning(res <- .compute_norm(big), "exceeds")
  expect_s3_class(res, "data.frame")
  expect_true(is.finite(res$W))
})

test_that(".ferx_compute_eta_normality flags non-normal etas", {
  set.seed(2)
  df <- data.frame(ID = 1:200, ETA_SKEW = rexp(200))  # clearly non-normal
  res <- .compute_norm(df)
  expect_true(res$p_val < 0.05)
  expect_identical(res$flag, "[!]")
})
