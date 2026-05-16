# Tests for parameter transform display: CI helpers, ferx_estimates(), print.ferx_fit
# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# .ferx_inv_logit helper ------------------------------------------------

test_that(".ferx_inv_logit gives expected values", {
  expect_equal(ferx:::.ferx_inv_logit(0),    0.5,     tolerance = 1e-10)
  expect_equal(ferx:::.ferx_inv_logit(Inf),  1,       tolerance = 1e-10)
  expect_equal(ferx:::.ferx_inv_logit(-Inf), 0,       tolerance = 1e-10)
  # logit(0.7) ≈ 0.8473; inv_logit(0.8473) ≈ 0.7
  expect_equal(ferx:::.ferx_inv_logit(log(0.7 / 0.3)), 0.7, tolerance = 1e-6)
})

# ferx_estimates() — identity transform ---------------------------------

test_that("ferx_estimates() identity theta: symmetric CI, no natural-scale columns", {
  fit <- make_fake_fit(
    theta         = c(CL = 0.134),
    se_theta      = c(CL = 0.02),
    theta_transforms = "identity",
    omega         = matrix(0.07, 1, 1),
    sigma         = 0.01,
    sigma_types   = "proportional"
  )
  est <- ferx_estimates(fit)
  theta_row <- est[est$param == "CL", ]
  expect_equal(theta_row$transform,   "identity")
  expect_equal(theta_row$lower_95,    0.134 - 1.96 * 0.02, tolerance = 1e-8)
  expect_equal(theta_row$upper_95,    0.134 + 1.96 * 0.02, tolerance = 1e-8)
  expect_true(is.na(theta_row$estimate_natural))
  expect_true(is.na(theta_row$lower_95_natural))
  expect_true(is.na(theta_row$upper_95_natural))
})

# ferx_estimates() — log transform --------------------------------------

test_that("ferx_estimates() log theta: asymmetric CI and natural back-transform", {
  log_est <- log(0.134)
  se      <- 0.15
  fit <- make_fake_fit(
    theta            = c(TVCL = log_est),
    se_theta         = c(TVCL = se),
    theta_transforms = "log",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01,
    sigma_types      = "proportional"
  )
  est <- ferx_estimates(fit)
  theta_row <- est[est$param == "TVCL", ]
  expect_equal(theta_row$transform,        "log")
  expect_equal(theta_row$estimate_natural, exp(log_est), tolerance = 1e-8)
  expect_equal(theta_row$lower_95_natural, exp(log_est - 1.96 * se), tolerance = 1e-8)
  expect_equal(theta_row$upper_95_natural, exp(log_est + 1.96 * se), tolerance = 1e-8)
  # Log-scale CI on estimate column (not back-transformed)
  expect_equal(theta_row$lower_95, log_est - 1.96 * se, tolerance = 1e-8)
  expect_equal(theta_row$upper_95, log_est + 1.96 * se, tolerance = 1e-8)
})

# ferx_estimates() — logit transform ------------------------------------

test_that("ferx_estimates() logit theta: asymmetric natural CI via inv_logit", {
  logit_est <- log(0.7 / 0.3)  # logit(0.7)
  se        <- 0.3
  fit <- make_fake_fit(
    theta            = c(THETA_F = logit_est),
    se_theta         = c(THETA_F = se),
    theta_transforms = "logit",
    omega            = matrix(0.10, 1, 1),
    sigma            = 0.01,
    sigma_types      = "proportional"
  )
  est <- ferx_estimates(fit)
  theta_row <- est[est$param == "THETA_F", ]
  inv_logit <- function(x) 1 / (1 + exp(-x))
  expect_equal(theta_row$transform,        "logit")
  expect_equal(theta_row$estimate_natural, inv_logit(logit_est), tolerance = 1e-6)
  expect_equal(theta_row$lower_95_natural, inv_logit(logit_est - 1.96 * se), tolerance = 1e-6)
  expect_equal(theta_row$upper_95_natural, inv_logit(logit_est + 1.96 * se), tolerance = 1e-6)
})

# ferx_estimates() — sigma types ----------------------------------------

test_that("ferx_estimates() sigma has correct transform column per type", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = c(0.01, 0.05),
    sigma_types      = c("proportional", "additive")
  )
  est <- ferx_estimates(fit)
  sig1 <- est[est$param == "SIGMA(1)", ]
  sig2 <- est[est$param == "SIGMA(2)", ]
  expect_equal(sig1$transform, "proportional")
  expect_equal(sig2$transform, "additive")
})

# ferx_estimates() — omega has variance transform -----------------------

test_that("ferx_estimates() omega rows carry transform = 'variance'", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- ferx_estimates(fit)
  omega_row <- est[grepl("OMEGA", est$param), ]
  expect_equal(omega_row$transform, "variance")
})

# ferx_estimates() — bare eta/sigma names when declared ------------------

test_that("ferx_estimates() uses bare eta_names for omega rows, not OMEGA() wrapper", {
  fit <- make_fake_fit(
    theta            = c(TVCL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.09, 1, 1),
    eta_names        = "ETA_CL",
    sigma            = 0.1,
    sigma_names      = "EPS_PROP",
    sigma_types      = "proportional"
  )
  est <- ferx_estimates(fit)
  expect_true("ETA_CL"  %in% est$param)
  expect_true("EPS_PROP" %in% est$param)
  expect_false(any(grepl("OMEGA\\(ETA", est$param)))
})

test_that("ferx_estimates() falls back to OMEGA(i,i) when eta_names absent", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- ferx_estimates(fit)
  expect_true("OMEGA(1,1)" %in% est$param)
})

test_that("print.ferx_fit uses bare eta_names in OMEGA section", {
  fit <- make_fake_fit(
    omega         = matrix(0.09, 1, 1),
    eta_names     = "ETA_CL",
    eta_param_types = "log_normal"
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("ETA_CL", out)))
  expect_false(any(grepl("OMEGA\\(ETA_CL\\)", out)))
})

# ferx_estimates() — no SE → NA columns ---------------------------------

test_that("ferx_estimates() returns NA for SE-derived columns when se_theta absent", {
  fit <- make_fake_fit(
    theta            = c(CL = 0.134),
    se_theta         = NULL,
    theta_transforms = "log",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- ferx_estimates(fit)
  theta_row <- est[est$param == "CL", ]
  expect_true(is.na(theta_row$se))
  expect_true(is.na(theta_row$lower_95))
  expect_true(is.na(theta_row$estimate_natural))
})

# print.ferx_fit — THETA section transform tags -------------------------

test_that("print.ferx_fit emits [log scale] tag for log-transformed theta", {
  fit <- make_fake_fit(
    theta            = c(TVCL = log(0.134)),
    se_theta         = c(TVCL = 0.15),
    theta_transforms = "log",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  theta_line <- out[grepl("TVCL", out)]
  expect_true(any(grepl("log scale", theta_line)))
  typical_line <- out[grepl("typical", out)]
  expect_true(length(typical_line) >= 1L)
  expect_true(any(grepl("95%", typical_line)))
})

test_that("print.ferx_fit emits [logit scale] tag for logit-transformed theta", {
  fit <- make_fake_fit(
    theta            = c(THETA_F = log(0.7 / 0.3)),
    se_theta         = c(THETA_F = 0.3),
    theta_transforms = "logit",
    omega            = matrix(0.10, 1, 1),
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  theta_line <- out[grepl("THETA_F", out)]
  expect_true(any(grepl("logit scale", theta_line)))
})

# print.ferx_fit — OMEGA section ETA type labels ------------------------

test_that("print.ferx_fit shows [log-normal] label and CV% for log_normal ETA", {
  fit <- make_fake_fit(
    theta           = c(CL = 1),
    omega           = matrix(0.40, 1, 1),
    eta_param_types = "log_normal",
    sigma           = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("log-normal", omega_line)))
  expect_true(any(grepl("CV%", omega_line)))
  expect_true(any(grepl("70\\.1", omega_line)))
})

test_that("print.ferx_fit shows [additive] label and SD for additive ETA", {
  fit <- make_fake_fit(
    theta            = c(TVTLAG = 0.5),
    theta_transforms = "identity",
    omega            = matrix(0.04, 1, 1),
    eta_param_types  = "additive",
    eta_linked_theta = "TVTLAG",
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("additive", omega_line)))
  expect_true(any(grepl("SD = 0\\.2", omega_line)))
})

test_that("print.ferx_fit shows [logit] label and +/-1SD range for logit ETA", {
  logit_theta <- log(0.7 / 0.3)
  fit <- make_fake_fit(
    theta            = c(THETA_F = logit_theta),
    theta_transforms = "logit",
    omega            = matrix(0.10, 1, 1),
    eta_param_types  = "logit",
    eta_linked_theta = "THETA_F",
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("logit", omega_line)))
  expect_true(any(grepl("+/-1SD", omega_line, fixed = TRUE)))
})

# print.ferx_fit — SIGMA section type labels ----------------------------

test_that("print.ferx_fit shows [proportional] label and CV% for proportional sigma", {
  fit <- make_fake_fit(
    theta       = c(CL = 1),
    omega       = matrix(0.07, 1, 1),
    sigma       = 0.01,
    sigma_types = "proportional"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("SIGMA\\(1\\)", out)]
  expect_true(any(grepl("proportional", sigma_line)))
  expect_true(any(grepl("CV%", sigma_line)))
})

test_that("print.ferx_fit shows [additive] label and SD for additive sigma", {
  fit <- make_fake_fit(
    theta       = c(CL = 1),
    omega       = matrix(0.07, 1, 1),
    sigma       = 1.0,
    sigma_types = "additive"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("SIGMA\\(1\\)", out)]
  expect_true(any(grepl("additive", sigma_line)))
  expect_true(any(grepl("SD = 1\\.0", sigma_line)))
})

# backward compatibility — absent transform fields ----------------------

test_that("print.ferx_fit works without eta_param_types (falls back to log_normal)", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), sigma = 0.01)
  # Remove transform fields to simulate old ferx-core binary
  fit$eta_param_types  <- NULL
  fit$theta_transforms <- NULL
  fit$sigma_types      <- NULL
  expect_no_error(capture.output(print(fit)))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("70\\.1", omega_line)))  # still shows log-normal CV%
})
