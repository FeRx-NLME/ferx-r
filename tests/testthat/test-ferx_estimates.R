
# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   .compute_estimates(), .ferx_est_row(), `.ferx_compute_eta_cov()`, `.ferx_compute_cor_matrix()`,
#   ferx_get_warnings(), and .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row           <- getFromNamespace(".ferx_est_row",           "ferx")
.compute_estimates <- getFromNamespace(".ferx_compute_estimates", "ferx")
.compute_norm  <- getFromNamespace(".ferx_compute_eta_normality", "ferx")

# ---------------------------------------------------------------------------
# `.ferx_compute_estimates()` — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# `.ferx_compute_eta_cov()` — message branches plus the correlation path
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# `.ferx_compute_cor_matrix()` — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that(".ferx_compute_estimates() builds a row per theta / omega / sigma / kappa", {
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
  tab <- .compute_estimates(fit)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("THETA1", "THETA2") %in% tab$param))
  expect_true(any(grepl("OMEGA", tab$param)))
  expect_true(any(grepl("SIGMA", tab$param)))
  expect_true(any(grepl("KAPPA", tab$param)))
})
test_that(".ferx_compute_estimates() carries a weighted kappa's weight expression", {
  # Sample-size-weighted IOV (ferx-core #1031). `estimate` on a weighted kappa
  # row is the *unweighted* gamma^2; without the weight beside it, a reader of
  # this table would take sqrt(estimate) for the between-occasion SD.
  fit <- make_fake_fit(
    theta       = c(TVCL = 1.0),
    omega       = 0.09,
    sigma       = 0.05,
    sigma_types = "proportional",
    omega_iov   = diag(c(0.04, 1.84)),
    kappa_names = c("KAPPA_CL", "KAPPA_EMAX"),
    kappa_weights = c(KAPPA_CL = NA_character_, KAPPA_EMAX = "NARM")
  )
  tab <- .compute_estimates(fit)
  expect_true("weight" %in% names(tab))
  expect_equal(tab$weight[tab$param == "KAPPA_EMAX"], "NARM")
  # Unweighted kappa, and every non-kappa row, stay NA.
  expect_true(is.na(tab$weight[tab$param == "KAPPA_CL"]))
  expect_true(all(is.na(tab$weight[!grepl("^KAPPA", tab$param)])))
})

test_that(".ferx_compute_estimates() leaves weight NA for an unweighted IOV model", {
  fit <- make_fake_fit(
    theta       = c(TVCL = 1.0),
    omega       = 0.09,
    sigma       = 0.05,
    sigma_types = "proportional",
    omega_iov   = 0.04,
    kappa_names = "KAPPA_CL"
  )
  tab <- .compute_estimates(fit)
  expect_true("weight" %in% names(tab))
  expect_true(all(is.na(tab$weight)))
})

test_that(".ferx_est_row back-transforms logit parameters", {
  row <- .est_row("LOGIT_P", estimate = 0.0, se = 0.2, transform = "logit",
                  init_as_sd = FALSE)
  expect_s3_class(row, "data.frame")
  expect_equal(row$estimate_natural, 0.5)            # inv_logit(0) == 0.5
  expect_true(row$lower_95_natural < row$upper_95_natural)
})

test_that(".ferx_est_row does not back-transform a logit_probability estimate (#371)", {
  # The engine reports a `logit_probability` theta already on (0, 1); a second
  # inv_logit() would report 0.6225 for a bioavailability of 0.5.
  row <- .est_row("F1", estimate = 0.5, se = 0.02,
                  transform = "logit_probability", init_as_sd = FALSE)
  expect_equal(row$estimate_natural, 0.5)
  expect_false(isTRUE(all.equal(row$estimate_natural, 1 / (1 + exp(-0.5)))))
  # `lower_95`/`upper_95` stay the symmetric Wald on the reported scale, as for
  # every other transform.
  expect_equal(row$lower_95, 0.5 - 1.96 * 0.02)
  expect_equal(row$upper_95, 0.5 + 1.96 * 0.02)
})

test_that(".ferx_est_row forms the logit_probability natural CI on the logit scale (#371)", {
  # A symmetric Wald on the probability scale escapes (0, 1) whenever the SE is
  # large relative to the estimate; the natural-scale interval is formed where
  # the parameter is unbounded and brought back.
  p <- 0.7923592; se <- 0.2485106
  row <- .est_row("F1", estimate = p, se = se,
                  transform = "logit_probability", init_as_sd = FALSE)
  se_logit <- se / (p * (1 - p))
  lg       <- log(p / (1 - p))
  inv      <- function(x) 1 / (1 + exp(-x))
  expect_equal(row$lower_95_natural, inv(lg - 1.96 * se_logit))
  expect_equal(row$upper_95_natural, inv(lg + 1.96 * se_logit))
  expect_true(row$lower_95_natural > 0 && row$upper_95_natural < 1)
  # the interval it replaces did not fit inside the parameter's own support
  expect_gt(row$upper_95, 1)
})

test_that("the two parameterisations of one probability agree on the natural columns (#371)", {
  # `logit_probability` reports p; `logit` reports logit(p) with the SE carried
  # by the delta method. The natural-scale columns must not be able to tell
  # which way the model was written.
  p  <- 0.7923592
  se <- 0.2485106
  a  <- .est_row("F_prob",  estimate = p,               se = se,
                 transform = "logit_probability")
  b  <- .est_row("F_logit", estimate = log(p / (1 - p)), se = se / (p * (1 - p)),
                 transform = "logit")
  expect_equal(a$estimate_natural, b$estimate_natural, tolerance = 1e-12)
  expect_equal(a$lower_95_natural, b$lower_95_natural, tolerance = 1e-12)
  expect_equal(a$upper_95_natural, b$upper_95_natural, tolerance = 1e-12)
})

test_that(".ferx_est_row suppresses the logit_probability natural CI when it is undefined", {
  # logit(1) is not finite and the Jacobian p(1-p) collapses; report no
  # interval rather than an infinite bound.
  row <- .est_row("F1", estimate = 1, se = 0.02, transform = "logit_probability")
  expect_equal(row$estimate_natural, 1)
  expect_true(is.na(row$lower_95_natural))
  expect_true(is.na(row$upper_95_natural))
  # and with no SE at all
  row2 <- .est_row("F1", estimate = 0.5, se = NA_real_, transform = "logit_probability")
  expect_true(is.na(row2$lower_95_natural))
  expect_true(is.na(row2$upper_95_natural))
})

test_that(".ferx_logit is vectorised and NA outside the open interval", {
  .logit <- getFromNamespace(".ferx_logit", "ferx")
  expect_equal(.logit(c(0.5, 0.25)), c(0, log(0.25 / 0.75)))
  expect_equal(.logit(c(0, 1, -1, 2, NA, NaN, Inf, -Inf)), rep(NA_real_, 8))
  expect_equal(.logit(numeric(0)), numeric(0))
  expect_equal(.logit(NULL), numeric(0))
  expect_equal(.logit(c(F1 = 0.5)), 0)   # names dropped, value kept
})

# ---- header from test-diagnostics.R ----
# check_diagnostics() — Tier 1











# --- .compute_estimates() init_as_sd column ---




test_that(".ferx_compute_estimates() returns init_as_sd column, TRUE for annotated omega/sigma", {
  fit <- structure(list(
    theta         = c(TVCL = 1.0),
    theta_names   = "TVCL",
    se_theta      = NULL,
    theta_transforms = "log",
    omega         = matrix(0.09, 1L, 1L),
    eta_names     = "ETA_CL",
    se_omega      = NULL,
    omega_init_as_sd = TRUE,
    sigma         = c(PROP_ERR = 0.01),
    sigma_names   = "PROP_ERR",
    se_sigma      = NULL,
    sigma_types   = "proportional",
    sigma_init_as_sd = FALSE,
    omega_iov     = NULL
  ), class = "ferx_fit")

  est <- .compute_estimates(fit)
  expect_true("init_as_sd" %in% names(est))
  expect_false(est[est$param == "TVCL", "init_as_sd"])
  expect_true(est[est$param == "ETA_CL", "init_as_sd"])
  expect_false(est[est$param == "PROP_ERR", "init_as_sd"])
})
test_that(".ferx_compute_estimates() init_as_sd is FALSE for all rows when flags absent (old fit)", {
  fit <- structure(list(
    theta         = c(TVCL = 1.0),
    theta_names   = "TVCL",
    se_theta      = NULL,
    theta_transforms = "identity",
    omega         = matrix(0.09, 1L, 1L),
    eta_names     = "ETA_CL",
    se_omega      = NULL,
    sigma         = c(PROP_ERR = 0.01),
    sigma_names   = "PROP_ERR",
    se_sigma      = NULL,
    sigma_types   = "proportional",
    omega_iov     = NULL
    # omega_init_as_sd, sigma_init_as_sd intentionally absent
  ), class = "ferx_fit")

  est <- .compute_estimates(fit)
  expect_true("init_as_sd" %in% names(est))
  expect_true(all(!est$init_as_sd))
})

# ---- header from test-transform-output.R ----
# Tests for parameter transform display: CI helpers, .compute_estimates(), print.ferx_fit
# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# .ferx_inv_logit helper ------------------------------------------------


# .compute_estimates() — identity transform ---------------------------------


# .compute_estimates() — log transform --------------------------------------


# .compute_estimates() — logit transform ------------------------------------


# .compute_estimates() — sigma types ----------------------------------------


# .compute_estimates() — omega has variance transform -----------------------


# .compute_estimates() — bare eta/sigma names when declared ------------------



# .compute_estimates() — no SE → NA columns ---------------------------------


# print.ferx_fit — THETA section transform tags -------------------------



# print.ferx_fit — OMEGA section ETA type labels ------------------------




# print.ferx_fit — SIGMA section type labels ----------------------------



# backward compatibility — absent transform fields ----------------------



test_that(".compute_estimates() identity theta: symmetric CI, no natural-scale columns", {
  fit <- make_fake_fit(
    theta         = c(CL = 0.134),
    se_theta      = c(CL = 0.02),
    theta_transforms = "identity",
    omega         = matrix(0.07, 1, 1),
    sigma         = 0.01,
    sigma_types   = "proportional"
  )
  est <- .compute_estimates(fit)
  theta_row <- est[est$param == "CL", ]
  expect_equal(theta_row$transform,   "identity")
  expect_equal(theta_row$lower_95,    0.134 - 1.96 * 0.02, tolerance = 1e-8)
  expect_equal(theta_row$upper_95,    0.134 + 1.96 * 0.02, tolerance = 1e-8)
  expect_true(is.na(theta_row$estimate_natural))
  expect_true(is.na(theta_row$lower_95_natural))
  expect_true(is.na(theta_row$upper_95_natural))
})
test_that(".compute_estimates() log theta: asymmetric CI and natural back-transform", {
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
  est <- .compute_estimates(fit)
  theta_row <- est[est$param == "TVCL", ]
  expect_equal(theta_row$transform,        "log")
  expect_equal(theta_row$estimate_natural, exp(log_est), tolerance = 1e-8)
  expect_equal(theta_row$lower_95_natural, exp(log_est - 1.96 * se), tolerance = 1e-8)
  expect_equal(theta_row$upper_95_natural, exp(log_est + 1.96 * se), tolerance = 1e-8)
  # Log-scale CI on estimate column (not back-transformed)
  expect_equal(theta_row$lower_95, log_est - 1.96 * se, tolerance = 1e-8)
  expect_equal(theta_row$upper_95, log_est + 1.96 * se, tolerance = 1e-8)
})
test_that(".compute_estimates() logit theta: asymmetric natural CI via inv_logit", {
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
  est <- .compute_estimates(fit)
  theta_row <- est[est$param == "THETA_F", ]
  inv_logit <- function(x) 1 / (1 + exp(-x))
  expect_equal(theta_row$transform,        "logit")
  expect_equal(theta_row$estimate_natural, inv_logit(logit_est), tolerance = 1e-6)
  expect_equal(theta_row$lower_95_natural, inv_logit(logit_est - 1.96 * se), tolerance = 1e-6)
  expect_equal(theta_row$upper_95_natural, inv_logit(logit_est + 1.96 * se), tolerance = 1e-6)
})
test_that(".compute_estimates() sigma has correct transform column per type", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = c(0.01, 0.05),
    sigma_types      = c("proportional", "additive")
  )
  est <- .compute_estimates(fit)
  sig1 <- est[est$param == "SIGMA(1)", ]
  sig2 <- est[est$param == "SIGMA(2)", ]
  expect_equal(sig1$transform, "proportional")
  expect_equal(sig2$transform, "additive")
})
test_that(".compute_estimates() omega rows carry transform = 'variance'", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- .compute_estimates(fit)
  omega_row <- est[grepl("OMEGA", est$param), ]
  expect_equal(omega_row$transform, "variance")
})
test_that(".compute_estimates() uses bare eta_names for omega rows, not OMEGA() wrapper", {
  fit <- make_fake_fit(
    theta            = c(TVCL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.09, 1, 1),
    eta_names        = "ETA_CL",
    sigma            = 0.1,
    sigma_names      = "EPS_PROP",
    sigma_types      = "proportional"
  )
  est <- .compute_estimates(fit)
  expect_true("ETA_CL"  %in% est$param)
  expect_true("EPS_PROP" %in% est$param)
  expect_false(any(grepl("OMEGA\\(ETA", est$param)))
})
test_that(".compute_estimates() falls back to OMEGA(i,i) when eta_names absent", {
  fit <- make_fake_fit(
    theta            = c(CL = 1),
    theta_transforms = "identity",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- .compute_estimates(fit)
  expect_true("OMEGA(1,1)" %in% est$param)
})
test_that(".compute_estimates() returns NA for SE-derived columns when se_theta absent", {
  fit <- make_fake_fit(
    theta            = c(CL = 0.134),
    se_theta         = NULL,
    theta_transforms = "log",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  est <- .compute_estimates(fit)
  theta_row <- est[est$param == "CL", ]
  expect_true(is.na(theta_row$se))
  expect_true(is.na(theta_row$lower_95))
  expect_true(is.na(theta_row$estimate_natural))
})

# ---------------------------------------------------------------------------
# `.ferx_compute_estimates()` - row names (#299)
# ---------------------------------------------------------------------------

test_that(".ferx_compute_estimates() sets row names to the parameter names", {
  fit <- make_fake_fit(
    theta       = c(TVCL = 1.0, TVV = 10.0),
    se_theta    = c(0.1, 0.5),
    omega       = 0.09,
    eta_names   = "ETA_CL",
    se_omega    = 0.01,
    sigma       = 0.05,
    sigma_names = "EPS_PROP",
    sigma_types = "proportional",
    se_sigma    = 0.005
  )
  est <- .compute_estimates(fit)
  expect_identical(rownames(est), c("TVCL", "TVV", "ETA_CL", "EPS_PROP"))
  # The papercut this fixes: name-indexing used to return a row of NA.
  expect_equal(est["TVCL", "estimate"], 1.0)
  expect_equal(est["ETA_CL", "estimate"], 0.09)
  # `param` still carries the name as declared.
  expect_identical(est$param, c("TVCL", "TVV", "ETA_CL", "EPS_PROP"))
})

test_that(".ferx_compute_estimates() qualifies a name declared in two blocks", {
  # A theta and an eta may collide. Both get block-qualified, so neither owns
  # the bare name and both stay addressable - the bare name addressing only the
  # first is the silent-wrong-coefficient failure #299 is about.
  fit <- make_fake_fit(
    theta     = c(CL = 1.0),
    omega     = 0.09,
    eta_names = "CL",
    sigma     = NULL
  )
  est <- .compute_estimates(fit)
  expect_identical(rownames(est), c("CL.theta", "CL.omega"))
  expect_identical(est$param, c("CL", "CL"))
  expect_equal(est["CL.theta", "estimate"], 1.0)
  expect_equal(est["CL.omega", "estimate"], 0.09)
})

test_that(".ferx_compute_estimates() leaves non-colliding names bare", {
  # The qualification is only for collisions; the ordinary table is untouched.
  fit <- make_fake_fit(
    theta       = c(TVCL = 1.0),
    omega       = 0.09,
    eta_names   = "ETA_CL",
    sigma       = 0.05,
    sigma_names = "EPS_PROP",
    sigma_types = "proportional"
  )
  expect_identical(rownames(.compute_estimates(fit)),
                   c("TVCL", "ETA_CL", "EPS_PROP"))
})
