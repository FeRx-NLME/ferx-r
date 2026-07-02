
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

# ---- header from test-map-estimates.R ----
# Per-subject EBE / individual-parameter outputs.
#
# Pins the contract that ferx_fit() returns three per-subject diagnostic
# tables — `ebe_etas`, `individual_estimates`, and (for IOV models)
# `ebe_kappas` — and that `sdtab` no longer carries ETA columns. Also
# checks the downstream consumers (`eta_normality`, `ferx_eta_cov`) read
# from `ebe_etas`.
#
# Uses the shared `warfarin_fit()` helper (helper-warfarin-fit.R) so the
# real FOCEI fit is reused across test files.







# Warfarin data has no non-NONMEM covariate columns, so a synthetic WT column
# is added to exercise the ferx_eta_cov() happy path. One value per subject
# (constant within subject) so it qualifies as a covariate.
warfarin_dat_with_cov <- function() {
  ex  <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  ids <- unique(dat$ID)
  set.seed(42L)
  wt_map <- setNames(rnorm(length(ids), mean = 70, sd = 10), ids)
  dat$WT <- wt_map[as.character(dat$ID)]
  dat
}





test_that("ferx_eta_cov() consumes ebe_etas (errors when absent)", {
  fit <- warfarin_fit()
  bad <- fit
  bad$ebe_etas <- NULL
  expect_error(ferx_eta_cov(bad, data.frame(ID = 1L)), regexp = "ebe_etas")
})
test_that("ferx_eta_cov() returns a data frame with eta and covariate columns", {
  fit <- warfarin_fit()
  res <- ferx_eta_cov(fit, warfarin_dat_with_cov())
  expect_s3_class(res, "data.frame")
  expect_true("eta" %in% names(res))
  expect_true("covariate" %in% names(res))
  eta_names <- setdiff(names(fit$ebe_etas), "ID")
  expect_true(all(eta_names %in% unique(res$eta)))
})
test_that("ferx_eta_cov() values are finite numerics", {
  fit <- warfarin_fit()
  res <- ferx_eta_cov(fit, warfarin_dat_with_cov())
  expect_true(is.numeric(res$r))
  expect_true(all(is.finite(res$r)))
})
test_that("ferx_eta_cov() has one row per ETA-covariate pair", {
  fit     <- warfarin_fit()
  res     <- ferx_eta_cov(fit, warfarin_dat_with_cov())
  n_etas  <- length(setdiff(names(fit$ebe_etas), "ID"))
  # One covariate (WT) × n_etas rows expected
  expect_equal(nrow(res), n_etas)
})
