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

test_that("fit$ebe_etas has one row per subject with named ETA columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$ebe_etas, "data.frame")
  expect_equal(nrow(fit$ebe_etas), fit$n_subjects)
  expect_true("ID" %in% names(fit$ebe_etas))
  # Warfarin declares ETA_CL, ETA_V, ETA_KA in [parameters].
  expect_setequal(
    setdiff(names(fit$ebe_etas), "ID"),
    c("ETA_CL", "ETA_V", "ETA_KA")
  )
  expect_true(all(vapply(fit$ebe_etas[, -1], is.numeric, logical(1L))))
})

test_that("fit$individual_estimates has one row per subject with named parameter columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$individual_estimates, "data.frame")
  expect_equal(nrow(fit$individual_estimates), fit$n_subjects)
  expect_true("ID" %in% names(fit$individual_estimates))
  # Warfarin declares CL, V, KA in [individual_parameters].
  expect_setequal(
    setdiff(names(fit$individual_estimates), "ID"),
    c("CL", "V", "KA")
  )
  # Values come out of pk_param_fn — must be finite and positive for a
  # log-normal parameterisation.
  for (p in c("CL", "V", "KA")) {
    vals <- fit$individual_estimates[[p]]
    expect_true(all(is.finite(vals)), info = sprintf("%s has non-finite values", p))
    expect_true(all(vals > 0), info = sprintf("%s has non-positive values", p))
  }
})

test_that("ebe_etas and individual_estimates align on ID and row count", {
  fit <- warfarin_fit()
  expect_equal(fit$ebe_etas$ID, fit$individual_estimates$ID)
})

test_that("sdtab no longer contains ETA columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$sdtab, "data.frame")
  eta_cols <- grep("^ETA", names(fit$sdtab), value = TRUE)
  expect_length(eta_cols, 0L)
})

test_that("eta_normality reads from ebe_etas, one row per ETA", {
  fit <- warfarin_fit()
  expect_s3_class(fit$eta_normality, "data.frame")
  expect_setequal(fit$eta_normality$eta, c("ETA_CL", "ETA_V", "ETA_KA"))
  expect_true(all(c("eta", "W", "p_val", "flag") %in% names(fit$eta_normality)))
})

test_that("ferx_eta_cov() consumes ebe_etas (errors when absent)", {
  fit <- warfarin_fit()
  bad <- fit
  bad$ebe_etas <- NULL
  expect_error(ferx_eta_cov(bad, data.frame(ID = 1L)), regexp = "ebe_etas")
})

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
