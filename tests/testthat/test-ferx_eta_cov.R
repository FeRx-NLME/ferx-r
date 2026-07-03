
# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   .ferx_compute_estimates(), .ferx_est_row(), .ferx_compute_eta_cov(),
#   .ferx_compute_cor_matrix(), ferx_get_warnings(), and
#   .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row           <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm      <- getFromNamespace(".ferx_compute_eta_normality", "ferx")
.compute_eta_cov   <- getFromNamespace(".ferx_compute_eta_cov", "ferx")

# ---------------------------------------------------------------------------
# .ferx_compute_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_cov — NULL branches plus the correlation path
# ---------------------------------------------------------------------------




# write a data.frame to a temp CSV and return its path, mirroring what
# .ferx_compute_eta_cov reads via fit$data_path.
.write_csv_path <- function(df) {
  path <- tempfile(fileext = ".csv")
  write.csv(df, path, row.names = FALSE)
  path
}

# ---------------------------------------------------------------------------
# .ferx_compute_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that(".ferx_compute_eta_cov returns NULL when there are no ETA columns", {
  path <- .write_csv_path(data.frame(ID = 1:3, WT = c(70, 80, 90)))
  res <- .compute_eta_cov(data.frame(ID = 1:3), path)
  expect_null(res)
})
test_that(".ferx_compute_eta_cov returns NULL when data has no numeric covariates", {
  path <- .write_csv_path(data.frame(ID = 1:3, SEX = c("M", "F", "M")))
  res <- .compute_eta_cov(data.frame(ID = 1:3, ETA_CL = c(.1, .2, .3)), path)
  expect_null(res)
})
test_that(".ferx_compute_eta_cov returns NULL when no covariate is constant per subject", {
  path <- .write_csv_path(data.frame(ID = c(1, 1, 2, 2), FOO = c(1, 2, 3, 4)))
  res <- .compute_eta_cov(data.frame(ID = c(1, 2), ETA_CL = c(.1, .2)), path)
  expect_null(res)
})
test_that(".ferx_compute_eta_cov returns a correlation table for constant covariates", {
  path <- .write_csv_path(data.frame(ID = 1:5, WT = c(70, 80, 90, 100, 110)))
  res <- .compute_eta_cov(data.frame(ID = 1:5, ETA_CL = c(.1, .2, .3, .4, .5)), path)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("eta", "covariate", "r", "p_val", "flag") %in% names(res)))
})
test_that(".ferx_compute_eta_cov yields NA stats when fewer than 3 pairs are usable", {
  path <- .write_csv_path(data.frame(ID = c(1, 2), WT = c(70, 80)))
  res <- .compute_eta_cov(data.frame(ID = c(1, 2), ETA_CL = c(.1, .2)), path)
  expect_true(is.na(res$r))
})
test_that(".ferx_compute_eta_cov returns NULL for missing ebe_etas or an unreadable data path", {
  expect_null(.compute_eta_cov(NULL, .write_csv_path(data.frame(x = 1))))
  expect_warning(
    res <- .compute_eta_cov(data.frame(ID = 1), file.path(tempdir(), "does-not-exist.csv")),
    "could not read data"
  )
  expect_null(res)
})
test_that(".ferx_compute_eta_cov stays silent when data_path is NA (no data bundled)", {
  expect_no_warning(res <- .compute_eta_cov(data.frame(ID = 1), NA_character_))
  expect_null(res)
})

# ---- header from test-map-estimates.R ----
# Per-subject EBE / individual-parameter outputs.
#
# Pins the contract that ferx_fit() returns three per-subject diagnostic
# tables — `ebe_etas`, `individual_estimates`, and (for IOV models)
# `ebe_kappas` — and that `sdtab` no longer carries ETA columns. Also
# checks the downstream consumer (`eta_normality`) reads from `ebe_etas`,
# and that `fit$eta_cov` is populated automatically from `ebe_etas` and
# `data_path` at the end of ferx_fit().
#
# Uses the shared `warfarin_fit()` helper (helper-warfarin-fit.R) so the
# real FOCEI fit is reused across test files.






# Warfarin data has no non-NONMEM covariate columns, so `fit$eta_cov` is NULL
# for a plain warfarin_fit(). A synthetic WT column (constant per subject) is
# written to its own CSV so .ferx_compute_eta_cov() has a real covariate to
# correlate against ebe_etas from a real fit.
warfarin_dat_with_cov_path <- function() {
  ex  <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  ids <- unique(dat$ID)
  set.seed(42L)
  wt_map <- setNames(rnorm(length(ids), mean = 70, sd = 10), ids)
  dat$WT <- wt_map[as.character(dat$ID)]
  .write_csv_path(dat)
}




test_that("fit$eta_cov is NULL for the warfarin example (no covariate columns)", {
  fit <- warfarin_fit()
  expect_null(fit$eta_cov)
})
test_that(".ferx_compute_eta_cov() returns a data frame with eta and covariate columns", {
  fit <- warfarin_fit()
  res <- .compute_eta_cov(fit$ebe_etas, warfarin_dat_with_cov_path())
  expect_s3_class(res, "data.frame")
  expect_true("eta" %in% names(res))
  expect_true("covariate" %in% names(res))
  eta_names <- setdiff(names(fit$ebe_etas), "ID")
  expect_true(all(eta_names %in% unique(res$eta)))
})
test_that(".ferx_compute_eta_cov() values are finite numerics", {
  fit <- warfarin_fit()
  res <- .compute_eta_cov(fit$ebe_etas, warfarin_dat_with_cov_path())
  expect_true(is.numeric(res$r))
  expect_true(all(is.finite(res$r)))
})
test_that(".ferx_compute_eta_cov() has one row per ETA-covariate pair", {
  fit    <- warfarin_fit()
  res    <- .compute_eta_cov(fit$ebe_etas, warfarin_dat_with_cov_path())
  n_etas <- length(setdiff(names(fit$ebe_etas), "ID"))
  # One covariate (WT) × n_etas rows expected
  expect_equal(nrow(res), n_etas)
})
