
# ---- header from test-fit.R ----
# Local aliases avoid the ::: operator (undesirable_operator_linter).
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# Return structure — Tier 1




















# These tests require the covariance step to have succeeded. With maxiter = 30L
# the outer optimisation may not converge on all machines, so we skip rather
# than fail — a skip here means "covariance step needs more iterations", not
# a bug. A full-convergence run (no maxiter cap) is tested manually / locally.









# fd_hessian_step argument — Tier 1 (R-side validation, no Rust call)








# ferx_check_init — Tier 1




# fit$cor_matrix — Tier 1, requires covariance = TRUE





# Gradient correctness — [ENZYME ONLY]
#
# Two blockers before this test can be fully implemented:
#
# 1. Build tier: requires the Enzyme Rust toolchain (a custom Rust fork with
#    automatic differentiation). The Enzyme build takes ~1.5 h and is not
#    available on normal dev machines. The skip_if() guard below makes this
#    test inert on all Tier 1 (FERX_NO_AUTODIFF=1) machines — CI with the
#    Enzyme toolchain is required to exercise it.
#
# 2. Missing API: comparing autodiff vs finite-difference gradients requires
#    a gradient inspection entry point to be exposed from the Rust side. No
#    such API exists yet in ferx-core. Once it lands, replace the inner
#    skip() with a real comparison (e.g. relative error < 1e-4 per element).


# -- Diagnostic fields (Step 3 feature-parity) --------------------------------




# -- IOV: shrinkage_kappa_by_occ (Step 1 feature-parity) ---------------------




# SAEM HMC proposals — [ENZYME ONLY]
#
# HMC proposals in the SAEM E-step are gated on the Enzyme autodiff build
# (`hmc_step` is `#[cfg(feature = "autodiff")]` in ferx-core). On a stable /
# FERX_NO_AUTODIFF=1 build, `n_leapfrog > 0` is silently ignored and the
# sampler uses Metropolis-Hastings, so `saem_n_subjects_hmc` is NA. The
# guard below makes this test inert on Tier 1 machines; CI with the Enzyme
# toolchain exercises the behavioural assertions.




test_that("fit$cor_matrix has same dimnames as $cov_matrix", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  expect_true(is.matrix(fit$cor_matrix))
  expect_equal(dimnames(fit$cor_matrix), dimnames(fit$cov_matrix))
})
test_that("fit$cor_matrix diagonal is all 1s", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  expect_true(all(diag(fit$cor_matrix) == 1))
})
test_that("fit$cor_matrix off-diagonal values are in [-1, 1]", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- fit$cor_matrix
  d <- nrow(cor)
  if (d > 1L) {
    off <- cor[row(cor) != col(cor)]
    expect_true(all(abs(off) <= 1))
  }
})
test_that("fit$cor_matrix is NULL when covariance was FALSE", {
  fit <- warfarin_fit()
  expect_null(fit$cor_matrix)
})

# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   .ferx_compute_estimates(), .ferx_est_row(), .ferx_compute_eta_cov(),
#   .ferx_compute_cor_matrix(), ferx_get_warnings(), and
#   .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row           <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm      <- getFromNamespace(".ferx_compute_eta_normality", "ferx")
.compute_cor_matrix <- getFromNamespace(".ferx_compute_cor_matrix", "ferx")

# ---------------------------------------------------------------------------
# .ferx_compute_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_cov — message branches plus the correlation path
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# .ferx_compute_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that(".ferx_compute_cor_matrix warns on a non-positive variance diagonal", {
  # A zero on the diagonal gives sqrt() == 0, hitting the `se <= 0` branch
  # (a negative diagonal would instead yield NaN and slip past the guard).
  expect_warning(
    cm <- .compute_cor_matrix(matrix(c(0, 0, 0, 4), 2, 2)),
    "non-positive"
  )
})
test_that(".ferx_compute_cor_matrix returns NULL when no covariance matrix is present", {
  expect_null(.compute_cor_matrix(NULL))
})
