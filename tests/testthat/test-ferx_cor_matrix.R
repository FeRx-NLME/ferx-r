
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




# -- Held (FIX) coordinates in cov_matrix (#424) ------------------------------
#
# The covariance step gives a coordinate it held an all-zero row and column.
# `fit$cov_fixed` is the engine's mask of those; cor_matrix is NA there
# without a warning, and still warns for a non-positive variance on a
# parameter that was estimated.

test_that("cor_matrix warns, naming it, for a negative variance on an estimated parameter", {
  cm <- diag(c(-1, 0, 4))
  dimnames(cm) <- list(c("A", "B", "C"), c("A", "B", "C"))
  # B is held; A was estimated and came back negative. Exactly one warning:
  # sqrt() of the negative variance would add R's own "NaNs produced".
  w <- testthat::capture_warnings(
    cor <- .compute_cor_matrix(cm, fixed = c(FALSE, TRUE, FALSE))
  )
  expect_identical(w, paste0("Non-positive variance for estimated parameter(s) A ",
                             "in cov_matrix; their correlations are NA."))
  expected <- outer(c(TRUE, TRUE, FALSE), c(TRUE, TRUE, FALSE), "|")
  dimnames(expected) <- dimnames(cm)
  expect_identical(is.na(cor), expected)
  expect_identical(cor[["C", "C"]], 1)
})

test_that("cor_matrix keys on the mask: an all-zero row the engine says it estimated warns", {
  # Without a mask this row would read as held; with one saying "estimated",
  # a zero variance is degenerate and the user is told.
  expect_warning(
    cor <- .compute_cor_matrix(diag(c(1, 0, 4)), fixed = c(FALSE, FALSE, FALSE)),
    "estimated parameter\\(s\\) #2 in cov_matrix"
  )
  expect_identical(is.na(cor), outer(1:3 == 2L, 1:3 == 2L, "|"))
  # ...and a row the mask holds is NA even when it is not zero, but the
  # disagreement between mask and matrix is reported. (A held row that IS
  # zero stays quiet: B in the test above.)
  cm <- diag(c(1, 2, 4))
  dimnames(cm) <- list(c("A", "B", "C"), c("A", "B", "C"))
  w <- testthat::capture_warnings(
    cor <- .compute_cor_matrix(cm, fixed = c(FALSE, TRUE, FALSE))
  )
  expect_identical(w, paste0("cov_matrix has non-zero entries for held parameter(s) B; ",
                             "their correlations are NA."))
  expect_identical(unname(is.na(cor)), outer(1:3 == 2L, 1:3 == 2L, "|"))
})

test_that("the derived-fields step hands fit$cov_fixed to cor_matrix", {
  # The step ferx_fit(), ferx_covariance() and ferx_load_fit() share. The
  # mask calls an all-zero row estimated, so it warns; reading zeros instead
  # would take the row as held and stay quiet.
  populate <- getFromNamespace(".ferx_populate_derived_fields", "ferx")
  fit <- make_fake_fit(cov_matrix = diag(c(1, 0)), cov_fixed = c(FALSE, FALSE),
                       theta = c(CL = 1), omega = matrix(0.1, 1, 1),
                       eta_names = "ETA_CL", data_path = NA_character_)
  expect_warning(out <- populate(fit), "estimated parameter\\(s\\) #2 in cov_matrix")
  expect_identical(is.na(out$cor_matrix), outer(1:2 == 2L, 1:2 == 2L, "|"))
})

test_that("cor_matrix without a mask reads only an all-zero row and column as held", {
  # A bundle written before #424 carries no mask.
  expect_no_warning(cor <- .compute_cor_matrix(matrix(c(0, 0, 0, 4), 2, 2)))
  expect_identical(is.na(cor), matrix(c(TRUE, TRUE, TRUE, FALSE), 2, 2))
  # A zero variance with a non-zero covariance is not how a held row looks.
  expect_warning(
    .compute_cor_matrix(matrix(c(0, 0.1, 0.1, 4), 2, 2)),
    "estimated parameter\\(s\\) #1 in cov_matrix"
  )
})

# A warfarin fit holding a theta in the middle of the packed layout and the
# sigma at its end, so a mask shifted by one coordinate lands on an estimated
# row.
warfarin_fix_fit <- function(covariance = TRUE) {
  ex <- ferx_example("warfarin")
  txt <- readLines(ex$model)
  txt <- sub("^(\\s*theta TVV\\(.*\\))\\s*$", "\\1 FIX", txt)
  txt <- sub("^(\\s*sigma PROP_ERR ~ .*\\(sd\\))\\s*$", "\\1 FIX", txt)
  stopifnot(sum(grepl("FIX$", txt)) == 2L)
  mod <- tempfile(fileext = ".ferx")
  writeLines(txt, mod)
  ferx_fit(mod, ex$data, method = "focei", verbose = FALSE,
           covariance = covariance, settings = list(maxiter = 30L))
}

test_that("a fit with FIX parameters is quiet, and cor_matrix is NA exactly on them", {
  expect_no_warning(fit <- warfarin_fix_fit(), message = "[Nn]on-positive|held parameter")
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  held <- c(TVCL = FALSE, TVV = TRUE, TVKA = FALSE, ETA_CL = FALSE,
            ETA_V = FALSE, ETA_KA = FALSE, PROP_ERR = TRUE)
  expect_identical(fit$cov_fixed, held)
  expected <- outer(held, held, "|")
  expect_identical(is.na(fit$cor_matrix), expected)
  expect_identical(unname(diag(fit$cor_matrix)[!held]), rep(1, sum(!held)))

  # The mask, and with it the matrix, survives a save / load, quietly.
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  ferx_save_fit(fit, f)
  expect_no_warning(fit2 <- ferx_load_fit(f), message = "[Nn]on-positive|held parameter")
  # unname(): ferx_load_fit() does not restore cov_matrix's dimnames, a gap
  # that predates #424 (#417). The values have to match.
  expect_identical(unname(fit2$cov_fixed), unname(fit$cov_fixed))
  expect_identical(unname(fit2$cor_matrix), unname(fit$cor_matrix))
})

test_that("ferx_covariance() on a FIX fit carries the mask and is quiet", {
  fit <- warfarin_fix_fit(covariance = FALSE)
  expect_null(fit$cov_fixed)
  expect_no_warning(out <- ferx_covariance(fit), message = "[Nn]on-positive|held parameter")
  skip_if(is.null(out$cov_matrix), "covariance step did not converge - skipping")
  expect_identical(names(which(out$cov_fixed)), c("TVV", "PROP_ERR"))
  expect_identical(which(is.na(diag(out$cor_matrix))), c(TVV = 2L, PROP_ERR = 7L))
})

test_that("the bundled mbma_placebo fit (sigma FIX, weight = SE) is quiet, in fit and reload", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  expect_no_warning(fit <- ferx_fit(ex$model, ex$data, verbose = FALSE),
                    message = "[Nn]on-positive|held parameter")
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  # 22 thetas, ETA_E0, then ADD_ERR (FIX), then KAPPA_ARM: the packed order
  # puts sigma before kappa.
  expect_identical(which(fit$cov_fixed), 24L)
  expect_identical(is.na(fit$cor_matrix), outer(1:25 == 24L, 1:25 == 24L, "|"))
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  ferx_save_fit(fit, f)
  expect_no_warning(fit2 <- ferx_load_fit(f), message = "[Nn]on-positive|held parameter")
  expect_identical(fit2$cov_fixed, fit$cov_fixed)
  expect_identical(fit2$cor_matrix, fit$cor_matrix)
})
test_that(".ferx_compute_cor_matrix returns NULL when no covariance matrix is present", {
  expect_null(.compute_cor_matrix(NULL))
})
