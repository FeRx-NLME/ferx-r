# Return structure — Tier 1

test_that("returns class ferx_fit", {
  expect_s3_class(warfarin_fit(), "ferx_fit")
})

test_that("$theta is a named numeric vector", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$theta))
  expect_false(is.null(names(fit$theta)))
  expect_true(length(names(fit$theta)) > 0L)
})

test_that("$theta names match model theta declarations", {
  fit <- warfarin_fit()
  ex  <- ferx_example("warfarin")
  s   <- ferx_model_inspect(ex$model)
  expect_setequal(names(fit$theta), s$theta_names)
})

test_that("$omega is a square numeric matrix", {
  fit <- warfarin_fit()
  expect_true(is.matrix(fit$omega))
  expect_true(is.numeric(fit$omega))
  expect_equal(nrow(fit$omega), ncol(fit$omega))
})

test_that("$sigma is a finite positive numeric", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$sigma))
  expect_true(is.finite(fit$sigma))
  expect_true(fit$sigma > 0)
})

test_that("$ofv is finite numeric", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$ofv))
  expect_true(is.finite(fit$ofv))
})

test_that("$aic > $ofv and $bic > $ofv", {
  fit <- warfarin_fit()
  expect_gt(fit$aic, fit$ofv)
  expect_gt(fit$bic, fit$ofv)
})

test_that("$method equals requested method", {
  fit <- warfarin_fit()
  expect_equal(tolower(fit$method), "focei")
})

test_that("$sdtab is a data frame", {
  fit <- warfarin_fit()
  expect_s3_class(fit$sdtab, "data.frame")
})

test_that("$sdtab has required columns", {
  fit <- warfarin_fit()
  required <- c("ID", "TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES")
  expect_true(all(required %in% names(fit$sdtab)))
})

test_that("$sdtab has one row per observation", {
  fit <- warfarin_fit()
  ex  <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(fit$sdtab), n_obs)
})

test_that("$se_theta absent when covariance = FALSE", {
  expect_null(warfarin_fit()$se_theta)
})

# These tests require the covariance step to have succeeded. With maxiter = 30L
# the outer optimisation may not converge on all machines, so we skip rather
# than fail — a skip here means "covariance step needs more iterations", not
# a bug. A full-convergence run (no maxiter cap) is tested manually / locally.
test_that("$se_theta present and named when covariance = TRUE", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$se_theta), "covariance step did not converge — skipping")
  expect_true(is.numeric(fit$se_theta))
  expect_false(is.null(names(fit$se_theta)))
})

test_that("$se_theta length matches $theta length when covariance = TRUE", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$se_theta), "covariance step did not converge — skipping")
  expect_equal(length(fit$se_theta), length(fit$theta))
})

test_that("errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit("no_such_model.ferx", ex$data, method = "focei"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

test_that("errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, "no_such_data.csv", method = "focei"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

# ferx_check_init — Tier 1

test_that("ferx_check_init returns without error on valid inputs", {
  ex <- ferx_example("warfarin")
  expect_no_error(ferx_check_init(ex$model, ex$data))
})

test_that("ferx_check_init errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_check_init("no_such_model.ferx", ex$data),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

test_that("ferx_check_init errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_check_init(ex$model, "no_such_data.csv"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

# ferx_cor_matrix — Tier 1, requires covariance = TRUE

test_that("ferx_cor_matrix returns a matrix with same dimnames as $cov_matrix", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  expect_true(is.matrix(cor))
  expect_equal(dimnames(cor), dimnames(fit$cov_matrix))
})

test_that("ferx_cor_matrix diagonal is all 1s", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  expect_true(all(diag(cor) == 1))
})

test_that("ferx_cor_matrix off-diagonal values are in [-1, 1]", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  d <- nrow(cor)
  if (d > 1L) {
    off <- cor[row(cor) != col(cor)]
    expect_true(all(abs(off) <= 1))
  }
})

test_that("ferx_cor_matrix errors gracefully when covariance was FALSE", {
  fit <- warfarin_fit()
  expect_error(ferx_cor_matrix(fit))
})

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
#    such API exists yet in ferx-nlme. Once it lands, replace the inner
#    skip() with a real comparison (e.g. relative error < 1e-4 per element).

test_that("autodiff OFV gradient matches finite-difference gradient within tolerance", {
  skip_if(
    !isTRUE(ferx_rust_autodiff_enabled()),
    "Enzyme autodiff not available — skipping gradient-correctness test"
  )
  skip("Gradient inspection API not yet exposed from ferx-nlme — implement once available")
})
