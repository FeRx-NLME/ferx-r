# Tests for the name-based parameter accessors ferx_coef() / ferx_se() (#299).
# Driven with a crafted fit - no live model fit needed. make_fake_fit() comes
# from helper-trace.R; the `estimates` table is the derived field the accessors
# read, so build it the way ferx_fit() does.

.compute_est <- getFromNamespace(".ferx_compute_estimates", "ferx")

fit_with_estimates <- function(...) {
  # modifyList() so an override of NULL *removes* the field (the "no covariance
  # step" case); passing it through make_fake_fit()'s `...` twice would keep the
  # first value instead.
  args <- modifyList(
    list(
      theta       = c(TVCL = 1.0, TVV = 10.0),
      se_theta    = c(0.1, 0.5),
      omega       = 0.09,
      eta_names   = "ETA_CL",
      se_omega    = 0.01,
      sigma       = 0.05,
      sigma_names = "EPS_PROP",
      sigma_types = "proportional",
      se_sigma    = 0.005
    ),
    list(...)
  )
  fit <- do.call(make_fake_fit, args)
  fit$estimates <- .compute_est(fit)
  fit
}

test_that("ferx_coef() returns every parameter when `param` is NULL", {
  fit <- fit_with_estimates()
  co  <- ferx_coef(fit)
  expect_identical(names(co), c("TVCL", "TVV", "ETA_CL", "EPS_PROP"))
  expect_equal(unname(co), c(1.0, 10.0, 0.09, 0.05))
})

test_that("ferx_coef() pulls named parameters in the requested order", {
  fit <- fit_with_estimates()
  expect_equal(ferx_coef(fit, "TVCL"), c(TVCL = 1.0))
  expect_equal(ferx_coef(fit, c("TVV", "TVCL")), c(TVV = 10.0, TVCL = 1.0))
  expect_equal(ferx_coef(fit, "ETA_CL"), c(ETA_CL = 0.09))
  expect_equal(ferx_coef(fit, "EPS_PROP"), c(EPS_PROP = 0.05))
})

test_that("ferx_se() pulls the standard errors", {
  fit <- fit_with_estimates()
  expect_equal(ferx_se(fit, c("TVCL", "TVV")), c(TVCL = 0.1, TVV = 0.5))
  expect_equal(ferx_se(fit, "ETA_CL"), c(ETA_CL = 0.01))
})

test_that("an unknown parameter errors instead of returning NA", {
  fit <- fit_with_estimates()
  expect_error(ferx_coef(fit, "TVET50"), "unknown parameter")
  expect_error(ferx_se(fit, "TVET50"), "unknown parameter")
  # Plural form names every offender.
  expect_error(ferx_coef(fit, c("TVCL", "NOPE1", "NOPE2")), "NOPE1.*NOPE2")
})

test_that("a near-miss name is answered with a suggestion", {
  fit <- fit_with_estimates()
  expect_error(ferx_coef(fit, "TVC"), "Did you mean.*TVCL")
  expect_error(ferx_coef(fit, "eta_cl"), "Did you mean.*ETA_CL")
  # Nothing close -> the available names instead of a misleading guess.
  expect_error(ferx_coef(fit, "ZZZZZZZZ"), "Available: ")
})

test_that("ferx_se() warns when the fit carries no standard errors", {
  fit <- fit_with_estimates(se_theta = NULL, se_omega = NULL, se_sigma = NULL)
  expect_warning(se <- ferx_se(fit, "TVCL"), "no standard errors")
  expect_true(is.na(se[["TVCL"]]))
  # ferx_coef() on the same fit is silent - the estimates are there.
  expect_silent(ferx_coef(fit, "TVCL"))
})

test_that("the accessors reject a non-fit and a fit with no estimates table", {
  expect_error(ferx_coef(list(theta = c(TVCL = 1))), "needs a `ferx_fit` object")
  expect_error(ferx_se("not a fit"), "needs a `ferx_fit` object")
  expect_error(ferx_coef(make_fake_fit()), "no `estimates` table")
})

test_that("`param` must be character", {
  fit <- fit_with_estimates()
  expect_error(ferx_coef(fit, 1L), "must be a character vector")
})
