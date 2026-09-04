# Synthetic fit carrying only what ferx_bic() reads, so the penalty arithmetic
# is tested without a live estimation.
bic_fit <- function(ofv = 100, n_parameters = 5L, n_subjects = 20L,
                    n_obs = 200L, theta_random = 2L, theta_fixed = 1L,
                    omega = 1L, kappa = 0L, sigma = 1L,
                    sigma_random = FALSE, bic_inputs = NULL) {
  inputs <- if (is.null(bic_inputs)) {
    list(n_obs = n_obs, theta_random = theta_random,
         theta_fixed = theta_fixed, omega = omega, kappa = kappa,
         sigma = sigma, sigma_random = sigma_random)
  } else {
    bic_inputs
  }
  structure(
    list(ofv = ofv, n_parameters = n_parameters, n_subjects = n_subjects,
         n_obs = n_obs, bic_inputs = inputs),
    class = "ferx_fit"
  )
}

test_that("the four variants apply the documented penalties", {
  fit <- bic_fit()
  # random class = theta_random + omega + kappa = 3; fixed class =
  # theta_fixed + sigma = 2 (sigma_random is FALSE).
  expect_equal(ferx_bic(fit, "mixed"),
               100 + 3 * log(20) + 2 * log(200))
  expect_equal(ferx_bic(fit, "fixed"), 100 + 5 * log(200))
  expect_equal(ferx_bic(fit, "random"), 100 + 5 * log(20))
  expect_equal(ferx_bic(fit, "iiv"), 100 + 1 * log(20))
})

test_that("mixed defaults and sigma_random moves sigma into the random class", {
  fit <- bic_fit()
  expect_equal(ferx_bic(fit), ferx_bic(fit, "mixed"))

  iiv_on_ruv <- bic_fit(sigma_random = TRUE)
  expect_equal(ferx_bic(iiv_on_ruv, "mixed"),
               100 + 4 * log(20) + 1 * log(200))
})

test_that("the fixed variant reproduces fit$bic on a real fit", {
  skip_on_cran()
  fit <- warfarin_fit()
  expect_equal(ferx_bic(fit, "fixed"), fit$bic)
  expect_false(is.na(ferx_bic(fit, "mixed")))
})

test_that("a tally that disagrees with n_parameters yields NA", {
  # 2+1+1+0+1 = 5 free, but the fit claims 6.
  expect_true(is.na(ferx_bic(bic_fit(n_parameters = 6L))))
})

test_that("a fit with no tally yields NA rather than a wrong penalty", {
  expect_true(is.na(ferx_bic(bic_fit(bic_inputs = list()))))
  no_tally <- bic_fit()
  no_tally$bic_inputs <- NULL
  expect_true(is.na(ferx_bic(no_tally)))
})

test_that("a count of zero never needs its logarithm", {
  # No free parameter: every variant is the OFV, whatever n_subjects /
  # n_obs say - including a bundle that recorded neither.
  none <- bic_fit(n_parameters = 0L, n_subjects = 0L, n_obs = 0L,
                  theta_random = 0L, theta_fixed = 0L, omega = 0L,
                  kappa = 0L, sigma = 0L)
  for (ty in c("mixed", "fixed", "iiv", "random")) {
    expect_equal(ferx_bic(none, ty), 100)
  }
})

test_that("a variant that needs a missing logarithm yields NA", {
  no_subjects <- bic_fit(n_subjects = 0L)
  expect_true(is.na(ferx_bic(no_subjects, "mixed")))
  expect_true(is.na(ferx_bic(no_subjects, "iiv")))
  expect_true(is.na(ferx_bic(no_subjects, "random")))
  # `fixed` only needs n_obs, which this fit has.
  expect_equal(ferx_bic(no_subjects, "fixed"), 100 + 5 * log(200))
})

test_that("input validation", {
  expect_error(ferx_bic(list(ofv = 1)), "ferx_fit")
  expect_error(ferx_bic(bic_fit(), "aic"), "'arg'")
})
