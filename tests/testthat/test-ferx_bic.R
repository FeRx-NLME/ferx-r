test_that("a fit carries the BIC tally and the three derived variants", {
  fit <- warfarin_fit()

  expect_type(fit$bic_inputs, "list")
  expect_setequal(
    names(fit$bic_inputs),
    c("n_obs", "theta_random", "theta_fixed", "omega", "kappa", "sigma",
      "sigma_random")
  )
  # The five parameter counts sum to n_parameters; n_obs is the record count
  # the observation-level penalty uses, not one of them.
  counts <- unlist(fit$bic_inputs[c("theta_random", "theta_fixed", "omega",
                                    "kappa", "sigma")])
  expect_equal(sum(counts), fit$n_parameters)
  expect_gt(fit$bic_inputs$n_obs, 0)
  expect_type(fit$bic_inputs$sigma_random, "logical")

  for (nm in c("bic_mixed", "bic_iiv", "bic_random")) {
    expect_true(is.numeric(fit[[nm]]) && length(fit[[nm]]) == 1L, info = nm)
    expect_false(is.na(fit[[nm]]), info = nm)
  }
})

test_that("ferx_bic() returns each convention", {
  fit <- warfarin_fit()

  expect_identical(ferx_bic(fit, "fixed"), as.numeric(fit$bic))
  expect_identical(ferx_bic(fit, "mixed"), as.numeric(fit$bic_mixed))
  expect_identical(ferx_bic(fit, "iiv"), as.numeric(fit$bic_iiv))
  expect_identical(ferx_bic(fit, "random"), as.numeric(fit$bic_random))
  # "mixed" is the default - the Delattre BIC Pharmpy ranks on.
  expect_identical(ferx_bic(fit), ferx_bic(fit, "mixed"))
})

test_that("the R penalty arithmetic agrees with the engine's", {
  fit <- warfarin_fit()

  # `ferx_bic(fit, "fixed")` short-circuits to the engine's own `bic`; the
  # recomputation from `bic_inputs` has to land on the same number, which
  # anchors the ln(n_obs) path against ferx-core.
  expect_equal(ferx:::.ferx_bic_variant(fit, "fixed"), as.numeric(fit$bic))
  # ... and these anchor the ln(n_subjects) paths against the values
  # ferx_core::bic() shipped on the fit.
  for (type in c("mixed", "iiv", "random")) {
    expect_equal(
      ferx:::.ferx_bic_variant(fit, type),
      as.numeric(fit[[paste0("bic_", type)]]),
      info = type
    )
  }
})

test_that("the penalties order as their sample sizes do", {
  fit <- warfarin_fit()
  skip_if_not(fit$n_subjects < fit$bic_inputs$n_obs)
  skip_if_not(fit$n_parameters > 0L)

  # Same count, smaller log: "random" penalises every parameter on
  # ln(n_subjects) where "fixed" uses ln(n_obs).
  expect_lt(ferx_bic(fit, "random"), ferx_bic(fit, "fixed"))
  # "iiv" penalises only the free omega elements, a subset of all parameters.
  expect_lte(ferx_bic(fit, "iiv"), ferx_bic(fit, "random"))
  expect_equal(ferx_bic(fit, "iiv") - fit$ofv,
               fit$bic_inputs$omega * log(fit$n_subjects))
})

test_that("a tally that cannot support the penalty gives NaN", {
  # What a .fitrx bundle saved before ferx-core #1177 loads as: an all-zero
  # tally next to a non-zero n_parameters.
  stale <- structure(
    list(
      ofv = 100, bic = 123, n_parameters = 4L, n_subjects = 32L,
      bic_inputs = list(n_obs = 0L, theta_random = 0L, theta_fixed = 0L,
                        omega = 0L, kappa = 0L, sigma = 0L,
                        sigma_random = FALSE)
    ),
    class = "ferx_fit"
  )
  expect_true(is.nan(ferx_bic(stale, "mixed")))
  expect_true(is.nan(ferx_bic(stale, "iiv")))
  expect_true(is.nan(ferx_bic(stale, "random")))
  # "fixed" is the engine's own value and stays usable.
  expect_identical(ferx_bic(stale, "fixed"), 123)

  no_tally <- structure(
    list(ofv = 100, bic = 123, n_parameters = 4L, n_subjects = 32L),
    class = "ferx_fit"
  )
  expect_true(is.nan(ferx_bic(no_tally, "mixed")))
})

test_that("a fit with no free parameter has one BIC", {
  # Every variant is the OFV, whatever n_subjects / n_obs say, so the four
  # agree and none needs a log.
  empty <- structure(
    list(
      ofv = 42, bic = 42, n_parameters = 0L, n_subjects = 0L,
      bic_inputs = list(n_obs = 0L, theta_random = 0L, theta_fixed = 0L,
                        omega = 0L, kappa = 0L, sigma = 0L,
                        sigma_random = FALSE)
    ),
    class = "ferx_fit"
  )
  for (type in c("mixed", "iiv", "random", "fixed")) {
    expect_equal(ferx_bic(empty, type), 42, info = type)
  }
})

test_that("ferx_bic() validates its arguments", {
  expect_error(ferx_bic(list(ofv = 1)), "must be a ferx_fit")
  expect_error(ferx_bic(warfarin_fit(), "delattre"), "should be one of")
})

test_that("the tally and the variants survive a .fitrx round trip", {
  fit <- warfarin_fit()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)
  ferx_save_fit(fit, tmp)
  back <- ferx_load_fit(tmp)

  expect_equal(back$bic_inputs, fit$bic_inputs)
  # The loaded fit carries no engine-computed variants; they are recomputed on
  # load from the tally, and must land on the same numbers.
  for (type in c("mixed", "iiv", "random", "fixed")) {
    expect_equal(ferx_bic(back, type), ferx_bic(fit, type), info = type)
  }
})
