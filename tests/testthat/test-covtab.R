# Covariate table (`fit$covtab`) + types (`fit$covariate_types`), surfaced from
# a `[covariates]` block. See ferx-core #182/#184 for the underlying DSL + table.

# Cached fit on the covariate example so the multiple tests below share one
# estimation (mirrors helper-warfarin-fit.R's memoised pattern).
cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("two_cpt_oral_cov")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "focei", verbose = FALSE,
        covariance = FALSE, settings = list(maxiter = 5L)
      )
    }
    fit
  }
})

test_that("covtab is a data frame: ID/TIME/EVID + declared covariates in order, one row per input record", {
  fit <- cov_fit()
  ex  <- ferx_example("two_cpt_oral_cov")

  expect_s3_class(fit$covtab, "data.frame")
  # Exact columns + order: ID, TIME, EVID, then declared covariates in
  # declaration order.
  expect_identical(names(fit$covtab), c("ID", "TIME", "EVID", "WT", "CRCL"))

  # One row per input dataset record (incl. dose / EVID rows), unlike sdtab.
  expect_equal(nrow(fit$covtab), nrow(utils::read.csv(ex$data)))
  expect_gt(nrow(fit$covtab), nrow(fit$sdtab)) # sdtab is observation rows only

  # Column types: ID character, TIME/covariates numeric, EVID integer.
  expect_type(fit$covtab$ID, "character")
  expect_type(fit$covtab$TIME, "double")
  expect_type(fit$covtab$WT, "double")
  expect_type(fit$covtab$EVID, "integer")

  # Dose rows are present (EVID == 1).
  expect_true(any(fit$covtab$EVID == 1L))
})

test_that("covariate_types carries the declared continuous/categorical tag", {
  fit <- cov_fit()
  expect_type(fit$covariate_types, "character")
  expect_identical(
    fit$covariate_types,
    c(WT = "continuous", CRCL = "continuous")
  )
})

test_that("covtab and covariate_types survive a .fitrx save/load round-trip", {
  fit <- cov_fit()
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)

  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)

  expect_s3_class(fit2$covtab, "data.frame")
  expect_identical(names(fit2$covtab), names(fit$covtab))
  expect_equal(nrow(fit2$covtab), nrow(fit$covtab))
  expect_type(fit2$covtab$ID, "character")
  expect_identical(fit2$covariate_types, fit$covariate_types)
})

test_that("covtab and covariate_types are NULL when the model declares no [covariates] block", {
  fit <- warfarin_fit() # warfarin example has no [covariates] block
  expect_null(fit$covtab)
  expect_null(fit$covariate_types)
})
