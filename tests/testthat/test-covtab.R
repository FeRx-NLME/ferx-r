# Covariate table (`fit$covtab`) — surfaced from a `[covariates]` block.
# See ferx-core #182/#184 for the underlying DSL + table.

test_that("covtab is a data frame: ID/TIME/EVID + declared covariates, one row per input record", {
  ex  <- ferx_example("two_cpt_oral_cov")
  fit <- ferx_fit(
    ex$model, ex$data,
    method     = "focei",
    verbose    = FALSE,
    covariance = FALSE,
    settings   = list(maxiter = 5L)
  )

  expect_s3_class(fit$covtab, "data.frame")
  # Leading metadata columns then the declared covariates in declaration order.
  expect_identical(names(fit$covtab)[1:3], c("ID", "TIME", "EVID"))
  expect_true(all(c("WT", "CRCL") %in% names(fit$covtab)))

  # One row per input dataset record (incl. dose / EVID rows), unlike sdtab.
  n_input_rows <- nrow(utils::read.csv(ex$data))
  expect_equal(nrow(fit$covtab), n_input_rows)
  expect_gt(nrow(fit$covtab), nrow(fit$sdtab)) # sdtab is observation rows only

  # Column types: ID character, TIME/covariates numeric, EVID integer.
  expect_type(fit$covtab$ID, "character")
  expect_type(fit$covtab$TIME, "double")
  expect_type(fit$covtab$WT, "double")
  expect_type(fit$covtab$EVID, "integer")

  # Dose rows are present (EVID == 1).
  expect_true(any(fit$covtab$EVID == 1L))
})

test_that("covtab is NULL when the model declares no [covariates] block", {
  fit <- warfarin_fit() # warfarin example has no [covariates] block
  expect_null(fit$covtab)
})
