# ferx-core #1620: when two covariate mu-referenced typical values share one
# estimated theta (THETA_WT on TVCL and TVV1 in two_cpt_oral_cov), SAEM and
# IMP/IMPMAP take one joint M-step for it and the fit carries an info note
# naming the etas and the shared theta. These tests pin the estimate and the
# note from R, so a pin bump or a glue change that drops either goes red here
# and not only in ferx-core (#422).

joint_mstep_rows <- function(fit) {
  df <- ferx_get_warnings(fit, as_df = TRUE)
  expect_true("message" %in% names(df))
  df[grepl("joint M-step", df$message, fixed = TRUE), , drop = FALSE]
}

test_that("SAEM moves THETA_WT off its pre-#1620 value on a shared theta", {
  skip_on_cran()
  fit <- two_cpt_oral_cov_saem_fit()
  # 0.0109 before #1620, 0.618 after; NONMEM SAEM gives 0.6618.
  expect_gt(fit$theta[["THETA_WT"]], 0.3)
})

test_that("the joint M-step note reaches R as one info row per method", {
  skip_on_cran()
  fits <- list(
    SAEM   = two_cpt_oral_cov_saem_fit(),
    IMPMAP = two_cpt_oral_cov_impmap_fit()
  )
  for (method in names(fits)) {
    rows <- joint_mstep_rows(fits[[method]])
    expect_equal(nrow(rows), 1L, label = paste(method, "joint M-step rows"))
    expect_identical(rows$severity, "info")
    expect_identical(rows$category, "mu_referencing")
    msg <- rows$message
    expect_true(startsWith(msg, paste0(method, ":")), label = msg)
    for (name in c("ETA_CL", "ETA_V1", "THETA_WT")) {
      expect_match(msg, name, fixed = TRUE)
    }
  }
})

test_that("no joint M-step note without a shared theta", {
  skip_on_cran()
  # warfarin_saem has no covariate mu-references, so no theta is shared.
  fit <- warfarin_saem_conddist_fit()
  expect_equal(nrow(joint_mstep_rows(fit)), 0L)
})

test_that("the legacy mu-ref detector ignores the joint M-step note", {
  skip_on_cran()
  fit <- two_cpt_oral_cov_saem_fit()
  expect_true(any(grepl("joint M-step", fit$warnings, fixed = TRUE)))
  expect_false(any(grepl(
    "joint M-step", .ferx_mu_ref_detection_names(fit$warnings), fixed = TRUE
  )))
})
