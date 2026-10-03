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
  # 0.0109 before #1620, 0.618 after; NONMEM SAEM gives 0.6618. The upper
  # side catches a joint step that overshoots towards the 5.0 bound.
  expect_gt(fit$theta[["THETA_WT"]], 0.3)
  expect_lt(fit$theta[["THETA_WT"]], 1.2)
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
    expect_identical(rows$severity, "info",
                     label = paste(method, "severity"))
    expect_identical(rows$category, "mu_referencing",
                     label = paste(method, "category"))
    msg <- rows$message
    expect_true(startsWith(msg, paste0(method, ":")), label = msg)
    for (eta in c("ETA_CL", "ETA_V1")) {
      expect_match(msg, eta, fixed = TRUE, label = paste(method, "note"))
    }
    # THETA_CRCL is also in the note's `reads` lists; this pins which theta
    # the note says is shared.
    expect_match(msg, "share THETA_WT and", fixed = TRUE,
                 label = paste(method, "note"))
  }
})

test_that("no joint M-step note when the covariate mu-refs share no theta", {
  skip_on_cran()
  # Same model, V1 on its own exponent THETA_WTV: still two covariate
  # mu-references, but nothing shared between them.
  fit <- two_cpt_oral_cov_separate_saem_fit()
  expect_true("THETA_WTV" %in% names(fit$theta))
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
