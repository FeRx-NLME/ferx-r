# Tests for [derived] and [output] blocks in ferx model files.
#
# These are integration tests that require a real fit; they are skipped on
# CRAN. The `warfarin_derived_fit()` helper in helper-warfarin-fit.R caches
# the result so the fit runs at most once per test session.
#
# Coverage:
#   - [derived] columns (per-row and aggregate) appear in fit$sdtab
#   - [output] columns (individual PK parameters) appear in fit$sdtab
#   - ETA columns are absent from sdtab (regression: ferx-core#188)
#   - fit$covariate_names is populated when the model uses covariates
#   - Derived column values are numerically plausible

# -- [derived] columns present in sdtab ------------------------------------

test_that("[derived] per-row columns appear in sdtab", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # KE and T_HALF are per-row expressions declared in warfarin_derived.ferx
  expect_true("KE"     %in% names(fit$sdtab),
              label = "KE (per-row [derived]) must be in sdtab")
  expect_true("T_HALF" %in% names(fit$sdtab),
              label = "T_HALF (per-row [derived]) must be in sdtab")
  expect_true(all(is.finite(fit$sdtab$KE)),     label = "KE must be finite")
  expect_true(all(is.finite(fit$sdtab$T_HALF)), label = "T_HALF must be finite")
})

test_that("[derived] aggregate columns appear in sdtab and are subject-constant", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # CMAX and AUC_0_72 are aggregate [derived] columns
  expect_true("CMAX"     %in% names(fit$sdtab))
  expect_true("AUC_0_72" %in% names(fit$sdtab))
  # Aggregates are constant within a subject (repeat the same value)
  for (col in c("CMAX", "AUC_0_72")) {
    by_id <- tapply(fit$sdtab[[col]], fit$sdtab$ID, function(x) length(unique(x)))
    expect_true(all(by_id == 1L),
                label = sprintf("%s must be constant within each subject", col))
  }
})

test_that("[derived] TAFD/TAD auxiliary columns are present and positive", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_true("DAY"      %in% names(fit$sdtab))
  expect_true("TAU_TIME" %in% names(fit$sdtab))
  expect_true(all(fit$sdtab$DAY >= 1L))
})

# -- [output] columns present in sdtab ------------------------------------

test_that("[output] individual PK parameters appear in sdtab", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # warfarin_derived.ferx has [output] CL V KA
  for (col in c("CL", "V", "KA")) {
    expect_true(col %in% names(fit$sdtab),
                label = sprintf("[output] column %s must be in sdtab", col))
    expect_true(all(is.finite(fit$sdtab[[col]])),
                label = sprintf("[output] column %s must be finite", col))
    expect_true(all(fit$sdtab[[col]] > 0),
                label = sprintf("[output] column %s must be positive", col))
  }
})

test_that("[output] PK params are constant within subject (no random noise)", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  for (col in c("CL", "V", "KA")) {
    by_id <- tapply(fit$sdtab[[col]], fit$sdtab$ID, function(x) length(unique(x)))
    expect_true(all(by_id == 1L),
                label = sprintf("%s must be constant within each subject", col))
  }
})

# -- ETAs absent from sdtab (regression guard) ----------------------------

test_that("[derived] fit: sdtab does not contain ETA columns", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  eta_cols <- grep("^ETA_|^ETA[0-9]", names(fit$sdtab), value = TRUE)
  expect_length(eta_cols, 0L)
})

# -- ETA values are in ebe_etas, accessible for merging -------------------

test_that("fit$ebe_etas has ETA columns matching model declarations", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_s3_class(fit$ebe_etas, "data.frame")
  expect_true("ID" %in% names(fit$ebe_etas))
  # warfarin_derived declares ETA_CL, ETA_V, ETA_KA
  for (eta in c("ETA_CL", "ETA_V", "ETA_KA")) {
    expect_true(eta %in% names(fit$ebe_etas),
                label = sprintf("%s must be in ebe_etas", eta))
  }
  expect_equal(nrow(fit$ebe_etas), fit$n_subjects)
})

test_that("ebe_etas can be merged onto sdtab by ID for downstream analysis", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  subj <- fit$sdtab[!duplicated(fit$sdtab$ID), ]
  merged <- merge(subj, fit$ebe_etas, by = "ID")
  expect_equal(nrow(merged), fit$n_subjects)
  expect_true("ETA_CL" %in% names(merged))
})

# -- Numerical plausibility -----------------------------------------------

test_that("T_HALF and KE are numerically consistent: T_HALF == 0.693 / KE", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_equal(fit$sdtab$T_HALF, 0.6931472 / fit$sdtab$KE,
               tolerance = 1e-6,
               label = "T_HALF must equal 0.6931472 / KE within tolerance")
})

test_that("CL / V == KE within tolerance", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_equal(fit$sdtab$CL / fit$sdtab$V, fit$sdtab$KE,
               tolerance = 1e-6,
               label = "KE must equal CL / V")
})

test_that("CMAX is non-negative and >= all IPRED values per subject", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_true(all(fit$sdtab$CMAX >= 0))
  max_ipred <- tapply(fit$sdtab$IPRED, fit$sdtab$ID, max)
  cmax_val  <- tapply(fit$sdtab$CMAX, fit$sdtab$ID, unique)
  expect_true(all(unlist(cmax_val) >= unlist(max_ipred) - 1e-6),
              label = "CMAX must be >= max(IPRED) for every subject")
})
