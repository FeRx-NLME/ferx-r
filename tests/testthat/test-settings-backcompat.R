# Back-compat shim: ferx-core #422 renamed the importance-sampling settings
# is_* -> imp_*. ferx_fit() accepts the legacy spellings, translates them, and
# warns. See .ferx_migrate_legacy_is_settings() in R/fit.R.

migrate   <- getFromNamespace(".ferx_migrate_legacy_is_settings", "ferx")
to_strings <- getFromNamespace(".ferx_settings_to_strings", "ferx")

test_that("legacy is_* IS settings translate to imp_* with a deprecation warning", {
  expect_warning(
    out <- migrate(list(
      is_samples           = 2000L,
      is_proposal_df       = 4,
      is_seed              = 1L,
      is_low_ess_threshold = 0.1
    )),
    "deprecated"
  )
  expect_named(out, c("imp_samples", "imp_proposal_df", "imp_seed",
                      "imp_low_ess_threshold"))
  expect_equal(out$imp_samples, 2000L)
  expect_equal(out$imp_seed, 1L)
})

test_that("the modern imp_* spelling wins when both are supplied", {
  expect_warning(
    out <- migrate(list(is_samples = 1L, imp_samples = 2L)),
    "deprecated"
  )
  expect_named(out, "imp_samples")
  expect_equal(out$imp_samples, 2L)
})

test_that("non-legacy settings pass through unchanged and silently", {
  expect_silent(out <- migrate(list(optimizer = "bobyqa", n_starts = 4L)))
  expect_identical(out, list(optimizer = "bobyqa", n_starts = 4L))
})

test_that("NULL and empty settings are returned as-is", {
  expect_identical(migrate(NULL), NULL)
  expect_identical(migrate(list()), list())
})

test_that("the shim runs end-to-end through .ferx_settings_to_strings", {
  expect_warning(parts <- to_strings(list(is_seed = 7L)), "deprecated")
  expect_identical(parts$keys, "imp_seed")
  expect_identical(parts$values, "7")
})
