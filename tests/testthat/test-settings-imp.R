to_strings <- getFromNamespace(".ferx_settings_to_strings", "ferx")

test_that("importance-sampling settings use imp_* names", {
  parts <- to_strings(list(
    imp_samples = 2000L,
    imp_proposal_df = 4,
    imp_seed = 1L,
    imp_low_ess_threshold = 0.1
  ))

  expect_identical(parts$keys, c(
    "imp_samples",
    "imp_proposal_df",
    "imp_seed",
    "imp_low_ess_threshold"
  ))
  expect_identical(parts$values, c("2000", "4", "1", "0.10000000000000001"))
})

test_that("short-lived is_* settings are not translated by R", {
  parts <- to_strings(list(is_seed = 7L))

  expect_identical(parts$keys, "is_seed")
  expect_identical(parts$values, "7")
})

test_that("duplicated settings are rejected before serialization", {
  expect_error(
    to_strings(list(imp_samples = 1L, imp_samples = 2L)),
    "uniquely-named"
  )
})
