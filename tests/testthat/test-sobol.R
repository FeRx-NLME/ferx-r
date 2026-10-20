# Sobol quasi-random importance sampling (#574, ferx-core #1965). `impmap_sobol`
# used to be ignored, with a warning, under the default Student-t proposal; it
# now covers every proposal form, and `imp_sobol` is its `imp` counterpart.

sobol_warnings <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("impmap_sobol runs on the default Student-t proposal without being ignored", {
  ex <- ferx_example("warfarin")
  res <- sobol_warnings(ferx_fit(
    ex$model, ex$data, method = "impmap", verbose = FALSE,
    settings = list(impmap_sobol = TRUE, impmap_iterations = 3L,
                    impmap_averaging = 1L, impmap_seed = 42L)
  ))
  expect_s3_class(res$value, "ferx_fit")
  expect_true(is.finite(res$value$ofv))
  expect_false(any(grepl("sobol", res$warnings, ignore.case = TRUE)))
  expect_false(any(grepl("ignored", res$warnings, ignore.case = TRUE)))
})

test_that("impmap_sobol combines with the defensive mixture", {
  ex <- ferx_example("warfarin")
  res <- sobol_warnings(ferx_fit(
    ex$model, ex$data, method = "impmap", verbose = FALSE,
    settings = list(impmap_sobol = TRUE, imp_defensive_alpha = 0.1,
                    impmap_iterations = 3L, impmap_averaging = 1L,
                    impmap_seed = 42L)
  ))
  expect_true(is.finite(res$value$ofv))
  expect_false(any(grepl("sobol", res$warnings, ignore.case = TRUE)))
})

test_that("imp_sobol is accepted by the estimating imp method", {
  ex <- ferx_example("warfarin")
  res <- sobol_warnings(ferx_fit(
    ex$model, ex$data, method = "imp", verbose = FALSE,
    settings = list(imp_sobol = TRUE, imp_iterations = 3L,
                    imp_averaging = 1L, imp_seed = 42L)
  ))
  expect_s3_class(res$value, "ferx_fit")
  expect_true(is.finite(res$value$ofv))
  expect_false(any(grepl("imp_sobol|ignored", res$warnings, ignore.case = TRUE)))
})

test_that("impmap_sobol refuses a Student-t proposal wider than 1000 df", {
  ex <- ferx_example("warfarin")
  expect_error(
    suppressWarnings(ferx_fit(
      ex$model, ex$data, method = "impmap", verbose = FALSE,
      settings = list(impmap_sobol = TRUE, impmap_proposal_df = 5000,
                      impmap_iterations = 2L, impmap_averaging = 1L)
    )),
    regexp = "at most 1000 degrees of freedom"
  )
})
