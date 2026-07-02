
# ---- header from test-runlog.R ----
# test-runlog.R
# Tests for ferx_runlog() -- shape, content, and graceful degradation.
# Uses warfarin_fit() (no covariance) and warfarin_fit_cov() (with covariance)
# from helper-warfarin-fit.R, which cap iterations at 30 for CI speed.

# -- Tier 1: basic contract ---------------------------------------------------




# -- Tier 2: header fields ----------------------------------------------------



# -- Tier 2: data summary -----------------------------------------------------







# -- Tier 2: parameter table --------------------------------------------------







# -- Tier 2: covariance step --------------------------------------------------




# -- Tier 2: estimation settings ----------------------------------------------






# -- Tier 2: diagnostics ------------------------------------------------------






# -- Tier 2: gradient ---------------------------------------------------------







# -- Tier 2: runtime ----------------------------------------------------------


# -- Tier 3: graceful degradation ---------------------------------------------



# -- Tier 3: NA-guard regression tests ----------------------------------------
# extendr maps Option<String>::None to NA_character_, not NULL.
# Verify that ferx_runlog() never emits NA tokens for these optional fields,
# even when they arrive as NA and another field (covariate_names) triggers
# has_settings = TRUE.





# -- Tier 2: iteration history -----------------------------------------------





# Local alias avoids the ::: operator which lintr flags.
.iter_table <- getFromNamespace(".runlog_iter_table", "ferx")













test_that("ferx_runlog_iters returns invisible character string", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  out <- ferx_runlog_iters(path, verbose = FALSE)
  expect_type(out, "character")
  expect_length(out, 1L)
  expect_match(out, "ITERATION HISTORY", fixed = TRUE)
})
test_that("ferx_runlog_iters verbose=TRUE prints and returns invisibly", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  expect_output(ferx_runlog_iters(path, verbose = TRUE), regexp = "ITERATION HISTORY")
})
test_that("ferx_runlog_iters errors on fit without trace", {
  fit <- warfarin_fit()
  fit$trace_path <- NULL
  expect_error(ferx_runlog_iters(fit), regexp = "No trace file found")
})
test_that("ferx_runlog_iters errors on non-ferx_fit non-string input", {
  expect_error(ferx_runlog_iters(42L), regexp = "ferx_fit object or a path")
})
test_that("ferx_runlog_iters output is pure ASCII", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  out <- ferx_runlog_iters(path, verbose = FALSE)
  expect_false(nchar(out, type = "bytes") != nchar(out, type = "chars"))
})
