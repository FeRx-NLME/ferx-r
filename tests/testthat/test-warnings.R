test_that("fit carries a structured warnings data frame", {
  fit <- warfarin_fit_cov()
  expect_true(is.data.frame(fit$warnings_structured))
  expect_named(fit$warnings_structured, c("severity", "category", "message", "source_method"))
  # Severity values are always from the fixed vocabulary
  if (nrow(fit$warnings_structured) > 0L) {
    expect_true(all(fit$warnings_structured$severity %in%
      c("critical", "warning", "info")))
  }
})

test_that("ferx_warnings(as_df = TRUE) returns the underlying data frame", {
  fit <- warfarin_fit_cov()
  df <- ferx_warnings(fit, as_df = TRUE)
  expect_identical(df, fit$warnings_structured)
})

test_that("ferx_warnings() prints a grouped summary", {
  fit <- warfarin_fit_cov()
  out <- capture.output(ferx_warnings(fit))
  expect_true(any(grepl("ferx fit warnings", out)))
  # Footer tallies are always present
  expect_true(any(grepl("CRITICAL", out)))
  expect_true(any(grepl("WARNING", out)))
  expect_true(any(grepl("INFO", out)))
})

test_that("ferx_warnings() rejects non-fit input", {
  expect_error(ferx_warnings(list()), "ferx_fit")
})

test_that("ferx_warnings() falls back to flat warnings when structured is absent", {
  fake <- structure(
    list(
      model_name = "legacy",
      warnings = c("Outer optimization did not converge", "something else"),
      warnings_structured = NULL
    ),
    class = "ferx_fit"
  )
  df <- ferx_warnings(fake, as_df = TRUE)
  expect_equal(nrow(df), 2L)
  expect_true(all(df$severity == "warning"))
  expect_true(all(df$category == "general"))
})

test_that("structured warnings include R-side condition number flag", {
  fake <- structure(
    list(
      model_name = "m",
      condition_number = 5000,
      eta_normality = NULL
    ),
    class = "ferx_fit"
  )
  raw <- list(
    warnings_severity = character(0),
    warnings_category = character(0),
    warnings_message  = character(0)
  )
  df <- ferx:::.ferx_assemble_structured_warnings(raw, fake)
  expect_true(any(df$category == "condition_number"))
  expect_equal(df$severity[df$category == "condition_number"], "critical")
})
