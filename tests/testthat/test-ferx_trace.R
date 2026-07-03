
# ---- header from test-trace.R ----








# plot.ferx_fit — Tier 1, no Rust needed





test_that("ferx_trace(path) reads a trace CSV into a data frame", {
  path <- write_fake_trace(n = 4)
  on.exit(unlink(path))

  tr <- ferx_trace(path)
  expect_s3_class(tr, "data.frame")
  expect_equal(nrow(tr), 4L)
  expect_true(all(c("iter", "method", "ofv", "grad_norm") %in% names(tr)))
  expect_type(tr$grad_norm, "double")
})
test_that("ferx_trace(fit) reads the path stored on a ferx_fit object", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  fit <- structure(list(trace_path = path), class = "ferx_fit")

  tr <- ferx_trace(fit)
  expect_equal(nrow(tr), 5L)
})
test_that("ferx_trace() with no argument reads the last-run trace", {
  reset_ferx_state()
  on.exit(reset_ferx_state(), add = TRUE)

  path <- write_fake_trace(n = 3)
  on.exit(unlink(path), add = TRUE)
  st <- ferx_state()
  st$last_trace_path  <- path
  st$last_trace_time  <- as.POSIXct("2026-04-27 12:34:56", tz = "UTC")
  st$last_trace_model <- "/some/dir/warfarin.ferx"

  expect_message(tr <- ferx_trace(), regexp = "Last ferx_fit\\(\\) started at")
  expect_message(ferx_trace(), regexp = "warfarin\\.ferx")
  expect_equal(nrow(tr), 3L)
})
test_that("ferx_trace() errors when no last-run state is present", {
  reset_ferx_state()
  on.exit(reset_ferx_state(), add = TRUE)

  expect_error(ferx_trace(), regexp = "No trace from a previous")
})
test_that("ferx_trace() rejects a ferx_fit object without a trace_path", {
  fit <- structure(list(trace_path = NULL), class = "ferx_fit")
  expect_error(ferx_trace(fit), regexp = "No trace path found")
})
test_that("ferx_trace() rejects unsupported argument types", {
  expect_error(ferx_trace(42),    regexp = "must be a ferx_fit object")
  expect_error(ferx_trace(list()), regexp = "must be a ferx_fit object")
})
test_that("ferx_trace(path) errors when the file does not exist", {
  expect_error(
    ferx_trace(file.path(tempdir(), "definitely-not-a-real-trace.csv")),
    regexp = "Trace file not found"
  )
})
test_that("last-run state is per-session: starts NULL on a fresh load", {
  # We can't truly spawn a second R session here cheaply, but we can verify
  # the contract: the state lives in a package env and resets to NULL when
  # cleared. A second session would re-run the .ferx_state initializer in
  # zzz.R and start in this same NULL state.
  reset_ferx_state()
  expect_null(ferx_state()$last_trace_path)
  expect_null(ferx_state()$last_trace_time)
  expect_null(ferx_state()$last_trace_model)
})
