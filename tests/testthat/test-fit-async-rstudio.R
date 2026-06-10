# Tests for the RStudio-jobs async backend in fit.R. RStudio is never present
# in CI, so rstudioapi is mocked with mockery::stub and the surrounding R logic
# (script generation, result collection, handle printing) is exercised for real.
#
# make_fake_fit() and write_fake_trace() come from helper-trace.R.

.async_rstudio   <- getFromNamespace(".ferx_async_rstudio",   "ferx")
.collect_rstudio <- getFromNamespace(".ferx_collect_rstudio", "ferx")
.print_job       <- getFromNamespace("print.ferx_job",        "ferx")

# ---------------------------------------------------------------------------
# .ferx_async_rstudio
# ---------------------------------------------------------------------------
test_that(".ferx_async_rstudio writes a runnable job script and returns a handle", {
  skip_if_not_installed("mockery")
  async <- .async_rstudio
  mockery::stub(async, "rstudioapi::jobRunScript", function(...) "job_123")

  h <- async("warfarin.ferx", "warfarin.csv", list(method = "focei"),
             "warfarin.ferx", 6L, 0.5, FALSE)
  on.exit(unlink(c(h$script_path, h$rds_path, h$sidecar_path)), add = TRUE)

  expect_s3_class(h, "ferx_job")
  expect_identical(h$backend, "rstudio")
  expect_identical(h$job_id, "job_123")
  expect_true(file.exists(h$script_path))

  # The generated script must drive a real ferx_fit in a background process.
  script <- readLines(h$script_path)
  args_line <- script[grepl("saved <- readRDS", script, fixed = TRUE)][1L]
  args_path <- sub('.*readRDS\\("([^"]+)"\\).*', "\\1", args_line)
  if (!is.na(args_line) && !identical(args_path, args_line)) on.exit(unlink(args_path), add = TRUE)
  expect_true(any(grepl("ferx::ferx_fit", script, fixed = TRUE)))
  expect_true(any(grepl("optimizer_trace <- TRUE", script, fixed = TRUE)))

test_that(".ferx_async_rstudio cleans up and rethrows when jobRunScript fails", {
  skip_if_not_installed("mockery")
  async <- .async_rstudio
  mockery::stub(async, "rstudioapi::jobRunScript", function(...) stop("no IDE here"))
  expect_error(
    async("m.ferx", "d.csv", list(), "m.ferx", 6L, 0.5, FALSE),
    "no IDE here"
  )
})

# ---------------------------------------------------------------------------
# .ferx_collect_rstudio
# ---------------------------------------------------------------------------
make_rstudio_handle <- function(rds_path = tempfile(fileext = ".rds"),
                                sidecar_path = tempfile(fileext = ".txt")) {
  structure(
    list(backend = "rstudio", job_id = "job_1", rds_path = rds_path,
         sidecar_path = sidecar_path, script_path = tempfile(fileext = ".R"),
         tail_n = 6L, poll_interval = 0.01, user_wanted_trace = FALSE),
    class = "ferx_job"
  )
}

test_that(".ferx_collect_rstudio returns the saved fit when the job succeeds", {
  skip_if_not_installed("mockery")
  h <- make_rstudio_handle()
  saveRDS(make_fake_fit(), h$rds_path)
  collect <- .collect_rstudio
  mockery::stub(collect, "rstudioapi::jobGetState", function(id) "succeeded")

  fit <- collect(h, verbose = FALSE)
  expect_s3_class(fit, "ferx_fit")
  expect_false(file.exists(h$rds_path))   # result rds cleaned up
})

test_that(".ferx_collect_rstudio streams the trace tail while the job is running (verbose)", {
  skip_if_not_installed("mockery")
  h <- make_rstudio_handle()
  saveRDS(make_fake_fit(), h$rds_path)
  writeLines(write_fake_trace(n = 4L), h$sidecar_path)  # sidecar -> real trace CSV
  collect <- .collect_rstudio
  calls <- 0L
  mockery::stub(collect, "rstudioapi::jobGetState", function(id) {
    calls <<- calls + 1L
    if (calls < 2L) "running" else "succeeded"
  })

  out <- capture.output(fit <- collect(h, verbose = TRUE))
  expect_s3_class(fit, "ferx_fit")
  expect_true(calls >= 2L)               # loop body ran at least once
})

test_that(".ferx_collect_rstudio errors when the job fails", {
  skip_if_not_installed("mockery")
  h <- make_rstudio_handle()
  collect <- .collect_rstudio
  mockery::stub(collect, "rstudioapi::jobGetState", function(id) "failed")
  expect_error(collect(h, verbose = FALSE), "Background fit failed")
})

test_that("ferx_collect dispatches to the rstudio backend", {
  skip_if_not_installed("mockery")
  fn_collect <- ferx::ferx_collect
  mockery::stub(fn_collect, ".ferx_collect_rstudio",
                function(handle, verbose) make_fake_fit(trace_path = NULL))
  handle <- structure(list(backend = "rstudio", user_wanted_trace = FALSE),
                      class = "ferx_job")
  fit <- fn_collect(handle)
  expect_s3_class(fit, "ferx_fit")
})

# ---------------------------------------------------------------------------
# print.ferx_job
# ---------------------------------------------------------------------------
test_that("print.ferx_job reports rstudio status and any trace progress", {
  skip_if_not_installed("mockery")
  sidecar <- tempfile(fileext = ".txt")
  writeLines(write_fake_trace(n = 3L), sidecar)
  handle <- structure(
    list(backend = "rstudio", job_id = "job_1", model_label = "warfarin.ferx",
         sidecar_path = sidecar),
    class = "ferx_job"
  )
  pj <- .print_job
  mockery::stub(pj, "rstudioapi::jobGetState", function(id) "running")

  out <- capture.output(res <- pj(handle))
  expect_identical(res, handle)
  expect_true(any(grepl("<ferx_job>", out)))
  expect_true(any(grepl("rstudio", out)))
  expect_true(any(grepl("iteration", out)))   # trace table footer
})

test_that("print.ferx_job reports callr status and skips trace when none present", {
  fake_bg <- list(is_alive = function() FALSE)
  handle <- structure(
    list(backend = "callr", bg = fake_bg, model_label = "m.ferx",
         sidecar_path = tempfile(fileext = ".txt")),  # nonexistent -> no trace
    class = "ferx_job"
  )
  out <- capture.output(res <- print(handle))
  expect_identical(res, handle)
  expect_true(any(grepl("done", out)))
  expect_true(any(grepl("callr", out)))
})
