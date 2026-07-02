
# ---- header from test-fit-async.R ----
# Tests for ferx_fit_async(), ferx_collect(), and their helpers.
# The async function is exercised via .ferx_print_trace_tail() against a
# synthetic trace CSV, plus a smoke test on the bundled warfarin example.
# We do not test the non-interactive fallback path against a live fit (which
# would just call ferx_fit()); the unit tests below cover the helper logic,
# and the smoke test confirms a real async run returns a valid ferx_fit object.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.print_trace_tail    <- getFromNamespace(".ferx_print_trace_tail",    "ferx")
.find_trace_lines    <- getFromNamespace(".ferx_find_trace_from_lines", "ferx")
.find_trace_bg       <- getFromNamespace(".ferx_find_trace_from_bg",  "ferx")
.async_callr         <- getFromNamespace(".ferx_async_callr",         "ferx")
.collect_callr       <- getFromNamespace(".ferx_collect_callr",       "ferx")














# -- callr sidecar tests ------------------------------------------------------






test_that("ferx_fit_async rejects invalid tail_n / poll_interval", {
  fn <- ferx::ferx_fit_async
  mockery::stub(fn, "interactive", function() TRUE)
  mockery::stub(fn, "requireNamespace", function(...) TRUE)
  mockery::stub(fn, "rstudioapi::isAvailable", function() FALSE)
  expect_error(fn("m", "d", tail_n = 0L), "tail_n")
  expect_error(fn("m", "d", tail_n = c(1L, 2L)), "tail_n")
  expect_error(fn("m", "d", poll_interval = 0), "poll_interval")
  expect_error(fn("m", "d", poll_interval = -1), "poll_interval")
})
test_that("ferx_fit_async returns ferx_job handle and ferx_collect retrieves result", {
  skip_if_not_installed("callr")
  fn_async    <- ferx::ferx_fit_async
  fn_collect  <- ferx::ferx_collect

  mockery::stub(fn_async, "interactive", function() TRUE)
  # Force callr backend regardless of IDE
  mockery::stub(fn_async, "requireNamespace",
                function(pkg, ...) if (pkg == "rstudioapi") FALSE else TRUE)

  ex <- ferx::ferx_example("warfarin")
  handle <- fn_async(ex$model, ex$data,
                     method = "focei", verbose = FALSE, covariance = FALSE,
                     settings = list(maxiter = 10L),
                     poll_interval = 0.2)

  expect_s3_class(handle, "ferx_job")
  expect_equal(handle$backend, "callr")

  fit <- fn_collect(handle)
  expect_s3_class(fit, "ferx_fit")
  expect_true(is.numeric(fit$ofv))
  expect_null(fit$trace_path)
})
test_that(".ferx_async_callr sets sidecar_path on handle", {
  skip_if_not_installed("callr")
  ex <- ferx::ferx_example("warfarin")
  h  <- .async_callr(ex$model, ex$data,
                     list(settings = list(maxiter = 5L)),
                     "warfarin", 6L, 0.5, FALSE)
  on.exit({
    if (h$bg$is_alive()) h$bg$kill()
    if (file.exists(h$sidecar_path)) unlink(h$sidecar_path)
  }, add = TRUE)
  expect_s3_class(h, "ferx_job")
  expect_equal(h$backend, "callr")
  expect_type(h$sidecar_path, "character")
  expect_true(nzchar(h$sidecar_path))
})
test_that("callr sidecar file is written and contains a valid trace path", {
  skip_if_not_installed("callr")
  ex <- ferx::ferx_example("warfarin")
  h  <- .async_callr(ex$model, ex$data,
                     list(settings = list(maxiter = 10L)),
                     "warfarin", 6L, 0.5, FALSE)
  on.exit({
    if (h$bg$is_alive()) h$bg$kill()
    if (file.exists(h$sidecar_path)) unlink(h$sidecar_path)
  }, add = TRUE)

  # Wait for the outer wrapper to finish (it finishes when the inner fit does).
  for (i in seq_len(60)) {
    if (!h$bg$is_alive()) break
    Sys.sleep(0.5)
  }
  expect_false(h$bg$is_alive(), label = "outer process should have finished")
  expect_true(file.exists(h$sidecar_path), label = "sidecar file should exist")

  trace_path <- trimws(readLines(h$sidecar_path, warn = FALSE)[1L])
  expect_true(nzchar(trace_path), label = "sidecar should contain a non-empty path")
  expect_true(file.exists(trace_path), label = "trace CSV should exist on disk")
})

# ---- header from test-fit-async-rstudio.R ----
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





# ---------------------------------------------------------------------------
# print.ferx_job
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
})
test_that(".ferx_async_rstudio cleans up and rethrows when jobRunScript fails", {
  skip_if_not_installed("mockery")
  async <- .async_rstudio
  mockery::stub(async, "rstudioapi::jobRunScript", function(...) stop("no IDE here"))
  expect_error(
    async("m.ferx", "d.csv", list(), "m.ferx", 6L, 0.5, FALSE),
    "no IDE here"
  )
})
