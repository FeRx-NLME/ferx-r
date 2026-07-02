
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






test_that(".ferx_print_trace_tail prints header and rows on first call", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(
    data.frame(
      iter = 1:3, method = "focei", phase = "", ofv = c(100, 90, 85),
      wall_ms = 0L, grad_norm = c(1.0, 0.5, 0.1), step_norm = NA,
      inner_iter_count = NA, optimizer = "slsqp",
      lm_lambda = NA, ofv_delta = NA, step_accepted = NA, cond_nll = NA,
      gamma = NA, mh_accept_rate = NA, n_ebe_unconverged = 0L,
      n_ebe_fallback = 0L
    ),
    tmp, row.names = FALSE
  )
  state <- list(last_iter = -1L, lines_printed = 0L)
  out <- capture.output(state <- .print_trace_tail(tmp, 6L, state))
  expect_true(any(grepl("iter", out)))
  expect_true(any(grepl("^\\s*1\\s", out)))
  expect_true(any(grepl("^\\s*3\\s", out)))
  expect_equal(state$last_iter, 3L)
})
test_that(".ferx_print_trace_tail only prints new rows on subsequent calls", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(
    data.frame(
      iter = 1:3, method = "focei", phase = "", ofv = c(100, 90, 85),
      wall_ms = 0L, grad_norm = c(1.0, 0.5, 0.1), step_norm = NA,
      inner_iter_count = NA, optimizer = "slsqp",
      lm_lambda = NA, ofv_delta = NA, step_accepted = NA, cond_nll = NA,
      gamma = NA, mh_accept_rate = NA, n_ebe_unconverged = 0L,
      n_ebe_fallback = 0L
    ),
    tmp, row.names = FALSE
  )
  state <- list(last_iter = 3L, lines_printed = 3L)
  out <- capture.output(state <- .print_trace_tail(tmp, 6L, state))
  expect_equal(state$last_iter, 3L)
  expect_length(out, 0L)
})
test_that(".ferx_print_trace_tail caps output at tail_n rows per call", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  n <- 20L
  utils::write.csv(
    data.frame(
      iter = seq_len(n), method = "focei", phase = "",
      ofv = seq.int(100, length.out = n, by = -1),
      wall_ms = 0L, grad_norm = NA, step_norm = NA,
      inner_iter_count = NA, optimizer = "slsqp",
      lm_lambda = NA, ofv_delta = NA, step_accepted = NA, cond_nll = NA,
      gamma = NA, mh_accept_rate = NA, n_ebe_unconverged = 0L,
      n_ebe_fallback = 0L
    ),
    tmp, row.names = FALSE
  )
  state <- list(last_iter = -1L, lines_printed = 0L)
  out <- capture.output(state <- .print_trace_tail(tmp, 6L, state))
  iter_lines <- out[grepl("^\\s+\\d+\\s", out)]
  expect_lte(length(iter_lines), 6L)
  expect_equal(state$last_iter, n)
})
test_that("ferx_collect rejects non-ferx_job input", {
  expect_error(ferx::ferx_collect(list()), "ferx_job")
  expect_error(ferx::ferx_collect("not a handle"), "ferx_job")
})
test_that("callr sidecar is removed after ferx_collect", {
  skip_if_not_installed("callr")
  ex <- ferx::ferx_example("warfarin")
  h  <- .async_callr(ex$model, ex$data,
                     list(settings = list(maxiter = 10L)),
                     "warfarin", 6L, 0.5, FALSE)

  fit <- .collect_callr(h, verbose = FALSE)
  expect_s3_class(fit, "ferx_fit")
  expect_false(file.exists(h$sidecar_path), label = "sidecar should be removed after collect")
})
test_that(".ferx_collect_callr polling loop runs when bg is alive at collect time", {
  skip_if_not_installed("callr")
  ex <- ferx::ferx_example("warfarin")
  # Use a longer fit so the outer wrapper is still alive when collect starts.
  h  <- .async_callr(ex$model, ex$data,
                     list(settings = list(maxiter = 200L)),
                     "warfarin", 6L, 0.5, FALSE)
  # Call collect immediately -- the outer wrapper should still be alive for
  # at least one poll iteration, exercising the while-loop body.
  fit <- .collect_callr(h, verbose = TRUE)
  expect_s3_class(fit, "ferx_fit")
  expect_false(file.exists(h$sidecar_path))
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
