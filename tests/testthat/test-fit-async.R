# Tests for ferx_fit_async() and its tail-print helper. The async function
# itself is exercised via .ferx_print_trace_tail() against a synthetic trace
# CSV, plus a smoke test on the bundled warfarin example. We do not test the
# non-interactive fallback path against a live fit (which would just call
# ferx_fit()); the unit tests below cover the helper logic, and the smoke
# test confirms a real async run returns a valid ferx_fit object.

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
  out <- capture.output(state <- ferx:::.ferx_print_trace_tail(tmp, 6L, state))
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
  # No new rows: helper must return state unchanged and produce no output
  out <- capture.output(state <- ferx:::.ferx_print_trace_tail(tmp, 6L, state))
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
  out <- capture.output(state <- ferx:::.ferx_print_trace_tail(tmp, 6L, state))
  # Header + 6 rows = 7 lines max
  iter_lines <- out[grepl("^\\s+\\d+\\s", out)]
  expect_lte(length(iter_lines), 6L)
  expect_equal(state$last_iter, n)
})

test_that(".ferx_find_trace_from_bg parses ferx-core's trace path message", {
  # Build a fake bg-like object whose read_output_lines / read_error_lines
  # return what callr would have captured. The trace announcement is on
  # stderr in ferx-core (eprintln!), so we put it in error_lines.
  fake_bg <- structure(
    list(
      stdout = c("some other line"),
      stderr = c(
        "Mu-referencing detected for: ETA_CL, ETA_KA",
        "[ferx] optimizer trace -> /tmp/ferx_trace_42_1700000000.csv",
        "another line"
      )
    ),
    class = "fake_bg"
  )
  fake_bg$read_output_lines <- function() fake_bg$stdout
  fake_bg$read_error_lines  <- function() fake_bg$stderr
  path <- ferx:::.ferx_find_trace_from_bg(fake_bg)
  expect_equal(path, "/tmp/ferx_trace_42_1700000000.csv")
})

test_that(".ferx_find_trace_from_bg returns NULL when no trace line present", {
  fake_bg <- structure(list(), class = "fake_bg")
  fake_bg$read_output_lines <- function() c("foo", "bar")
  fake_bg$read_error_lines  <- function() character(0)
  expect_null(ferx:::.ferx_find_trace_from_bg(fake_bg))
})

test_that("ferx_fit_async rejects invalid tail_n / poll_interval", {
  # The interactive() check returns early in non-interactive mode (Rscript /
  # R CMD check). To exercise the argument validators we stub interactive()
  # to TRUE via mockery (already in Suggests).
  fn <- ferx::ferx_fit_async
  mockery::stub(fn, "interactive", function() TRUE)
  expect_error(fn("m", "d", tail_n = 0L), "tail_n")
  expect_error(fn("m", "d", tail_n = c(1L, 2L)), "tail_n")
  expect_error(fn("m", "d", poll_interval = 0), "poll_interval")
  expect_error(fn("m", "d", poll_interval = -1), "poll_interval")
})

test_that("ferx_fit_async returns a valid ferx_fit object end-to-end", {
  skip_if_not_installed("callr")
  # The async path requires interactive() == TRUE; we mock it so the test
  # exercises the real code path even under Rscript / R CMD check.
  fn <- ferx::ferx_fit_async
  mockery::stub(fn, "interactive", function() TRUE)

  ex <- ferx_example("warfarin")
  fit <- fn(ex$model, ex$data,
            method = "focei", verbose = FALSE, covariance = FALSE,
            settings = list(maxiter = 10L),
            poll_interval = 0.2)

  expect_s3_class(fit, "ferx_fit")
  expect_true(is.numeric(fit$ofv))
  # We did not request the trace, so it should have been cleaned up.
  expect_null(fit$trace_path)
})
