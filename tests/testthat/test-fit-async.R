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

test_that(".ferx_find_trace_from_lines parses ferx-core's trace path message", {
  lines <- c(
    "Mu-referencing detected for: ETA_CL, ETA_KA",
    "[ferx] optimizer trace -> /tmp/ferx_trace_42_1700000000.csv",
    "another line"
  )
  path <- .find_trace_lines(lines)
  expect_equal(path, "/tmp/ferx_trace_42_1700000000.csv")
})

test_that(".ferx_find_trace_from_lines returns NULL when no trace line present", {
  expect_null(.find_trace_lines(c("foo", "bar")))
  expect_null(.find_trace_lines(character(0)))
})

test_that(".ferx_find_trace_from_bg delegates to .ferx_find_trace_from_lines", {
  fake_bg <- structure(list(), class = "fake_bg")
  fake_bg$read_output_lines <- function() character(0)
  fake_bg$read_error_lines  <- function() {
    "[ferx] optimizer trace -> /tmp/ferx_trace_99.csv"
  }
  path <- .find_trace_bg(fake_bg)
  expect_equal(path, "/tmp/ferx_trace_99.csv")
})

test_that(".ferx_find_trace_from_lines tolerates spaces in the path", {
  lines <- c(
    "Some other line",
    "[ferx] optimizer trace -> /Users/Some User/tmp/ferx_trace_1.csv"
  )
  expect_equal(
    .find_trace_lines(lines),
    "/Users/Some User/tmp/ferx_trace_1.csv"
  )
})

test_that(".ferx_find_trace_from_lines strips trailing CR (Windows line endings)", {
  lines <- c(
    "[ferx] optimizer trace -> /tmp/ferx_trace_42.csv\r",
    "next line\r"
  )
  expect_equal(
    .find_trace_lines(lines),
    "/tmp/ferx_trace_42.csv"
  )
})

test_that(".ferx_find_trace_from_lines picks the announcement line over noise", {
  # Earlier line mentions "trace" but doesn't end in .csv; later line is the
  # real announcement.
  lines <- c(
    "WARNING: trace optimizer setting changed",
    "[ferx] optimizer trace -> /var/folders/xx/zz/T/ferx_trace_7.csv",
    "post-fit chatter"
  )
  expect_equal(
    .find_trace_lines(lines),
    "/var/folders/xx/zz/T/ferx_trace_7.csv"
  )
})

test_that(".ferx_find_trace_from_lines tolerates spaces via fallback path", {
  # No "optimizer trace" keyword, so the fallback branch runs.
  lines <- c(
    "Starting optimizer",
    "trace output -> /Users/Some User/tmp/ferx_trace_5.csv"
  )
  expect_equal(
    .find_trace_lines(lines),
    "/Users/Some User/tmp/ferx_trace_5.csv"
  )
})

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

test_that("ferx_collect rejects non-ferx_job input", {
  expect_error(ferx::ferx_collect(list()), "ferx_job")
  expect_error(ferx::ferx_collect("not a handle"), "ferx_job")
})

# -- callr sidecar tests ------------------------------------------------------

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
