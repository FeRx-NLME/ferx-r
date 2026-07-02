
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
