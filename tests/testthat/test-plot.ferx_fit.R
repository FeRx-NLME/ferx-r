# Tests for plot.ferx_fit(), plot.ferx_job(), and .add_phase_lines(). Plots
# are drawn to a null PDF device so nothing is rendered; we assert the call
# returns the trace invisibly and exercises each method-specific metric
# branch. write_fake_trace() / make_fake_fit() come from helper-trace.R.

.add_phase_lines <- getFromNamespace(".add_phase_lines", "ferx")

# Write a trace CSV with the columns plot.ferx_fit() reads.
write_trace_csv <- function(method = "focei", phase = "",
                            mh = NA_real_, lm = NA_real_, grad = NA_real_,
                            ofv = NULL, n = 5L) {
  if (is.null(ofv)) ofv <- seq.int(100L, length.out = n, by = -1L)
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      iter = seq_len(n), method = method, phase = phase,
      ofv = ofv,
      wall_ms = 0L, grad_norm = grad, step_norm = NA_real_,
      inner_iter_count = NA_integer_, optimizer = "slsqp",
      lm_lambda = lm, ofv_delta = NA_real_, step_accepted = NA_integer_,
      cond_nll = NA_real_, gamma = NA_real_, mh_accept_rate = mh,
      stringsAsFactors = FALSE
    ),
    path, row.names = FALSE
  )
  path
}

test_that("plot.ferx_fit draws the SAEM accept-rate panel and phase lines", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  path <- write_trace_csv(method = "saem",
                       phase = rep(c("burn-in", "stochastic"), length.out = 5L),
                       mh = seq(0.2, 0.6, length.out = 5L))
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  res <- withVisible(plot(fit))
  expect_false(res$visible)
  tr <- res$value
  expect_s3_class(tr, "data.frame")
})
test_that("plot.ferx_fit draws the GN LM-lambda panel", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  path <- write_trace_csv(method = "gn", lm = seq(1, 0.1, length.out = 5L))
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_invisible(plot(fit))
})
test_that("plot.ferx_fit shows a placeholder when the metric is unavailable", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # focei with an all-NA gradient column -> 'metric not available' branch.
  path <- write_trace_csv(method = "focei", grad = NA_real_)
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_invisible(plot(fit))
})
test_that("plot.ferx_fit honours log_ofv", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  path <- write_trace_csv(method = "focei", grad = seq(1, 0.1, length.out = 5L))
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_invisible(plot(fit, log_ofv = TRUE))
})
test_that("plot.ferx_fit(fit) returns invisibly without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  path <- write_fake_trace()
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_no_error({
    res <- withVisible(plot(fit))
  })
  expect_false(res$visible)
})
test_that("plot.ferx_fit errors on object with no trace", {
  fit <- structure(list(trace_path = NULL), class = "ferx_fit")
  expect_error(plot(fit))
})

test_that("plot.ferx_fit applies cummin to FOCEI OFV by default (monotonic)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # Non-monotonic raw OFV: dips then rises then dips again.
  raw_ofv <- c(100, 95, 98, 90, 92)
  path <- write_trace_csv(method = "focei", grad = seq(1, 0.1, length.out = 5L),
                           ofv = raw_ofv)
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")

  core <- getFromNamespace(".ferx_plot_trace_core", "ferx")
  tr <- ferx_trace(fit)
  expect_false(identical(tr$ofv, cummin(tr$ofv)))  # raw trace is non-monotonic

  expect_invisible(core(tr, monotonic = TRUE))
  expect_invisible(core(tr, monotonic = FALSE))
  # the returned trace itself stays raw; only the drawn panel is transformed
  expect_equal(withVisible(core(tr, monotonic = TRUE))$value$ofv, raw_ofv)
})
test_that("plot.ferx_fit leaves GN/SAEM OFV untouched regardless of monotonic=", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  raw_ofv <- c(100, 95, 98, 90, 92)
  path <- write_trace_csv(method = "gn", lm = seq(1, 0.1, length.out = 5L),
                           ofv = raw_ofv)
  on.exit(unlink(path), add = TRUE)
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_invisible(plot(fit, monotonic = TRUE))
  expect_invisible(plot(fit, monotonic = FALSE))
})
test_that(".ferx_monotonic_ofv only applies cummin to FOCE/FOCEI rows in a chained trace", {
  monotonic_ofv <- getFromNamespace(".ferx_monotonic_ofv", "ferx")
  # SAEM segment dips-rises-dips, then a FOCEI segment does the same.
  method <- c("saem", "saem", "saem", "focei", "focei", "focei")
  ofv    <- c(200,     180,    185,    160,     155,     158)

  out <- monotonic_ofv(ofv, method, monotonic = TRUE)
  expect_equal(out[method == "saem"], ofv[method == "saem"])       # untouched
  expect_equal(out[method == "focei"], cummin(ofv[method == "focei"]))
})
test_that(".ferx_monotonic_ofv is a no-op when monotonic = FALSE or no FOCE rows exist", {
  monotonic_ofv <- getFromNamespace(".ferx_monotonic_ofv", "ferx")
  ofv <- c(200, 180, 185, 160, 155, 158)
  expect_equal(monotonic_ofv(ofv, rep("focei", 6L), monotonic = FALSE), ofv)
  expect_equal(monotonic_ofv(ofv, rep("gn", 6L),    monotonic = TRUE),  ofv)
})
test_that("plot.ferx_fit(y = ...) errors instead of silently mis-binding to log_ofv", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_error(plot(fit, y = 1), "does not use `y`")
})
test_that("plot.ferx_job(y = ...) errors instead of silently mis-binding to log_ofv", {
  handle <- structure(
    list(backend = "callr", model_label = "m.ferx",
         sidecar_path = tempfile(fileext = ".txt")),
    class = "ferx_job"
  )
  expect_error(plot(handle, y = 1), "does not use `y`")
})

test_that("plot.ferx_job plots the trace accumulated so far", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  sidecar <- tempfile(fileext = ".txt")
  writeLines(write_fake_trace(n = 3L), sidecar)
  handle <- structure(
    list(backend = "callr", model_label = "warfarin.ferx",
         sidecar_path = sidecar),
    class = "ferx_job"
  )
  res <- withVisible(plot(handle))
  expect_false(res$visible)
  expect_s3_class(res$value, "data.frame")
})
test_that("plot.ferx_job shows a placeholder when no trace exists yet", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  handle <- structure(
    list(backend = "callr", model_label = "warfarin.ferx",
         sidecar_path = tempfile(fileext = ".txt")),  # never written
    class = "ferx_job"
  )
  res <- withVisible(plot(handle))
  expect_false(res$visible)
  expect_null(res$value)
})

test_that(".add_phase_lines draws boundaries only when phases change", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(1:4, 1:4)  # need an active plot for abline()
  expect_invisible(.add_phase_lines(data.frame(iter = 1:4,
                                               phase = c("a", "a", "b", "b"))))
  expect_invisible(.add_phase_lines(data.frame(iter = 1:2, phase = c("", ""))))
  expect_invisible(.add_phase_lines(data.frame(iter = 1:2)))  # no phase column
})
