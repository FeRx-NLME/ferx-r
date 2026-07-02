
# ---- header from test-trace-plot.R ----
# Tests for ferx_plot_trace() and .add_phase_lines(). Plots are drawn to a
# null PDF device so nothing is rendered; we assert the call returns the trace
# invisibly and exercises each method-specific metric branch.

.add_phase_lines <- getFromNamespace(".add_phase_lines", "ferx")

# Write a trace CSV with the columns ferx_plot_trace() reads.
write_trace_csv <- function(method = "focei", phase = "",
                            mh = NA_real_, lm = NA_real_, grad = NA_real_,
                            n = 5L) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      iter = seq_len(n), method = method, phase = phase,
      ofv = seq.int(100L, length.out = n, by = -1L),
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







test_that("ferx_plot_trace draws the SAEM accept-rate panel and phase lines", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- write_trace_csv(method = "saem",
                       phase = rep(c("burn-in", "stochastic"), length.out = 5L),
                       mh = seq(0.2, 0.6, length.out = 5L))
  on.exit(unlink(p), add = TRUE)
  res <- withVisible(ferx_plot_trace(p))
  expect_false(res$visible)
  tr <- res$value
  expect_s3_class(tr, "data.frame")
})
test_that("ferx_plot_trace draws the GN LM-lambda panel", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- write_trace_csv(method = "gn", lm = seq(1, 0.1, length.out = 5L))
  on.exit(unlink(p), add = TRUE)
  expect_invisible(ferx_plot_trace(p))
})
test_that("ferx_plot_trace shows a placeholder when the metric is unavailable", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # focei with an all-NA gradient column -> 'metric not available' branch.
  p <- write_trace_csv(method = "focei", grad = NA_real_)
  on.exit(unlink(p), add = TRUE)
  expect_invisible(ferx_plot_trace(p))
})
test_that("ferx_plot_trace honours log_ofv", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- write_trace_csv(method = "focei", grad = seq(1, 0.1, length.out = 5L))
  on.exit(unlink(p), add = TRUE)
  expect_invisible(ferx_plot_trace(p, log_ofv = TRUE))
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

# ---- header from test-trace.R ----








# ferx_plot_trace — Tier 1, no Rust needed





test_that("ferx_plot_trace(fit) returns invisibly without error", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_no_error({
    res <- withVisible(ferx_plot_trace(fit))
  })
  expect_false(res$visible)
})
test_that("ferx_plot_trace(fit, log_ofv = TRUE) runs without error", {
  path <- write_fake_trace()
  on.exit(unlink(path))
  fit <- structure(list(trace_path = path), class = "ferx_fit")
  expect_no_error(ferx_plot_trace(fit, log_ofv = TRUE))
})
test_that("ferx_plot_trace errors on object with no trace", {
  fit <- structure(list(trace_path = NULL), class = "ferx_fit")
  expect_error(ferx_plot_trace(fit))
})
