# ferx_covariance() / ferx_sir() run on the fit's worker-pool width unless
# told otherwise, and report the width they ran on as `n_threads_used`
# (#527). Before, both ran on the engine default (`FitOptions::default()`,
# threads = None) and kept reporting the fit's width.

test_that(".ferx_step_threads: NULL takes the fit's width, 0 / n pass through (#527)", {
  step_threads <- getFromNamespace(".ferx_step_threads", "ferx")
  expect_identical(step_threads(NULL, list(n_threads_used = 2L)), 2L)
  expect_identical(step_threads(NULL, list(n_threads_used = 3)), 3L)
  # A fit that records no usable width gets the engine default.
  expect_identical(step_threads(NULL, list()), 0L)
  expect_identical(step_threads(NULL, list(n_threads_used = NA_integer_)), 0L)
  expect_identical(step_threads(NULL, list(n_threads_used = 0L)), 0L)
  # An explicit value wins over the fit's.
  expect_identical(step_threads(0, list(n_threads_used = 2L)), 0L)
  expect_identical(step_threads(5L, list(n_threads_used = 2L)), 5L)

  for (bad in list(-1L, 1.5, NA_real_, Inf, c(1L, 2L), "2", 2^40)) {
    expect_error(step_threads(bad, list()), "`threads` must be NULL",
                 info = format(bad))
  }
})

test_that("ferx_covariance / ferx_sir hand the binding the step's width (#527)", {
  skip_if_not_installed("mockery")
  fit <- warfarin_fit()
  fit$n_threads_used <- 2L
  # ferx_sir() refuses a fit without a covariance matrix before the binding;
  # any square matrix of the right size gets it there.
  n <- length(fit$theta) + nrow(fit$omega) + length(fit$sigma)
  sir_fit <- fit
  sir_fit$cov_matrix <- diag(n)

  seen <- function(f, binding, x, ...) {
    cap <- capture_binding_args()
    mockery::stub(f, binding, cap$fake)
    expect_error(f(x, ...), "captured")
    cap$seen()$threads
  }
  for (case in list(
    list(f = ferx_covariance, binding = "ferx_rust_covariance", x = fit),
    list(f = ferx_sir, binding = "ferx_rust_sir", x = sir_fit)
  )) {
    expect_identical(seen(case$f, case$binding, case$x), 2L)
    expect_identical(seen(case$f, case$binding, case$x, threads = 0), 0L)
    expect_identical(seen(case$f, case$binding, case$x, threads = 1L), 1L)
    expect_error(case$f(case$x, threads = -1L), "`threads` must be NULL")
  }
})

test_that("ferx_covariance reports the width it ran on (#527)", {
  fit <- warfarin_fit()
  # The fit's width is what a NULL `threads` runs on.
  fit$n_threads_used <- 2L
  out <- ferx_covariance(fit)
  expect_identical(out$n_threads_used, 2L)

  pinned <- ferx_covariance(fit, threads = 1L)
  expect_identical(pinned$n_threads_used, 1L)
  # The width moves resources, not the bits.
  expect_identical(pinned$cov_matrix, out$cov_matrix)

  default <- ferx_covariance(fit, threads = 0)
  expect_true(is.integer(default$n_threads_used) && default$n_threads_used >= 1L)
})

test_that("ferx_sir reports the width it ran on (#527)", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  fit$n_threads_used <- 2L
  run <- function(...) {
    ferx_sir(fit, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L, ...)
  }
  out <- run()
  expect_identical(out$n_threads_used, 2L)

  pinned <- run(threads = 1L)
  expect_identical(pinned$n_threads_used, 1L)
  expect_identical(pinned$sir_ci_theta, out$sir_ci_theta)

  default <- run(threads = 0)
  expect_true(is.integer(default$n_threads_used) && default$n_threads_used >= 1L)
})
