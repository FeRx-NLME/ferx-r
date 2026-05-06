# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# print.ferx_fit CV% formula --------------------------------------------

test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is absent (fallback)", {
  # No eta_param_types field -> defaults to log_normal.
  # omega = 0.40: exact CV% = sqrt(exp(0.40) - 1) * 100 ≈ 70.1 (not 63.2 from sqrt(0.40)*100)
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(length(omega_line) == 1L)
  expect_false(grepl("63\\.2", omega_line))  # old approximate value
  expect_true(grepl("70\\.1", omega_line))   # exact log-normal CV%
})

test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is 'log_normal'", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "log_normal")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("70\\.1", omega_line))
})

test_that("print.ferx_fit shows CV% = N/A for non-log-normal eta_param_types (handled in #53)", {
  # logit/additive display is deferred to #53; show N/A rather than a misleading 0
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "logit")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("CV% = N/A", omega_line))
})

test_that("print.ferx_fit CV% equals zero when omega diagonal is zero", {
  fit <- make_fake_fit(omega = matrix(0, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("CV% = 0\\.0", omega_line))
})

test_that("print.ferx_fit uses exact log-normal CV% for OMEGA_IOV when kappa_param_types absent", {
  # kappa = 0.20 -> exact CV% = sqrt(exp(0.20) - 1) * 100 ≈ 47.1 (not 44.7 from sqrt(0.20)*100)
  fit <- make_fake_fit(
    omega           = matrix(0.10, 1, 1),
    omega_iov       = matrix(0.20, 1, 1),
    kappa_names     = "KAPPA1",
    se_kappa        = NULL,
    shrinkage_kappa = NULL
  )
  out <- capture.output(print(fit))
  kappa_line <- out[grepl("KAPPA1", out)]
  expect_true(length(kappa_line) >= 1L)
  expect_false(grepl("44\\.7", kappa_line[1]))  # old approximate value
  expect_true(grepl("47\\.1", kappa_line[1]))   # exact log-normal CV%
})
