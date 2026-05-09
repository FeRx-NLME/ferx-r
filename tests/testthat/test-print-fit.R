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

test_that("print.ferx_fit shows SD_logit for logit eta_param_types", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "logit")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("logit", omega_line))
  expect_true(grepl("SD_logit", omega_line))
  # CV% is not meaningful for logit ETAs — must not appear
  expect_false(grepl("CV% =", omega_line))
})

test_that("print.ferx_fit CV% equals zero when omega diagonal is zero", {
  fit <- make_fake_fit(omega = matrix(0, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("CV% = 0\\.0", omega_line))
})

# SIGMA display (ferx#59 / ferx-nlme#57): sigma is on the SD scale, so
# print() must show variance = sigma^2 always and CV% = sigma * 100 only
# for proportional components.

test_that("print.ferx_fit shows variance + CV% for proportional sigma using its declared name", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = "PROP_ERR",
    sigma_types = "proportional"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("PROP_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.1^2 = 0.01; CV% = 0.1 * 100 = 10.0.
  expect_true(grepl("var = 0\\.010000", sigma_line))
  expect_true(grepl("CV% = 10\\.0", sigma_line))
  # No legacy SIGMA(1) label when a name is supplied.
  expect_false(any(grepl("SIGMA\\(1\\)", out)))
})

test_that("print.ferx_fit shows variance but no CV% for additive sigma", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.5,
    sigma_names = "ADD_ERR",
    sigma_types = "additive"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("ADD_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.5^2 = 0.25; no CV% on observation-unit scale.
  expect_true(grepl("var = 0\\.250000", sigma_line))
  expect_false(grepl("CV%", sigma_line))
})

test_that("print.ferx_fit handles combined error: CV% on prop component only", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = c(0.2, 0.5),
    sigma_names = c("PROP_ERR", "ADD_ERR"),
    sigma_types = c("proportional", "additive")
  )
  out <- capture.output(print(fit))
  prop_line <- out[grepl("PROP_ERR", out)]
  add_line  <- out[grepl("ADD_ERR",  out)]
  expect_length(prop_line, 1L)
  expect_length(add_line,  1L)
  expect_true(grepl("CV% = 20\\.0", prop_line))
  expect_false(grepl("CV%", add_line))
  expect_true(grepl("var = 0\\.040000", prop_line))
  expect_true(grepl("var = 0\\.250000", add_line))
})

test_that("print.ferx_fit falls back to SIGMA(i) when sigma_names is missing", {
  # sigma_names absent (older Rust binary or unit test that didn't pass them).
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = NULL,
    sigma_types = NULL
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("SIGMA\\(1\\)", out)))
})

# Block-omega correlations (#60): print uses omega_param_corr when the engine
# provides it (bivariate-lognormal formula), and falls back to the eta-level
# Pearson formula when it is absent.

test_that("print.ferx_fit uses omega_param_corr value and 'param corr' label when present", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = pc,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\)", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
})

test_that("print.ferx_fit falls back to eta-level corr label when omega_param_corr is NULL", {
  # cov = 0.025, vars = 0.10 -> Pearson corr = 0.025 / 0.10 = 0.25
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = NULL,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\)", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("\\(corr = 0\\.2500\\)", corr_line))
  expect_false(grepl("param corr", corr_line))
})

test_that("print.ferx_fit uses omega_iov_param_corr when present for IOV correlations", {
  om     <- matrix(0.10, 1, 1)
  iov    <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  iov_pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega                 = om,
    omega_iov             = iov,
    omega_iov_param_corr  = iov_pc,
    kappa_names           = c("KAPPA1", "KAPPA2"),
    se_kappa              = NULL,
    shrinkage_kappa       = NULL
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("KAPPA2 ~ KAPPA1", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
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
