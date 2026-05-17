# Tests for SDE (diffusion / EKF) surface on the R side.
# All tests use make_fake_fit() from helper-trace.R — no real fit needed.

# uses_sde flag in print.ferx_fit ----------------------------------------

test_that("print.ferx_fit shows '+ SDE (EKF)' when uses_sde is TRUE", {
  fit <- make_fake_fit(uses_sde = TRUE, omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  method_line <- out[grepl("Estimation method", out)]
  expect_length(method_line, 1L)
  expect_match(method_line, "SDE \\(EKF\\)")
})

test_that("print.ferx_fit does not show SDE tag when uses_sde is FALSE", {
  fit <- make_fake_fit(uses_sde = FALSE, omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("SDE", out)))
})

test_that("print.ferx_fit does not show SDE tag when uses_sde is absent", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("SDE", out)))
})

test_that("print.ferx_fit shows SDE tag after method chain when uses_sde TRUE", {
  fit <- make_fake_fit(
    uses_sde     = TRUE,
    method       = "focei",
    method_chain = c("saem", "focei"),
    omega        = matrix(0.09, 1, 1)
  )
  out <- capture.output(print(fit))
  chain_line <- out[grepl("Estimation chain", out)]
  expect_length(chain_line, 1L)
  expect_match(chain_line, "SDE \\(EKF\\)")
})

# uses_sde flag in summary / print.ferx_summary ---------------------------

test_that("summary.ferx_fit carries uses_sde = TRUE into ferx_summary", {
  fit <- make_fake_fit(uses_sde = TRUE)
  s   <- summary(fit)
  expect_true(isTRUE(s$uses_sde))
})

test_that("summary.ferx_fit carries uses_sde = FALSE into ferx_summary", {
  fit <- make_fake_fit(uses_sde = FALSE)
  s   <- summary(fit)
  expect_false(s$uses_sde)
})

test_that("summary.ferx_fit defaults uses_sde to FALSE when field absent", {
  fit <- make_fake_fit()
  s   <- summary(fit)
  expect_false(s$uses_sde)
})

test_that("print.ferx_summary shows '+ SDE (EKF)' when uses_sde TRUE", {
  fit <- make_fake_fit(uses_sde = TRUE)
  out <- capture.output(print(summary(fit)))
  method_line <- out[grepl("Method:", out)]
  expect_length(method_line, 1L)
  expect_match(method_line, "SDE \\(EKF\\)")
})

test_that("print.ferx_summary omits SDE tag when uses_sde FALSE", {
  fit <- make_fake_fit(uses_sde = FALSE)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("SDE", out)))
})

# DIFF_* theta parameters flow through naturally --------------------------

test_that("DIFF_* theta parameters appear in print.ferx_fit theta table", {
  fit <- make_fake_fit(
    uses_sde  = TRUE,
    theta     = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    se_theta  = NULL,
    omega     = matrix(0.09, 1, 1)
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("DIFF_CENTRAL", out)))
})

# ferx_estimates() with DIFF_* theta -------------------------------------

test_that("ferx_estimates() includes DIFF_* theta rows", {
  fit <- make_fake_fit(
    uses_sde = TRUE,
    theta    = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    omega    = matrix(0.09, 1, 1)
  )
  est <- ferx_estimates(fit)
  expect_true(any(est$param == "DIFF_CENTRAL"))
  diff_row <- est[est$param == "DIFF_CENTRAL", ]
  expect_equal(diff_row$estimate, 0.5)
})

# ferx_save_fit / ferx_load_fit round-trip --------------------------------

test_that("uses_sde survives a ferx_save_fit / ferx_load_fit round-trip (TRUE)", {
  fit <- make_fake_fit(
    uses_sde = TRUE,
    theta    = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    omega    = matrix(0.09, 1, 1),
    sigma    = 1.0,
    sigma_names = "ADD"
  )
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp))
  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)
  expect_true(isTRUE(fit2$uses_sde))
})

test_that("uses_sde survives a ferx_save_fit / ferx_load_fit round-trip (FALSE)", {
  fit <- make_fake_fit(
    uses_sde = FALSE,
    omega    = matrix(0.09, 1, 1),
    sigma    = 1.0,
    sigma_names = "ADD"
  )
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp))
  ferx_save_fit(fit, tmp)
  fit2 <- ferx_load_fit(tmp)
  expect_false(isTRUE(fit2$uses_sde))
})
