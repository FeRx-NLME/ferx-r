# print.ferx_fit() refuses a fit whose OMEGA or OMEGA_IOV row fields it cannot
# read (.ferx_check_row_fields(): #545, #562) before it prints anything, so a
# capture.output() or knitr chunk gets the whole report or none of it (#567).
# Before, the header and the THETA table were already written when the OMEGA
# rows refused.

# What print(fit) wrote before it stopped, and the condition it stopped with.
# capture.output() drops its lines on an error, so sink to a connection here.
print_until_refusal <- function(fit) {
  lines <- character()
  con <- textConnection("lines", "w", local = TRUE)
  sink(con)
  err <- tryCatch(print(fit), ferx_fit_field_error = identity,
                  finally = { sink(); close(con) })
  list(lines = lines, err = err)
}

expect_silent_refusal <- function(fit, field) {
  res <- print_until_refusal(fit)
  expect_s3_class(res$err, "ferx_fit_field_error")
  expect_identical(res$err$field, field)
  expect_identical(res$lines, character())
}

test_that("print refuses an unreadable OMEGA row field before the first line (#567)", {
  base <- make_fake_fit(omega = diag(c(0.1, 0.2)), eta_names = c("ETA_CL", "ETA_V"),
                        eta_param_types = c("additive", "log_normal"),
                        covariance_status = "computed")
  # The good fit prints the header and THETA table the refusals must not.
  out <- print_until_refusal(base)
  expect_false(inherits(out$err, "condition"))
  expect_true(any(grepl("THETA", out$lines, fixed = TRUE)))
  expect_true(any(grepl("OMEGA", out$lines, fixed = TRUE)))

  f <- base; f$eta_param_types <- "additive"
  expect_silent_refusal(f, "eta_param_types")
  f <- base; f$eta_param_types <- c("lognormal", "additive")
  expect_silent_refusal(f, "eta_param_types")
  f <- base; f$covariance_status <- "done"
  expect_silent_refusal(f, "covariance_status")
})

test_that("print refuses an unreadable OMEGA_IOV row field before the first line (#567)", {
  base <- make_fake_fit(
    omega = matrix(0.10, 1, 1),
    omega_iov = diag(c(0.04, 0.02)), kappa_names = c("KAPPA_CL", "KAPPA_V"),
    kappa_param_types = c("additive", "log_normal"), covariance_status = "computed"
  )
  out <- print_until_refusal(base)
  expect_false(inherits(out$err, "condition"))
  expect_true(any(grepl("OMEGA_IOV", out$lines, fixed = TRUE)))

  f <- base; f$kappa_param_types <- "additive"
  expect_silent_refusal(f, "kappa_param_types")
  f <- base; f$kappa_param_types <- c("additive", NA)
  expect_silent_refusal(f, "kappa_param_types")
})
