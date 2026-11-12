# Reporting a VI fit (#507). VI reports a bound, not an objective, so under
# the default `vi_final_ofv = none` a VI fit's OFV, AIC and BIC are NaN by
# design. print(), summary(), ferx_bic() and ferx_runlog() have to say so -
# and how to get a comparable OFV - instead of passing the NaN on silently.

nan_note <- "OFV is NaN because VI estimates a bound"

test_that("print() of a VI fit shows the VI block and explains the NaN OFV (#507)", {
  fit <- vi_fit()
  expect_true(is.nan(fit$ofv))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "VARIATIONAL INFERENCE", fixed = TRUE)
  expect_match(out, sprintf("-2*ELBO: %.4f  (a bound, not an OFV", fit$vi$neg_two_elbo),
               fixed = TRUE)
  expect_match(out, nan_note, fixed = TRUE)
  expect_match(out, "vi_final_ofv = \"laplace\"", fixed = TRUE)
  expect_match(out, "agq_eval_only = TRUE", fixed = TRUE)
  expect_match(out, sprintf("Iterations: %d", as.integer(fit$vi$n_iterations)),
               fixed = TRUE)
  expect_match(out, sprintf("Family: %s", fit$vi$family), fixed = TRUE)
  expect_match(out, "ELBO tightness ratio:", fixed = TRUE)
  expect_match(out, "posterior means, not modes", fixed = TRUE)

  sout <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(sout, "Variational inference:", fixed = TRUE)
  expect_match(sout, nan_note, fixed = TRUE)
})

test_that("ferx_bic() and ferx_runlog() give NA for a VI fit with no OFV (#507)", {
  fit <- vi_fit()
  for (type in c("mixed", "fixed", "iiv", "random")) {
    expect_warning(b <- ferx_bic(fit, type), "vi_final_ofv", fixed = TRUE)
    expect_identical(b, NA_real_)
  }
  log <- ferx_runlog(fit, verbose = FALSE)
  expect_match(log, "OFV   = NA\n", fixed = TRUE)
  expect_match(log, "BIC   = NA\n", fixed = TRUE)
  # The only NaN left is the note saying why.
  expect_false(grepl("= NaN", log, fixed = TRUE))
  expect_match(log, "vi_final_ofv", fixed = TRUE)
  expect_match(log, "-2*ELBO:", fixed = TRUE)
})

test_that("a VI fit with a Laplace OFV shows the bound but no NaN note (#507)", {
  fit <- vi_fit("laplace")
  expect_true(is.finite(fit$ofv))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "VARIATIONAL INFERENCE", fixed = TRUE)
  expect_false(grepl(nan_note, out, fixed = TRUE))
  expect_silent(b <- ferx_bic(fit, "fixed"))
  expect_true(is.finite(b))
})

test_that("the VI block flags a loose bound, FD subjects and a superseding stage (#507)", {
  fake <- list(
    ofv = 123.4,
    vi = list(
      neg_two_elbo = 130, data_term = 50, kl_term = 15,
      n_iterations = 400, converged = FALSE, family = "mean_field",
      n_mc_samples = 32, kl = "analytic", n_kl_fallback_subjects = 0,
      n_fd_subjects = 3, elbo_tightness_ratio = 78.2,
      superseded_by = "FOCEI"
    )
  )
  lines <- ferx:::.ferx_vi_report_lines(fake)
  out <- paste(lines, collapse = "\n")
  expect_match(out, "ELBO tightness ratio: 78.20  [!] above 25", fixed = TRUE)
  expect_match(out, "Converged: NO", fixed = TRUE)
  expect_match(out, "Finite-difference eta gradients for 3 subject(s)", fixed = TRUE)
  expect_match(out, "Superseded by FOCEI", fixed = TRUE)
  # Superseded: the reported EBEs are the later stage's, not VI's means.
  expect_false(grepl("posterior means, not modes", out, fixed = TRUE))
  expect_false(grepl(nan_note, out, fixed = TRUE))

  # At the limit itself the bound is still usable, as in the engine (`>`).
  fake$vi$elbo_tightness_ratio <- 25
  fake$vi$n_fd_subjects <- 0
  fake$vi$superseded_by <- NA_character_
  out <- paste(ferx:::.ferx_vi_report_lines(fake), collapse = "\n")
  expect_false(grepl("[!]", out, fixed = TRUE))
  expect_false(grepl("Finite-difference", out, fixed = TRUE))
  expect_false(grepl("Superseded", out, fixed = TRUE))

  # No VI result, no block, no NA from ferx_bic's VI arm.
  expect_null(ferx:::.ferx_vi_report_lines(list(ofv = 1)))
  expect_false(ferx:::.ferx_vi_ofv_missing(list(ofv = NaN)))
})
