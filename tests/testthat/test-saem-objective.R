# What a SAEM fit's OFV is (ferx-core #2973): NONMEM's SAEM objective by
# default, with no AIC/BIC; `saem_ofv` picks a likelihood instead.

saem_obj_fit <- function(...) {
  ex <- ferx_example("warfarin_saem")
  ferx_fit(
    ex$model, ex$data,
    method = "saem", verbose = FALSE, covariance = FALSE,
    settings = c(list(n_exploration = 20L, n_convergence = 20L, seed = 42L),
                 list(...))
  )
}

test_that("a default SAEM fit reports the SAEM objective and no AIC/BIC", {
  fit <- saem_obj_fit()
  so <- fit$saem_objective
  expect_type(so, "list")
  expect_identical(so$kind, "saem")
  expect_false(so$is_likelihood)
  expect_true(is.finite(fit$ofv))
  expect_true(is.finite(so$foce_ofv))
  expect_true(is.na(so$mc_standard_error))
  expect_true(is.na(fit$aic))
  expect_true(is.na(fit$bic))

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "AIC: n/a    BIC: n/a", fixed = TRUE)
  expect_match(out, "SAEM objective (NONMEM SAEMOBJ convention; not a likelihood)",
               fixed = TRUE)
  expect_match(out, "saem_ofv = \"imp\"", fixed = TRUE)
  expect_match(out, sprintf("FOCE/Laplace OFV at these estimates: %.4f", so$foce_ofv),
               fixed = TRUE)

  sout <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(sout, "AIC: n/a  BIC: n/a", fixed = TRUE)
  expect_match(sout, "not a likelihood", fixed = TRUE)

  expect_warning(b <- ferx_bic(fit), "saem_ofv", fixed = TRUE)
  expect_identical(b, NA_real_)
})

test_that("saem_ofv = 'focei' keeps the FOCE/Laplace OFV and its AIC/BIC", {
  fit <- saem_obj_fit(saem_ofv = "focei")
  expect_null(fit$saem_objective)
  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
  expect_gt(fit$aic, fit$ofv)
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_false(grepl("SAEM objective", out, fixed = TRUE))
})

test_that("saem_ofv = 'imp' reports an importance-sampled -2 log L", {
  fit <- saem_obj_fit(saem_ofv = "imp", saem_ofv_samples = 200L)
  so <- fit$saem_objective
  expect_identical(so$kind, "imp")
  expect_true(so$is_likelihood)
  expect_true(is.finite(so$mc_standard_error))
  expect_equal(so$n_samples, 200)
  expect_true(is.finite(fit$aic))
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "importance-sampled -2 log L", fixed = TRUE)
})

test_that("non-SAEM fits carry no saem_objective", {
  expect_null(warfarin_fit()$saem_objective)
})
