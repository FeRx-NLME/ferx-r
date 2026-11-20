# End-to-end R verification suite for SAEM regular-vine copula random effects
# (`omega_dist = "vine"`, ferx-core #225).

test_that("ferx_example registers vine_tail_iv, vine_tail_iv_gaussian, and warfarin_vine", {
  avail <- ferx_example()
  expect_true(all(c("vine_tail_iv", "vine_tail_iv_gaussian", "warfarin_vine") %in% avail))

  ex_vine <- ferx_example("vine_tail_iv")
  ex_gauss <- ferx_example("vine_tail_iv_gaussian")
  ex_warf <- ferx_example("warfarin_vine")
  expect_true(file.exists(ex_vine$model))
  expect_true(file.exists(ex_vine$data))
  expect_identical(ex_vine$data, ex_gauss$data)
  expect_true(file.exists(ex_warf$model))
  expect_identical(ex_warf$data, ferx_example("warfarin")$data)
})

test_that("invalid omega_dist = 'vine' combinations raise coded ferx_engine_error refusals", {
  ex <- ferx_example("vine_tail_iv")

  # 1. Non-SAEM estimator -> E_OMEGA_DIST_NO_SAEM
  err_focei <- expect_error(
    suppressWarnings(ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE)),
    class = "ferx_engine_error"
  )
  expect_identical(err_focei$code, "E_OMEGA_DIST_NO_SAEM")

  # 2. Method chain c("saem", "focei") -> E_OMEGA_DIST_CHAIN
  err_chain <- expect_error(
    suppressWarnings(ferx_fit(ex$model, ex$data, method = c("saem", "focei"), verbose = FALSE)),
    class = "ferx_engine_error"
  )
  expect_identical(err_chain$code, "E_OMEGA_DIST_CHAIN")

  # 3. sir = TRUE -> E_OMEGA_DIST_SIR
  err_sir <- expect_error(
    ferx_fit(ex$model, ex$data, method = "saem", covariance = TRUE, sir = TRUE, verbose = FALSE),
    class = "ferx_engine_error"
  )
  expect_identical(err_sir$code, "E_OMEGA_DIST_SIR")
})

test_that("vine_tail_iv recovers Clayton lower-tail dependence, beats Gaussian, round-trips .fitrx, and simulates from the vine", {
  ex_vine <- ferx_example("vine_tail_iv")
  ex_gauss <- ferx_example("vine_tail_iv_gaussian")

  fit_vine <- ferx_fit(ex_vine$model, ex_vine$data, verbose = FALSE)
  fit_gauss <- ferx_fit(ex_gauss$model, ex_gauss$data, verbose = FALSE)

  expect_true(fit_vine$converged)
  expect_null(fit_gauss$vine_copula)
  expect_null(fit_gauss$vine_corrected_ofv)

  # 1. Structured `fit$vine_copula` surface
  vc <- fit_vine$vine_copula
  expect_type(vc, "list")
  expect_s3_class(vc$marginals, "data.frame")
  expect_identical(vc$marginals$eta, c("ETA_CL", "ETA_V"))
  expect_equal(vc$marginals$mean, c(0, 0))
  expect_true(all(vc$marginals$sd > 0))
  expect_length(vc$trees, 1L)
  expect_identical(vc$trees[[1L]]$tree, 1L)
  expect_s3_class(vc$pairs, "data.frame")
  expect_equal(nrow(vc$pairs), 1L)
  expect_gt(vc$pairs$kendall_tau[[1L]], 0.25)
  expect_gt(vc$pairs$tail_dep_lower[[1L]], 0.30)

  # 2. Marginal QMC IS OFV, AIC, BIC, and parameter counts
  expect_true(is.finite(fit_vine$ofv))
  expect_equal(fit_vine$vine_corrected_ofv, fit_vine$ofv)
  expect_true(is.finite(fit_vine$aic) && is.finite(fit_vine$bic))
  expect_lt(fit_vine$ofv, fit_gauss$ofv - 3.0)
  expect_equal(fit_vine$bic_inputs$omega, vc$n_free_params)

  # 3. `print()` and `summary()` include VINE COPULA section cleanly
  pr <- paste(capture.output(print(fit_vine)), collapse = "\n")
  expect_match(pr, "VINE COPULA")
  expect_match(pr, "lambda_L")
  sm <- paste(capture.output(print(summary(fit_vine))), collapse = "\n")
  expect_match(sm, "Vine copula:")

  # 4. `ferx_vine_density()`: 2-D joint quadrature == 1.0, 1-D marginal quadrature == 1.0, and `at =`
  grid2 <- ferx_vine_density(fit_vine, n = 80L, range_sd = 5.0)
  expect_s3_class(grid2, "data.frame")
  expect_identical(colnames(grid2), c("ETA_CL", "ETA_V", "density", "log_density"))
  u1 <- sort(unique(grid2$ETA_CL))
  u2 <- sort(unique(grid2$ETA_V))
  dx1 <- u1[[2L]] - u1[[1L]]
  dx2 <- u2[[2L]] - u2[[1L]]
  integral_2d <- sum(grid2$density) * dx1 * dx2
  expect_equal(integral_2d, 1.0, tolerance = 0.015)

  grid1 <- ferx_vine_density(fit_vine, etas = "ETA_CL", n = 200L, range_sd = 6.0)
  dx_m <- grid1$ETA_CL[[2L]] - grid1$ETA_CL[[1L]]
  expect_equal(sum(grid1$density) * dx_m, 1.0, tolerance = 0.005)

  pts <- ferx_vine_density(fit_vine, at = rbind(c(0, 0), c(-0.4, -0.4), c(0.4, 0.4)))
  expect_equal(nrow(pts), 3L)
  # Lower-tail density (-0.4, -0.4) exceeds symmetric upper-tail (+0.4, +0.4) under Clayton!
  expect_gt(pts$density[[2L]], pts$density[[3L]])

  # 5. `.fitrx` save/load round-trip preserves `vine_copula` and `vine_corrected_ofv`
  tmp <- tempfile(fileext = ".fitrx")
  on.exit(unlink(tmp), add = TRUE)
  ferx_save_fit(fit_vine, tmp)
  reloaded <- ferx_load_fit(tmp)
  expect_equal(reloaded$vine_corrected_ofv, fit_vine$vine_corrected_ofv)
  expect_equal(reloaded$vine_copula$pairs, fit_vine$vine_copula$pairs)
  expect_equal(reloaded$vine_copula$params, fit_vine$vine_copula$params)
  expect_equal(
    ferx_vine_density(reloaded, at = rbind(c(-0.3, -0.3)))$log_density,
    ferx_vine_density(fit_vine, at = rbind(c(-0.3, -0.3)))$log_density
  )

  # 6. `ferx_simulate(..., fit = fit_vine)` draws ETAs from the fitted Clayton R-vine
  sim <- ferx_simulate(ex_vine$model, ex_vine$data, n_sim = 25L, seed = 2026L, fit = fit_vine)
  sim_p <- ferx_sim_parameters(sim)
  e1 <- sim_p$ETA_CL
  e2 <- sim_p$ETA_V
  q10_1 <- stats::quantile(e1, 0.10)
  q10_2 <- stats::quantile(e2, 0.10)
  q90_1 <- stats::quantile(e1, 0.90)
  q90_2 <- stats::quantile(e2, 0.90)
  p_low <- mean(e1[e2 <= q10_2] <= q10_1)
  p_high <- mean(e1[e2 >= q90_2] >= q90_1)
  expect_gt(p_low, p_high + 0.12)

  # 7. `ferx_calc_npde()` evaluates predictive distribution errors under the fitted R-vine
  with_npde <- ferx_calc_npde(fit_vine, nsim = 100L, seed = 42L)
  expect_true("NPDE" %in% colnames(with_npde$sdtab))
  expect_true(all(is.finite(with_npde$sdtab$NPDE)))
})

test_that("warfarin_vine fits a 3-ETA regular vine cleanly", {
  ex <- ferx_example("warfarin_vine")
  fit <- suppressWarnings(ferx_fit(
    ex$model, ex$data,
    covariance = FALSE, verbose = FALSE,
    settings = list(n_exploration = 120L, n_convergence = 120L)
  ))
  expect_true(fit$converged)
  expect_type(fit$vine_copula, "list")
  expect_equal(nrow(fit$vine_copula$marginals), 3L)
  expect_length(fit$vine_copula$trees, 2L)
  expect_equal(nrow(fit$vine_copula$pairs), 3L)
  expect_true(is.finite(fit$ofv))
})
