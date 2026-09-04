# A plain (non-FIX) `block_sigma` estimates its off-diagonal correlation
# (ferx-core #847). Two things follow that the R layer used to get wrong:
#
#   * the fitted rho has to reach R and travel with the fit - everything that
#     reconstructs parameters from a fit (predict / simulate / NPDE / SIR /
#     covariance) otherwise falls back to the model file's declared value;
#   * the engine packs those correlations *last*, after sigma, so a covariance
#     matrix that counts every non-theta/non-sigma coordinate as omega labels
#     its trailing rows at the wrong offset.

# Model + data written to a temp dir: the package bundles no `block_sigma`
# example, and the data is the two-point-per-subject set the engine's own
# correlated-residual tests use.
rho_case <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    dir <- tempfile("ferx-rho-")
    dir.create(dir)
    model <- file.path(dir, "rho.ferx")
    data  <- file.path(dir, "rho.csv")
    writeLines(c(
      "[parameters]",
      "  theta TVCL(1.0, 0.01, 10.0)",
      "  theta TVV(10.0, 0.1, 100.0)",
      "  omega ETA_CL ~ 0.04",
      "  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30]",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV",
      "[structural_model]",
      "  pk one_cpt_iv(cl=CL, v=V)",
      "[error_model]",
      "  DV ~ combined(PROP_ERR, ADD_ERR)",
      "[fit_options]",
      "  method = focei"
    ), model)
    write.csv(data.frame(
      ID   = rep(1:3, each = 4),
      TIME = rep(c(0, 1, 3, 8), 3),
      DV   = c(0, 9.1, 6.8, 3.0, 0, 8.4, 6.2, 2.6, 0, 9.7, 7.3, 3.4),
      EVID = rep(c(1, 0, 0, 0), 3),
      AMT  = rep(c(100, 0, 0, 0), 3),
      CMT  = 1,
      MDV  = rep(c(1, 0, 0, 0), 3)
    ), data, row.names = FALSE, quote = FALSE)
    cached <<- list(model = model, data = data,
                    fit = ferx_fit(model, data, verbose = FALSE,
                                   covariance = TRUE))
    cached
  }
})

test_that("a free block_sigma correlation reaches the fit", {
  skip_on_cran()
  fit <- rho_case()$fit

  rc <- fit$residual_correlations
  expect_s3_class(rc, "data.frame")
  expect_equal(nrow(rc), 1L)
  expect_named(rc, c("sigma_i", "sigma_j", "name", "rho", "fixed", "se"))
  # 1-based indices into fit$sigma, and the off-diagonal label convention.
  expect_true(all(rc$sigma_i %in% seq_along(fit$sigma)))
  expect_true(all(rc$sigma_j %in% seq_along(fit$sigma)))
  expect_match(rc$name, " ~ ")
  expect_false(rc$fixed)
  # The declared init is 0.30; a fitted value equal to it would mean the
  # estimate never reached R.
  expect_true(is.finite(rc$rho))
  expect_true(abs(rc$rho) < 1)
  expect_false(isTRUE(all.equal(rc$rho, 0.30)))
})

test_that("covariance labels account for the trailing rho coordinates", {
  skip_on_cran()
  fit <- rho_case()$fit
  skip_if(is.null(fit$cov_matrix), "covariance step produced no matrix")

  nms <- rownames(fit$cov_matrix)
  expect_length(nms, nrow(fit$cov_matrix))
  # Every coordinate is named - the off-by-n_rho bug left an empty string
  # where an omega label had been shifted onto a sigma coordinate.
  expect_true(all(nzchar(nms)))
  expect_equal(nms[seq_along(fit$theta)], names(fit$theta))
  # The correlations are packed last, so they label the final rows.
  expect_equal(tail(nms, nrow(fit$residual_correlations)),
               fit$residual_correlations$name)
  # ...and the sigma names sit immediately before them, not shifted.
  n_rho <- nrow(fit$residual_correlations)
  sig_slots <- seq.int(length(nms) - n_rho - length(fit$sigma) + 1L,
                       length(nms) - n_rho)
  expect_equal(nms[sig_slots], fit$sigma_names)
})

test_that("the fitted correlation survives a .fitrx round-trip", {
  skip_on_cran()
  fit <- rho_case()$fit

  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)

  expect_equal(loaded$residual_correlations, fit$residual_correlations,
               tolerance = 1e-12)
})

test_that("reconstruction paths run against a free block_sigma fit", {
  skip_on_cran()
  case <- rho_case()

  pred <- ferx_predict(case$fit, model = case$model, data = case$data)
  expect_s3_class(pred, "data.frame")
  expect_true(all(is.finite(pred$PRED)))

  sim <- ferx_simulate(case$fit, model = case$model, data = case$data,
                       n_sim = 2, seed = 1)
  expect_s3_class(sim, "data.frame")
  expect_gt(nrow(sim), 0L)

  # ferx_covariance() refreshes the correlation SEs and labels the matrix the
  # same way the fit does.
  refreshed <- ferx_covariance(case$fit)
  expect_equal(refreshed$residual_correlations$rho,
               case$fit$residual_correlations$rho, tolerance = 1e-12)
  expect_equal(rownames(refreshed$cov_matrix), rownames(case$fit$cov_matrix))
})
