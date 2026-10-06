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
#
# `rho_fix_case()` and `rho_diag_case()` are the same model with the block
# held `FIX` and with two independent sigmas: the controls for #480, where a
# fit without its fitted rho is refused only when the model estimates one.
rho_fixture <- function(sigma_lines) {
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
      sigma_lines,
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
}
rho_case <- rho_fixture(
  "  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30]"
)
rho_fix_case <- rho_fixture(
  "  block_sigma (PROP_ERR, ADD_ERR) = [0.04, 0.10, 0.30] FIX"
)
rho_diag_case <- rho_fixture(c(
  "  sigma PROP_ERR ~ 0.04",
  "  sigma ADD_ERR ~ 0.30"
))

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
  # The declared init is 0.10 / sqrt(0.04 * 0.30) = 0.913; a fitted value
  # equal to it would mean the estimate never reached R.
  expect_true(is.finite(rc$rho))
  expect_true(abs(rc$rho) < 1)
  expect_false(isTRUE(all.equal(rc$rho, 0.10 / sqrt(0.04 * 0.30))))
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

# -- #480: a free correlation without its fitted value is refused --
#
# A from-fit path rebuilds the model at the fit's own correlations. When the
# model estimates one and the fit carries no value for it - a fit made before
# the field existed (ferx 0.3.x), or one whose field was dropped - the engine
# used to fill in the declared initial value, which moved ferx_covariance()'s
# se_theta by 23 % on this fixture with no error. Every from-fit binding now
# refuses, with the prefix it already uses for its other fit errors.

# The fit as an older ferx (or a caller's list surgery) leaves it.
without_rho <- function(fit) {
  fit$residual_correlations <- NULL
  fit
}

declared_rho <- 0.10 / sqrt(0.04 * 0.30)

test_that("covariance and SIR refuse a free fit without its fitted rho", {
  skip_on_cran()
  fit <- without_rho(rho_case()$fit)

  expect_error(ferx_covariance(fit), "^ferx_covariance: the model estimates")
  expect_error(ferx_sir(fit, sir_samples = 50L, sir_resamples = 20L,
                        sir_seed = 7L),
               "^ferx_sir: the model estimates")
})

test_that("predict, simulate and NPDE refuse a free fit without its fitted rho", {
  skip_on_cran()
  case <- rho_case()
  fit <- without_rho(case$fit)

  # predict and NPDE do not read rho numerically, so the refusal is the only
  # thing that shows they go through the same check.
  expect_error(ferx_predict(case$model, case$data, fit = fit),
               "^Fit error: the model estimates")
  expect_error(ferx_simulate(case$model, case$data, n_sim = 1L, seed = 1L,
                             fit = fit),
               "^Fit error: the model estimates")
  expect_error(ferx_calc_npde(fit, nsim = 20L, seed = 1L,
                              model = case$model, data = case$data),
               "^Fit error: the model estimates")
})

test_that("simulate_with_uncertainty refuses a free fit without its fitted rho", {
  skip_on_cran()
  case <- rho_case()
  fit <- without_rho(case$fit)

  # The refusal fires while the fit is rebuilt, before any draw is sampled.
  expect_error(ferx_simulate_with_uncertainty(case$model, case$data, fit,
                                              n_uncertainty_draws = 5L),
               "^Uncertainty error: the model estimates")
})

test_that("a zero-row residual_correlations frame is refused like NULL", {
  skip_on_cran()
  fit <- rho_case()$fit
  fit$residual_correlations <- fit$residual_correlations[0, ]

  expect_error(ferx_covariance(fit), "^ferx_covariance: the model estimates")
})

test_that("an intact free fit is rebuilt at its fitted rho", {
  skip_on_cran()
  case <- rho_case()
  fit <- case$fit

  # The inline and the standalone covariance step see the same parameters
  # only if the fitted rho reaches the engine: at the declared 0.913 the
  # standalone se_theta moves by 23 %.
  expect_identical(ferx_covariance(fit)$se_theta, fit$se_theta)

  # simulate reads rho through the residual draw, so moving it to the
  # declared value has to change DV_SIM.
  at_declared <- fit
  at_declared$residual_correlations$rho <- declared_rho
  sim_fit <- ferx_simulate(case$model, case$data, n_sim = 2L, seed = 1L,
                           fit = fit)
  sim_declared <- ferx_simulate(case$model, case$data, n_sim = 2L, seed = 1L,
                                fit = at_declared)
  expect_false(identical(sim_fit$DV_SIM, sim_declared$DV_SIM))
})

test_that("a FIX correlation needs no fitted value", {
  skip_on_cran()
  case <- rho_fix_case()
  fit <- case$fit
  stripped <- without_rho(fit)
  expect_true(all(fit$residual_correlations$fixed))

  expect_identical(ferx_covariance(stripped)$se_theta,
                   ferx_covariance(fit)$se_theta)
  expect_identical(
    ferx_simulate(case$model, case$data, n_sim = 2L, seed = 1L,
                  fit = stripped)$DV_SIM,
    ferx_simulate(case$model, case$data, n_sim = 2L, seed = 1L,
                  fit = fit)$DV_SIM
  )
  expect_identical(ferx_predict(case$model, case$data, fit = stripped)$PRED,
                   ferx_predict(case$model, case$data, fit = fit)$PRED)
})

test_that("a diagonal-sigma fit runs without a residual_correlations field", {
  skip_on_cran()
  case <- rho_diag_case()
  fit <- case$fit
  expect_null(fit$residual_correlations)

  expect_length(ferx_covariance(fit)$se_theta, length(fit$theta))
  sim <- ferx_simulate(case$model, case$data, n_sim = 1L, seed = 1L, fit = fit)
  expect_true(all(is.finite(sim$DV_SIM)))
})

test_that("a non-finite or out-of-range fitted rho is refused", {
  skip_on_cran()
  case <- rho_case()
  label <- case$fit$residual_correlations$name

  na_fit <- case$fit
  na_fit$residual_correlations$rho <- NA_real_
  msg <- tryCatch(ferx_covariance(na_fit), error = conditionMessage)
  expect_match(msg, "^ferx_covariance: the block_sigma correlation ")
  expect_match(msg, label, fixed = TRUE)
  expect_match(msg, "must be finite and strictly between -1 and 1",
               fixed = TRUE)
  expect_match(msg, "fit$residual_correlations$rho", fixed = TRUE)
  expect_no_match(msg, "0.3.x", fixed = TRUE)

  one_fit <- case$fit
  one_fit$residual_correlations$rho <- 1
  msg <- tryCatch(
    ferx_simulate(case$model, case$data, n_sim = 1L, seed = 1L,
                  fit = one_fit),
    error = conditionMessage
  )
  expect_match(msg, "^Fit error: the block_sigma correlation ")
  expect_match(msg, "rho = 1 ", fixed = TRUE)
  expect_match(msg, "must be finite and strictly between -1 and 1",
               fixed = TRUE)
})

test_that("the refusal names the correlation and the remedy", {
  skip_on_cran()
  fit <- rho_case()$fit
  label <- fit$residual_correlations$name
  msg <- tryCatch(ferx_covariance(without_rho(fit)), error = conditionMessage)

  expect_match(msg, "^ferx_covariance: ")
  expect_match(
    msg,
    paste0("the model estimates the block_sigma correlation ", label, ","),
    fixed = TRUE
  )
  expect_match(
    msg,
    paste("but the fit carries no value for it",
          "(fit$residual_correlations is missing or empty), so it would be",
          "rebuilt at the declared initial correlation."),
    fixed = TRUE
  )
  expect_match(
    msg,
    paste("A fit made with ferx 0.3.x or earlier held this correlation at",
          "its declared value and predates the field."),
    fixed = TRUE
  )
  expect_match(msg, "Re-fit via ferx_fit(model, data).", fixed = TRUE)
  # Editing the model to FIX is no remedy: covariance and SIR check the
  # model hash against the fit.
  expect_no_match(msg, "FIX", fixed = TRUE)
})

test_that("a v0.3.x-shaped fit gets the rho refusal from SIR, not a dimension error", {
  skip_on_cran()
  fit <- without_rho(rho_case()$fit)
  skip_if(is.null(fit$cov_matrix), "covariance step produced no matrix")
  # A v0.3.x fit has no rho coordinate in its covariance matrix either.
  n <- nrow(fit$cov_matrix) - 1L
  fit$cov_matrix <- fit$cov_matrix[seq_len(n), seq_len(n)]

  expect_error(ferx_sir(fit, sir_samples = 50L, sir_resamples = 20L,
                        sir_seed = 7L),
               "^ferx_sir: .*estimates the block_sigma correlation")
})
