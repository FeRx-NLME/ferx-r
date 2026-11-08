# The omega/sigma SEs, the eta shrinkage and the omega parameter correlation
# carry the names `omega` and `sigma` carry (#722). `se_omega` and `se_sigma`
# are named exactly as their rows of `cov_matrix`, slot for slot, and keep the
# names through a ferx_save_fit()/ferx_load_fit() round trip.

.round_trip <- function(fit) {
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

.expect_field_names <- function(fit, se_omega_names) {
  expect_identical(names(fit$se_omega), se_omega_names)
  expect_identical(names(fit$se_sigma), fit$sigma_names)
  expect_identical(names(fit$shrinkage_eta), fit$eta_names)
  # By position: `cov_matrix` packs theta, then omega in `se_omega`'s order,
  # then sigma.
  rn <- rownames(fit$cov_matrix)
  n_th <- length(fit$theta)
  n_om <- length(fit$se_omega)
  expect_identical(rn[n_th + seq_len(n_om)], names(fit$se_omega))
  expect_identical(rn[n_th + n_om + seq_along(fit$se_sigma)], names(fit$se_sigma))
}

test_that("a diagonal fit names se_omega, se_sigma and shrinkage_eta, and a reload keeps them", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  expect_identical(fit$covariance_status, "computed")
  expect_identical(fit$eta_names, c("ETA_CL", "ETA_V", "ETA_KA"))
  .expect_field_names(fit, fit$eta_names)
  .expect_field_names(.round_trip(fit), fit$eta_names)
})

test_that("a block fit names its lower-triangle se_omega and omega_param_corr, and a reload keeps them", {
  skip_on_cran()
  ex <- ferx_example("warfarin_block_omega")
  fit <- ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                  covariance = TRUE, settings = list(maxiter = 30L))
  expect_identical(fit$covariance_status, "computed")
  lt <- c("ETA_CL", "ETA_V ~ ETA_CL", "ETA_KA ~ ETA_CL",
          "ETA_V", "ETA_KA ~ ETA_V", "ETA_KA")
  expect_identical(fit$eta_names, c("ETA_CL", "ETA_V", "ETA_KA"))
  .expect_field_names(fit, lt)
  expect_identical(dimnames(fit$omega_param_corr), dimnames(fit$omega))
  back <- .round_trip(fit)
  .expect_field_names(back, lt)
  expect_identical(dimnames(back$omega_param_corr), dimnames(fit$omega))
})
