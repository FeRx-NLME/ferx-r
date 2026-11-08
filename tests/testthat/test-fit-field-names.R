# The omega/sigma SEs, the eta shrinkage and the omega parameter correlation
# carry the names `omega` and `sigma` carry (#722), spelled as `cov_matrix`'s
# labels are, and keep them through a ferx_save_fit()/ferx_load_fit() round
# trip.

.omega_se_names <- getFromNamespace(".ferx_omega_se_names", "ferx")

.round_trip <- function(fit) {
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

.expect_field_names <- function(fit, se_omega_names) {
  eta <- fit$eta_names
  expect_identical(names(fit$se_omega), se_omega_names)
  expect_identical(names(fit$se_sigma), fit$sigma_names)
  expect_identical(names(fit$shrinkage_eta), eta)
  # The SE names are the omega/sigma rows of the covariance matrix.
  expect_true(all(c(se_omega_names, fit$sigma_names) %in% rownames(fit$cov_matrix)))
}

test_that("se_omega labels follow cov_matrix's spelling, with OMEGA(i,j) fallbacks", {
  expect_identical(.omega_se_names(c("A", "B"), 2L, 2L), c("A", "B"))
  expect_identical(.omega_se_names(c("A", "B"), 2L, 3L), c("A", "B ~ A", "B"))
  expect_identical(.omega_se_names(NULL, 2L, 3L),
                   c("OMEGA(1,1)", "OMEGA(2,1)", "OMEGA(2,2)"))
  expect_identical(.omega_se_names(c("A", ""), 2L, 3L),
                   c("A", "OMEGA(2,1)", "OMEGA(2,2)"))
  expect_null(.omega_se_names(c("A", "B"), 2L, 4L))
  expect_null(.omega_se_names(character(0), 0L, 0L))
})

test_that("a diagonal fit names se_omega, se_sigma and shrinkage_eta, and a reload keeps them", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  expect_identical(fit$covariance_status, "computed")
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
