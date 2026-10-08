# A non-Gaussian endpoint's per-subject score (ferx-core #1744).
#
# The score feeds `covariance_method = "s"` / `"rsr"` and the gradient-based
# outer optimizers. Before #1744 it covered only the Gaussian observations, so
# on the bundled exponential TTE model - which keeps placeholder PK blocks - it
# was zero without random effects: `rsr` reported SE = 0, and `lbfgs` stayed at
# the initial estimates and called that convergence. With random effects `rsr`
# reported an SE that looked plausible and was a third of the right one.
# Reddened by the pin before #1744 (e2f9f641).

tte_model_without_eta <- function(env = parent.frame()) {
  text <- readLines(ferx_example("tte_exponential")$model)
  text <- text[!grepl("omega ETA_LAMBDA", text, fixed = TRUE)]
  text <- gsub(" * exp(ETA_LAMBDA)", "", text, fixed = TRUE)
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(text, path)
  path
}

tte_fit <- function(model, ...) {
  suppressWarnings(ferx_fit(model, ferx_example("tte_exponential")$data,
                            verbose = FALSE, ...))
}

test_that("S1: covariance_method = rsr gives the TTE hazard a nonzero SE", {
  ex <- ferx_example("tte_exponential")
  for (model in c(ex$model, tte_model_without_eta())) {
    fit <- tte_fit(model, covariance = TRUE,
                   settings = list(covariance_method = "rsr"))
    expect_identical(fit$covariance_status, "computed")
    expect_gt(fit$se_theta[["TVLAMBDA"]], 1e-3)
  }
  # With ETA_LAMBDA: the RSR SE core measured on this model (0.0086, #1744),
  # not the 0.0028 it reported before.
  fit <- tte_fit(ex$model, covariance = TRUE,
                 settings = list(covariance_method = "rsr"))
  expect_equal(fit$se_theta[["TVLAMBDA"]], 0.0086, tolerance = 0.02)
})

test_that("S2: optimizer = lbfgs moves a TTE fit off its initial estimates", {
  model <- tte_model_without_eta()
  init <- 0.05
  lbfgs <- tte_fit(model, covariance = FALSE, settings = list(optimizer = "lbfgs"))
  bobyqa <- tte_fit(model, covariance = FALSE)
  expect_gt(abs(lbfgs$theta[["TVLAMBDA"]] - init), 0.01)
  # Without random effects it reaches BOBYQA's optimum (#1744).
  expect_equal(lbfgs$theta[["TVLAMBDA"]], bobyqa$theta[["TVLAMBDA"]], tolerance = 1e-4)
  expect_equal(lbfgs$ofv, bobyqa$ofv, tolerance = 1e-6)
})
