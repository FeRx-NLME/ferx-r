# A residual `correlation cs|ar1|car1(RHO) by ...` line next to kappa and on
# nested levels (ferx-core #1082 PR5, #2102). Before it the engine refused
# both. The oracles are the uncorrelated fit of the same model and data (a
# free RHO with the uncorrelated fit nested at RHO = 0 cannot score worse),
# and the fit's own in-fit covariance step for ferx_covariance().

quiet <- function(x) suppressWarnings(suppressMessages(x))

# mbma_placebo with RHO in [parameters] and `line` under the `DV ~` line.
mbma_with_line <- function(line, rho = "  theta RHO(0.3, 0.0, 0.99)") {
  ex <- ferx_example("mbma_placebo")
  txt <- readLines(ex$model)
  at_theta <- grep("^  theta PLACEBO\\[", txt)
  at_dv <- grep("^  DV ~ additive", txt)
  stopifnot(length(at_theta) == 1L, length(at_dv) == 1L)
  txt <- append(txt, paste0("  ", line), after = at_dv)
  txt <- append(txt, rho, after = at_theta)
  model <- tempfile(fileext = ".ferx")
  writeLines(txt, model)
  list(model = model, data = ex$data)
}

mbma_car1_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      m <- mbma_with_line("correlation car1(RHO) by OCC")
      fit <<- quiet(ferx_fit(m$model, m$data, verbose = FALSE))
    }
    fit
  }
})

cov_parts <- function(fit) {
  lapply(fit[c("se_theta", "se_omega", "se_sigma", "se_kappa", "cov_matrix")], unname)
}

test_that("car1 by OCC fits next to mbma_placebo's arm kappa", {
  skip_on_cran()
  fit <- mbma_car1_fit()
  ex <- ferx_example("mbma_placebo")
  base <- quiet(ferx_fit(ex$model, ex$data, verbose = FALSE))

  expect_true("RHO" %in% names(fit$theta))
  rho <- fit$theta[["RHO"]]
  expect_true(is.finite(rho))
  expect_gte(rho, 0)
  expect_lt(rho, 0.99)
  # The kappa is still estimated alongside it.
  expect_true(all(is.finite(diag(fit$omega_iov))))
  expect_lte(fit$ofv, base$ofv)

  out <- capture.output(quiet(print(fit)))
  expect_true(any(grepl("^RHO\\s", out)))
  expect_true(any(grepl("Residual:    additive, correlation car1(RHO) by OCC",
                        out, fixed = TRUE)))
})

test_that("ferx_covariance() repeats the car1 + kappa fit's SEs, live and reloaded", {
  skip_on_cran()
  fit <- mbma_car1_fit()
  expect_identical(fit$covariance_status, "computed")
  expect_true(is.finite(fit$se_theta[["RHO"]]))
  expect_identical(cov_parts(quiet(ferx_covariance(fit))), cov_parts(fit))

  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$theta, fit$theta)
  expect_identical(cov_parts(quiet(ferx_covariance(back))), cov_parts(fit))
})

test_that("cs over the kappa's own occasions is refused as confounded", {
  skip_on_cran()
  m <- mbma_with_line("correlation cs(RHO) by OCC")
  e <- tryCatch(quiet(ferx_fit(m$model, m$data, verbose = FALSE)),
                error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_RESIDUAL_CORRELATION_CONFOUNDED_KAPPA")
  expect_match(conditionMessage(e), "KAPPA_ARM", fixed = TRUE)
})

test_that("ar1 by a nested level fits on the three-level fixture", {
  skip_on_cran()
  model <- test_path("fixtures", "nested_levels_3.ferx")
  data <- test_path("fixtures", "nested_levels_3.csv")
  txt <- readLines(model)
  txt <- append(txt, "  correlation ar1(RHO) by ARM", after = grep("^  DV ~", txt))
  txt <- append(txt, "  theta RHO(0.3, 0.0, 0.99)", after = grep("^  theta ET50", txt))
  corr_model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(txt, corr_model)

  fit <- quiet(ferx_fit(corr_model, data, verbose = FALSE))
  base <- quiet(ferx_fit(model, data, verbose = FALSE))
  expect_true(fit$converged)
  rho <- fit$theta[["RHO"]]
  expect_true(is.finite(rho))
  expect_gte(rho, 0)
  expect_lt(rho, 0.99)
  expect_lte(fit$ofv, base$ofv)
  expect_identical(fit$kappa_levels, c(K_ST = "STUDY", K_ARM = "ARM"))
  expect_identical(fit$model_structure$residual,
                   "additive, correlation ar1(RHO) by ARM")
})
