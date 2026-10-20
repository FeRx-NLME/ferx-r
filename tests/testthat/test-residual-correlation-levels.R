# A residual `correlation cs|ar1|car1(RHO) by ...` line next to kappa and on
# nested levels (ferx-core #1082 PR5, #2102). Before it the engine refused
# both. The oracles are the uncorrelated fit of the same model and data, and
# the fit's own in-fit covariance step for ferx_covariance(). The uncorrelated
# model is the correlated one at RHO = 0, so a converged base bounds the
# correlated OFV from above; the nested_levels_3 base stops unconverged, so
# there it is a stopping point, not an optimum. An engine that parsed the line
# but left it out of the likelihood would tie the base OFV and leave RHO at
# its start (measured on mbma_placebo: OFV 97.0700 both, RHO 0.3, SE 0), so
# each fit must beat the base by a margin and move RHO.
#
# With the line next to kappa the gradients used to be finite differences, and
# on that route `optimizer = "auto"` picked BOBYQA, which stopped on its own
# tolerance short of the optimum at a point that moved with the platform and
# the start (mbma_placebo car1: 89.82 on macOS, 99.18 on Linux from RHO = 0.3;
# 100.01 / 118.43 from 0.01), each reported as converged. Since ferx-core
# #2127 (#2122) `auto` takes L-BFGS on this route, and since ferx-core #2093
# (#2101) its outer gradient and the covariance step's R-matrix are exact.
# With an explicit nlopt_lbfgs it reached 89.451278 (RHO 0.7598) from RHO =
# 0.01 and 0.3 on macOS and Linux (#592); under `auto` it was measured from
# 0.01 and 0.6 on macOS, on finite differences and again on the exact gradient
# (89.4512784 both, 34 s -> 0.3 s per fit). On nested_levels_3 ar1 it is
# 105.639739 (RHO 0.11512). The fits below run on the default optimizer, and
# the optima are pinned to those values.

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

car1_rho <- function(start) sprintf("  theta RHO(%s, 0.0, 0.95)", start)

mbma_car1_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      m <- mbma_with_line("correlation car1(RHO) by OCC", car1_rho("0.6"))
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
  expect_lt(rho, 0.95)
  # The kappa is still estimated alongside it.
  expect_true(all(is.finite(diag(fit$omega_iov))))
  # Measured gap 7.62 (97.070036 -> 89.451278).
  expect_lt(fit$ofv, base$ofv - 1)
  expect_gt(abs(rho - 0.6), 1e-3)
  expect_lt(abs(fit$ofv - 89.451278), 1e-4)
  expect_lt(abs(rho - 0.7598), 1e-3)

  out <- capture.output(quiet(print(fit)))
  expect_true(any(grepl("^RHO\\s", out)))
  expect_true(any(grepl("Residual:    additive, correlation car1(RHO) by OCC",
                        out, fixed = TRUE)))
})

test_that("auto takes L-BFGS on car1 + kappa and reaches one optimum from every start", {
  skip_on_cran()
  # ferx-core #2127 (#2122). From 0.01 the old BOBYQA default stopped at 100.01
  # on macOS and 118.43 on Linux, both above the kappa-only 97.07; from 0.6 at
  # 98.31 on Linux. L-BFGS lands within 1.2e-7 of 89.4512784 from both.
  ex <- ferx_example("mbma_placebo")
  base <- quiet(ferx_fit(ex$model, ex$data, verbose = FALSE))
  m <- mbma_with_line("correlation car1(RHO) by OCC", car1_rho("0.01"))
  low <- quiet(ferx_fit(m$model, m$data, verbose = FALSE))
  high <- mbma_car1_fit()
  for (fit in list(low, high)) {
    expect_true(fit$converged)
    expect_identical(fit$optimizer, "auto (nlopt_lbfgs)")
    expect_lte(fit$ofv, base$ofv)
  }
  expect_lt(abs(low$ofv - high$ofv), 1e-5)
  expect_lt(abs(low$theta[["RHO"]] - high$theta[["RHO"]]), 1e-3)

  # ferx-core #2093 (#2101): the route takes the exact adjoint gradient. Before
  # it the outer gradient was finite differences and the fit warned twice: the
  # correlation line had no analytic gradient, and the model's 25 theta + eta +
  # kappa axes were past the analytic path's 24-axis cap. The adjoint has no
  # theta cap, so neither warning applies.
  for (fit in list(low, high)) {
    expect_identical(fit$gradient_method_outer, "analytic (Dual2)")
    w <- fit$warnings
    expect_false(any(grepl("together with kappa", w, fixed = TRUE)))
    expect_false(any(grepl("24-axis cap", w, fixed = TRUE)))
  }
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
  expect_gt(abs(rho - 0.3), 1e-3)
  # Measured gap 1.40 (107.043571 -> 105.639739), against an unconverged base.
  expect_lt(fit$ofv, base$ofv - 0.5)
  expect_lt(abs(fit$ofv - 105.639739), 1e-4)
  expect_lt(abs(rho - 0.11512), 1e-3)
  expect_identical(fit$kappa_levels, c(K_ST = "STUDY", K_ARM = "ARM"))
  expect_identical(fit$model_structure$residual,
                   "additive, correlation ar1(RHO) by ARM")
})

test_that("ferx_simulate() and ferx_calc_npde() run on the car1 + kappa fit", {
  skip_on_cran()
  fit <- mbma_car1_fit()
  m <- mbma_with_line("correlation car1(RHO) by OCC", car1_rho("0.6"))
  sim <- quiet(ferx_simulate(m$model, fit = fit, n_sim = 2L))
  expect_s3_class(sim, "data.frame")
  expect_true(all(is.finite(sim$DV_SIM)))
  expect_no_error(quiet(ferx_calc_npde(fit, nsim = 50L)))
})
