# fit$ode_linear_exact and ode_method = "linear" (#578, ferx-core #2021).
#
# ferx-core solves a linear `[odes]` system exactly by the matrix exponential
# under the default `ode_method = "auto"`, and reports it on the fit. A named
# integrator opts a linear model out; a nonlinear system never qualifies.

linear_exact_fit <- function(name, ...) {
  ex <- ferx_example(name)
  # The overrides of the model files' [fit_options] each warn; not under test.
  suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                            verbose = FALSE,
                            settings = list(maxiter = 5L, ...)))
}

test_that("a linear [odes] model is solved exactly under ode_method = auto", {
  skip_on_cran()
  fit <- linear_exact_fit("one_cpt_iv_ode")
  expect_true(fit$uses_ode_solver)
  expect_identical(fit$ode_linear_exact, TRUE)

  out <- capture.output(print(fit))
  expect_true(any(grepl("ODE solution: exact (linear system, matrix exponential)",
                        out, fixed = TRUE)))
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("exact (linear system, matrix exponential)",
                        out, fixed = TRUE)))
})

test_that("ode_method = linear takes the exact path, rk45 integrates", {
  skip_on_cran()
  expect_identical(
    linear_exact_fit("one_cpt_iv_ode", ode_method = "linear")$ode_linear_exact,
    TRUE
  )
  fit <- linear_exact_fit("one_cpt_iv_ode", ode_method = "rk45")
  expect_true(fit$uses_ode_solver)
  expect_identical(fit$ode_linear_exact, FALSE)
  expect_false(any(grepl("matrix exponential", capture.output(print(fit)),
                         fixed = TRUE)))
})

test_that("a nonlinear or analytical model is not on the exact path", {
  skip_on_cran()
  mm <- linear_exact_fit("mm_oral")
  expect_true(mm$uses_ode_solver)
  expect_identical(mm$ode_linear_exact, FALSE)

  an <- linear_exact_fit("one_cpt_iv")
  expect_false(an$uses_ode_solver)
  expect_identical(an$ode_linear_exact, FALSE)
})

test_that("ode_method = linear on a nonlinear model is refused", {
  skip_on_cran()
  err <- expect_error(linear_exact_fit("mm_oral", ode_method = "linear"),
                      class = "ferx_engine_error")
  expect_identical(err$code, "E_ODE_METHOD_LINEAR")
  expect_match(conditionMessage(err), "not linear in the states", fixed = TRUE)

  # From the model file the parser refuses it.
  ex <- ferx_example("mm_oral")
  src <- sub("[fit_options]", "[fit_options]\n  ode_method = linear",
             paste(readLines(ex$model), collapse = "\n"), fixed = TRUE)
  model <- tempfile(fileext = ".ferx")
  writeLines(src, model)
  expect_error(ferx_fit(model, ex$data, verbose = FALSE),
               "not linear in the states", fixed = TRUE)
})

test_that("ode_linear_exact survives ferx_save_fit() / ferx_load_fit()", {
  skip_on_cran()
  fit <- linear_exact_fit("one_cpt_iv_ode")
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  expect_identical(fitrx_fit_json(path)$ode_linear_exact, TRUE)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$ode_linear_exact, TRUE)
  expect_identical(loaded$uses_ode_solver, TRUE)

  # A bundle written before the field existed has no key, and reads back
  # FALSE, as ferx-core's `#[serde(default)]` does.
  old <- fitrx_edit_bundle(path, function(w) {
    w$ode_linear_exact <- NULL
    w
  })
  expect_identical(ferx_load_fit(old)$ode_linear_exact, FALSE)
})
