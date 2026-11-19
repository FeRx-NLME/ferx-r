# `mstep_optimizer` (ferx-core #3001): the numerical theta/sigma M-step
# optimizer of SAEM / IMP / IMPMAP, and `fit$mstep_optimizer_used`.

mstep_fit <- function(..., method = "saem") {
  skip_on_cran()
  ex <- ferx_example("warfarin_saem")
  ferx_fit(
    ex$model, ex$data,
    method = method, verbose = FALSE, covariance = FALSE,
    settings = c(list(n_exploration = 15L, n_convergence = 15L, seed = 42L),
                 list(...))
  )
}

test_that("an explicit mstep_optimizer is run and reported", {
  for (opt in c("bobyqa", "lbfgs", "newton")) {
    fit <- mstep_fit(method = "imp", mstep_optimizer = opt)
    expect_identical(fit$mstep_optimizer_used, opt)
    expect_true(is.finite(fit$ofv), info = opt)
  }
})

test_that("auto resolves to newton for a Gaussian model inside its scope", {
  expect_identical(mstep_fit(method = "imp")$mstep_optimizer_used, "newton")
  expect_identical(
    mstep_fit(method = "imp", mstep_optimizer = "auto")$mstep_optimizer_used,
    "newton"
  )
})

test_that("a closed-form SAEM M-step runs no numerical optimizer", {
  # warfarin_saem's M-step is closed form, so there is nothing to report.
  expect_null(mstep_fit()$mstep_optimizer_used)
})

test_that("an unknown mstep_optimizer is refused", {
  skip_on_cran()
  ex <- ferx_example("warfarin_saem")
  expect_error(
    ferx_fit(ex$model, ex$data, method = "saem", verbose = FALSE,
             settings = list(mstep_optimizer = "nelder")),
    "mstep_optimizer"
  )
})

test_that("score_sa and non-SAEM methods report no mstep_optimizer_used", {
  expect_null(mstep_fit(mstep_solver = "score_sa")$mstep_optimizer_used)
  expect_null(warfarin_fit()$mstep_optimizer_used)
})
