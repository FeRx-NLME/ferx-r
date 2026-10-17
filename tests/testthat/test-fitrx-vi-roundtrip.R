# A VI fit through ferx_save_fit() / ferx_load_fit() (#508). VI reports a
# bound, not an objective, so under the default `vi_final_ofv = none` its
# OFV, AIC and BIC are NaN by design; the bundle has to bring them back as
# NaN, and the method as the engine's "VI" label. `fit$vi` itself (#507) is
# not on the fit yet, so its round trip is not tested here.

vi_fit <- function(final_ofv = NULL) {
  ex <- ferx_example("warfarin")
  src <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", src)
  src <- c(
    src[seq_len(at - 1L)],
    "[fit_options]",
    "  method     = vi",
    "  vi_iters   = 30",
    "  vi_seed    = 42",
    if (!is.null(final_ofv)) paste0("  vi_final_ofv = ", final_ofv),
    "  covariance = false"
  )
  model <- tempfile(fileext = ".ferx")
  writeLines(src, model)
  suppressWarnings(ferx_fit(model, ex$data, verbose = FALSE))
}

test_that("a VI fit reloads with method \"VI\" and its NaN OFV as NaN (#508)", {
  fit <- vi_fit()
  # The premise: a live VI fit is labelled by the engine and has no OFV.
  expect_identical(fit$method, "VI")
  expect_true(is.nan(fit$ofv))

  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)

  # NaN is written as JSON null (JSON has no NaN), not as a string.
  wire <- fitrx_fit_json(path)
  expect_identical(wire$method, "vi")
  expect_null(wire$ofv)

  back <- ferx_load_fit(path)
  expect_identical(back$method, "VI")
  expect_identical(back$method_chain, fit$method_chain)
  for (field in c("ofv", "ofv_data", "aic", "bic")) {
    expect_length(back[[field]], 1L)
    expect_identical(is.nan(back[[field]]), is.nan(fit[[field]]), info = field)
  }
  expect_true(is.nan(back$ofv))

  # What a VI fit does carry today: the estimates, and the posterior means as
  # the EBEs.
  expect_equal(back$theta, fit$theta, tolerance = 1e-12)
  expect_equal(back$omega, fit$omega, tolerance = 1e-12)
  # (A reloaded sigma carries names a live one does not; not VI-specific.)
  expect_equal(unname(back$sigma), unname(fit$sigma), tolerance = 1e-12)
  eta_cols <- setdiff(names(fit$ebe_etas), "ID")
  expect_equal(back$ebe_etas[eta_cols], fit$ebe_etas[eta_cols], tolerance = 1e-12)

  # And it survives a second trip unchanged.
  path2 <- tempfile(fileext = ".fitrx")
  ferx_save_fit(back, path2)
  again <- ferx_load_fit(path2)
  expect_identical(again$method, "VI")
  expect_true(is.nan(again$ofv))
})

test_that("a VI fit with a Laplace OFV keeps the finite OFV (#508)", {
  fit <- vi_fit("laplace")
  expect_identical(fit$method, "VI")
  expect_true(is.finite(fit$ofv))
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  back <- ferx_load_fit(path)
  expect_identical(back$method, "VI")
  expect_identical(back$ofv, fit$ofv)
  expect_identical(back$aic, fit$aic)
})
