# Allometric scaling from R.
#
# The convention (0.75 on a clearance, 1.0 on a volume, `power(center = 70)`)
# is the engine's; what is asserted here is that the R surface exposes it both
# as a transform (`fit = FALSE`, a model you can go on to edit) and as a tool
# (`fit = TRUE`, two fits side by side), and that the numbers it reports are the
# engine's own.

test_that("the entry forms are mutually exclusive", {
  cfg <- tempfile(fileext = ".ferxsearch")
  writeLines(c('base = "model.ferx"', "[space]", 'mfl = "ALLOMETRY(WT, 70)"'), cfg)

  expect_error(ferx_allometry(config = cfg, covariate = "WT"), "covariate")
  expect_error(ferx_allometry(config = cfg, model = "model.ferx"), "model")
  expect_error(ferx_allometry(), "`config`.*or `model`")
})

test_that("argument validation happens in R", {
  ex <- ferx_example("two_cpt_oral_base")

  expect_error(ferx_allometry(ex$model, ex$data, reference = 0), "positive")
  expect_error(ferx_allometry(ex$model, ex$data, parameters = 1), "character vector")
  expect_error(ferx_allometry(ex$model, ex$data, exponents = "0.75"), "numeric vector")
  expect_error(ferx_allometry(ex$model, ex$data, estimate = NA), "TRUE or FALSE")
  # The engine owns the pairing rule, and says so by name.
  expect_error(
    ferx_allometry(ex$model, ex$data, parameters = c("CL", "V"), exponents = 0.75,
                   fit = FALSE),
    "exponent"
  )
})

test_that("fit = FALSE is a model transform, not a run", {
  ex <- ferx_example("two_cpt_oral_base")
  res <- ferx_allometry(ex$model, ex$data, fit = FALSE)

  expect_s3_class(res, "ferx_allometry")
  expect_false(res$fitted)
  expect_null(res$fit)

  # Every scaled parameter is a row, with the convention's exponent.
  expect_true(nrow(res$scalings) > 0L)
  expect_true(all(res$scalings$fixed))
  expect_true(all(res$scalings$exponent %in% c(0.75, 1)))

  # The transform's product is a model you can go on to fit or edit.
  expect_true(file.exists(res$model_path))
  expect_match(res$model, "[covariate_model]", fixed = TRUE)
  expect_match(res$model, "power", fixed = TRUE)
  expect_equal(res$covariate, "WT")
  expect_equal(res$reference, 70)
})

test_that("estimate = TRUE declares a theta instead of fixing the exponent", {
  ex <- ferx_example("two_cpt_oral_base")
  res <- ferx_allometry(ex$model, ex$data, estimate = TRUE, fit = FALSE)

  expect_false(any(res$scalings$fixed))
  expect_true(all(!is.na(res$scalings$theta)))
})

test_that("the covariate and reference reach the relations", {
  ex <- ferx_example("two_cpt_oral_base")
  res <- ferx_allometry(ex$model, ex$data, covariate = "WT", reference = 55,
                        fit = FALSE)
  expect_equal(res$reference, 55)
  expect_match(res$model, "55", fixed = TRUE)
})

test_that("print() shows the scaling and says nothing was fitted", {
  ex <- ferx_example("two_cpt_oral_base")
  res <- ferx_allometry(ex$model, ex$data, fit = FALSE)
  out <- paste(capture.output(print(res)), collapse = "\n")

  expect_match(out, "ferx allometric scaling")
  expect_match(out, "power(center = 70", fixed = TRUE)
  expect_match(out, "Not fitted", fixed = TRUE)
})

test_that("fit = TRUE fits the base and the scaled model side by side", {
  skip_on_cran()
  ex <- ferx_example("two_cpt_oral_base")
  res <- ferx_allometry(ex$model, ex$data, retries = 0,
                        directory = file.path(tempdir(), "allometry-run"))

  expect_true(res$fitted)
  expect_equal(res$comparison$model, c("base", "scaled"))
  expect_true(all(is.finite(res$comparison$ofv)))
  # The verdict and the termination status are columns, not a filter.
  expect_true(is.logical(res$comparison$converged))
  expect_true(is.logical(res$comparison$passed))

  # dOFV is the engine's own subtraction, reported rather than recomputed.
  expect_equal(res$dofv, res$comparison$ofv[1] - res$comparison$ofv[2],
               tolerance = 1e-8)

  expect_s3_class(res$fit, "ferx_fit")
  expect_s3_class(res$base_fit, "ferx_fit")
  expect_equal(unname(res$fit$ofv), res$comparison$ofv[2], tolerance = 1e-8)
})
