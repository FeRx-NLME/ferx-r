# Random-effect shape transforms from R (#806, ferx-core #1721 / #1716 part A):
# `[eta_shape]` and the inline boxcox() reach the fit, its print, the saved
# bundle and the written-back model. The shaped ETA is reported as `custom`
# with no CV% until core reports a shape-aware spread (ferx-core #1716 part B).

eta_shape_sections <- function(individual_parameters, eta_shape = NULL) {
  sections <- list(
    parameters = c(
      "  theta TVCL(0.134, 0.001, 10.0)",
      "  theta TVV(8.1, 0.1, 500.0)",
      "  theta TVKA(1.0, 0.01, 50.0)",
      "  omega ETA_CL ~ 0.07",
      "  omega ETA_V  ~ 0.02",
      "  omega ETA_KA ~ 0.40",
      "  sigma PROP_ERR ~ 0.01 (sd)"
    ),
    individual_parameters = individual_parameters
  )
  if (!is.null(eta_shape)) sections$eta_shape <- eta_shape
  c(sections, list(
    structural_model = "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    error_model = "  DV ~ proportional(PROP_ERR)"
  ))
}

unshaped_individual_parameters <- c(
  "  CL = TVCL * exp(ETA_CL)",
  "  V  = TVV  * exp(ETA_V)",
  "  KA = TVKA * exp(ETA_KA)"
)

eta_shape_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      model <- write_test_model(eta_shape_sections(
        unshaped_individual_parameters,
        eta_shape = "  ETA_CL ~ boxcox()"
      ))
      fit <<- ferx_fit(model, ferx_example("warfarin")$data, method = "focei",
                       verbose = FALSE, covariance = FALSE,
                       settings = list(maxiter = 30L))
    }
    fit
  }
})

test_that("[eta_shape] declares the shape theta and marks the ETA custom (#806)", {
  skip_on_cran()
  fit <- eta_shape_fit()
  expect_identical(names(fit$theta), c("TVCL", "TVV", "TVKA", "LAMBDA_ETA_CL"))
  expect_true(is.finite(fit$theta[["LAMBDA_ETA_CL"]]))
  expect_identical(unname(fit$eta_param_types),
                   c("custom", "log_normal", "log_normal"))
  scale <- fit$estimates$scale
  names(scale) <- fit$estimates$param
  expect_identical(unname(scale[c("ETA_CL", "ETA_V", "ETA_KA")]),
                   c("custom", "log_normal", "log_normal"))
})

test_that("print() shows the shaped ETA as [custom] with no CV% (#806)", {
  skip_on_cran()
  out <- capture.output(print(eta_shape_fit()))
  row_of <- function(eta) {
    row <- grep(paste0("^\\s+", eta, "\\s+\\["), out, value = TRUE)
    expect_length(row, 1L)
    row
  }
  shaped <- row_of("ETA_CL")
  expect_match(shaped, "[custom]", fixed = TRUE)
  expect_false(grepl("CV%", shaped, fixed = TRUE))
  expect_match(row_of("ETA_V"), "CV% =", fixed = TRUE)
})

test_that("a saved shaped fit loads with its shape theta and ETA types (#806)", {
  skip_on_cran()
  fit <- eta_shape_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$theta, fit$theta)
  expect_identical(loaded$eta_param_types, fit$eta_param_types)
  expect_identical(loaded$estimates$scale, fit$estimates$scale)
  expect_equal(loaded$ofv, fit$ofv)
})

test_that("ferx_update_model() writes the shape back and refits to the same OFV (#806)", {
  skip_on_cran()
  fit <- eta_shape_fit()
  out <- tempfile(fileext = ".ferx")
  on.exit(unlink(out), add = TRUE)
  ferx_update_model(fit, out)
  text <- readLines(out)
  expect_true(any(grepl("^\\s*theta LAMBDA_ETA_CL\\(", text)))
  expect_true(any(grepl("^\\s*ETA_CL ~ boxcox\\(LAMBDA_ETA_CL\\)", text)))

  refit <- ferx_fit(out, ferx_example("warfarin")$data, method = "focei",
                    verbose = FALSE, covariance = FALSE,
                    settings = list(maxiter = 0L))
  expect_equal(refit$theta, fit$theta, tolerance = 1e-10)
  expect_equal(refit$ofv, fit$ofv, tolerance = 1e-6)
})

test_that("inline boxcox() and [eta_shape] give the same objective (#806)", {
  skip_on_cran()
  fit <- eta_shape_fit()
  block_text <- strsplit(ferx_update_model(fit), "\n", fixed = TRUE)[[1]]

  # The same model at the same estimates, with the shape written inline
  # instead of in [eta_shape].
  header <- grep("^\\s*\\[", block_text)
  shape_start <- grep("^\\s*\\[eta_shape\\]", block_text)
  expect_length(shape_start, 1L)
  shape_end <- min(c(header[header > shape_start], length(block_text) + 1L)) - 1L
  inline_text <- block_text[-(shape_start:shape_end)]
  cl <- grep("^\\s*CL = TVCL \\* exp\\(ETA_CL\\)\\s*$", inline_text)
  expect_length(cl, 1L)
  inline_text[cl] <- "  CL = TVCL * exp(boxcox(ETA_CL, LAMBDA_ETA_CL))"

  evaluate <- function(lines) {
    path <- tempfile(fileext = ".ferx")
    on.exit(unlink(path), add = TRUE)
    writeLines(lines, path)
    ferx_fit(path, ferx_example("warfarin")$data, method = "focei",
             verbose = FALSE, covariance = FALSE,
             settings = list(maxiter = 0L))
  }
  block <- evaluate(block_text)
  inline <- evaluate(inline_text)
  expect_equal(inline$ofv, block$ofv, tolerance = 1e-8)
  expect_identical(inline$eta_param_types, block$eta_param_types)
  expect_identical(names(inline$theta), names(block$theta))
})
