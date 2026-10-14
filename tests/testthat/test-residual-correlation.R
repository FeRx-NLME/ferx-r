# Within-group residual correlation, `correlation cs|ar1(RHO) by COL` in
# [error_model] (ferx-core #1081). The `by` column is an ordinary data column
# no model expression reads, so these check that it reaches the engine from a
# CSV and from a data.frame (#419), and that RHO is reported as a theta.

.rescorr_case <- function(kind = "cs") {
  ex <- ferx_example("warfarin")
  mtxt <- readLines(ex$model)
  mtxt <- sub("theta TVKA(1.0, 0.01, 50.0)",
              "theta TVKA(1.0, 0.01, 50.0)\n  theta RHO(0.2, -0.1, 0.95)",
              mtxt, fixed = TRUE)
  mtxt <- sub("DV ~ proportional(PROP_ERR)",
              sprintf("DV ~ proportional(PROP_ERR)\n  correlation %s(RHO) by GRP", kind),
              mtxt, fixed = TRUE)
  mtxt <- sub("covariance = true", "covariance = false", mtxt, fixed = TRUE)
  model <- tempfile(fileext = ".ferx")
  writeLines(mtxt, model)
  d <- utils::read.csv(ex$data, na.strings = ".")
  # Two blocks per subject: the absorption phase and the rest.
  d$GRP <- ifelse(d$TIME < 24, 1L, 2L)
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(d, csv, row.names = FALSE, na = ".")
  list(model = model, data = d, csv = csv)
}

test_that("a cs residual correlation fits from a CSV and from a data.frame alike", {
  skip_on_cran()
  cs <- .rescorr_case("cs")
  f_csv <- ferx_fit(cs$model, cs$csv, verbose = FALSE)
  f_df  <- ferx_fit(cs$model, cs$data, verbose = FALSE)

  expect_true("RHO" %in% names(f_csv$theta))
  expect_true(is.finite(f_csv$theta[["RHO"]]))
  # RHO is free and the blocks hold several rows, so it moves off its init.
  expect_false(isTRUE(all.equal(f_csv$theta[["RHO"]], 0.2)))
  expect_identical(f_df$theta, f_csv$theta)
  expect_identical(f_df$ofv, f_csv$ofv)

  expect_true("RHO" %in% rownames(f_csv$estimates))
  expect_equal(f_csv$estimates["RHO", "estimate"], f_csv$theta[["RHO"]],
               ignore_attr = TRUE)

  out <- capture.output(suppressWarnings(print(f_csv)))
  expect_true(any(grepl("^RHO\\s", out)))
  expect_true(any(grepl("Residual:    proportional, correlation cs(RHO) by GRP",
                        out, fixed = TRUE)))
})

test_that("the residual correlation changes the objective from the independent model", {
  skip_on_cran()
  cs <- .rescorr_case("cs")
  fixed0 <- readLines(cs$model)
  fixed0 <- sub("theta RHO(0.2, -0.1, 0.95)", "theta RHO(0, FIX)", fixed0, fixed = TRUE)
  m0 <- tempfile(fileext = ".ferx")
  writeLines(fixed0, m0)
  f1 <- ferx_fit(cs$model, cs$data, verbose = FALSE)
  f0 <- ferx_fit(m0, cs$data, verbose = FALSE)
  expect_false(isTRUE(all.equal(f1$ofv, f0$ofv)))
})

test_that("an ar1 residual correlation fits from a data.frame", {
  skip_on_cran()
  ar <- .rescorr_case("ar1")
  f <- ferx_fit(ar$model, ar$data, verbose = FALSE)
  expect_true(is.finite(f$theta[["RHO"]]))
  expect_identical(f$model_structure$residual,
                   "proportional, correlation ar1(RHO) by GRP")
})

test_that("a data.frame without the by column is refused as a missing column", {
  skip_on_cran()
  cs <- .rescorr_case("cs")
  d0 <- cs$data
  d0$GRP <- NULL
  expect_error(ferx_fit(cs$model, d0, verbose = FALSE), "GRP")
})

test_that("the pre-fit structure names the residual correlation as the fit does", {
  cs <- .rescorr_case("cs")
  ms <- ferx:::.ferx_parse_structure(cs$model)
  expect_identical(ms$residual, "proportional, correlation cs(RHO) by GRP")
})
