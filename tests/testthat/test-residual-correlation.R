# Within-group residual correlation, `correlation cs|ar1(RHO) by COL` in
# [error_model] (ferx-core #1081). The `by` column is an ordinary data column
# no model expression reads, so these check that it reaches the engine from a
# CSV and from a data.frame (#419), and that RHO is reported as a theta.

# warfarin with RHO and a `correlation` line. `line` is the line as written;
# `above = TRUE` puts it before the `DV ~` line instead of after it.
.rescorr_case <- function(kind = "cs", line = NULL, above = FALSE,
                          rho = "theta RHO(0.2, -0.1, 0.95)") {
  if (is.null(line)) line <- sprintf("correlation %s(RHO) by GRP", kind)
  ex <- ferx_example("warfarin")
  mtxt <- readLines(ex$model)
  mtxt <- sub("theta TVKA(1.0, 0.01, 50.0)",
              paste0("theta TVKA(1.0, 0.01, 50.0)\n  ", rho),
              mtxt, fixed = TRUE)
  err <- "DV ~ proportional(PROP_ERR)"
  mtxt <- sub(err,
              if (above) paste0(line, "\n  ", err) else paste0(err, "\n  ", line),
              mtxt, fixed = TRUE)
  mtxt <- sub("covariance = true", "covariance = false", mtxt, fixed = TRUE)
  # Left to the default, so `settings = list(maxiter = 0L)` overrides nothing.
  mtxt <- mtxt[!grepl("^\\s*maxiter\\s*=", mtxt)]
  model <- tempfile(fileext = ".ferx")
  writeLines(mtxt, model)
  d <- utils::read.csv(ex$data, na.strings = ".")
  # Two blocks per subject: the absorption phase and the rest. GRP2 splits
  # them once more, for the two-column `by`.
  d$GRP <- ifelse(d$TIME < 24, 1L, 2L)
  d$GRP2 <- ifelse(d$TIME < 72, 1L, 2L)
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

test_that("RHO at its init scores a different objective than RHO FIX at 0", {
  skip_on_cran()
  # One evaluation each at the same estimates (maxiter = 0), so the gap is
  # the correlation's alone, not two optimizer runs' noise.
  f1 <- ferx_fit(.rescorr_case("cs", rho = "theta RHO(0.2, FIX)")$model,
                 .rescorr_case("cs")$data, settings = list(maxiter = 0L), verbose = FALSE)
  f0 <- ferx_fit(.rescorr_case("cs", rho = "theta RHO(0, FIX)")$model,
                 .rescorr_case("cs")$data, settings = list(maxiter = 0L), verbose = FALSE)
  expect_gt(abs(f1$ofv - f0$ofv), 1)
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

# The pre-fit structure is the engine's own (ferx-r #573 review), so every
# spelling labels the model as the fit does - including the line above `DV ~`,
# a RHO FIX at 0 (the engine drops the line), and a parenthesised `by`.
test_that("the pre-fit residual label is the fit's for every spelling", {
  skip_on_cran()
  cases <- list(
    "after DV"      = list(),
    "above DV"      = list(above = TRUE),
    "RHO FIX 0"     = list(rho = "theta RHO(0, FIX)"),
    "by (GRP)"      = list(line = "correlation cs(RHO) by (GRP)"),
    "by (GRP,GRP2)" = list(line = "correlation cs(RHO) by (GRP,GRP2)"),
    "upper case"    = list(line = "CORRELATION CS(RHO) BY GRP")
  )
  expected <- c(
    "proportional, correlation cs(RHO) by GRP",
    "proportional, correlation cs(RHO) by GRP",
    "proportional",
    "proportional, correlation cs(RHO) by GRP",
    "proportional, correlation cs(RHO) by (GRP, GRP2)",
    "proportional, correlation cs(RHO) by GRP"
  )
  for (i in seq_along(cases)) {
    case <- do.call(.rescorr_case, cases[[i]])
    pre <- ferx:::.ferx_parse_structure(case$model)
    fit <- ferx_fit(case$model, case$data, settings = list(maxiter = 0L), verbose = FALSE)
    what <- names(cases)[[i]]
    expect_identical(pre$residual, fit$model_structure$residual, info = what)
    expect_identical(pre$residual, expected[[i]], info = what)
    expect_identical(pre$theta_names, as.character(fit$model_structure$theta_names), info = what)
  }
})
