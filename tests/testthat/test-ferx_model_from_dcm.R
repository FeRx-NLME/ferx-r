# ferx_model_from_dcm() (#344): a fitted DCM written out as a plain model.
#
# Two halves. The fixture tests run in every build: the export never evaluates
# a network, so a build without the `nn` cargo feature - which cannot parse a
# [covariate_nn] block, let alone fit one - still exports. The fit tests need
# `nn` to fit the DCM, and skip, by name, without it.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ferx-core's own export fixture (tests/data/dcm_export at the pinned
# revision): a 2 -> 3 -> 2 -> 2 network, a theta vector, and the export the
# engine commits for them. On a pin bump that changes the export's text, copy
# model.ferx and expected_analytic.ferx again from the new revision.
dcm_fixture <- function(name) test_path("fixtures", "dcm_export", name)

# theta.csv is core's, plus a `hex` column holding each value's exact bits
# (Python's float.hex() of the decimal). The theta is read from `hex`: R's
# decimal reader is not correctly rounded, and reads three of these
# shortest-round-trip literals one ULP off - which the byte comparison below
# then reports as a different export, correctly.
dcm_fixture_theta <- function() {
  tab <- utils::read.csv(dcm_fixture("theta.csv"), stringsAsFactors = FALSE,
                         colClasses = "character")
  stats::setNames(as.numeric(tab$hex), tab$name)
}

# The warfarin DCM fitted for two outer iterations: the weights need not be
# converged, only fitted. NULL in a build without `nn`, where the model does
# not parse.
fit_warfarin_dcm <- function() {
  ex <- ferx_example("warfarin_dcm")
  fit <- NULL
  utils::capture.output(fit <- tryCatch(
    suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                              settings = list(maxiter = 2L))),
    error = function(e) NULL
  ))
  if (is.null(fit) || length(fit$neural_networks) == 0L) return(NULL)
  fit
}

skip_without_nn <- function(fit) {
  skip_if(is.null(fit), "ferx built without the `nn` feature: cannot fit a [covariate_nn] model")
}

# `theta NAME(value, ...)` declarations of a model file: each value as the
# literal written (a character vector - see dcm_fixture_theta() on why not
# as.numeric()), and whether each line ends in FIX.
declared_thetas <- function(path) {
  lines <- readLines(path, warn = FALSE)
  m <- regmatches(lines, regexec("^\\s*theta\\s+(\\w+)\\s*\\(\\s*([^,)]+)", lines))
  hit <- lengths(m) == 3L
  list(
    value = stats::setNames(trimws(vapply(m[hit], `[`, "", 3L)),
                            vapply(m[hit], `[`, "", 2L)),
    fix = stats::setNames(grepl("\\bFIX\\s*(#.*)?$", lines[hit]),
                          vapply(m[hit], `[`, "", 2L))
  )
}

# ---------------------------------------------------------------------------
# Every build: export from a committed model + theta vector, no fit
# ---------------------------------------------------------------------------

test_that("the fixture exports to ferx-core's committed text, byte for byte", {
  out <- tempfile(fileext = ".ferx")
  res <- ferx_model_from_dcm(dcm_fixture("model.ferx"), theta = dcm_fixture_theta(),
                             output_model = out)
  expect_s3_class(res, "ferx_model")
  expect_identical(normalizePath(res$model), normalizePath(out))
  # Bytes, not lines: a weight literal one ULP off, or one dropped from the
  # theta vector on its way to the engine, is a different file.
  expect_identical(
    readBin(out, "raw", file.size(out)),
    readBin(dcm_fixture("expected_analytic.ferx"), "raw",
            file.size(dcm_fixture("expected_analytic.ferx")))
  )
})

test_that("the export has no [covariate_nn] block and parses in this build", {
  out <- tempfile(fileext = ".ferx")
  ferx_model_from_dcm(dcm_fixture("model.ferx"), theta = dcm_fixture_theta(),
                      output_model = out)
  expect_false(any(grepl("^\\s*\\[covariate_nn", readLines(out, warn = FALSE))))
  # Parsed by the engine - in a build without `nn` too.
  utils::capture.output(info <- ferx_model_inspect(out))
  expect_identical(info$theta_names, "TVV")
  expect_identical(info$iiv, c("ETA_CL", "ETA_V", "ETA_KA"))
})

test_that("weights = 'fixed_theta' declares every weight as a FIXed theta at its value", {
  theta <- dcm_fixture_theta()
  weights <- theta[grepl("^[WB]_TYPICAL_PK_", names(theta))]
  out <- tempfile(fileext = ".ferx")
  ferx_model_from_dcm(dcm_fixture("model.ferx"), theta = theta,
                      output_model = out, weights = "fixed_theta")

  # The engine reads every weight back as a theta of the exported model ...
  utils::capture.output(info <- ferx_model_inspect(out))
  expect_setequal(info$theta_names, c("TVV", names(weights)))
  # ... declared FIX at the value it was given, to the bit: the engine writes
  # the shortest literal that reads back as the same double, which is unique,
  # so it must be ferx-core's own literal for that weight.
  tab <- utils::read.csv(dcm_fixture("theta.csv"), colClasses = "character")
  decl <- declared_thetas(out)
  expect_identical(decl$value[names(weights)],
                   stats::setNames(tab$value, tab$name)[names(weights)])
  expect_true(all(decl$fix[names(weights)]))
  expect_false(decl$fix[["TVV"]])
  # Literal is the default and declares none of them.
  lit <- tempfile(fileext = ".ferx")
  ferx_model_from_dcm(dcm_fixture("model.ferx"), theta = theta, output_model = lit)
  expect_identical(names(declared_thetas(lit)$value), "TVV")
})

test_that("the default output path is <stem>_expanded.ferx next to the model", {
  dir <- tempfile("dcm")
  dir.create(dir)
  model <- file.path(dir, "net.ferx")
  file.copy(dcm_fixture("model.ferx"), model)
  res <- ferx_model_from_dcm(model, theta = dcm_fixture_theta())
  expect_identical(res$model, file.path(normalizePath(dir), "net_expanded.ferx"))
  expect_true(file.exists(res$model))
})

test_that("ferx_model_from_dcm() refuses bad inputs before the engine", {
  model <- dcm_fixture("model.ferx")
  theta <- dcm_fixture_theta()
  expect_error(ferx_model_from_dcm(model), "exactly one of `fit`")
  expect_error(ferx_model_from_dcm(model, fit = structure(list(), class = "ferx_fit"),
                                   theta = theta),
               "exactly one of `fit`")
  expect_error(ferx_model_from_dcm(model, fit = list(theta = 1)), "ferx_fit object")
  expect_error(ferx_model_from_dcm(model, theta = unname(theta)), "named numeric")
  expect_error(ferx_model_from_dcm(model, theta = theta, weights = "nope"),
               "should be one of")
  expect_error(ferx_model_from_dcm(model, theta = theta, network_only = NA),
               "TRUE or FALSE")
  expect_error(ferx_model_from_dcm(file.path(tempdir(), "absent.ferx"), theta = theta),
               "Model file not found")

  # The model file itself is never the output.
  dir <- tempfile("dcm")
  dir.create(dir)
  copy <- file.path(dir, "net.ferx")
  file.copy(model, copy)
  before <- readLines(copy)
  expect_error(ferx_model_from_dcm(copy, theta = theta, output_model = copy),
               "the model file itself")
  expect_identical(readLines(copy), before)
})

test_that("the engine refuses a model with no network, or a weight it is not given", {
  out <- tempfile(fileext = ".ferx")
  expect_error(
    ferx_model_from_dcm(ferx_example("warfarin")$model, theta = c(TVCL = 1),
                        output_model = out),
    "no \\[covariate_nn\\] block"
  )
  theta <- dcm_fixture_theta()
  expect_error(
    ferx_model_from_dcm(dcm_fixture("model.ferx"),
                        theta = theta[names(theta) != "W_TYPICAL_PK_1_1_1"],
                        output_model = out),
    "W_TYPICAL_PK_1_1_1"
  )
  expect_false(file.exists(out))
})

# ---------------------------------------------------------------------------
# `nn` builds: export from a fit, and predict identically
# ---------------------------------------------------------------------------

test_that("the export of a fitted DCM predicts identically to the DCM", {
  fit <- fit_warfarin_dcm()
  skip_without_nn(fit)
  ex <- ferx_example("warfarin_dcm")

  for (w in c("literal", "fixed_theta")) {
    # A ferx_model carries its dataset into the result.
    plain <- ferx_model_from_dcm(ferx_model(data = ex$data, model = ex$model), fit,
                                 output_model = tempfile(fileext = ".ferx"),
                                 weights = w)
    expect_identical(plain$data, ex$data)
    expect_false(any(grepl("^\\s*\\[covariate_nn", readLines(plain$model, warn = FALSE))))

    # The DCM at the fit's estimates; the export at its own declared values,
    # which are the same estimates. `expect_identical`, not `expect_equal`: the
    # target is bit equality, and a tolerance would hide a truncated weight.
    dcm_pred <- suppressWarnings(ferx_predict(ex$model, ex$data, fit = fit))
    plain_pred <- ferx_predict(plain$model, ex$data)
    expect_identical(plain_pred$PRED, dcm_pred$PRED, info = w)
    expect_identical(plain_pred[c("ID", "TIME")], dcm_pred[c("ID", "TIME")], info = w)

    # Simulation draws from the fitted omega / sigma, which the export writes
    # into its declarations: the same seed draws the same subjects.
    dcm_sim <- suppressWarnings(ferx_simulate(ex$model, ex$data, n_sim = 2L,
                                              seed = 7L, fit = fit))
    plain_sim <- ferx_simulate(plain$model, ex$data, n_sim = 2L, seed = 7L)
    expect_identical(plain_sim$IPRED, dcm_sim$IPRED, info = w)
    expect_identical(plain_sim$DV_SIM, dcm_sim$DV_SIM, info = w)
  }
})

test_that("a fit saved and loaded again exports the same model", {
  fit <- fit_warfarin_dcm()
  skip_without_nn(fit)
  ex <- ferx_example("warfarin_dcm")
  bundle <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, bundle)
  loaded <- suppressWarnings(ferx_load_fit(bundle))
  a <- ferx_model_from_dcm(ex$model, fit, output_model = tempfile(fileext = ".ferx"))
  b <- ferx_model_from_dcm(ex$model, loaded, output_model = tempfile(fileext = ".ferx"))
  expect_identical(readLines(b$model), readLines(a$model))
})

test_that("network_only = TRUE keeps the model's declared non-network values", {
  fit <- fit_warfarin_dcm()
  skip_without_nn(fit)
  ex <- ferx_example("warfarin_dcm")
  plain <- ferx_model_from_dcm(ex$model, fit, network_only = TRUE,
                               output_model = tempfile(fileext = ".ferx"))
  decl <- declared_thetas(plain$model)$value
  expect_identical(decl[["TVCL"]], "4.0")
  fitted <- ferx_model_from_dcm(ex$model, fit,
                                output_model = tempfile(fileext = ".ferx"))
  # The fit's estimate. To R's decimal reader only (one ULP); the prediction
  # test above holds the export to the bit.
  expect_equal(as.numeric(declared_thetas(fitted$model)$value[["TVCL"]]),
               unname(fit$theta[["TVCL"]]), tolerance = 4 * .Machine$double.eps)
  expect_false(isTRUE(all.equal(unname(fit$theta[["TVCL"]]), 4.0)))
})

test_that("a fit of a different model is refused", {
  fit <- fit_warfarin_dcm()
  skip_without_nn(fit)
  ex <- ferx_example("warfarin_dcm")

  # Same network, one equation changed: every weight is found, so only the
  # engine's model check stands between the fit and a wrong export.
  edited <- tempfile(fileext = ".ferx")
  lines <- readLines(ex$model, warn = FALSE)
  lines <- sub("KA = TVKA * TYPICAL_PK.KA * exp(ETA_KA)",
               "KA = TVKA * TYPICAL_PK.KA", lines, fixed = TRUE)
  writeLines(lines, edited)
  out <- tempfile(fileext = ".ferx")
  expect_error(ferx_model_from_dcm(edited, fit, output_model = out),
               "not the one the fit was made with")
  expect_false(file.exists(out))

  # A fit of a model without the network has none of its weights.
  wf <- ferx_example("warfarin")
  utils::capture.output(plain_fit <- suppressWarnings(
    ferx_fit(wf$model, wf$data, covariance = FALSE, settings = list(maxiter = 2L))
  ))
  expect_error(ferx_model_from_dcm(ex$model, plain_fit, output_model = out),
               "W_TYPICAL_PK")
  expect_false(file.exists(out))
})
