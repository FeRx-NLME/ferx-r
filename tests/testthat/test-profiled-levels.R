# A theta level block declared `profile` (ferx-core #1065, #2130): its levels
# are solved exactly inside the FOCE / FOCEI objective at every outer step and
# the optimizer searches only the structural parameters. The oracle is the
# same model unprofiled, which has the same objective, so the two fits must
# agree on the OFV, every theta and the standard errors, within the bands
# #2130 measured against its own unprofiled restart (theta 1e-3, SE 2e-2;
# its A1 anchor puts the OFV within 1e-5). Measured here on mbma_placebo
# (macOS arm64): OFV 2.7e-9 below the unprofiled fit, structural theta within
# 8.6e-6 relative, the levels within 2.8e-5 of their SE, structural SEs within
# 0.7%, and ferx_covariance() within 4.2e-6 of the fit's SEs (#2130's replay
# gap is 5.5e-6, band 1e-4). An engine that ignored `profile` would fit the
# same numbers, so the fit must also say which blocks it profiled.

quiet <- function(x) suppressWarnings(suppressMessages(x))

structural <- c("TVE0", "EMAX", "ED50", "ET50")

# mbma_placebo with its PLACEBO block profiled. A profiled block is solved
# without constraints, so the bounds go (E_PROFILE_SCOPE otherwise).
mbma_profiled_model <- function(bounds = FALSE) {
  ex <- ferx_example("mbma_placebo")
  txt <- readLines(ex$model)
  at <- grep("^  theta PLACEBO\\[", txt)
  stopifnot(length(at) == 1L, grepl("(0.0, -50.0, 50.0)", txt[at], fixed = TRUE))
  if (!bounds) txt[at] <- sub("(0.0, -50.0, 50.0)", "(0.0)", txt[at], fixed = TRUE)
  txt[at] <- paste(txt[at], "profile")
  model <- tempfile(fileext = ".ferx")
  writeLines(txt, model)
  list(model = model, data = ex$data)
}

fits <- local({
  out <- NULL
  function() {
    if (is.null(out)) {
      ex <- ferx_example("mbma_placebo")
      m <- mbma_profiled_model()
      out <<- list(
        profiled = quiet(ferx_fit(m$model, m$data, verbose = FALSE)),
        base = quiet(ferx_fit(ex$model, ex$data, verbose = FALSE))
      )
    }
    out
  }
})

rel <- function(a, b) max(abs(a - b) / abs(b))

se_parts <- function(fit) {
  lapply(fit[c("se_theta", "se_omega", "se_kappa")], unname)
}

test_that("a profiled PLACEBO block fits to the unprofiled optimum", {
  skip_on_cran()
  f <- fits()
  p <- f$profiled
  b <- f$base
  expect_true(p$converged)
  expect_identical(p$covariance_status, "computed")
  expect_lt(abs(p$ofv - b$ofv), 1e-5)
  expect_identical(names(p$theta), names(b$theta))
  expect_lt(rel(p$theta[structural], b$theta[structural]), 1e-3)
  expect_lt(rel(p$se_theta[structural], b$se_theta[structural]), 2e-2)
  # Every level is still a theta, counted in n_parameters and BIC.
  expect_identical(p$n_parameters, b$n_parameters)
  expect_equal(p$bic, b$bic, tolerance = 1e-6)

  pl <- p$theta_levels
  bl <- b$theta_levels
  expect_identical(pl[c("block", "index", "label", "group", "contrast", "theta_name")],
                   bl[c("block", "index", "label", "group", "contrast", "theta_name")])
  expect_lt(max(abs(pl$value - bl$value) / bl$se), 1e-3)
  # Level SEs come from the profiled covariance, derived levels included.
  expect_true(all(is.finite(pl$se)))
  expect_lt(rel(pl$se, bl$se), 2e-2)
})

test_that("the fit says which blocks it profiled", {
  skip_on_cran()
  f <- fits()
  expect_identical(f$profiled$profiled_blocks, "PLACEBO")
  expect_true(all(f$profiled$theta_levels$profiled))
  expect_identical(f$base$profiled_blocks, character(0))
  expect_false(any(f$base$theta_levels$profiled))
  # `profiled` sits beside `contrast`: a property of the block.
  expect_identical(match("profiled", names(f$base$theta_levels)),
                   match("contrast", names(f$base$theta_levels)) + 1L)

  # mbma's 18-coefficient block is under the compact threshold, so it has no
  # THETA BLOCKS line; the Structural line marks it, in print and summary.
  structural_line <- function(out) grep("Structural:|Structure:", out, value = TRUE)
  out <- utils::capture.output(quiet(print(f$profiled)))
  expect_match(structural_line(out), "PLACEBO[STUDY, TIME] (18 coefficients, profiled)",
               fixed = TRUE)
  out <- utils::capture.output(quiet(print(summary(f$profiled))))
  expect_match(structural_line(out), "PLACEBO[STUDY, TIME] (18 coefficients, profiled)",
               fixed = TRUE)
  out <- utils::capture.output(quiet(print(f$base)))
  expect_match(structural_line(out), "PLACEBO[STUDY, TIME] (18 coefficients))",
               fixed = TRUE)
})

test_that("ferx_covariance() repeats the profiled covariance step, live and reloaded", {
  skip_on_cran()
  p <- fits()$profiled
  check <- function(cv) {
    a <- se_parts(cv)
    e <- se_parts(p)
    for (k in names(e)) expect_lt(rel(a[[k]], e[[k]]), 1e-4, label = k)
    expect_true(all(cv$theta_levels$profiled))
    expect_lt(rel(cv$theta_levels$se, p$theta_levels$se), 1e-4)
  }
  check(quiet(ferx_covariance(p)))
  # The glue rebuilds the rows' `profiled` from the fit's record, and the
  # refresh compares it: a record that disagrees with the fit's own column is
  # refused, never written onto the wrong rows.
  bad <- p
  bad$profiled_blocks <- character(0)
  expect_error(quiet(ferx_covariance(bad)), "laid the theta levels out differently",
               fixed = TRUE)

  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(p, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$profiled_blocks, "PLACEBO")
  expect_identical(back$theta_levels, p$theta_levels)
  expect_identical(back$theta, p$theta)
  check(quiet(ferx_covariance(back)))
})

test_that("profiled_blocks travels in fit.json the way ferx-core writes it", {
  skip_on_cran()
  f <- fits()
  read_fit_json <- function(fit) {
    path <- withr::local_tempfile(fileext = ".fitrx")
    quiet(ferx_save_fit(fit, path))
    fitrx_fit_json(path)
  }
  # An array, even for one block.
  expect_identical(read_fit_json(f$profiled)$profiled_blocks, list("PLACEBO"))
  # Omitted, not null, on a fit without one: ferx-core's loader reads a
  # missing key as none and refuses a null.
  expect_false("profiled_blocks" %in% names(read_fit_json(f$base)))

  # ferx-core's own loader reads the R-written key, and its writer writes it
  # back where ferx_load_fit() finds it.
  a <- withr::local_tempfile(fileext = ".fitrx")
  b <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(f$profiled, a, include_data = TRUE))
  quiet(getFromNamespace("ferx_rust_fitrx_engine_resave", "ferx")(a, b))
  back <- quiet(ferx_load_fit(b))
  expect_identical(back$profiled_blocks, "PLACEBO")
  expect_true(all(back$theta_levels$profiled))
})

test_that("a bounded profiled block is refused with E_PROFILE_SCOPE", {
  skip_on_cran()
  m <- mbma_profiled_model(bounds = TRUE)
  e <- tryCatch(quiet(ferx_fit(m$model, m$data, verbose = FALSE)),
                error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_PROFILE_SCOPE")
  expect_match(conditionMessage(e), "takes none", fixed = TRUE)
})

test_that("a profiled fit's failed covariance step is reported in the covariance class", {
  skip_on_cran()
  # ferx-core #2158 (#2176): every "Covariance step failed" warning is
  # `covariance_failed`, whatever its reason. The reason that used to be filed
  # under `convergence`, "the profile solve did not converge", comes from
  # replay_covariance() only when the level solve breaks at the estimates,
  # which no bundled model does; the classifier table in
  # test-ferx_get_warnings.R holds its text. This fit fails the step on a
  # reason any profiled model can reach: an unused theta gives the profile
  # objective's Hessian a zero row.
  m <- mbma_profiled_model()
  txt <- readLines(m$model)
  txt <- append(txt, "  theta TVUNUSED(1.0, 0.01, 100.0)",
                after = grep("^  theta PLACEBO\\[", txt))
  writeLines(txt, m$model)
  fit <- quiet(ferx_fit(m$model, m$data, verbose = FALSE))
  expect_identical(fit$profiled_blocks, "PLACEBO")
  expect_identical(fit$covariance_status, "failed")

  check <- function(ws) {
    cov <- grepl("^Covariance step failed", ws$message)
    expect_gte(sum(cov), 1L)
    expect_true(all(ws$category[cov] %in% c("covariance_failed", "covariance")))
    expect_match(ws$message[cov][1], "Hessian of the profile objective", fixed = TRUE)
    g <- ferx:::.ferx_warning_guidance(ws$category[cov][1], message = ws$message[cov][1])
    expect_match(g, "confounded with the level effects", fixed = TRUE)
  }
  check(fit$warnings_structured)
  check(ferx_get_warnings(quiet(ferx_covariance(fit)), as_df = TRUE))
})
