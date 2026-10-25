# A theta level block declared `profile` (ferx-core #1065, #2130): its levels
# are solved exactly inside the FOCE / FOCEI objective at every outer step and
# the optimizer searches only the structural parameters. The oracle is the
# same model unprofiled, which has the same objective, so the two fits must
# agree on the OFV, every theta and the standard errors, within the bands
# #2130 measured against its own unprofiled restart (theta 1e-3, SE 2e-2;
# its A1 anchor puts the OFV within 1e-5). Measured here on mbma_placebo
# (macOS arm64) at ferx-core d66c0c51: OFV 7.8e-9 below the unprofiled fit
# (6.6e-7 above at 5bd862b8, before #2285's verified stop; 1.8e-6 at
# 59a4fe3f), structural theta within 2.0e-6 relative, the levels within
# 7.1e-6 of their SE, structural SEs within 0.4%, and ferx_covariance()
# identical to the fit's SEs (#2173). An engine that ignored `profile` would fit the
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

test_that("a profiled fit writes its optimizer trace (ferx-core #2136)", {
  skip_on_cran()
  # The profiled outer loop skipped the trace, checkpoints, the flat-theta
  # freeze and the SIR fallback until ferx-core #2136.
  m <- mbma_profiled_model()
  p <- quiet(ferx_fit(m$model, m$data, verbose = FALSE, optimizer_trace = TRUE))
  expect_identical(p$profiled_blocks, "PLACEBO")
  tr <- ferx_trace(p)
  expect_true(is.data.frame(tr))
  expect_gt(nrow(tr), 0L)
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
  # Bit for bit since ferx-core #2173: the replay re-solves the profile as the
  # fit's own covariance step did (it used to differ at round-off, 5.5e-6).
  check <- function(cv) {
    expect_identical(se_parts(cv), se_parts(p))
    expect_true(all(cv$theta_levels$profiled))
    expect_identical(cv$theta_levels$se, p$theta_levels$se)
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

test_that("ferx_covariance() refuses s / rsr on a profiled fit with E_PROFILE_SCOPE (#597)", {
  skip_on_cran()
  # The profiled step computes `r` only (ferx-core #2108). The engine's replay
  # answers `s` / `rsr` with no matrix and a warning, which used to come back
  # as a fit with no standard errors and no R condition; ferx_fit() refuses
  # the same request up front, and so does ferx_covariance() now.
  p <- fits()$profiled
  for (m in c("s", "rsr", "cross_product", "SANDWICH")) {
    e <- tryCatch(quiet(ferx_covariance(p, covariance_method = m)),
                  error = function(e) e)
    expect_s3_class(e, "ferx_engine_error")
    expect_identical(e$code, "E_PROFILE_SCOPE")
    expect_identical(e$block, "parameters")
    expect_match(conditionMessage(e), "theta PLACEBO[...] is `profile`d",
                 fixed = TRUE)
    label <- c(s = "s", rsr = "rsr", cross_product = "s", SANDWICH = "rsr")[[m]]
    expect_match(conditionMessage(e),
                 sprintf("asks for `covariance_method = %s`.", label), fixed = TRUE)
  }
  # `r` and an unprofiled fit are untouched by the refusal.
  expect_identical(quiet(ferx_covariance(p, covariance_method = "hessian"))$covariance_status,
                   "computed")
  # The unprofiled fit runs its `s` step (its S is singular on mbma_placebo's
  # few studies, so the step reports failure rather than being refused).
  expect_s3_class(quiet(ferx_covariance(fits()$base, covariance_method = "s")),
                  "ferx_fit")
})

test_that("the profiled covariance refusal names every profiled block (#597)", {
  refuse <- ferx:::.ferx_refuse_profiled_covariance
  expect_null(refuse(list(profiled_blocks = character(0)), "s"))
  expect_null(refuse(list(), "rsr"))
  expect_null(refuse(list(profiled_blocks = "A"), "r"))
  e <- tryCatch(refuse(list(profiled_blocks = c("A", "B", "C")), "sandwich"),
                error = function(e) e)
  expect_identical(e$code, "E_PROFILE_SCOPE")
  expect_match(conditionMessage(e),
               "theta A[...], B[...] and C[...] are `profile`d", fixed = TRUE)
  expect_match(conditionMessage(e), "`covariance_method = rsr`", fixed = TRUE)
  expect_match(conditionMessage(e), "[E_PROFILE_SCOPE]", fixed = TRUE)
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

# A profiled fit whose covariance step fails on a reason any profiled model can
# reach: a theta added to TVE0 is confounded with it, so the profile
# objective's Hessian is singular. This used an unused theta until ferx-core
# #2136, which made a profiled fit freeze a theta the objective never reads
# (as an unprofiled fit does), so that fit's covariance step now succeeds.
confounded_theta_fit <- local({
  out <- list()
  function(covariance = TRUE) {
    key <- as.character(covariance)
    if (is.null(out[[key]])) {
      m <- mbma_profiled_model()
      txt <- readLines(m$model)
      txt <- append(txt, "  theta TVSHIFT(1.0, -100.0, 100.0)",
                    after = grep("^  theta PLACEBO\\[", txt))
      txt <- sub("BASE = TVE0 +", "BASE = TVE0 + TVSHIFT +", txt, fixed = TRUE)
      writeLines(txt, m$model)
      out[[key]] <<- quiet(ferx_fit(m$model, m$data, verbose = FALSE,
                                    covariance = covariance))
    }
    out[[key]]
  }
})

test_that("a profiled fit freezes a theta the objective never reads (ferx-core #2136)", {
  skip_on_cran()
  m <- mbma_profiled_model()
  txt <- readLines(m$model)
  txt <- append(txt, "  theta TVUNUSED(1.0, 0.01, 100.0)",
                after = grep("^  theta PLACEBO\\[", txt))
  writeLines(txt, m$model)
  fit <- quiet(ferx_fit(m$model, m$data, verbose = FALSE))
  expect_identical(fit$profiled_blocks, "PLACEBO")
  expect_identical(names(fit$theta)[fit$theta_fixed], "TVUNUSED")
  expect_identical(fit$covariance_status, "computed")
  # With the warning an unprofiled fit gives for the same theta.
  ws <- fit$warnings_structured
  flat <- ws[ws$category == "flat_parameter", , drop = FALSE]
  expect_identical(nrow(flat), 1L)
  expect_match(flat$message, "`TVUNUSED` has no effect on the objective", fixed = TRUE)
  expect_match(flat$message, "Freezing it at its initial value", fixed = TRUE)
})

test_that("a profiled fit's failed covariance step is reported in the covariance class", {
  skip_on_cran()
  # ferx-core #2158 (#2176): every "Covariance step failed" warning is
  # `covariance_failed`, whatever its reason. The reason that used to be filed
  # under `convergence`, "the profile solve did not converge", comes from
  # replay_covariance() only when the level solve breaks at the estimates,
  # which no bundled model does; the classifier table in
  # test-ferx_get_warnings.R holds its text.
  fit <- confounded_theta_fit()
  expect_identical(fit$profiled_blocks, "PLACEBO")
  expect_identical(fit$covariance_status, "failed")

  # `categories` is exact, one per matching row in table order. The replay's
  # row is classified by the engine like the fit's own (#611), so it is the
  # same row and ferx_covariance() leaves one.
  check <- function(ws, categories) {
    cov <- grepl("^Covariance step failed", ws$message)
    expect_identical(ws$category[cov], categories)
    for (i in which(cov)) {
      expect_match(ws$message[i], "Hessian of the profile objective", fixed = TRUE)
      g <- ferx:::.ferx_warning_guidance(ws$category[i], message = ws$message[i])
      expect_match(g, "confounded with the level effects", fixed = TRUE)
    }
  }
  check(fit$warnings_structured, "covariance_failed")
  check(ferx_get_warnings(quiet(ferx_covariance(fit)), as_df = TRUE),
        "covariance_failed")
})

test_that("a replayed covariance warning keeps the in-fit row's class and severity (#611)", {
  skip_on_cran()
  # ferx_covariance() used to fold the replay's flat warnings in as
  # `warning`/`covariance` rows: the fit's critical `covariance_failed` row
  # then sat beside a second, milder copy of itself, and a fit whose own step
  # never ran printed the replay's failure as a plain warning, missing from the
  # STATUS line.
  cols <- c("severity", "category", "message", "source_method")
  in_fit <- confounded_theta_fit()$warnings_structured
  in_fit <- in_fit[grepl("^Covariance step failed", in_fit$message), cols]
  expect_identical(nrow(in_fit), 1L)
  expect_identical(in_fit$severity, "critical")

  # Replayed on top of the failed in-fit step: one row, the fit's.
  ws <- quiet(ferx_covariance(confounded_theta_fit()))$warnings_structured
  rows <- ws[grepl("^Covariance step failed", ws$message), cols]
  rownames(rows) <- rownames(in_fit) <- NULL
  expect_identical(rows, in_fit)

  # Replayed on a fit whose step never ran: the replay's row alone, classed as
  # the fit's would have been, and it reaches the STATUS line.
  recov <- quiet(ferx_covariance(confounded_theta_fit(covariance = FALSE)))
  ws <- recov$warnings_structured
  rows <- ws[grepl("^Covariance step failed", ws$message), cols]
  rownames(rows) <- NULL
  expect_identical(rows, in_fit)
  status <- grep("STATUS:", capture.output(print(recov)), value = TRUE)
  expect_match(status, "(covariance_failed)", fixed = TRUE)
})
