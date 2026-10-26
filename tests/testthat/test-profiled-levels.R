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

# mbma_placebo with its PLACEBO block profiled, without its bounds (the
# unbounded profile is the oracle the tests below compare with; a bounded one
# is tested further down).
mbma_profiled_model <- function() {
  ex <- ferx_example("mbma_placebo")
  txt <- readLines(ex$model)
  at <- grep("^  theta PLACEBO\\[", txt)
  stopifnot(length(at) == 1L, grepl("(0.0, -50.0, 50.0)", txt[at], fixed = TRUE))
  txt[at] <- sub("(0.0, -50.0, 50.0)", "(0.0)", txt[at], fixed = TRUE)
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

test_that("covariance_wall_time_secs times the profiled covariance step (ferx-core #2214)", {
  skip_on_cran()
  # The step's time includes the re-solve of the levels at the estimates
  # since ferx-core #2214. The re-solve cannot be timed apart from R, so this
  # checks the field is the step's, live, replayed and reloaded.
  p <- fits()$profiled
  expect_true(is.numeric(p$covariance_wall_time_secs))
  expect_gt(p$covariance_wall_time_secs, 0)
  expect_lt(p$covariance_wall_time_secs, p$wall_time_secs)
  # No step, no time; ferx_covariance() then reports its own step's.
  m <- mbma_profiled_model()
  none <- quiet(ferx_fit(m$model, m$data, verbose = FALSE, covariance = FALSE))
  expect_identical(none$covariance_wall_time_secs, 0)
  expect_gt(quiet(ferx_covariance(none))$covariance_wall_time_secs, 0)

  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(p, path))
  expect_identical(quiet(ferx_load_fit(path))$covariance_wall_time_secs,
                   p$covariance_wall_time_secs)
  # Under ferx-core's own key: an R writer and reader that agreed on another
  # name would pass the round trip above while the engine's loader read 0.
  trip <- quiet(fitrx_engine_trip(p))
  expect_identical(fitrx_fit_json(trip$core_path)$covariance_wall_time_secs,
                   p$covariance_wall_time_secs)
})

# mbma_placebo has six subjects, fewer than its free parameters, so its `s`
# step cannot be formed (ferx-core: "still needs more subjects than free
# parameters"), and each of its [STUDY, TIME] levels is read by one subject,
# which W_COV_RSR_ONE_SUBJECT_LEVELS flags under `rsr`. This fixture has
# neither: mbma_placebo's studies five times over with new IDs and noise on
# DV (30 subjects), and the levels on TIME alone, so every level is read by
# 30 subjects and the profiled block has four free coefficients.
visit_profiled_model <- function() {
  ex <- ferx_example("mbma_placebo")
  d <- utils::read.csv(ex$data)
  set.seed(2108)
  d <- do.call(rbind, lapply(0:4, function(r) {
    x <- d
    x$ID <- x$ID + 6L * r
    x$STUDY <- x$STUDY + 6L * r
    obs <- x$MDV == 0
    x$DV[obs] <- round(x$DV[obs] + stats::rnorm(sum(obs)), 2)
    x
  }))
  data <- tempfile(fileext = ".csv")
  utils::write.csv(d, data, row.names = FALSE, quote = FALSE)
  txt <- readLines(ex$model)
  at <- grep("^  theta PLACEBO\\[", txt)
  stopifnot(length(at) == 1L)
  txt[at] <- "  theta PLACEBO[TIME, contrast = sum_to_zero](0.0) profile"
  model <- tempfile(fileext = ".ferx")
  writeLines(txt, model)
  list(model = model, data = data)
}

test_that("s and rsr run on a profiled fit, and ferx_covariance() repeats them (ferx-core #2108)", {
  skip_on_cran()
  # Both were refused with E_PROFILE_SCOPE until ferx-core #2108, by
  # ferx_fit() and by ferx_covariance() (#597). Now the in-fit step computes
  # them and the replay repeats it bit for bit, live and after a reload.
  m <- visit_profiled_model()
  for (meth in c("s", "rsr")) {
    f <- quiet(ferx_fit(m$model, m$data, verbose = FALSE,
                        settings = list(covariance_method = meth)))
    expect_identical(f$profiled_blocks, "PLACEBO")
    expect_identical(f$covariance_status, "computed", info = meth)
    expect_true(all(is.finite(f$se_theta)) && all(f$se_theta > 0), info = meth)
    expect_true(all(is.finite(f$theta_levels$se)), info = meth)
    expect_false(any(grepl("W_COV_RSR_ONE_SUBJECT_LEVELS", f$warnings, fixed = TRUE)))

    check <- function(cv) {
      expect_identical(cv$covariance_status, "computed")
      expect_identical(se_parts(cv), se_parts(f))
      expect_identical(cv$cov_matrix, f$cov_matrix)
      expect_identical(cv$theta_levels$se, f$theta_levels$se)
    }
    check(quiet(ferx_covariance(f, covariance_method = meth)))
    path <- withr::local_tempfile(fileext = ".fitrx")
    quiet(ferx_save_fit(f, path))
    check(quiet(ferx_covariance(quiet(ferx_load_fit(path)), covariance_method = meth)))
  }
})

test_that("s and rsr on mbma_placebo's profiled block are diagnosed, not refused", {
  skip_on_cran()
  p <- fits()$profiled
  # Six subjects: S is singular, which the step reports as a failure with
  # its own advice (#308), not as E_PROFILE_SCOPE.
  s <- quiet(ferx_covariance(p, covariance_method = "cross_product"))
  expect_identical(s$covariance_status, "failed")
  w <- ferx_get_warnings(s, as_df = TRUE)
  row <- grepl("score cross-product matrix S", w$message, fixed = TRUE)
  expect_identical(w$category[row], "covariance_failed")
  # rsr runs, and says its level SEs are not usable: one subject per study
  # reads each [STUDY, TIME] level (ferx-core #2415).
  r <- quiet(ferx_covariance(p, covariance_method = "SANDWICH"))
  expect_identical(r$covariance_status, "computed")
  w <- ferx_get_warnings(r, as_df = TRUE)
  row <- grepl("W_COV_RSR_ONE_SUBJECT_LEVELS", w$message, fixed = TRUE)
  expect_identical(sum(row), 1L)
  expect_match(w$message[row], "24 of the 24 levels of theta PLACEBO[...]", fixed = TRUE)
  expect_identical(w$category[row], "covariance_step")
  g <- ferx:::.ferx_warning_guidance(w$category[row], message = w$message[row])
  expect_match(g, "were produced", fixed = TRUE)
  expect_match(g, "covariance_method = \"r\"", fixed = TRUE)
  expect_false(grepl("unavailable", g, fixed = TRUE))
  # `r`, the profiled fit's own step, carries no such warning.
  expect_false(any(grepl("W_COV_RSR_ONE_SUBJECT_LEVELS", p$warnings, fixed = TRUE)))
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

# mbma_placebo's PLACEBO block profiled with the box (-10, 4): the unbounded
# fit puts PLACEBO[STUDY=3,TIME=0] at 4.42 and every other free level inside
# the box (ferx-core #2112, #2402, #2416). Measured at ferx-core 32dbb309
# (macOS arm64): the profiled fit's OFV is 6.5e-9 below the same model
# unprofiled with the same bounds, whose optimizer stops the level at the
# bound too, and the structural thetas agree within 1e-5 relative. The
# standard errors are a different quantity - the unprofiled fit reports one
# for the level on its bound; the profiled step holds it - so they are not
# compared.
bounded_fits <- local({
  out <- NULL
  function() {
    if (is.null(out)) {
      ex <- ferx_example("mbma_placebo")
      mk <- function(profile) {
        txt <- readLines(ex$model)
        at <- grep("^  theta PLACEBO\\[", txt)
        txt[at] <- sub("(0.0, -50.0, 50.0)", "(0.0, -10.0, 4.0)", txt[at], fixed = TRUE)
        if (profile) txt[at] <- paste(txt[at], "profile")
        model <- tempfile(fileext = ".ferx")
        writeLines(txt, model)
        model
      }
      out <<- list(
        profiled = quiet(ferx_fit(mk(TRUE), ex$data, verbose = FALSE)),
        plain = quiet(ferx_fit(mk(FALSE), ex$data, verbose = FALSE))
      )
    }
    out
  }
})

held <- "PLACEBO[STUDY=3,TIME=0]"

test_that("a bounded profiled block fits, holding a level on its bound (ferx-core #2112)", {
  skip_on_cran()
  f <- bounded_fits()
  b <- f$profiled
  expect_true(b$converged)
  expect_identical(b$profiled_blocks, "PLACEBO")
  expect_identical(b$covariance_status, "computed")
  expect_lt(abs(b$ofv - f$plain$ofv), 1e-6)
  expect_lt(rel(b$theta[structural], f$plain$theta[structural]), 1e-4)
  expect_identical(unname(b$theta[held]), 4)

  # The held level, by name; a level inside the box is not held.
  expect_identical(b$profiled_levels_at_bound, held)
  free <- "PLACEBO[STUDY=3,TIME=2]"
  expect_gt(b$theta[[free]], -10)
  expect_lt(b$theta[[free]], 4)

  # Its SE is not available anywhere R reports one; fit$se_theta keeps the
  # engine's 0, as for a FIX theta. The level beside it keeps a finite SE.
  expect_identical(unname(b$se_theta[held]), 0)
  expect_gt(b$se_theta[[free]], 0)
  lv <- b$theta_levels
  expect_true(is.na(lv$se[which(lv$theta_name == held)]))
  expect_true(all(is.finite(lv$se[which(lv$theta_name != held | is.na(lv$theta_name))])))
  est <- b$estimates
  expect_true(all(is.na(est[est$param == held, c("se", "rse_pct", "lower_95", "upper_95")])))
  expect_false(est$fixed[est$param == held])
  expect_true(is.finite(est$se[est$param == free]))

  # print(): estimated, so no [FIX]; `N/A` as ferx-core's console prints it.
  out <- utils::capture.output(quiet(print(b)))
  row <- strsplit(trimws(grep(held, out, value = TRUE, fixed = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(row, c(held, "4.000000", "N/A", "N/A"))
  # ferx_runlog() too: it printed the engine's 0 as the SE and 0.0% as %RSE.
  rl <- utils::capture.output(quiet(ferx_runlog(b)))
  row <- strsplit(trimws(grep(held, rl, value = TRUE, fixed = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(tail(row, 2L), c("N/A", "N/A"))
  row <- strsplit(trimws(grep(free, rl, value = TRUE, fixed = TRUE)), "[[:space:]]+")[[1]]
  expect_match(tail(row, 1L), "%", fixed = TRUE)

  # One BoundaryEstimate warning names it.
  w <- ferx_get_warnings(b, as_df = TRUE)
  row <- grepl("Profiled level(s) on a declared bound", w$message, fixed = TRUE)
  expect_identical(sum(row), 1L)
  expect_identical(w$category[row], "boundary_estimate")
  expect_match(w$message[row], paste0(held, " (4.0000 at upper bound)"), fixed = TRUE)
  expect_match(w$message[row], "Their standard errors are not available", fixed = TRUE)

  # The unbounded profiled fit holds nothing.
  expect_identical(fits()$profiled$profiled_levels_at_bound, character(0))
})

test_that("a level held on its bound survives ferx_covariance() and .fitrx (ferx-core #2112)", {
  skip_on_cran()
  b <- bounded_fits()$profiled
  check <- function(cv) {
    expect_identical(cv$profiled_levels_at_bound, held)
    expect_identical(se_parts(cv), se_parts(b))
    expect_identical(cv$cov_matrix, b$cov_matrix)
    expect_identical(cv$theta_levels$se, b$theta_levels$se)
  }
  check(quiet(ferx_covariance(b)))
  # The step rebuilds the list from the estimates and the box, and the
  # replay takes the step's: a fit that lost it gets it back.
  lost <- b
  lost$profiled_levels_at_bound <- character(0)
  expect_identical(quiet(ferx_covariance(lost))$profiled_levels_at_bound, held)

  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(b, path))
  # ferx-core's key and shape: 0-based theta indices, an array.
  expect_identical(fitrx_fit_json(path)$profiled_levels_at_bound,
                   list(match(held, names(b$theta)) - 1L))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$profiled_levels_at_bound, held)
  expect_identical(back$theta_levels, b$theta_levels)
  expect_identical(back$estimates, b$estimates)
  check(quiet(ferx_covariance(back)))
  # The engine's own loader and writer round-trip the R-written key.
  trip <- quiet(fitrx_engine_trip(b))
  expect_identical(fitrx_fit_json(trip$core_path)$profiled_levels_at_bound,
                   list(match(held, names(b$theta)) - 1L))
  # Omitted, not empty, on a fit that holds none, as the engine writes it.
  quiet(ferx_save_fit(fits()$profiled, path))
  expect_false("profiled_levels_at_bound" %in% names(fitrx_fit_json(path)))

  # rsr holds it as well (ferx-core #2108 + #2112).
  r <- quiet(ferx_covariance(b, covariance_method = "rsr"))
  expect_identical(r$covariance_status, "computed")
  expect_identical(r$profiled_levels_at_bound, held)
  expect_true(is.na(r$theta_levels$se[which(r$theta_levels$theta_name == held)]))
})

test_that("a held level the fit does not have is refused, not dropped", {
  skip_on_cran()
  b <- bounded_fits()$profiled
  bad <- b
  bad$profiled_levels_at_bound <- "PLACEBO[STUDY=9,TIME=0]"
  expect_error(quiet(ferx_save_fit(bad, withr::local_tempfile(fileext = ".fitrx"))),
               "names a theta the fit does not have", fixed = TRUE)
  expect_error(quiet(ferx_sir(bad, sir_samples = 10L, sir_resamples = 5L)),
               "which is not one of the fit's thetas", fixed = TRUE)
  expect_error(ferx:::.fitrx_theta_names_at(list(3L), list("A", "B")),
               "outside the bundle's 2 thetas", fixed = TRUE)
  expect_identical(ferx:::.fitrx_theta_names_at(list(1L, 0L), list("A", "B")), c("B", "A"))
  expect_identical(ferx:::.fitrx_theta_names_at(NULL, list("A")), character(0))
})

# A profiled fit whose covariance step fails on a reason any profiled model can
# reach: a theta added to TVE0 is confounded with it, so the profile
# objective's Hessian is singular. This used an unused theta until ferx-core
# #2136, which made a profiled fit freeze a theta the objective never reads
# (as an unprofiled fit does), so that fit's covariance step now succeeds.
#
# At fd_hessian_step = 1e-2. Since ferx-core #2414 a profiled fit's
# structural stencil steps at 3e-3 by default, and at that step this
# confounded pair's Hessian comes out positive definite (smallest eigenvalue
# ~1e-5 against -0.0014 at 1e-2, measured at 32dbb309): the step reports
# standard errors, flagged only by a TVE0 ~ TVSHIFT correlation of -1.00 and
# a condition number of 1.6e5. The same model failed at either step before.
# ferx_covariance() replays the fit's step, so its replays fail too.
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
                                    covariance = covariance,
                                    fd_hessian_step = 0.01))
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
