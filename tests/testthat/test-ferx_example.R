test_that("known name returns list with $model and $data", {
  ex <- ferx_example("warfarin")
  expect_true(is.list(ex))
  expect_true(all(c("model", "data") %in% names(ex)))
})
test_that("$model path exists on disk", {
  ex <- ferx_example("warfarin")
  expect_true(file.exists(ex$model))
})
test_that("$data path exists on disk", {
  ex <- ferx_example("warfarin")
  expect_true(file.exists(ex$data))
})
test_that("$model has .ferx extension", {
  ex <- ferx_example("warfarin")
  expect_true(grepl("\\.ferx$", ex$model))
})
test_that("$data has .csv extension", {
  ex <- ferx_example("warfarin")
  expect_true(grepl("\\.csv$", ex$data))
})
test_that("no-arg call returns character vector of available names", {
  res <- ferx_example()
  expect_type(res, "character")
  expect_true(length(res) > 0L)
  expect_true("warfarin" %in% res)
})
test_that("unknown name errors with available names in message", {
  expect_error(ferx_example("not_a_real_example"), regexp = "warfarin")
})
test_that("$data path exists on disk for every bundled example", {
  # Every example -- whether it ships its own dataset or resolves to one via the
  # .data_aliases map in ferx_example() -- must point at a CSV that exists.
  # Driven off ferx_example() rather than a hand-maintained list so newly added
  # examples and aliases are covered automatically (a previous hardcoded list
  # silently missed every alias added after it was written).
  for (nm in ferx_example()) {
    ex <- ferx_example(nm)
    expect_true(
      is.character(ex$data) && nzchar(ex$data) && file.exists(ex$data),
      label = paste0("$data for '", nm, "' must exist on disk")
    )
  }
})
test_that("emax_timecourse bundles a compartment-free model and synthetic data", {
  ex <- ferx_example("emax_timecourse")

  structural <- ferx:::.ferx_extract_blocks(ex$model)$structural_model
  expect_true(any(grepl("^EFF\\s*=", structural)))
  expect_true(any(grepl("^y\\s*=", structural)))
  expect_false(any(grepl("^(pk|ode|ode_template)\\b", structural)))
  capture.output(s <- ferx_model_inspect(ex$model))
  expect_equal(s$model_type, "compartment-free")

  dat <- utils::read.csv(ex$data)
  expect_identical(names(dat), c("ID", "TIME", "DV", "MDV"))
  expect_equal(nrow(dat), 240L)
  expect_equal(length(unique(dat$ID)), 30L)
  expect_false("AMT" %in% names(dat))
})
test_that("warfarin_derived $data points to warfarin.csv (alias check)", {
  ex <- ferx_example("warfarin_derived")
  expect_true(grepl("warfarin\\.csv$", ex$data))
})

# -- mbma_placebo (#370) -------------------------------------------------------
# A model-based meta-analysis of synthetic arm-level data: 6 studies x 4 visits
# = 24 (STUDY, TIME) cells in a PLACEBO[STUDY, TIME] level block, 18 of them
# free under sum_to_zero_within. With TVE0, EMAX, ED50, ET50, omega and gamma^2
# that is 24 estimated parameters, the most the analytic outer gradient takes.

mbma_cols <- c("ID", "STUDY", "OCC", "TIME", "DV", "MDV", "NARM", "SE", "DOSE")

# The generator's truth (T2 ties every value back to it). TVE0 is the realised
# one: under sum_to_zero_within TVE0 + ETA_E0 is a study's mean placebo level
# over its own visits, 50 + eta_s + mean_t P(s, t), so the truth is the mean of
# that over the 6 studies. GAMMA2 is the between-arm variance of a 1-patient
# arm.
mbma_truth <- c(TVE0 = 48.1949, EMAX = 12, ED50 = 40, ET50 = 3, GAMMA2 = 100)

# Fit the example once for T3-T5. Engine warnings land in fit$warnings; R-level
# warnings raised on the way are kept beside the fit.
mbma_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      ex <- ferx_example("mbma_placebo")
      r_warnings <- character(0)
      fit <- withCallingHandlers(
        ferx_fit(ex$model, ex$data, verbose = FALSE),
        warning = function(w) {
          r_warnings <<- c(r_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      cached <<- list(fit = fit, r_warnings = r_warnings)
    }
    cached
  }
})

test_that("T1: mbma_placebo bundles a level-block MBMA model and its arm-level data", {
  ex <- ferx_example("mbma_placebo")
  expect_true(file.exists(system.file("examples", "ex_mbma_placebo.R", package = "ferx")))
  expect_true(file.exists(
    system.file("examples", "mbma_placebo", "simulate_dataset.R", package = "ferx")
  ))

  blocks <- ferx:::.ferx_extract_blocks(ex$model)
  expect_true(any(grepl("^theta PLACEBO\\[STUDY, TIME", blocks$parameters)))
  expect_true(any(grepl("^kappa KAPPA_ARM .*weight = NARM$", blocks$parameters)))
  expect_true(any(grepl("^sigma ADD_ERR .*FIX$", blocks$parameters)))
  expect_true(any(grepl("^DV ~ additive\\(ADD_ERR\\) weight = SE$", blocks$error_model)))
  capture.output(s <- ferx_model_inspect(ex$model))
  expect_equal(s$model_type, "compartment-free")

  d <- utils::read.csv(ex$data)
  expect_identical(names(d), mbma_cols)
  expect_equal(nrow(d), 88L)
  expect_equal(length(unique(d$STUDY)), 6L)
  expect_equal(nrow(unique(d[c("STUDY", "TIME")])), 24L)
  expect_identical(d$ID, d$STUDY)
  expect_true(all(d$NARM > 0))
  expect_true(all(d$SE > 0))
})

test_that("T2: the bundled CSV is what mbma_placebo/simulate_dataset.R writes", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  gen <- system.file("examples", "mbma_placebo", "simulate_dataset.R", package = "ferx")
  env <- new.env()
  env$mbma_out_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(env$mbma_out_csv), add = TRUE)
  capture.output(suppressMessages(sys.source(gen, envir = env)))

  # The draws were measured bit-identical (17 significant digits) on aarch64
  # macOS (R 4.6.1) and on rocker/r-ver 4.4.2 under linux/amd64 and
  # linux/arm64, so this holds on CI too.
  expect_equal(utils::read.csv(env$mbma_out_csv), utils::read.csv(ex$data),
               tolerance = 1e-12)
  # mbma_truth is what the generator used ...
  expect_equal(env$tve0_realised, mbma_truth[["TVE0"]], tolerance = 1e-5)
  expect_identical(
    c(EMAX = env$EMAX, ED50 = env$ED50, ET50 = env$ET50, GAMMA2 = env$GAMMA2),
    mbma_truth[c("EMAX", "ED50", "ET50", "GAMMA2")]
  )
  # ... and so is the truth table the example script prints.
  script <- readLines(system.file("examples", "ex_mbma_placebo.R", package = "ferx"))
  truth_line <- grep("^truth <- c\\(", script, value = TRUE)
  expect_length(truth_line, 1L)
  expect_identical(eval(parse(text = sub("^truth <- ", "", truth_line))),
                   mbma_truth[c("TVE0", "EMAX", "ED50", "ET50")])
})

test_that("T3: the mbma_placebo fit recovers its simulation truth", {
  skip_on_cran()
  fit <- mbma_fit()$fit

  expect_true(fit$converged)
  expect_identical(fit$covariance_status, "computed")
  expect_false(any(grepl("regulari[sz]", unlist(fit$warnings), ignore.case = TRUE)))

  tl <- fit$theta_levels
  expect_equal(nrow(tl), 24L)
  expect_identical(tl$contrast, rep("sum_to_zero_within", 24L))
  # One derived level per study; every other level is a named fitted theta.
  expect_identical(sort(tl$group[is.na(tl$theta_name)]), 0:5)
  expect_equal(sum(!is.na(tl$theta_name)), 18L)
  expect_true(all(tl$theta_name[!is.na(tl$theta_name)] %in% names(fit$theta)))

  for (p in c("TVE0", "EMAX", "ED50", "ET50")) {
    z <- (fit$theta[[p]] - mbma_truth[[p]]) / fit$se_theta[[p]]
    expect_lt(abs(z), 1.96, label = paste0("|z| for ", p))
  }
  # weight = NARM: the kappa estimate is the unweighted gamma^2.
  z_gamma2 <- (fit$omega_iov[1, 1] - mbma_truth[["GAMMA2"]]) / fit$se_kappa[[1]]
  expect_lt(abs(z_gamma2), 1.96, label = "|z| for gamma^2")
  expect_lt(fit$shrinkage_kappa[[1]], 0.30)
})

test_that("T4: the mbma_placebo fit raises no false 'not referenced' warning", {
  skip_on_cran()
  # Every free PLACEBO level is read through the block, so none is unused.
  # Before ferx-core#1640 (#1628) each of the 18 drew this warning.
  res <- mbma_fit()
  expect_false(any(grepl("not referenced", unlist(res$fit$warnings))))
  expect_false(any(grepl("not referenced", res$r_warnings)))
})

test_that("T5: PLACEBO on its own line fits like the shipped form (to 1e-8)", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  fit <- mbma_fit()$fit

  txt <- readLines(ex$model)
  shipped <- "  BASE = TVE0 + PLACEBO + ETA_E0 + KAPPA_ARM"
  expect_equal(sum(txt == shipped), 1L)
  txt[txt == shipped] <- "  PL   = PLACEBO\n  BASE = TVE0 + PL + ETA_E0 + KAPPA_ARM"
  split <- tempfile(fileext = ".ferx")
  on.exit(unlink(split), add = TRUE)
  writeLines(txt, split)

  # Before ferx-core#1640 (#1628) the split form left every level at its init:
  # OFV 284.76 against 97.07. The two texts may compile to differently ordered
  # expressions, so they are compared to 1e-8 rather than bit for bit.
  fit_split <- suppressWarnings(ferx_fit(split, ex$data, verbose = FALSE))
  expect_equal(fit_split$ofv, fit$ofv, tolerance = 1e-8)
  expect_equal(fit_split$theta, fit$theta, tolerance = 1e-8)
})
