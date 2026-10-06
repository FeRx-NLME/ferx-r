# Symbolic [covariate_model] statistics from R (#412, #487).
#
# A relation may state its centre as a statistic of the data
# (`center = median`). ferx_fit() resolves it on the fitted data and records
# the values in `fit$covariate_stats`; every later use of the fit centres the
# relations on those values, never on the design's. The oracle is the twin: the
# same model with the literal centre the fit resolved, which ferx-core's
# desugar makes the same model to the last bit
# (tests/covariate_model_equivalence.rs). The fixture is that test's data,
# `two_cpt_oral_cov.csv`, with one relation, `CL ~ WT power(center = C)`; C9
# adds a categorical `GRP` for the `levels` and `mode` statistics.
# The fits stop after two outer iterations: the assertions are about binding,
# not about where the optimizer ends up.

# `covariates` is the [covariates] line and `relation` the [covariate_model]
# line; by default the continuous `WT` power relation centred on `center`.
cs_model <- function(center, covariates = "WT continuous",
                     relation = sprintf(
                       "CL ~ WT power(center = %s) => THETA_CL_WT(0.6, 0.01, 5.0)",
                       center
                     )) {
  path <- tempfile(fileext = ".ferx")
  writeLines(sprintf("
[parameters]
  theta TVCL(4.0, 0.1, 100.0)
  theta TVV1(40.0, 1.0, 500.0)
  theta TVQ(8.0, 0.1, 100.0)
  theta TVV2(80.0, 1.0, 500.0)
  theta TVKA(1.0, 0.01, 10.0)
  omega ETA_CL ~ 0.15
  omega ETA_V1 ~ 0.15
  sigma PROP_ERR ~ 0.04 (sd)

[individual_parameters]
  CL = TVCL * exp(ETA_CL)
  V1 = TVV1 * exp(ETA_V1)
  Q  = TVQ
  V2 = TVV2
  KA = TVKA

[covariates]
  %s

[covariate_model]
  %s

[structural_model]
  pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2, ka=KA)

[error_model]
  DV ~ proportional(PROP_ERR)

[fit_options]
  method   = focei
  maxiter  = 2
  covariance = false
", covariates, relation), path)
  path
}

# Per-subject WT, one value per subject: how core summarises a static
# covariate (PsN's weighting).
cs_subject_wt <- function(rows) {
  vapply(split(rows$WT, rows$ID), function(x) unique(x)[1L], numeric(1L))
}

# The symbolic fit, its literal twin, and a design of the heavier half of the
# subjects, whose median weight is not the fit's. Shared by the tests that only
# read them.
cs_cache <- new.env(parent = emptyenv())
cs_base <- function() {
  if (is.null(cs_cache$fit)) {
    data <- ferx_example("two_cpt_oral_cov")$data
    rows <- utils::read.csv(data)
    wt <- cs_subject_wt(rows)
    heavy <- rows[rows$ID %in% as.numeric(names(wt)[wt > stats::median(wt)]), ]
    design <- tempfile(fileext = ".csv")
    utils::write.csv(heavy, design, row.names = FALSE, quote = FALSE, na = ".")
    sym <- cs_model("median")
    fit <- ferx_fit(sym, data, verbose = FALSE)
    twin <- cs_model(sprintf("%.17g", fit$covariate_stats$median))
    cs_cache$data <- data
    cs_cache$wt <- wt
    cs_cache$design <- design
    cs_cache$design_wt <- cs_subject_wt(heavy)
    cs_cache$sym <- sym
    cs_cache$twin <- twin
    cs_cache$fit <- fit
    cs_cache$twin_fit <- ferx_fit(twin, data, verbose = FALSE)
  }
  as.list(cs_cache)
}

# A copy of `fit` carrying a hand-made positive-definite covariance matrix in
# the packed space, so SIR and the covariance step reach the engine without
# relying on a two-iteration covariance step.
cs_with_cov <- function(fit) {
  n <- length(fit$theta) + nrow(fit$omega) + length(fit$sigma)
  fit$cov_matrix <- diag(1e-4, n)
  fit
}

cs_msg <- function(expr) {
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))
}

# --- C1: the fit binds -------------------------------------------------------

test_that("C1: ferx_fit() binds center = median and reports covariate_stats", {
  b <- cs_base()
  cs <- b$fit$covariate_stats
  expect_s3_class(cs, "data.frame")
  expect_identical(
    names(cs), c("covariate", "median", "mean", "min", "max", "mode", "levels")
  )
  expect_identical(cs$covariate, "WT")
  # One value per subject, as core summarises a static covariate.
  expect_equal(cs$median, stats::median(b$wt), tolerance = 1e-14)
  expect_equal(cs$mean, mean(b$wt), tolerance = 1e-12)
  expect_identical(cs$min, min(b$wt))
  expect_identical(cs$max, max(b$wt))
  expect_identical(cs$levels[[1L]], sort(unique(unname(b$wt))))
  # The median differs from the design's, or C3 would prove nothing.
  expect_gt(abs(stats::median(b$design_wt) - cs$median), 1)
})

test_that("C1b: the symbolic fit is its literal twin, to the last bit", {
  b <- cs_base()
  expect_identical(b$fit$ofv, b$twin_fit$ofv)
  expect_identical(b$fit$theta, b$twin_fit$theta)
  # A model with no symbolic relation records no statistics.
  expect_identical(nrow(b$twin_fit$covariate_stats), 0L)
})

test_that("C2: predict and simulate without a fit centre on the design", {
  b <- cs_base()
  at_design <- cs_model(sprintf("%.17g", stats::median(b$design_wt)))
  expect_identical(
    ferx_predict(b$sym, b$design)$PRED,
    ferx_predict(at_design, b$design)$PRED
  )
  expect_identical(
    ferx_simulate(b$sym, b$design, n_sim = 1L, seed = 3L)$IPRED,
    ferx_simulate(at_design, b$design, n_sim = 1L, seed = 3L)$IPRED
  )
})

# --- C3: a fit centres the design on the fit's statistics (T4 of #487) -------

test_that("C3: from-fit predict / simulate centre on the fit's median, not the design's", {
  b <- cs_base()
  sym <- ferx_predict(b$sym, b$design, fit = b$fit)$PRED
  expect_identical(sym, ferx_predict(b$twin, b$design, fit = b$fit)$PRED)
  # The shift this prevents: centring on the design's own median moves PRED.
  at_design <- cs_model(sprintf("%.17g", stats::median(b$design_wt)))
  expect_gt(
    max(abs(sym - ferx_predict(at_design, b$design, fit = b$fit)$PRED)), 1e-3
  )
  expect_identical(
    ferx_simulate(b$sym, b$design, fit = b$fit, n_sim = 1L, seed = 3L)$IPRED,
    ferx_simulate(b$twin, b$design, fit = b$fit, n_sim = 1L, seed = 3L)$IPRED
  )
})

test_that("C4: ferx_covariance / ferx_sir run on a symbolic fit like its twin", {
  b <- cs_base()
  c_sym <- ferx_covariance(cs_with_cov(b$fit))
  c_twin <- ferx_covariance(cs_with_cov(b$twin_fit))
  expect_identical(c_sym$ofv, c_twin$ofv)
  expect_identical(unname(c_sym$se_theta), unname(c_twin$se_theta))
  # SIR re-reads the fit's own data and refuses any other (`data hash
  # mismatch`), so no design can make a wrong centre differ from the right
  # one here: the twin's intervals are the strongest check SIR admits (#494).
  sir <- function(fit) {
    ferx_sir(cs_with_cov(fit), sir_samples = 50L, sir_resamples = 20L,
             sir_seed = 1L)
  }
  s_sym <- sir(b$fit)
  s_twin <- sir(b$twin_fit)
  expect_s3_class(s_sym, "ferx_fit")
  expect_identical(s_sym$sir_ci_theta, s_twin$sir_ci_theta)
  expect_identical(s_sym$sir_ci_omega, s_twin$sir_ci_omega)
  expect_identical(s_sym$sir_ci_sigma, s_twin$sir_ci_sigma)
})

# --- C5: a fit without statistics (T3 of #487) -------------------------------

test_that("C5: a symbolic model with a fit lacking statistics gets core's refusal on every path", {
  b <- cs_base()
  fit <- cs_with_cov(b$fit)
  # The shape of a fit saved before ferx recorded the statistics.
  fit$covariate_stats <- NULL
  msgs <- c(
    predict = cs_msg(ferx_predict(b$sym, b$design, fit = fit)),
    simulate = cs_msg(ferx_simulate(b$sym, b$design, fit = fit)),
    covariance = cs_msg(ferx_covariance(fit)),
    sir = cs_msg(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L))
  )
  for (k in names(msgs)) {
    expect_match(msgs[[k]], "this fit carries no data-derived bindings",
                 fixed = TRUE, info = k)
    expect_match(msgs[[k]], "a statistic of `WT` symbolically", fixed = TRUE,
                 info = k)
    # Nothing an R user cannot act on, and not the level-block refusal.
    expect_no_match(msgs[[k]], "ferx_core::api", fixed = TRUE, info = k)
    expect_no_match(msgs[[k]], "bind_covariate_stats", fixed = TRUE, info = k)
    expect_no_match(msgs[[k]], "theta level", fixed = TRUE, info = k)
  }
  # One engine writes it for every path.
  expect_length(unique(unname(msgs)), 1L)
})

test_that("C6: statistics for a covariate the model does not read are refused", {
  b <- cs_base()
  fit <- b$fit
  fit$covariate_stats$covariate <- "AGE"
  msg <- cs_msg(ferx_predict(b$sym, b$design, fit = fit))
  expect_match(msg, "carry no entry for it", fixed = TRUE)
})

test_that("C7: malformed statistics are refused in R's terms", {
  b <- cs_base()
  dup <- b$fit
  dup$covariate_stats <- rbind(dup$covariate_stats, dup$covariate_stats)
  expect_match(cs_msg(ferx_predict(b$sym, b$design, fit = dup)),
               "covariate `WT` is listed twice", fixed = TRUE)
  na <- b$fit
  na$covariate_stats$median <- NA_real_
  expect_match(cs_msg(ferx_predict(b$sym, b$design, fit = na)),
               "covariate `WT` has a `median` that is not a finite number",
               fixed = TRUE)
})

# --- C8: persistence ---------------------------------------------------------

test_that("C8: ferx_save_fit / ferx_load_fit keep covariate_stats", {
  b <- cs_base()
  bundle <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(b$fit, bundle)
  re <- ferx_load_fit(bundle)
  expect_identical(re$covariate_stats, b$fit$covariate_stats)
  expect_identical(
    ferx_predict(b$sym, b$design, fit = re)$PRED,
    ferx_predict(b$sym, b$design, fit = b$fit)$PRED
  )
  # In ferx-core's own slot, so the engine reads them too.
  fit_json <- file.path(withr::local_tempdir(), "x")
  utils::unzip(bundle, files = "fit.json", exdir = dirname(fit_json))
  wire <- jsonlite::read_json(file.path(dirname(fit_json), "fit.json"))
  expect_identical(wire$data_bindings$covariate_stats$WT$median,
                   b$fit$covariate_stats$median)
  expect_null(wire$r_extras$covariate_stats)
})

test_that("C8b: no statistic round-trips as zero rows; an absent slot loads NULL", {
  b <- cs_base()
  bundle <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(b$twin_fit, bundle)
  expect_identical(ferx_load_fit(bundle)$covariate_stats,
                   b$twin_fit$covariate_stats)
  # A bundle that does not record them (before #412): unknown, not empty.
  expect_null(ferx:::.fitrx_covariate_stats_from_wire(NULL))
})

test_that("C8c: a malformed statistics entry in a bundle is refused by name", {
  w <- list(covariate_stats = list(WT = list(
    median = 70, mean = 70, min = 45, max = 93, levels = list(45, 93)
  )))
  expect_error(ferx:::.fitrx_covariate_stats_from_wire(w),
               "covariate `WT` lacks one of median, mean, min, max, mode",
               fixed = TRUE)
  w$covariate_stats$WT$mode <- list(45)
  expect_error(ferx:::.fitrx_covariate_stats_from_wire(w),
               "covariate `WT` has a `mode` that is not a number", fixed = TRUE)
})

test_that("C8d: the engine's own .fitrx loader reads the statistics slot", {
  # `[priors] from_fit` is read by ferx-core's loader while the model file is
  # parsed, so a slot it cannot deserialise (a `median` written as an array, a
  # `levels` written as a scalar) fails validation with a parse error.
  # The covariance step's result: an import needs standard errors, and the
  # step keeps the fit's statistics.
  b <- cs_base()
  fit <- ferx_covariance(cs_with_cov(b$fit))
  expect_identical(fit$covariate_stats, b$fit$covariate_stats)
  bundle <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, bundle)
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(
    c(readLines(b$sym), "", "[priors]",
      paste0("  from_fit = ", gsub("/+", "/", normalizePath(bundle)))),
    model
  )
  utils::capture.output(res <- ferx_model_validate(model, b$data))
  expect_true(
    isTRUE(res$ok),
    info = paste(utils::capture.output(print(res$diagnostics)), collapse = "\n")
  )
})

# --- C9: levels = auto / ref = mode (review of #493, row 1) -----------------

# A categorical relation reads the statistics a continuous centre never does:
# `levels = auto` sets the relation's theta count from `levels`, and
# `ref = mode` its reference level from `mode`. `GRP` cycles 1, 2, 2, 3 over
# the subjects, so the mode (2) is neither the smallest nor the largest level.
cs_categorical <- function() {
  if (is.null(cs_cache$cat_fit)) {
    rows <- utils::read.csv(ferx_example("two_cpt_oral_cov")$data)
    ids <- sort(unique(rows$ID))
    grp <- stats::setNames(rep(c(1, 2, 2, 3), length.out = length(ids)), ids)
    rows$GRP <- grp[as.character(rows$ID)]
    data <- tempfile(fileext = ".csv")
    utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
    # Levels 1 and 2 only: the fit's level set and reference must be kept.
    design <- tempfile(fileext = ".csv")
    utils::write.csv(rows[rows$GRP %in% c(1, 2), ], design, row.names = FALSE,
                     quote = FALSE, na = ".")
    model <- function(levels, ref) {
      cs_model(covariates = sprintf("GRP categorical(levels = %s)", levels),
               relation = sprintf("CL ~ GRP categorical(ref = %s)", ref))
    }
    cs_cache$cat_data <- data
    cs_cache$cat_design <- design
    cs_cache$cat_sym <- model("auto", "mode")
    cs_cache$cat_twin <- model("[1, 2, 3]", "2")
    cs_cache$cat_fit <- ferx_fit(cs_cache$cat_sym, data, verbose = FALSE)
    cs_cache$cat_twin_fit <- ferx_fit(cs_cache$cat_twin, data, verbose = FALSE)
  }
  list(data = cs_cache$cat_data, design = cs_cache$cat_design,
       sym = cs_cache$cat_sym, twin = cs_cache$cat_twin,
       fit = cs_cache$cat_fit, twin_fit = cs_cache$cat_twin_fit)
}

test_that("C9: levels = auto / ref = mode fit and bind from the fit like their twin", {
  b <- cs_categorical()
  cs <- b$fit$covariate_stats
  expect_identical(cs$covariate, "GRP")
  expect_identical(cs$levels, list(c(1, 2, 3)))
  expect_identical(cs$mode, 2)
  # Reference 2: one theta per other level.
  expect_identical(
    names(b$fit$theta)[6:7], c("THETA_CL_GRP_1", "THETA_CL_GRP_3")
  )
  expect_identical(b$fit$ofv, b$twin_fit$ofv)
  expect_identical(b$fit$theta, b$twin_fit$theta)
  # The design holds levels 1 and 2 only. Its own `levels = auto` would drop
  # level 3 and change the theta count, so this reads the fit's `levels`.
  expect_identical(
    ferx_predict(b$sym, b$design, fit = b$fit)$PRED,
    ferx_predict(b$twin, b$design, fit = b$fit)$PRED
  )
  # And through a bundle, which carries `levels` as a JSON array.
  bundle <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(b$fit, bundle)
  re <- ferx_load_fit(bundle)
  expect_identical(re$covariate_stats, cs)
  expect_identical(
    ferx_predict(b$sym, b$design, fit = re)$PRED,
    ferx_predict(b$twin, b$design, fit = b$fit)$PRED
  )
})

# --- C10: uncertainty simulation and NPDE from a symbolic fit (#494) ---------

# The fit's data with every weight raised by 15: the same ID / TIME rows, so
# its NPDE still aligns with `fit$sdtab`, but a median 15 above the fit's. On the
# fit's own data an unbound symbolic model centres on the fit's median anyway,
# so NPDE there could not tell a dropped bind from a working one.
cs_shifted <- function() {
  if (is.null(cs_cache$shifted)) {
    rows <- utils::read.csv(cs_base()$data)
    rows$WT <- rows$WT + 15
    path <- tempfile(fileext = ".csv")
    utils::write.csv(rows, path, row.names = FALSE, quote = FALSE, na = ".")
    cs_cache$shifted <- path
    cs_cache$shifted_wt <- cs_subject_wt(rows)
  }
  list(data = cs_cache$shifted, wt = cs_cache$shifted_wt)
}

test_that("C10a: ferx_simulate_with_uncertainty() centres on the fit's median", {
  b <- cs_base()
  fit <- cs_with_cov(b$fit)
  swu <- function(model) {
    ferx_simulate_with_uncertainty(model, b$design, fit,
                                   n_uncertainty_draws = 3L, seed = 7L)
  }
  sym <- swu(b$sym)
  expect_gt(nrow(sym), 0L)
  expect_identical(sym, swu(b$twin))
  # Centring on the design's own median moves IPRED, so the identity above is
  # not true of a bind from the wrong source.
  at_design <- cs_model(sprintf("%.17g", stats::median(b$design_wt)))
  expect_gt(max(abs(sym$IPRED - swu(at_design)$IPRED)), 1e-3)
})

test_that("C10b: ferx_calc_npde() centres on the fit's median, not the data's", {
  b <- cs_base()
  s <- cs_shifted()
  npde <- function(model) {
    ferx_calc_npde(b$fit, nsim = 50L, seed = 5L, model = model,
                   data = s$data)$sdtab$NPDE
  }
  sym <- npde(b$sym)
  expect_true(all(is.finite(sym)))
  expect_identical(sym, npde(b$twin))
  expect_gt(abs(stats::median(s$wt) - b$fit$covariate_stats$median), 1)
  at_data <- cs_model(sprintf("%.17g", stats::median(s$wt)))
  expect_gt(max(abs(sym - npde(at_data))), 1e-3)
})

# --- C11: survival from a symbolic-centre TTE fit (#494) ---------------------

# The bundled exponential TTE data with a per-subject WT and the hazard scaled
# on it. `[individual_parameters]` needs the structural and error blocks, so
# the model keeps the bundled example's FIX dummy one-compartment triple.
cs_tte_model <- function(center) {
  path <- tempfile(fileext = ".ferx")
  writeLines(sprintf("
[parameters]
  theta TVLAMBDA(0.05, 0.001, 10.0)
  theta DUMMY_CL(1.0, FIX)
  theta DUMMY_V(1.0, FIX)
  omega ETA_LAMBDA ~ 0.09
  sigma SIGMA_DV ~ 0.01 FIX

[individual_parameters]
  LAMBDA = TVLAMBDA * exp(ETA_LAMBDA)
  CL     = DUMMY_CL
  V      = DUMMY_V

[covariates]
  WT continuous

[covariate_model]
  LAMBDA ~ WT power(center = %s) => THETA_LAMBDA_WT(0.8, 0.01, 5.0)

[structural_model]
  pk one_cpt_iv(cl=CL, v=V)

[error_model]
  DV ~ additive(SIGMA_DV)

[event_model]
  cmt    = 2
  family = exponential
  scale  = LAMBDA

[fit_options]
  method   = focei
  maxiter  = 2
  covariance = false
", center), path)
  path
}

# The symbolic TTE fit, its twin, and a design of the heavier half.
cs_tte <- function() {
  if (is.null(cs_cache$tte_fit)) {
    rows <- utils::read.csv(ferx_example("tte_exponential")$data)
    # Heavier subjects tend to have their events earlier (a fixed scramble of
    # the event order), so the fit estimates an exponent inside its bounds: at
    # the lower one the design's median and the fit's would give almost the
    # same survival.
    ids <- rows$ID[order(rows$TIME, rows$ID)]
    rank <- seq_along(ids)
    wt <- stats::setNames(100 - 2 * rank + 12 * ((rank * 7) %% 5 - 2), ids)
    rows$WT <- wt[as.character(rows$ID)]
    data <- tempfile(fileext = ".csv")
    utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
    heavy <- names(wt)[wt > stats::median(wt)]
    design <- tempfile(fileext = ".csv")
    utils::write.csv(rows[as.character(rows$ID) %in% heavy, ], design,
                     row.names = FALSE, quote = FALSE, na = ".")
    sym <- cs_tte_model("median")
    fit <- ferx_fit(sym, data, verbose = FALSE)
    twin <- cs_tte_model(sprintf("%.17g", fit$covariate_stats$median))
    cs_cache$tte_design <- design
    cs_cache$tte_design_wt <- unname(wt[heavy])
    cs_cache$tte_sym <- sym
    cs_cache$tte_twin <- twin
    cs_cache$tte_fit <- fit
    cs_cache$tte_twin_fit <- ferx_fit(twin, data, verbose = FALSE)
  }
  list(design = cs_cache$tte_design, design_wt = cs_cache$tte_design_wt,
       sym = cs_cache$tte_sym, twin = cs_cache$tte_twin,
       fit = cs_cache$tte_fit, twin_fit = cs_cache$tte_twin_fit)
}

test_that("C11: ferx_predict_survival() centres on a TTE fit's median", {
  b <- cs_tte()
  expect_identical(b$fit$covariate_stats$covariate, "WT")
  expect_identical(b$fit$ofv, b$twin_fit$ofv)
  expect_identical(b$fit$theta, b$twin_fit$theta)
  surv <- function(model) {
    ferx_predict_survival(model, b$design, times = c(5, 10, 20),
                          fit = b$fit)$survival
  }
  sym <- surv(b$sym)
  expect_gt(length(sym), 0L)
  expect_identical(sym, surv(b$twin))
  at_design <- cs_tte_model(sprintf("%.17g", stats::median(b$design_wt)))
  expect_gt(max(abs(sym - surv(at_design))), 1e-3)
})
