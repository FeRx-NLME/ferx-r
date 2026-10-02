# theta NAME[COL, ...] level blocks from R (#370).
#
# A level block declares one theta per observed combination of its columns, so
# its theta layout is a property of the data it is bound to. ferx_fit() binds
# it to the fitted data; every later use of the fit binds the design against
# the fit's layout (`fit$theta_levels`), so a theta is never read at a position
# the fit did not give it. The fixtures are ferx-core's own
# (tests/theta_level_blocks.rs): two studies x TIME {1, 4, 12}, one subject per
# study, plus a PLA_IDX column holding the same design in the counted form.
# The fits stop after two outer iterations: the assertions are about binding,
# not about where the optimizer ends up.

tl_data <- "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY,PLA_IDX
1,0,.,1,100,1,0,1,1,1
1,1,8.1,0,.,1,0,0,1,1
1,4,6.2,0,.,1,0,0,1,2
1,12,3.1,0,.,1,0,0,1,3
2,0,.,1,100,1,0,1,2,4
2,1,7.4,0,.,1,0,0,2,4
2,4,5.5,0,.,1,0,0,2,5
2,12,2.8,0,.,1,0,0,2,6"

# Three studies x three times, one subject per study, for the nested
# (sum_to_zero_within) contrast.
tl_data3 <- "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY
1,0,.,1,100,1,0,1,1
1,1,8.1,0,.,1,0,0,1
1,4,6.2,0,.,1,0,0,1
1,12,3.1,0,.,1,0,0,1
2,0,.,1,100,1,0,1,2
2,1,7.4,0,.,1,0,0,2
2,4,5.5,0,.,1,0,0,2
2,12,2.8,0,.,1,0,0,2
3,0,.,1,100,1,0,1,3
3,1,7.9,0,.,1,0,0,3
3,4,5.9,0,.,1,0,0,3
3,12,2.5,0,.,1,0,0,3"

tl_fit_options <- "
[fit_options]
  maxiter = 2
  inner_maxiter = 3
  covariance = false
"

# A one-compartment IV model whose clearance reads `cl`, with the theta lines
# `thetas` and the structural block `structure` (analytical by default).
tl_model <- function(thetas, cl,
                     structure = "[structural_model]\n  pk one_cpt_iv(cl=CL, v=V)\n",
                     eta = "  omega ETA_V ~ 0.04", v = "TVV * exp(ETA_V)") {
  paste0(
    "[parameters]\n", thetas, "\n",
    "  theta TVV(10.0, 0.1, 500.0)\n",
    eta, "\n",
    "  sigma PROP_ERR ~ 0.05\n\n",
    "[individual_parameters]\n",
    "  CL = ", cl, "\n",
    "  V  = ", v, "\n\n",
    structure, "\n",
    "[error_model]\n  DV ~ proportional(PROP_ERR)\n",
    tl_fit_options
  )
}

tl_write <- function(text, ext) {
  path <- tempfile(fileext = ext)
  writeLines(text, path)
  path
}

# The level-block model with an intercept: global sum-to-zero (auto).
tl_col_model <- function() {
  tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO"
  ), ".ferx")
}

tl_fit <- function(model, data, ...) {
  ferx_fit(model, data, verbose = FALSE, ...)
}

# The T1 fit, shared by the tests that only read it.
tl_cache <- new.env(parent = emptyenv())
tl_base <- function() {
  if (is.null(tl_cache$fit)) {
    tl_cache$model <- tl_col_model()
    tl_cache$data <- tl_write(tl_data, ".csv")
    tl_cache$fit <- tl_fit(tl_cache$model, tl_cache$data)
  }
  list(model = tl_cache$model, data = tl_cache$data, fit = tl_cache$fit)
}

# A copy of `fit` carrying a hand-made positive-definite covariance matrix in
# the packed space (theta, omega, sigma), so the uncertainty and SIR legs reach
# the engine without relying on a two-iteration covariance step.
tl_with_cov <- function(fit) {
  n <- length(fit$theta) + nrow(fit$omega) + length(fit$sigma)
  fit$cov_matrix <- diag(1e-4, n)
  fit
}

# --- T1: the fit binds ---------------------------------------------------------

test_that("T1: ferx_fit() binds a level block and reports theta_levels", {
  b <- tl_base()
  fit <- b$fit
  expect_identical(names(fit$theta), c(
    "TVCL",
    "PLACEBO[STUDY=1,TIME=1]", "PLACEBO[STUDY=1,TIME=4]",
    "PLACEBO[STUDY=1,TIME=12]", "PLACEBO[STUDY=2,TIME=1]",
    "PLACEBO[STUDY=2,TIME=4]",
    "TVV"
  ))
  tl <- fit$theta_levels
  expect_s3_class(tl, "data.frame")
  expect_identical(
    names(tl), c("block", "index", "label", "group", "contrast", "theta_name")
  )
  expect_identical(nrow(tl), 6L)
  expect_identical(tl$block, rep("PLACEBO", 6L))
  expect_identical(tl$index, 1:6)
  expect_identical(tl$label, c(
    "STUDY=1,TIME=1", "STUDY=1,TIME=4", "STUDY=1,TIME=12",
    "STUDY=2,TIME=1", "STUDY=2,TIME=4", "STUDY=2,TIME=12"
  ))
  expect_identical(tl$contrast, rep("sum_to_zero", 6L))
  expect_identical(tl$group, rep(0L, 6L))
  # Only the level the contrast derives from the others has no theta.
  expect_identical(which(is.na(tl$theta_name)), 6L)
  expect_identical(tl$theta_name[1:5], names(fit$theta)[2:6])
})

# --- T2: the named form is the counted form -----------------------------------

test_that("T2: contrast = none fits exactly like the counted form", {
  data <- tl_write(tl_data, ".csv")
  # No intercept: unconstrained levels plus an intercept would be
  # rank-deficient.
  col <- tl_write(tl_model(
    "  theta PLACEBO[STUDY, TIME, contrast = none](2.0, 0.001, 20.0)", "PLACEBO"
  ), ".ferx")
  counted <- tl_write(tl_model(
    "  theta PLACEBO[6](2.0, 0.001, 20.0)", "PLACEBO[PLA_IDX]"
  ), ".ferx")
  f_col <- tl_fit(col, data)
  f_cnt <- tl_fit(counted, data)
  expect_identical(f_col$theta_levels$contrast, rep("none", 6L))
  expect_false(anyNA(f_col$theta_levels$theta_name))
  # Measured bit-identical (OFV 24.578051826322419 for both, theta max abs
  # difference 0): the same records read the same theta through the same
  # gather. The tolerance only absorbs a platform's floating-point noise.
  expect_equal(f_col$ofv, f_cnt$ofv, tolerance = 1e-10)
  expect_equal(unname(f_col$theta), unname(f_cnt$theta), tolerance = 1e-10)
})

# --- T3: stamps on the model survive the bind's re-parse -----------------------

test_that("T3: a settings ODE override survives the level bind", {
  data <- tl_write(tl_data, ".csv")
  ode <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO",
    structure = paste0(
      "[structural_model]\n  ode(obs_cmt=central, states=[central])\n\n",
      "[odes]\n  d/dt(central) = -(CL/V) * central\n\n",
      "[scaling]\n  obs_scale = V\n"
    )
  ), ".ferx")
  has_ode_diag <- function(fit) {
    any(startsWith(fit$warnings, "W_ODE_SOLVER_DIAGNOSTICS"))
  }
  # Control: the default solver budget is ample for this model.
  expect_false(has_ode_diag(tl_fit(ode, data, gradient = "fd")))
  starved <- tl_fit(ode, data, gradient = "fd",
                    settings = list(ode_max_steps = 2))
  expect_true(has_ode_diag(starved))
  expect_identical(starved$gradient_used, "fd")
  # The engine carries a non-default ODE option to the integrator itself, so
  # the stamp on the model decides the outcome only when the call restores the
  # *default* over a value the model file pinned: lose the stamp to the bind's
  # re-parse and the file's starved budget comes back.
  pinned <- tl_write(sub(
    "[fit_options]\n", "[fit_options]\n  ode_max_steps = 2\n",
    paste(readLines(ode), collapse = "\n"), fixed = TRUE
  ), ".ferx")
  expect_true(has_ode_diag(tl_fit(pinned, data, gradient = "fd")))
  expect_warning(
    restored <- tl_fit(pinned, data, gradient = "fd",
                       settings = list(ode_max_steps = 10000)),
    "overrides it with `10000`", fixed = TRUE
  )
  expect_false(has_ode_diag(restored))
})

test_that("T3b: a bloq_method argument survives the level bind", {
  # STUDY 2's TIME 12 record is below the quantification limit: M3 scores it
  # as a censored likelihood, `drop` leaves it out, so the OFVs differ only if
  # the method stamped on the model reaches the fit.
  cens <- sub("ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY,PLA_IDX",
              "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY,PLA_IDX,CENS", tl_data)
  rows <- strsplit(cens, "\n")[[1]]
  rows[-1] <- paste0(rows[-1], ",0")
  rows[9] <- sub(",0$", ",1", rows[9])
  data <- tl_write(paste(rows, collapse = "\n"), ".csv")
  model <- tl_col_model()
  m3 <- tl_fit(model, data, bloq_method = "m3")
  dropped <- tl_fit(model, data, bloq_method = "drop")
  expect_false(isTRUE(all.equal(m3$ofv, dropped$ofv)))
})

test_that("T3c: the block binds to the data after ignore =", {
  data <- tl_write(tl_data, ".csv")
  col <- tl_write(tl_model(
    "  theta PLACEBO[STUDY, TIME, contrast = none](2.0, 0.001, 20.0)", "PLACEBO"
  ), ".ferx")
  fit <- tl_fit(col, data, ignore = "STUDY == 2")
  expect_identical(fit$theta_levels$label, c(
    "STUDY=1,TIME=1", "STUDY=1,TIME=4", "STUDY=1,TIME=12"
  ))
  expect_identical(names(fit$theta), c(
    "PLACEBO[STUDY=1,TIME=1]", "PLACEBO[STUDY=1,TIME=4]",
    "PLACEBO[STUDY=1,TIME=12]", "TVV"
  ))
})

# --- T4: a design is read at the fit's positions -------------------------------

test_that("T4: a subset design reads the fit's theta, checked in closed form", {
  b <- tl_base()
  fit <- b$fit
  fit$theta[] <- c(2, 0.1, 0.2, 0.3, 0.4, 0.5, 10)
  # Study 2 only: its levels are the fit's 4, 5 and 6.
  rows <- strsplit(tl_data, "\n")[[1]]
  design <- tl_write(paste(rows[c(1, 6:9)], collapse = "\n"), ".csv")
  pred <- ferx_predict(b$model, design, fit = fit)
  # Global sum-to-zero: the last level is minus the sum of the five free ones.
  p <- c(0.4, 0.5, -(0.1 + 0.2 + 0.3 + 0.4 + 0.5))
  v <- 10
  dt <- c(1, 3, 8)
  # IV bolus of 100 at TIME 0; CL = TVCL + P_k on the interval ending at the
  # k-th observation.
  expected <- 100 / v * exp(-cumsum((2 + p) / v * dt))
  expect_equal(pred$TIME, c(1, 4, 12))
  expect_equal(pred$PRED, expected, tolerance = 1e-6)
})

# --- T5: a design level the fit never saw is refused ---------------------------

test_that("T5: every from-fit entry point refuses labels the fit never saw", {
  b <- tl_base()
  # Same level count as the fit, different labels: TIME 12 becomes 24. The
  # design keeps real DVs, because npde reads it with the fitting policy.
  design <- tl_write(gsub(",12,", ",24,", tl_data, fixed = TRUE), ".csv")
  check <- function(expr) {
    err <- tryCatch({
      expr
      NULL
    }, error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    expect_match(msg, "theta PLACEBO[STUDY, TIME]", fixed = TRUE)
    expect_match(msg, "`STUDY=1,TIME=24`", fixed = TRUE)
    expect_match(msg, "`STUDY=2,TIME=24`", fixed = TRUE)
    expect_match(msg, "can only be simulated at the fit's observation times",
                 fixed = TRUE)
    # No stage prefix, which would let the single-error fallback attach an
    # unrelated diagnostic code.
    expect_false(startsWith(msg, "Error reading data:"))
    expect_false(startsWith(msg, "Error parsing model:"))
  }
  check(ferx_simulate(b$model, design, fit = b$fit))
  check(ferx_predict(b$model, design, fit = b$fit))
  check(ferx_simulate_with_uncertainty(b$model, design, tl_with_cov(b$fit),
                                       n_uncertainty_draws = 2L))
  check(ferx_calc_npde(b$fit, nsim = 10L, model = b$model, data = design))
  check(ferx_predict_survival(b$model, design, times = c(1, 2), fit = b$fit))
})

# --- T6: design order does not matter -----------------------------------------

test_that("T6: a reordered design predicts the same rows", {
  b <- tl_base()
  rows <- strsplit(tl_data, "\n")[[1]]
  reordered <- tl_write(paste(rows[c(1, 6:9, 2:5)], collapse = "\n"), ".csv")
  original <- ferx_predict(b$model, b$data, fit = b$fit)
  swapped <- ferx_predict(b$model, reordered, fit = b$fit)
  key <- function(p) p[order(p$ID, p$TIME), c("ID", "TIME", "PRED")]
  o <- key(original)
  s <- key(swapped)
  rownames(o) <- rownames(s) <- NULL
  expect_equal(s, o)
})

# --- T7: the resolved contrast travels with the fit ---------------------------

tl_nested_model <- function() {
  tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 10.0)\n  theta PLACEBO[STUDY, TIME](0.0, -10.0, 10.0)",
    "TVCL * exp(ETA_CL) + PLACEBO",
    eta = "  omega ETA_CL ~ 0.09", v = "TVV"
  ), ".ferx")
}

test_that("T7: sum_to_zero_within survives a design with more subjects", {
  model <- tl_nested_model()
  data <- tl_write(tl_data3, ".csv")
  fit <- tl_fit(model, data)
  tl <- fit$theta_levels
  # One subject per study and an eta on the parameter that reads the block:
  # auto resolves to sum-to-zero within each study.
  expect_identical(tl$contrast, rep("sum_to_zero_within", 9L))
  expect_identical(tl$group, rep(0:2, each = 3L))
  expect_identical(length(fit$theta), 8L)
  # Two subjects per study. Re-discovered on this design the contrast would
  # resolve differently (study no longer identifies a subject); bound from the
  # fit, it keeps the fit's.
  rows <- strsplit(tl_data3, "\n")[[1]]
  dup <- vapply(rows[-1], function(r) {
    id <- as.integer(sub(",.*", "", r))
    sub("^[0-9]+,", paste0(id + 3L, ","), r)
  }, character(1), USE.NAMES = FALSE)
  design <- tl_write(paste(c(rows, dup), collapse = "\n"), ".csv")
  original <- ferx_predict(model, data, fit = fit)
  doubled <- ferx_predict(model, design, fit = fit)
  expect_true(all(is.finite(doubled$PRED)))
  expect_identical(nrow(doubled), 2L * nrow(original))
  first <- doubled[doubled$ID %in% c("1", "2", "3"), ]
  copies <- doubled[doubled$ID %in% c("4", "5", "6"), ]
  expect_equal(first$PRED, original$PRED)
  expect_equal(copies$PRED, original$PRED)
})

test_that("T7b: from-fit predictions reproduce the fit's own PRED, per contrast", {
  # The contrast picks which level of a group carries no theta (the first for
  # ref, the last for sum-to-zero, none for none), so it has to come back from
  # R exactly: read as sum_to_zero, a ref fit keeps its theta count and reads
  # every level at the wrong position. The fit's sdtab PRED was computed by the
  # engine on the fit's own binding, so it is the oracle.
  data <- tl_write(tl_data, ".csv")
  data3 <- tl_write(tl_data3, ".csv")
  cases <- list(
    sum_to_zero = list(tl_col_model(), data),
    none = list(tl_write(tl_model(
      "  theta PLACEBO[STUDY, TIME, contrast = none](2.0, 0.001, 20.0)",
      "PLACEBO"
    ), ".ferx"), data),
    ref = list(tl_write(tl_model(
      "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[STUDY, TIME, contrast = ref](0.0, -5.0, 5.0)",
      "TVCL + PLACEBO"
    ), ".ferx"), data),
    sum_to_zero_within = list(tl_nested_model(), data3)
  )
  for (contrast in names(cases)) {
    model <- cases[[contrast]][[1]]
    d <- cases[[contrast]][[2]]
    fit <- tl_fit(model, d)
    expect_identical(unique(fit$theta_levels$contrast), contrast)
    # A refusal is caught per case, so every contrast reports on its own.
    pred <- tryCatch(ferx_predict(model, d, fit = fit), error = function(e) e)
    expect_false(inherits(pred, "error"), info = contrast)
    if (!inherits(pred, "error")) {
      expect_equal(pred$PRED, fit$sdtab$PRED, tolerance = 1e-6, info = contrast)
    }
  }
  # A two-iteration fit leaves the ref levels near 0, where the first level
  # pinned and the last one derived agree, so ref is also checked in closed
  # form at distinct values: level 1 is held at exactly 0 and levels 2..6
  # carry the five free thetas.
  ref_model <- cases$ref[[1]]
  fit <- tl_fit(ref_model, data)
  fit$theta[] <- c(2, 0.1, 0.2, 0.3, 0.4, 0.5, 10)
  level <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
  dt <- c(1, 3, 8)
  expected <- c(
    100 / 10 * exp(-cumsum((2 + level[1:3]) / 10 * dt)),
    100 / 10 * exp(-cumsum((2 + level[4:6]) / 10 * dt))
  )
  expect_equal(ferx_predict(ref_model, data, fit = fit)$PRED, expected,
               tolerance = 1e-6)
})

# --- T8: persistence ---------------------------------------------------------

tl_roundtrip <- function(fit) {
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

tl_expect_roundtrip <- function(model, data, fit) {
  fit2 <- tl_roundtrip(fit)
  expect_identical(fit2$theta_levels, fit$theta_levels)
  # Not `identical()`: a .fitrx round trip already moves theta, omega and sigma
  # in their last bit, level block or not (measured at ecccbac: the warfarin
  # fit's predictions move by 1.4e-14), so the reloaded fit is held to that.
  expect_equal(
    ferx_predict(model, data, fit = fit2),
    ferx_predict(model, data, fit = fit),
    tolerance = 1e-12
  )
}

test_that("T8: ferx_save_fit / ferx_load_fit keep theta_levels, per contrast", {
  data <- tl_write(tl_data, ".csv")
  # sum_to_zero
  b <- tl_base()
  tl_expect_roundtrip(b$model, b$data, b$fit)
  # none
  none <- tl_write(tl_model(
    "  theta PLACEBO[STUDY, TIME, contrast = none](2.0, 0.001, 20.0)", "PLACEBO"
  ), ".ferx")
  f_none <- tl_fit(none, data)
  expect_identical(unique(f_none$theta_levels$contrast), "none")
  tl_expect_roundtrip(none, data, f_none)
  # ref
  ref <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[STUDY, TIME, contrast = ref](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO"
  ), ".ferx")
  f_ref <- tl_fit(ref, data)
  expect_identical(unique(f_ref$theta_levels$contrast), "ref")
  expect_identical(which(is.na(f_ref$theta_levels$theta_name)), 1L)
  tl_expect_roundtrip(ref, data, f_ref)
  # sum_to_zero_within
  nested <- tl_nested_model()
  data3 <- tl_write(tl_data3, ".csv")
  f_within <- tl_fit(nested, data3)
  expect_identical(unique(f_within$theta_levels$contrast), "sum_to_zero_within")
  tl_expect_roundtrip(nested, data3, f_within)
})

test_that("T8: a one-level block survives the round trip as a frame", {
  # Every column has length 1, which JSON writing would unbox to a scalar.
  one <- "ID,TIME,DV,EVID,AMT,CMT,RATE,MDV,STUDY
1,0,.,1,100,1,0,1,1
1,1,8.1,0,.,1,0,0,1
2,0,.,1,100,1,0,1,1
2,1,7.4,0,.,1,0,0,1"
  data <- tl_write(one, ".csv")
  model <- tl_write(tl_model(
    "  theta PLACEBO[STUDY, TIME, contrast = none](2.0, 0.001, 20.0)", "PLACEBO"
  ), ".ferx")
  fit <- tl_fit(model, data)
  expect_identical(nrow(fit$theta_levels), 1L)
  tl_expect_roundtrip(model, data, fit)
})

test_that("T8: a model without a level block round-trips an empty frame", {
  # The counted form declares no level block.
  counted <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[6](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO[PLA_IDX]"
  ), ".ferx")
  fit <- tl_fit(counted, tl_write(tl_data, ".csv"))
  expect_identical(nrow(fit$theta_levels), 0L)
  fit2 <- tl_roundtrip(fit)
  expect_identical(fit2$theta_levels, fit$theta_levels)
})

# --- T9: a fit without bindings ------------------------------------------------

test_that("T9: a level-block model with a fit lacking bindings is refused", {
  b <- tl_base()
  fit <- b$fit
  # The shape of a .fitrx written by ferx-core, or saved before #370.
  fit$theta_levels <- NULL
  err <- tryCatch(ferx_simulate(b$model, b$data, fit = fit), error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "this fit carries no theta level bindings", fixed = TRUE)
  expect_match(msg, "written by ferx-core", fixed = TRUE)
  expect_match(msg, "Refit with `ferx_fit()`", fixed = TRUE)
  expect_no_match(msg, "was the model edited", fixed = TRUE)
})

# --- T10: uses without a fit ---------------------------------------------------

test_that("T10: predict and simulate without a fit bind the design's own levels", {
  b <- tl_base()
  pred <- ferx_predict(b$model, b$data)
  expect_identical(nrow(pred), 6L)
  expect_true(all(is.finite(pred$PRED)))
  # The engine's unused-theta check does not see a gathered theta as used and
  # warns about every level - the counted form `PLACEBO[6]` gets the same
  # warning, so it predates #370. Only that warning is muffled.
  sim <- withCallingHandlers(
    ferx_simulate(b$model, b$data, n_sim = 1L),
    warning = function(w) {
      if (grepl("is declared in [parameters] but not referenced",
                conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_identical(nrow(sim), 6L)
  expect_true(all(is.finite(sim$IPRED)))
})

# --- T11: SIR and the standalone covariance step -------------------------------

test_that("T11: ferx_sir and ferx_covariance refuse a level-block fit", {
  b <- tl_base()
  check <- function(expr, who) {
    err <- tryCatch({
      expr
      NULL
    }, error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    expect_match(msg, paste0(who, ": not supported yet"), fixed = TRUE)
    expect_match(msg, "`PLACEBO`", fixed = TRUE)
    expect_match(msg, "FeRx-NLME/ferx-core#1622", fixed = TRUE)
    expect_no_match(msg, "values but this model has")
  }
  check(ferx_sir(tl_with_cov(b$fit), sir_samples = 20L, sir_resamples = 10L),
        "ferx_sir")
  check(ferx_covariance(b$fit), "ferx_covariance")
})

# --- T12: a tampered bundle ----------------------------------------------------

test_that("T12: an unknown contrast token in a bundle is refused by name", {
  b <- tl_base()
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(b$fit, path)
  staging <- tempfile("fitrx")
  dir.create(staging)
  files <- utils::unzip(path, exdir = staging)
  json <- file.path(staging, "fit.json")
  w <- jsonlite::read_json(json, simplifyVector = FALSE)
  w$r_extras$theta_levels$contrast <-
    lapply(w$r_extras$theta_levels$contrast, function(x) "sum_to_zero_wthin")
  jsonlite::write_json(w, json, auto_unbox = TRUE, pretty = TRUE, digits = NA,
                       null = "null", na = "null")
  tampered <- tempfile(fileext = ".fitrx")
  old <- setwd(staging)
  on.exit(setwd(old), add = TRUE)
  utils::zip(tampered, basename(files), flags = "-q")
  setwd(old)
  fit2 <- ferx_load_fit(tampered)
  expect_identical(fit2$theta_levels$contrast, rep("sum_to_zero_wthin", 6L))
  expect_error(
    ferx_predict(b$model, b$data, fit = fit2),
    "unknown contrast `sum_to_zero_wthin`", fixed = TRUE
  )
})

# --- T13: the counted form is untouched ---------------------------------------

test_that("T13: the counted form keeps its names and an empty theta_levels", {
  data <- tl_write(tl_data, ".csv")
  counted <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[6](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO[PLA_IDX]"
  ), ".ferx")
  fit <- tl_fit(counted, data)
  expect_identical(
    names(fit$theta),
    c("TVCL", paste0("PLACEBO[", 1:6, "]"), "TVV")
  )
  expect_identical(nrow(fit$theta_levels), 0L)
  # The counted form drives a from-fit prediction with no bindings at all.
  expect_true(all(is.finite(ferx_predict(counted, data, fit = fit)$PRED)))
})
