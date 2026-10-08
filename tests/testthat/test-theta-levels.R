# theta NAME[COL, ...] level blocks from R (#370).
#
# A level block declares one theta per observed combination of its columns, so
# its theta layout is a property of the data it is bound to. ferx_fit() binds
# it to the fitted data; every later use of the fit binds the design against
# the fit's layout (`fit$theta_levels`), so a theta is never read at a position
# the fit did not give it. The fixtures are ferx-core's own
# (tests/theta_level_blocks.rs): two studies x TIME {1, 4, 12}, one subject per
# study, plus a PLA_IDX column holding the same design in the counted form. As
# in core since #1675, the default random effect sits on a parameter `y` never
# reads (see tl_model()).
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
#
# By default the random effect sits on `Z`, a parameter `y` never reads, as in
# ferx-core's fixture since #1675: with one subject per study the [STUDY, TIME]
# block takes a level at every observation, so an eta that reached `y` would be
# absorbed by it - `contrast = auto` then resolves to sum_to_zero_within and an
# explicit none / sum_to_zero / ref is refused (T14 pins that). `Z` is reported
# in [derived], which does not reach `y`, so the unused-parameter check stays
# quiet. ETA_V then has no effect on the OFV and its omega is not identified;
# it is kept only so the fit has an eta at all, since a fit without one cannot
# make the persistence round trips T8 and T12 take (FeRx-NLME/ferx-r#461).
# Pass `z = NULL` to drop both.
tl_model <- function(thetas, cl,
                     structure = "[structural_model]\n  pk one_cpt_iv(cl=CL, v=V)\n",
                     eta = "  omega ETA_V ~ 0.04", v = "TVV",
                     z = "TVV * exp(ETA_V)") {
  paste0(
    "[parameters]\n", thetas, "\n",
    "  theta TVV(10.0, 0.1, 500.0)\n",
    eta, "\n",
    "  sigma PROP_ERR ~ 0.05\n\n",
    "[individual_parameters]\n",
    "  CL = ", cl, "\n",
    "  V  = ", v, "\n",
    if (!is.null(z)) paste0("  Z  = ", z, "\n"), "\n",
    structure, "\n",
    if (!is.null(z)) "[derived]\n  Z_OUT = Z\n\n",
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
  # Measured bit-identical (OFV -1.792364950549102 for both, theta max abs
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

# Every from-fit level refusal carries the engine's code (ferx-core #1791, the
# leftover cell of ferx-r #498): before #1791 the binders returned a bare
# `String`, so these were plain errors - re-validating the model and design
# finds nothing wrong, since without the fit the design's own levels bind.
# SIR and the covariance step reach the binders through `layout_from_fit`
# (R8c's missing block); the adaptive path binds the design's own levels
# (T14b). Mutation that reddens this: a
# glue binder that formats the `EngineError` to text (`e.to_string()`)
# instead of handing it to `engine_refusal`.
# The refusal `expr` raises, checked to carry the level-block code; its message.
tl_coded <- function(expr, who) {
  e <- tryCatch({
    expr
    NULL
  }, error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_THETA_LEVEL_BINDING", info = who)
  expect_identical(e$block, "parameters", info = who)
  msg <- if (inherits(e, "condition")) conditionMessage(e) else ""
  expect_match(msg, "[E_THETA_LEVEL_BINDING]", fixed = TRUE, info = who)
  msg
}

test_that("T5b: from-fit level refusals carry E_THETA_LEVEL_BINDING on every path", {
  b <- tl_base()
  design <- tl_write(gsub(",12,", ",24,", tl_data, fixed = TRUE), ".csv")
  coded <- tl_coded
  coded(ferx_predict(b$model, design, fit = b$fit), "predict")
  coded(ferx_simulate(b$model, design, fit = b$fit), "simulate")
  coded(ferx_simulate_with_uncertainty(b$model, design, tl_with_cov(b$fit),
                                       n_uncertainty_draws = 2L),
        "simulate_with_uncertainty")
  coded(ferx_calc_npde(b$fit, nsim = 10L, model = b$model, data = design),
        "npde")
  coded(ferx_predict_survival(b$model, design, times = c(1, 2), fit = b$fit),
        "predict_survival")

  model <- tl_write(tl_model(
    paste0("  theta TVCL(2.0, 0.001, 20.0)\n",
           "  theta PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)\n",
           "  theta VSHIFT[STUDY](0.0, -5.0, 5.0)"),
    "TVCL + PLACEBO", v = "TVV * exp(VSHIFT)"
  ), ".ferx")
  fit <- tl_with_cov(tl_fit(model, tl_write(tl_data, ".csv")))
  fit$theta_levels <- fit$theta_levels[fit$theta_levels$block != "VSHIFT", ]
  for (who in c("ferx_covariance", "ferx_sir")) {
    msg <- coded(if (who == "ferx_sir") {
      ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L)
    } else {
      ferx_covariance(fit)
    }, who)
    expect_match(msg, "the fit's level bindings carry no `VSHIFT`", fixed = TRUE,
                 info = who)
  }
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
    eta = "  omega ETA_CL ~ 0.09", v = "TVV", z = NULL
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
  # A .fitrx round trip gives back every double bit for bit (#415). PRED reads
  # theta alone, so omega and sigma are checked directly; sigma unnamed, as a
  # reload names it and the fit does not (#417).
  expect_identical(fit2$theta, fit$theta)
  expect_identical(fit2$omega, fit$omega)
  expect_identical(unname(fit2$sigma), unname(fit$sigma))
  expect_identical(
    ferx_predict(model, data, fit = fit2),
    ferx_predict(model, data, fit = fit)
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
  data <- tl_write(tl_data, ".csv")
  fit <- tl_fit(counted, data)
  expect_identical(nrow(fit$theta_levels), 0L)
  tl_expect_roundtrip(counted, data, fit)
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

test_that("T9b: ferx_sir / ferx_covariance give the predict paths' refusal", {
  b <- tl_base()
  fit <- tl_with_cov(b$fit)
  fit$theta_levels <- NULL
  msg_of <- function(expr) {
    tryCatch({
      expr
      NA_character_
    }, error = function(e) conditionMessage(e))
  }
  # Byte for byte the text T9 anchors: one function writes it for every path.
  sim <- msg_of(ferx_simulate(b$model, b$data, fit = fit))
  expect_match(sim, "this fit carries no theta level bindings", fixed = TRUE)
  expect_identical(msg_of(ferx_predict(b$model, b$data, fit = fit)), sim)
  expect_identical(msg_of(ferx_covariance(fit)), sim)
  expect_identical(
    msg_of(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L)), sim
  )
  expect_no_match(sim, "theta length", fixed = TRUE)
  expect_no_match(sim, "expected)", fixed = TRUE)
})

test_that("T9c: bindings that no longer lay out the fit's theta are named", {
  b <- tl_base()
  fit <- tl_with_cov(b$fit)
  # An edited fit: one level row gone, so the layout is one theta short.
  fit$theta_levels <- fit$theta_levels[-nrow(fit$theta_levels), ]
  check <- function(expr, who) {
    err <- tryCatch({
      expr
      NULL
    }, error = function(e) e)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    # The count is the whole model's (TVCL, TVV and four free levels), not the
    # block's alone.
    expect_match(msg, paste0(
      who, ": the model file, laid out on the fit's theta level bindings ",
      "(`fit$theta_levels`) for its level block(s) `PLACEBO`, has 6 thetas, ",
      "but the fit carries 7 (`fit$theta`). Either the model file was edited ",
      "since the fit, or `fit$theta_levels` was"
    ), fixed = TRUE)
    expect_no_match(msg, "does not match model", fixed = TRUE)
  }
  check(ferx_covariance(fit), "ferx_covariance")
  check(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L), "ferx_sir")
})

test_that("T9d: a model file edited after the fit is named as a possible cause", {
  # The model file is checked against `fit$model_hash` before the theta
  # count (#492), so with a hash the edit is reported as an edit. Only a fit
  # that carries no hash reaches the count check, and there an edited model
  # file must not be blamed on `fit$theta_levels` alone.
  b <- tl_base()
  model <- tl_col_model()
  fit <- tl_with_cov(tl_fit(model, b$data))
  text <- readLines(model)
  at <- which(text == "[parameters]")
  writeLines(append(text, "  theta EXTRA(1.0, 0.1, 10.0)", after = at), model)
  msg_of <- function(expr) {
    err <- tryCatch({
      expr
      NULL
    }, error = function(e) e)
    expect_s3_class(err, "error")
    conditionMessage(err)
  }
  calls <- list(
    ferx_covariance = function(f) ferx_covariance(f),
    ferx_sir = function(f) ferx_sir(f, sir_samples = 20L, sir_resamples = 10L)
  )
  no_hash <- fit
  no_hash$model_hash <- NA_character_
  for (who in names(calls)) {
    msg <- msg_of(calls[[who]](fit))
    expect_match(msg, paste0(who, ": model hash mismatch for "), fixed = TRUE)
    expect_no_match(msg, "laid out", fixed = TRUE)

    msg <- msg_of(suppressWarnings(calls[[who]](no_hash)))
    expect_match(msg, paste0(who, ": the model file, laid out"), fixed = TRUE)
    expect_match(msg, "has 8 thetas, but the fit carries 7", fixed = TRUE)
    expect_match(msg, "Either the model file was edited since the fit", fixed = TRUE)
  }
})

# --- T10: uses without a fit ---------------------------------------------------

test_that("T10: predict and simulate without a fit bind the design's own levels", {
  b <- tl_base()
  pred <- ferx_predict(b$model, b$data)
  expect_identical(nrow(pred), 6L)
  expect_true(all(is.finite(pred$PRED)))
  # A gathered theta counts as used (FeRx-NLME/ferx-core#1628): before it,
  # the unused-theta check warned about every level, here and for the
  # counted form `PLACEBO[6]` alike.
  expect_no_warning(sim <- ferx_simulate(b$model, b$data, n_sim = 1L))
  expect_identical(nrow(sim), 6L)
  expect_true(all(is.finite(sim$IPRED)))
})

# --- T11: SIR and the standalone covariance step -------------------------------

#
# Both run on the fit's own level layout (#463): the glue lays the model out on
# `fit$theta_levels` and hands the bindings to the engine, which binds the
# re-read data to them. The oracle is a twin: the bundled mbma_placebo data
# with no random effect reaching `y`, once as
# `PLACEBO[STUDY, TIME, contrast = none]` and once as the counted `PLACEBO[24]`
# read through a `PLA_IDX` column. The counted form binds no level block, so it
# ran through both entry points before #463; the two are the same model, so
# every number must agree bit for bit. An eta or kappa reaching `y` would be
# absorbed by a [STUDY, TIME] block under contrast = none
# (FeRx-NLME/ferx-core#1675), and tl_data's six observations are too few for
# an identified covariance step. `eta = TRUE` adds a FIX eta on a parameter `y`
# never reads, as tl_model() does, so the fit can take a .fitrx round trip
# (FeRx-NLME/ferx-r#461) without moving any number.

tl_twin_model <- function(placebo, base, extra_covariate = "", eta = FALSE) {
  paste0(
    "[parameters]\n",
    "  theta EMAX(8.0, 0.0, 100.0)\n",
    "  theta ED50(20.0, 0.1, 1000.0)\n",
    "  theta ET50(2.0, 0.01, 100.0)\n",
    "  theta ", placebo, "(45.0, 0.0, 200.0)\n",
    if (eta) "  omega ETA_Z ~ 0.04 FIX\n",
    "  sigma ADD_ERR ~ 1.0 (variance) FIX\n\n",
    "[covariates]\n",
    "  STUDY categorical\n",
    "  NARM  continuous\n",
    "  SE    continuous\n",
    "  DOSE  continuous\n",
    extra_covariate, "\n",
    "[individual_parameters]\n",
    "  BASE = ", base, "\n",
    if (eta) "  Z = EMAX * exp(ETA_Z)\n", "\n",
    "[structural_model]\n",
    "  DRUG = EMAX * DOSE / (ED50 + DOSE) * TIME / (TIME + ET50)\n",
    "  y    = BASE - DRUG\n\n",
    if (eta) "[derived]\n  Z_OUT = Z\n\n",
    "[error_model]\n",
    "  DV ~ additive(ADD_ERR) weight = SE\n\n",
    "[fit_options]\n",
    "  method     = focei\n",
    "  covariance = true\n"
  )
}

tl_twin_data <- function() {
  if (is.null(tl_cache$twin_data)) {
    d <- utils::read.csv(ferx_example("mbma_placebo")$data)
    lab <- paste0("STUDY=", d$STUDY, ",TIME=", d$TIME)
    d$PLA_IDX <- match(lab, unique(lab))
    path <- tempfile(fileext = ".csv")
    utils::write.csv(d, path, row.names = FALSE, quote = FALSE)
    tl_cache$twin_data <- path
  }
  tl_cache$twin_data
}

tl_twin_col_model <- function(eta = FALSE) {
  tl_write(tl_twin_model("PLACEBO[STUDY, TIME, contrast = none]", "PLACEBO",
                         eta = eta), ".ferx")
}

# Both twin fits, shared by T11a and T11b.
tl_twin <- function() {
  if (is.null(tl_cache$twin)) {
    data <- tl_twin_data()
    n_levels <- max(utils::read.csv(data)$PLA_IDX)
    counted <- tl_write(tl_twin_model(
      sprintf("PLACEBO[%d]", n_levels), "PLACEBO[PLA_IDX]",
      "  PLA_IDX continuous\n"
    ), ".ferx")
    tl_cache$twin <- list(
      col = tl_fit(tl_twin_col_model(), data),
      cnt = tl_fit(counted, data)
    )
  }
  tl_cache$twin
}

test_that("T11a: ferx_covariance on a level-block fit matches the counted form", {
  tw <- tl_twin()
  expect_identical(nrow(tw$col$theta_levels), 24L)
  expect_identical(nrow(tw$cnt$theta_levels), 0L)
  expect_identical(tw$col$ofv, tw$cnt$ofv)
  c_col <- ferx_covariance(tw$col)
  c_cnt <- ferx_covariance(tw$cnt)
  expect_true(all(is.finite(c_col$se_theta)))
  expect_identical(unname(c_col$cov_matrix), unname(c_cnt$cov_matrix))
  expect_identical(unname(c_col$se_theta), unname(c_cnt$se_theta))
  # Each also agrees with its own fit's in-fit covariance step. Bit for bit on
  # the macOS FD build, but 3e-12 apart (relative) on CI's Linux build: the
  # standalone step rebuilds Omega from `fit$omega`. A layout the fit never
  # had moves these SEs by up to 3593x.
  expect_equal(unname(c_col$se_theta), unname(tw$col$se_theta), tolerance = 1e-10)
  expect_equal(unname(c_cnt$se_theta), unname(tw$cnt$se_theta), tolerance = 1e-10)
  expect_identical(names(c_col$theta), names(tw$col$theta))
})

test_that("T11b: ferx_sir on a level-block fit matches the counted form", {
  tw <- tl_twin()
  sir <- function(fit) {
    ferx_sir(fit, sir_samples = 300L, sir_resamples = 100L, sir_seed = 5L)
  }
  s_col <- sir(tw$col)
  s_cnt <- sir(tw$cnt)
  expect_true(is.finite(s_col$sir_ess))
  expect_identical(s_col$sir_ess, s_cnt$sir_ess)
  expect_identical(unname(as.matrix(s_col$sir_ci_theta)),
                   unname(as.matrix(s_cnt$sir_ci_theta)))
})

test_that("T11c: ferx_covariance on the bundled mbma_placebo fit matches the in-fit step", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  fit <- ferx_fit(ex$model, ex$data, verbose = FALSE)
  expect_gt(nrow(fit$theta_levels), 0L)
  sa <- ferx_covariance(fit)
  # Not bit for bit: the standalone step rebuilds Omega from `fit$omega`, and
  # this fit's free placebo levels are ill-conditioned, which amplifies the
  # re-decomposition. Measured 1.5e-5 (theta), 2.4e-6 (omega), 4.0e-7 (kappa);
  # a layout the fit never had moves SEs by orders of magnitude.
  rel <- function(a, b) max(abs(unname(a) - unname(b)) / abs(unname(b)))
  expect_lt(rel(sa$se_theta, fit$se_theta), 1e-4)
  expect_lt(rel(sa$se_omega, fit$se_omega), 1e-4)
  expect_lt(rel(sa$se_kappa, fit$se_kappa), 1e-4)
})

test_that("T11d: ferx_sir runs on a block whose levels straddle 0", {
  # FeRx-NLME/ferx-core#1701: SIR used to reject every sample holding a
  # theta <= 0, whatever that theta's bounds, so a level block centred on 0
  # never passed ("All N SIR samples had invalid weights"). A theta now only
  # has to sit inside its declared bounds. With 20 samples the ESS is ~3, so
  # this pins that SIR runs and keeps negative levels negative, not the
  # interval widths (the bundled mbma_placebo run in the bump PR measures
  # those against the Wald intervals).
  b <- tl_base()
  s <- ferx_sir(tl_with_cov(b$fit), sir_samples = 20L, sir_resamples = 10L,
                sir_seed = 1L)
  # An all-rejected SIR is an engine error, so reaching this line is the
  # #1701 check; the intervals below are what the old engine could not give.
  ci <- s$sir_ci_theta
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[, "lower"] <= ci[, "upper"]))
  lvl <- grepl("^PLACEBO\\[", rownames(ci))
  neg <- lvl & b$fit$theta[rownames(ci)] < 0
  expect_gt(sum(neg), 0L)
  expect_true(all(ci[neg, "lower"] < 0))
  expect_true(all(ci[neg, "upper"] < 0))
})

test_that("T8b: a reloaded level-block fit runs the covariance step identically", {
  # R's .fitrx carries the bindings under r_extras (#466 moves them to the
  # native slot), and that is enough for the standalone step.
  fit <- tl_fit(tl_twin_col_model(eta = TRUE), tl_twin_data())
  fit2 <- tl_roundtrip(fit)
  expect_identical(fit2$theta_levels, fit$theta_levels)
  c1 <- ferx_covariance(fit)
  c2 <- ferx_covariance(fit2)
  expect_true(all(is.finite(c1$se_theta)))
  expect_identical(c2$cov_matrix, c1$cov_matrix)
  expect_identical(c2$se_theta, c1$se_theta)
  # The unread FIX eta moves nothing: the T11a twin's SEs (its in-fit step,
  # so the same 1e-10 band as T11a).
  expect_equal(unname(c1$se_theta), unname(tl_twin()$col$se_theta), tolerance = 1e-10)
})

test_that("R7: bindings on a model without a level block are refused", {
  b <- tl_base()
  counted <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[6](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO[PLA_IDX]"
  ), ".ferx")
  fit <- tl_with_cov(tl_fit(counted, b$data))
  # Stray rows from another fit: before #463 both steps ignored them.
  fit$theta_levels <- b$fit$theta_levels
  check <- function(expr) {
    err <- tryCatch({
      expr
      NULL
    }, error = function(e) e)
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err),
                 "carry the block(s) `PLACEBO`, which this model does not declare",
                 fixed = TRUE)
  }
  check(ferx_covariance(fit))
  check(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L))
})

# --- R8: the skeleton is laid out by core's layout_from_fit (ferx-r #469) ------

test_that("R8a: the SIR / covariance skeleton keeps the fit's theta count, names and FIX flags", {
  # A FIX theta declared *after* the level block, so its position in the FIX
  # mask depends on the layout: last of the bound thetas, but third of the
  # unbound parse's, which carries no PLACEBO theta at all. Declared before the
  # block, as the twin model has it, ET50 is third in both layouts and no mask
  # check can tell them apart (ferx-r #488).
  #
  # Only `cov_fixed` and the `cov_matrix` dimnames come from the glue's
  # layout. `names(cv$theta)`, `estimates$fixed` and the SIR row names are
  # built in R from `fit$theta` / `fit$theta_fixed`, so they can only check
  # that the glue's count agreed; they are kept for that.
  #
  # Mutations that redden this (ferx-r #488, each built into a scratch
  # library): drop the `layout_from_fit` call in bind_layout_from_fit() (the
  # count check refuses); return the unbound parse's FIX mask, padded with free
  # entries to the bound count (`cov_fixed` flags the third theta; green before
  # ET50 moved below the block); return the unbound parse's mask unpadded (no
  # mask comes back); label the theta rows from the unbound parse, padded with
  # `THETA<i>` (dimnames); reverse the level labels before the layout
  # (dimnames).
  placebo <- "PLACEBO[STUDY, TIME, contrast = none]"
  text <- sub(
    paste0("  theta ET50(2.0, 0.01, 100.0)\n  theta ", placebo, "(45.0, 0.0, 200.0)\n"),
    paste0("  theta ", placebo, "(45.0, 0.0, 200.0)\n  theta ET50(2.0, FIX)\n"),
    tl_twin_model(placebo, "PLACEBO"),
    fixed = TRUE
  )
  fit <- tl_fit(tl_write(text, ".ferx"), tl_twin_data())
  expect_gt(nrow(fit$theta_levels), 1L)
  n_theta <- length(fit$theta)
  expect_identical(names(fit$theta)[n_theta], "ET50")
  expect_identical(unname(which(fit$cov_fixed[seq_len(n_theta)])), n_theta)
  cv <- ferx_covariance(fit)
  expect_identical(names(cv$theta), names(fit$theta))
  expect_identical(rownames(cv$cov_matrix), rownames(fit$cov_matrix))
  expect_identical(cv$cov_fixed, fit$cov_fixed)
  expect_identical(cv$estimates$fixed, fit$estimates$fixed)
  s <- ferx_sir(fit, sir_samples = 50L, sir_resamples = 20L, sir_seed = 1L)
  expect_identical(rownames(s$sir_ci_theta), names(fit$theta))
})

test_that("R8b: malformed level bindings are refused before the skeleton is built", {
  # Before ferx-r #469 the glue re-parsed on the edited bindings itself: a split group
  # failed in the parser's words ("is not contiguous"), and an `auto` contrast
  # was laid out and refused only later, by the engine (`run_sir:` /
  # `run_covariance:`). Now core's layout_from_fit, the validation
  # bind_from_fit runs, refuses both up front. Mutation that reddens this:
  # restore the hand-copied `parse_full_model_with` re-parse.
  refused <- function(fit, text) {
    for (who in c("ferx_covariance", "ferx_sir")) {
      msg <- tryCatch({
        if (who == "ferx_sir") {
          ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L)
        } else {
          ferx_covariance(fit)
        }
        NA_character_
      }, error = function(e) conditionMessage(e))
      expect_match(msg, text, fixed = TRUE, info = who)
      expect_no_match(msg, "run_sir:", fixed = TRUE, info = who)
      expect_no_match(msg, "run_covariance:", fixed = TRUE, info = who)
      expect_no_match(msg, "is not contiguous", fixed = TRUE, info = who)
    }
  }
  # A split contrast group, on the nested (sum_to_zero_within) fit.
  nested <- tl_with_cov(tl_fit(tl_nested_model(), tl_write(tl_data3, ".csv")))
  expect_identical(nested$theta_levels$group, rep(0:2, each = 3L))
  split <- nested
  split$theta_levels$group <- c(0L, 0L, 1L, 0L, 1L, 1L, 2L, 2L, 2L)
  refused(split, paste0(
    "theta PLACEBO[STUDY, TIME]: the fit's level bindings are malformed: the ",
    "levels of contrast group 0 are split"
  ))
  # `auto` recorded as the contrast: a fit records the contrast it resolved to.
  auto <- tl_with_cov(tl_base()$fit)
  auto$theta_levels$contrast <- "auto"
  refused(auto, paste0(
    "theta PLACEBO[STUDY, TIME]: the fit's level bindings are malformed: they ",
    "record the contrast `auto`"
  ))
  # A repeated label never reaches core: the glue's own table check names it
  # (as for the predict paths, R1). Mutation that reddens this: delete the
  # label loop in level_bindings_from_r(), and core's wording comes back.
  dup <- tl_with_cov(tl_base()$fit)
  dup$theta_levels$label[2] <- dup$theta_levels$label[1]
  refused(dup, "gives the label `STUDY=1,TIME=1` to levels 1 and 2")
})

test_that("R8c: a fit that lost one block's rows names the missing block", {
  # The case ferx-r #469 was filed for. Mutation that reddens this: restore
  # the hand-copied re-parse, which never checks that every declared block is
  # bound.
  model <- tl_write(tl_model(
    paste0("  theta TVCL(2.0, 0.001, 20.0)\n",
           "  theta PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)\n",
           "  theta VSHIFT[STUDY](0.0, -5.0, 5.0)"),
    "TVCL + PLACEBO", v = "TVV * exp(VSHIFT)"
  ), ".ferx")
  fit <- tl_with_cov(tl_fit(model, tl_write(tl_data, ".csv")))
  expect_setequal(unique(fit$theta_levels$block), c("PLACEBO", "VSHIFT"))
  fit$theta_levels <- fit$theta_levels[fit$theta_levels$block != "VSHIFT", ]
  text <- paste0("theta VSHIFT[STUDY]: the fit's level bindings carry no `VSHIFT`, ",
                 "so there is no fitted layout to bind the design against ",
                 "(was the model edited since the fit?)")
  expect_error(ferx_covariance(fit), text, fixed = TRUE)
  expect_error(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L), text,
               fixed = TRUE)
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

test_that("T13b: a plain fit with no bindings drives every from-fit path", {
  # The control for the one cell R refuses itself (#487): empty level bindings
  # are refused only on a model that declares a level block. A model that
  # declares nothing data-derived reaches core, which binds nothing.
  data <- tl_write(tl_data, ".csv")
  counted <- tl_write(tl_model(
    "  theta TVCL(2.0, 0.001, 20.0)\n  theta PLACEBO[6](0.0, -5.0, 5.0)",
    "TVCL + PLACEBO[PLA_IDX]"
  ), ".ferx")
  fit <- tl_with_cov(tl_fit(counted, data))
  expect_identical(nrow(fit$theta_levels), 0L)
  expect_true(all(is.finite(ferx_predict(counted, data, fit = fit)$PRED)))
  sim <- ferx_simulate(counted, data, fit = fit, n_sim = 1L, seed = 1L)
  expect_gt(nrow(sim), 0L)
  expect_s3_class(ferx_covariance(fit), "ferx_fit")
  expect_s3_class(ferx_sir(fit, sir_samples = 20L, sir_resamples = 10L), "ferx_fit")
})

# --- T14: an eta the block can reproduce (FeRx-NLME/ferx-core#1675) -----------

# The fixture before #1675: ETA_V on `V`, which `y` reads as central / V. With
# one subject per study, PLACEBO[STUDY, TIME] takes a level at every
# observation, so it can reproduce any per-subject effect; only a contrast
# that sums to zero within each study keeps the two apart.
tl_absorbing_model <- function(block) {
  tl_write(tl_model(
    paste0("  theta TVCL(2.0, 0.001, 20.0)\n  theta ", block),
    "TVCL + PLACEBO", v = "TVV * exp(ETA_V)", z = NULL
  ), ".ferx")
}

test_that("T14: auto resolves to sum_to_zero_within next to an eta that reaches y", {
  data <- tl_write(tl_data, ".csv")
  fit <- tl_fit(tl_absorbing_model("PLACEBO[STUDY, TIME](0.0, -5.0, 5.0)"), data)
  tl <- fit$theta_levels
  expect_identical(tl$contrast, rep("sum_to_zero_within", 6L))
  expect_identical(tl$group, c(0L, 0L, 0L, 1L, 1L, 1L))
  # One level per study is derived from the others.
  expect_identical(which(is.na(tl$theta_name)), c(3L, 6L))
})

test_that("T14: an explicit global contrast next to such an eta is refused", {
  data <- tl_write(tl_data, ".csv")
  for (contrast in c("none", "sum_to_zero", "ref")) {
    model <- tl_absorbing_model(sprintf(
      "PLACEBO[STUDY, TIME, contrast = %s](0.0, -5.0, 5.0)", contrast
    ))
    e <- tryCatch(tl_fit(model, data), error = function(e) e)
    expect_s3_class(e, "ferx_engine_error")
    expect_identical(e$code, "E_THETA_LEVEL_BINDING", label = contrast)
    expect_match(conditionMessage(e), "ETA_V", fixed = TRUE, label = contrast)
    expect_match(conditionMessage(e), "sum_to_zero_within", fixed = TRUE,
                 label = contrast)
  }
})

# The adaptive path binds the design's own levels (`bind_design`) and is
# refused there with the code too (ferx-core #1791). Its glue puts its own name
# in front of the engine's text, which must not cost the code.
test_that("T14b: ferx_simulate_adaptive() refuses the global contrast with its code", {
  model <- tl_absorbing_model("PLACEBO[STUDY, TIME, contrast = none](0.0, -5.0, 5.0)")
  msg <- tl_coded(ferx_simulate_adaptive(with_adaptive_block(model),
                                         tl_write(tl_data, ".csv"),
                                         n_sim = 1L, seed = 1L),
                  "simulate_adaptive")
  expect_match(msg, "^ferx_simulate_adaptive: ")
  expect_match(msg, "sum_to_zero_within", fixed = TRUE)
})

# --- R1: a tampered theta_levels is refused by the glue -------------------------

test_that("R1: the glue refuses a theta_levels table it cannot trust", {
  b <- tl_base()
  refused <- function(tl, pattern) {
    fit <- b$fit
    fit$theta_levels <- tl
    expect_error(ferx_predict(b$model, b$data, fit = fit), pattern, fixed = TRUE)
  }
  tl <- b$fit$theta_levels
  # Columns that are not parallel (a list, since a data frame cannot be).
  ragged <- as.list(tl)
  ragged$label <- ragged$label[-1]
  refused(ragged, "columns have 6, 6, 5, 6 and 6 rows; they must be parallel")
  # A level missing from the block.
  refused(tl[-3, ], "has the level indices [1, 2, 4, 5, 6]; they must be exactly 1..5")
  # An NA index, shown as R shows it.
  na_index <- tl
  na_index$index[2] <- NA_integer_
  refused(na_index, "has the level indices [NA, 1, 3, 4, 5, 6]")
  # Two levels under one label: the engine would blame an unseen level.
  dup <- tl
  dup$label[2] <- dup$label[1]
  refused(dup, "gives the label `STUDY=1,TIME=1` to levels 1 and 2")
  # A contrast that changes inside the block.
  mixed <- tl
  mixed$contrast[2] <- "none"
  refused(mixed, "mixes the contrasts")
  # A missing or negative group.
  na_group <- tl
  na_group$group[2] <- NA_integer_
  refused(na_group, "level 2 has the group NA")
  neg_group <- tl
  neg_group$group[2] <- -1L
  refused(neg_group, "level 2 has the group -1")
})

# --- R6: a malformed bundle entry is refused by the loader ----------------------

test_that("R6: the loader refuses a malformed r_extras$theta_levels", {
  wire <- list(
    block = list("P", "P"), index = list(1L, 2L), label = list("A=1", "A=2"),
    group = list(0L, 0L), contrast = list("none", "none"),
    theta_name = list("P[A=1]", NULL)
  )
  from_wire <- ferx:::.fitrx_theta_levels_from_wire
  # The well-formed entry, with JSON null read back as NA.
  ok <- from_wire(wire)
  expect_identical(ok$theta_name, c("P[A=1]", NA))
  expect_identical(ok$index, 1:2)
  missing_col <- wire
  missing_col$contrast <- NULL
  expect_error(from_wire(missing_col), "expected the columns", fixed = TRUE)
  ragged <- wire
  ragged$label <- list("A=1")
  expect_error(from_wire(ragged), "the columns differ in length", fixed = TRUE)
  null_label <- wire
  null_label$label <- list("A=1", NULL)
  expect_error(from_wire(null_label), "column `label` has a null", fixed = TRUE)
  fractional <- wire
  fractional$index <- list(1, 2.5)
  expect_error(from_wire(fractional),
               "column `index` has an entry that is not a single integer",
               fixed = TRUE)
})

# --- R2: search tools' final fits carry theta_levels ---------------------------

test_that("R2: a search tool's final fit carries theta_levels and drives predict", {
  skip_on_cran()
  b <- tl_base()
  # The search gates on a converged input fit, which the two-iteration fixture
  # is not built to give. It used to pass the gate only because ETA_V sat on V,
  # where the block absorbed it and its omega ran off to 1 in two iterations
  # (FeRx-NLME/ferx-core#1649); with the eta off `y` (#1675) the fit needs a
  # real budget to converge.
  budget <- "  maxiter = 2\n  inner_maxiter = 3\n"
  text <- paste(readLines(b$model), collapse = "\n")
  # Fail loudly if tl_fit_options moves on, rather than search from the
  # two-iteration fit and report the gate.
  stopifnot(grepl(budget, text, fixed = TRUE))
  model <- tl_write(sub(
    budget, "  maxiter = 200\n  inner_maxiter = 50\n", text, fixed = TRUE
  ), ".ferx")
  res <- ferx_ruvsearch(model, b$data, progress = FALSE)
  expect_s3_class(res$fit, "ferx_fit")
  expect_identical(nrow(res$fit$theta_levels), 6L)
  expect_identical(res$fit$theta_levels$label, b$fit$theta_levels$label)
  pred <- ferx_predict(res$final_model_path, b$data, fit = res$fit)
  expect_identical(nrow(pred), 6L)
  expect_true(all(is.finite(pred$PRED)))
})
