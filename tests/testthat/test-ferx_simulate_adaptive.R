
# ---- header from test-adaptive.R ----
# Tests for ferx_simulate_adaptive() - the declarative [adaptive_dosing] /
# feedback-dosing simulation path (ferx-core #585, epic #391).













# Write a minimal 1-cpt IV adaptive model with a custom [adaptive_dosing] rule,
# reusing the bundled dose-free TDM dataset, so the hold/stop outcome branches of
# the result converter can be exercised (the TDM example only produces "dosed").
.write_adaptive_model <- function(rule) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(4.0, 0.5, 20.0)",
    "  theta TVV(50.0, 5.0, 200.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma PROP ~ 0.10 (sd)",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  ode(states=[central])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central",
    "[scaling]",
    "  y = central / V",
    "[error_model]",
    "  DV ~ proportional(PROP)",
    "[adaptive_dosing]",
    "  observe     = central / V",
    "  at          = every 12 from 0 to 48",
    "  start_dose  = 1000",
    "  route       = bolus(cmt=1)",
    "  dose_bounds = [0, 2000]",
    "  confirm     = 1",
    rule
  ), path)
  path
}




test_that("ferx_simulate_adaptive returns four named data frames", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_type(res, "list")
  expect_named(res, c("trajectories", "doses", "decisions", "metrics"))
  expect_s3_class(res$trajectories, "data.frame")
  expect_s3_class(res$doses, "data.frame")
  expect_s3_class(res$decisions, "data.frame")
  expect_s3_class(res$metrics, "data.frame")
})
test_that("the result tables carry the documented columns", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "IPRED", "DV_SIM") %in%
    names(res$trajectories)))
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "AMT", "CMT", "RATE", "DECISION",
    "SIGNAL", "RULE") %in% names(res$doses)))
  expect_true(all(c("DRAW", "SIM", "ID", "DECISION", "TIME", "SIGNAL", "OUTCOME",
    "N_DOSED") %in% names(res$decisions)))
  expect_true(all(c("DRAW", "SIM", "ID", "CUM_DOSE", "N_DOSES", "N_INCREASES",
    "N_DECREASES", "N_HOLDS", "DISCONTINUED", "TIME_TO_DISCONT", "SIGNAL_MIN",
    "SIGNAL_MAX", "SIGNAL_MEAN", "PCT_TIME_IN_WINDOW") %in% names(res$metrics)))
})
test_that("n_sim = 3 produces three replicates across all three tables", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 3L, seed = 1L)
  expect_equal(sort(unique(res$trajectories$SIM)), 1:3)
  # The ledger and decision log must carry the per-replicate SIM too - the
  # single-subject driver emits SIM = 0 and the orchestrator stamps the real
  # replicate index. Every subject doses here, so all three replicates appear.
  expect_equal(sort(unique(res$doses$SIM)), 1:3)
  expect_equal(sort(unique(res$decisions$SIM)), 1:3)
})
test_that("same seed is deterministic", {
  ex <- ferx_example("adaptive_tdm")
  a  <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 7L)
  b  <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 7L)
  expect_equal(a, b)
})
test_that("trajectories carry MIXNUM on a [mixture] model only (ferx-core #2372)", {
  data <- ferx_example("adaptive_tdm")$data
  plain <- ferx_simulate_adaptive(.write_adaptive_model("  when signal < 10 : increase 25%"),
                                  data, n_sim = 1L, seed = 1L)
  expect_false("MIXNUM" %in% names(plain$trajectories))

  # The same model with two classes, P(class 1) = 0.3, class 2 clearing faster.
  path <- tempfile(fileext = ".ferx")
  lines <- readLines(.write_adaptive_model("  when signal < 10 : increase 25%"))
  lines <- sub("  CL = TVCL * exp(ETA_CL)",
               "  CL = if (MIXNUM == 1) TVCL * exp(ETA_CL) else 2 * TVCL * exp(ETA_CL)",
               lines, fixed = TRUE)
  # The rewrite must have landed, or the model has no class-dependent CL.
  expect_true(any(grepl("if (MIXNUM == 1)", lines, fixed = TRUE)))
  lines <- append(lines, c("  theta MIXL(-0.8472979, -10.0, 10.0)"), after = 3L)
  lines <- append(lines, c("[mixture]", "  nsub = 2", "  logit(1) = MIXL"),
                  after = which(lines == "[individual_parameters]") - 1L)
  writeLines(lines, path)

  # The engine calls MIXL "not referenced" (ferx-core #2483); see
  # helper-mixture-lint.R.
  res <- muffle_mixl_lint(ferx_simulate_adaptive(path, data, n_sim = 50L, seed = 1L))
  traj <- res$trajectories
  expect_identical(tail(names(traj), 1L), "MIXNUM")
  expect_type(traj$MIXNUM, "integer")
  per_run <- unique(traj[c("SIM", "ID", "MIXNUM")])
  expect_equal(nrow(per_run), nrow(unique(traj[c("SIM", "ID")])))
  expect_true(all(per_run$MIXNUM %in% 1:2))
  # Five SDs of the class-1 share around 0.3 over these subject-replicates.
  expect_lt(abs(mean(per_run$MIXNUM == 1L) - 0.3), 5 * sqrt(0.21 / nrow(per_run)))
})
test_that("the controller titrates: a sub-therapeutic trough escalates the dose", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 5L, seed = 1L)
  # start_dose is 1000 mg; troughs below 10 mg/L trigger `increase 25%`, so the
  # realized doses must climb above the starting dose and the escalation rule
  # must be named in the ledger's RULE column.
  expect_gt(max(res$doses$AMT), 1000)
  expect_true(any(grepl("increase", res$doses$RULE)))
  # Every decision is logged (incl. holds), so there are at least as many
  # decision rows as dose rows.
  expect_gte(nrow(res$decisions), nrow(res$doses))
})
test_that("a model without an [adaptive_dosing] block is rejected", {
  ex <- ferx_example("warfarin")
  # The engine raises an R error (not a NULL return) when the model has no
  # [adaptive_dosing] block - this entry point requires one.
  expect_error(
    ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L),
    "adaptive_dosing"
  )
})
test_that("a non-positive n_sim is rejected", {
  ex <- ferx_example("adaptive_tdm")
  expect_error(ferx_simulate_adaptive(ex$model, ex$data, n_sim = 0L), "n_sim")
  expect_error(ferx_simulate_adaptive(ex$model, ex$data, n_sim = -1L), "n_sim")
})
test_that("verify = FALSE skips the replay verifier and still returns results", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L,
                                verify = FALSE)
  expect_named(res, c("trajectories", "doses", "decisions", "metrics"))
  expect_gt(nrow(res$doses), 0L)
})
test_that("a non-logical `verify` is rejected before touching the engine", {
  ex <- ferx_example("adaptive_tdm")
  expect_error(ferx_simulate_adaptive(ex$model, ex$data, verify = "yes"), "verify")
  expect_error(ferx_simulate_adaptive(ex$model, ex$data, verify = c(TRUE, FALSE)),
               "verify")
})
test_that("decisions tag every realized dose with a valid outcome (Dosed path)", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_true(all(res$decisions$OUTCOME %in% c("dosed", "hold", "stop")))
  # The TDM titration has no hold/stop rule, so every decision issues a dose;
  # the converter must mark each "dosed" with N_DOSED >= 1.
  expect_true(all(res$decisions$OUTCOME == "dosed"))
  expect_true(all(res$decisions$N_DOSED >= 1L))
  # Every ledger dose joins back to a logged decision.
  expect_true(all(res$doses$DECISION %in% res$decisions$DECISION))
})
test_that("metrics: one row per (subject, draw, replicate), a faithful ledger reduction", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 3L, seed = 1L)
  # Exactly one metrics row per realized (DRAW, SIM, ID) - the same key as the
  # other artifacts.
  key <- res$metrics[, c("DRAW", "SIM", "ID")]
  expect_false(any(duplicated(key)))
  expect_equal(nrow(res$metrics),
    nrow(unique(res$decisions[, c("DRAW", "SIM", "ID")])))
  # CUM_DOSE is the summed ledger AMT and N_DOSES the ledger row count for that
  # key - the metrics summarise the ledger rather than re-deriving anything.
  for (i in seq_len(nrow(res$metrics))) {
    r    <- res$metrics[i, ]
    rows <- res$doses$DRAW == r$DRAW & res$doses$SIM == r$SIM &
      res$doses$ID == r$ID
    expect_equal(r$CUM_DOSE, sum(res$doses$AMT[rows]))
    expect_equal(r$N_DOSES, sum(rows))
  }
})
test_that("metrics: target_window populates PCT_TIME_IN_WINDOW in [0, 1]", {
  # The bundled adaptive_tdm model declares target_window = [10, 20], so the
  # attainment metric is reported (not NA) and is a valid fraction.
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 5L, seed = 1L)
  expect_false(any(is.na(res$metrics$PCT_TIME_IN_WINDOW)))
  expect_true(all(res$metrics$PCT_TIME_IN_WINDOW >= 0 &
    res$metrics$PCT_TIME_IN_WINDOW <= 1))
  # This example only titrates (no hold/stop rule): nothing discontinues, so the
  # discontinuation time is NA throughout.
  expect_true(all(!res$metrics$DISCONTINUED))
  expect_true(all(is.na(res$metrics$TIME_TO_DISCONT)))
})
test_that("a hold rule yields hold decisions and an empty dose ledger", {
  model <- .write_adaptive_model("  when signal < 1000000 : hold")
  data  <- ferx_example("adaptive_tdm")$data
  res   <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L)
  # No drug is ever given, so the signal stays 0 and the controller always holds.
  expect_equal(nrow(res$doses), 0L)
  expect_true(all(res$decisions$OUTCOME == "hold"))
  expect_true(all(res$decisions$N_DOSED == 0L))
  # Metrics reflect the all-hold run: no doses or dose changes, every decision a
  # hold, no discontinuation. This model declares no target_window, so the
  # attainment metric is NA.
  expect_true(all(res$metrics$N_DOSES == 0L))
  expect_true(all(res$metrics$N_INCREASES == 0L & res$metrics$N_DECREASES == 0L))
  expect_true(all(res$metrics$N_HOLDS > 0L))
  expect_true(all(!res$metrics$DISCONTINUED))
  expect_true(all(is.na(res$metrics$PCT_TIME_IN_WINDOW)))
})
test_that("a stop rule discontinues after the first decision", {
  model <- .write_adaptive_model("  when signal < 1000000 : stop")
  data  <- ferx_example("adaptive_tdm")$data
  res   <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L)
  # The stop fires at the first decision (signal 0 < 1e6); nothing is logged after.
  expect_true(all(res$decisions$OUTCOME == "stop"))
  expect_true(all(res$decisions$N_DOSED == 0L))
  expect_equal(nrow(res$doses), 0L)
  # Exactly one decision per subject, then silence.
  expect_equal(nrow(res$decisions), length(unique(res$decisions$ID)))
  # Metrics flag the discontinuation: DISCONTINUED is TRUE and TIME_TO_DISCONT is
  # the stop time - the first decision, t = 0 here - with no dose ever given.
  expect_true(all(res$metrics$DISCONTINUED))
  expect_true(all(res$metrics$TIME_TO_DISCONT == 0))
  expect_true(all(res$metrics$N_DOSES == 0L))
})

# ---- Pre-scheduled base regimen (loading dose) + system resets (ferx-core #702 / #716 / #929) ----

# A vancomycin loading-dose model: 1-cpt IV, decisions every 24 h from 24 to 144,
# titrating the maintenance dose on the (un-noised) latent trough. The base
# regimen rides on the data, not the model. Optional [adaptive_dosing] extras and
# an IOV kappa let the same skeleton exercise the still-rejected combinations.
.write_vanco_model <- function(dosing_extra = character(0),
                               cl_expr = "TVCL * exp(ETA_CL)",
                               extra_params = character(0)) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(4.0, 0.1, 50.0)",
    "  theta TVV(80.0, 5.0, 300.0)",
    "  omega ETA_CL ~ 1e-10",
    extra_params,
    "  sigma ADD ~ 1.0",
    "[individual_parameters]",
    paste0("  CL = ", cl_expr),
    "  V  = TVV",
    "[structural_model]",
    "  ode(states=[CENT])",
    "[odes]",
    "  d/dt(CENT) = -(CL / V) * CENT",
    "[scaling]",
    "  y = CENT / V",
    "[error_model]",
    "  DV ~ additive(ADD)",
    "[adaptive_dosing]",
    "  observe       = CENT / V",
    "  at            = every 24 from 24 to 144",
    "  start_dose    = 750",
    "  route         = bolus(cmt=1)",
    "  dose_bounds   = [250, 4000]",
    dosing_extra
  ), path)
  path
}

# A one-subject dataset carrying a 1500 mg loading dose at t=0, an optional
# reset at t=12 (EVID=3 pure reset, or EVID=4 reset + 300 mg redose), and the
# daily observation grid.
.write_vanco_data <- function(reset = c("none", "evid3", "evid4")) {
  reset <- match.arg(reset)
  rows <- "1,0,0,1500,1,1,1"
  if (reset == "evid3") rows <- c(rows, "1,12,0,0,3,1,1")
  if (reset == "evid4") rows <- c(rows, "1,12,0,300,4,1,1")
  rows <- c(rows, sprintf("1,%d,0,0,0,1,0", seq(24, 144, by = 24)))
  path <- tempfile(fileext = ".csv")
  writeLines(c("ID,TIME,DV,AMT,EVID,CMT,MDV", rows), path)
  path
}

# The controller's observed signal at the first decision (t = 24 here).
.first_signal <- function(res) {
  d <- res$decisions
  d$SIGNAL[which.min(d$TIME)]
}

test_that("a dosed base subject (pre-scheduled loading dose) now simulates", {
  # This used to be rejected outright ("requires a dose-free base subject").
  ex  <- ferx_example("adaptive_vanco_loading")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_named(res, c("trajectories", "doses", "decisions", "metrics"))
  expect_gt(nrow(res$doses), 0L)
})

test_that("the pre-scheduled loading dose is integrated (non-zero day-1 trough)", {
  ex  <- ferx_example("adaptive_vanco_loading")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L)
  # 1500 mg into V = 80 L with ke = CL/V = 0.05/h decays to ~5.6 mg/L by t = 24.
  # A dose-free base subject would show 0 at the first decision; a materially
  # positive trough proves the base dose was integrated.
  s1 <- .first_signal(res)
  expect_gt(s1, 3)
  expect_lt(s1, 9)
})

test_that("the ledger and CUM_DOSE hold controller doses only (base dose excluded)", {
  ex  <- ferx_example("adaptive_vanco_loading")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  # The 1500 mg loading dose is a base dose, so it never appears in the ledger.
  expect_false(any(res$doses$AMT == 1500))
  expect_gt(nrow(res$doses), 0L)
  # CUM_DOSE sums controller doses only - it equals the ledger sum for the key
  # and therefore excludes the loading dose.
  for (i in seq_len(nrow(res$metrics))) {
    r    <- res$metrics[i, ]
    rows <- res$doses$DRAW == r$DRAW & res$doses$SIM == r$SIM &
      res$doses$ID == r$ID
    expect_equal(r$CUM_DOSE, sum(res$doses$AMT[rows]))
  }
})

test_that("an EVID=3 system reset zeroes the state", {
  model     <- ferx_example("adaptive_vanco_loading")$model
  s_none  <- .first_signal(
    ferx_simulate_adaptive(model, .write_vanco_data("none"),  n_sim = 1L, seed = 1L))
  s_reset <- .first_signal(
    ferx_simulate_adaptive(model, .write_vanco_data("evid3"), n_sim = 1L, seed = 1L))
  # The reset at t = 12 zeroes CENT; nothing is dosed between the reset and the
  # t = 24 decision, so the trough collapses to ~0, far below the no-reset ~5.6.
  expect_lt(s_reset, 1)
  expect_gt(s_none, 3)
  expect_lt(s_reset, s_none)
})

test_that("an EVID=4 row resets then doses (reset before dose)", {
  model     <- ferx_example("adaptive_vanco_loading")$model
  s_none  <- .first_signal(
    ferx_simulate_adaptive(model, .write_vanco_data("none"),  n_sim = 1L, seed = 1L))
  s_evid3 <- .first_signal(
    ferx_simulate_adaptive(model, .write_vanco_data("evid3"), n_sim = 1L, seed = 1L))
  s_evid4 <- .first_signal(
    ferx_simulate_adaptive(model, .write_vanco_data("evid4"), n_sim = 1L, seed = 1L))
  # EVID=4 zeroes at t = 12 (like EVID=3) and THEN gives 300 mg, so the t = 24
  # trough sits above the pure-reset case but below the persisting-loading case.
  expect_gt(s_evid4, s_evid3)
  expect_lt(s_evid4, s_none)
})

test_that("unsupported adaptive combinations raise a typed error", {
  # auc_target combined with a system reset is rejected (the per-window trapezoid
  # grid carries no node at the reset instant), not silently mis-scored.
  auc_model <- .write_vanco_model(dosing_extra = c(
    "  target_window = [10, 15]",
    "  auc_target    = [400, 600]",
    "  when signal < 10 : increase 25%"
  ))
  expect_error(
    ferx_simulate_adaptive(auc_model, .write_vanco_data("evid3"), n_sim = 1L, seed = 1L),
    "auc_target"
  )
  # A base regimen combined with a reset UNDER IOV is rejected (base x reset is
  # supported on the constant-covariate path only).
  iov_model <- .write_vanco_model(
    cl_expr      = "TVCL * exp(ETA_CL + KAPPA_CL)",
    extra_params = "  kappa KAPPA_CL ~ 0.09",
    dosing_extra = c(
      "  target_window = [10, 15]",
      "  when signal < 10 : increase 25%"
    ))
  expect_error(
    ferx_simulate_adaptive(iov_model, .write_vanco_data("evid3"), n_sim = 1L, seed = 1L),
    "base regimen combined with system"
  )
})

test_that("a decision on a diverged ODE state is an error naming the state (ferx-core #1578)", {
  # The right-hand side blows up within minutes of the first 1000 mg dose, so
  # the state is NaN at the t = 12 decision. Before #1578 the controller was
  # handed that NaN: with `verify = FALSE` the run returned results, every
  # later decision logged as "dosed" on a NaN signal, and with the replay
  # verifier on it failed on "prediction 1 is NaN" after the fact.
  model <- .write_adaptive_model("  when signal < 10 : increase 25%")
  writeLines(sub("-(CL / V) * central", "(CL / V) * central * central",
                 readLines(model), fixed = TRUE), model)
  data <- ferx_example("adaptive_tdm")$data
  for (verify in c(TRUE, FALSE)) {
    err <- expect_error(
      ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L, verify = verify),
      "decision 1 at t=12: the ODE state is non-finite (central)", fixed = TRUE
    )
    # A runtime failure of the model and its parameter draw, not a refusal of
    # the caller's input: ferx-core gives it no code.
    expect_null(err$code)
  }
})

# --- fit = (ferx-core #1594) --------------------------------------------------

# A one-compartment IV ODE model with initial values `init` (TVCL, TVV, the
# ETA_CL variance, the PROP sd), each written with 17 significant digits so a
# fitted double comes back from the file bit for bit. The readout is Form C,
# `y = central / V`: an `obs_scale` readout is wrong on the adaptive path until
# ferx-core #2142. The controller titrates on the concentration, so a change of
# parameters reaches the dose ledger as well as the trajectories.
.write_fit_adaptive_model <- function(init) {
  v <- sprintf("%.17g", init)
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    paste0("  theta TVCL(", v[1], ", 0.01, 50.0)"),
    paste0("  theta TVV(", v[2], ", 1.0, 500.0)"),
    paste0("  omega ETA_CL ~ ", v[3]),
    paste0("  sigma PROP ~ ", v[4], " (sd)"),
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  ode(states=[central])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central",
    "[scaling]",
    "  y = central / V",
    "[error_model]",
    "  DV ~ proportional(PROP)",
    "[fit_options]",
    "  maxiter = 2",
    "  covariance = false",
    "[adaptive_dosing]",
    "  observe     = central / V",
    "  at          = every 12 from 12 to 36",
    "  start_dose  = 100",
    "  route       = bolus(cmt=1)",
    "  dose_bounds = [0, 400]",
    "  confirm     = 1",
    "  when signal < 2 : increase 25%"
  ), path)
  path
}

# Three subjects, a 100 mg base dose each, observed out to 48 h. The DVs were
# made with CL near 3 and V near 15, away from the initial values the tests
# start from, so two outer iterations move every parameter.
.fit_adaptive_data <- function() {
  path <- tempfile(fileext = ".csv")
  writeLines(c(
    "ID,TIME,DV,EVID,AMT,CMT,MDV",
    "1,0,.,1,100,1,1", "1,1,5.9,0,.,1,0", "1,4,3.5,0,.,1,0",
    "1,12,0.62,0,.,1,0", "1,24,.,0,.,1,0", "1,36,.,0,.,1,0", "1,48,.,0,.,1,0",
    "2,0,.,1,100,1,1", "2,1,6.3,0,.,1,0", "2,4,4.4,0,.,1,0",
    "2,12,1.1,0,.,1,0", "2,24,.,0,.,1,0", "2,36,.,0,.,1,0", "2,48,.,0,.,1,0",
    "3,0,.,1,100,1,1", "3,1,6.1,0,.,1,0", "3,4,3.1,0,.,1,0",
    "3,12,0.48,0,.,1,0", "3,24,.,0,.,1,0", "3,36,.,0,.,1,0", "3,48,.,0,.,1,0"
  ), path)
  path
}

# The 24-48 h rows carry no DV, so a simulation keeps them as design points and
# says so (`W_DESIGN_DV`); that one warning is expected here, any other is not.
.quiet_design_dv <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("W_DESIGN_DV", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# Before ferx-core #1594 the adaptive glue simulated at the model file's
# initial values whatever the caller had fitted. `fit =` must run at the fitted
# ones: the same output as a model file that declares them literally, and not
# the output at the initial values. Mutation that reddens this: the from-fit
# glue passing `&parsed.model.default_params` instead of the fitted params.
test_that("fit = simulates at the fitted estimates, not the initial values", {
  init <- c(2, 10, 0.09, 0.1)
  model <- .write_fit_adaptive_model(init)
  data <- .fit_adaptive_data()
  fit <- ferx_fit(model, data, verbose = FALSE)
  # The fixture tells the two apart only if the fit moved off every initial.
  expect_true(all(fit$theta != init[1:2]))
  expect_true(fit$omega[1, 1] != init[3])
  expect_true(fit$sigma[[1]] != init[4])

  run <- function(model, fit = NULL) {
    .quiet_design_dv(
      ferx_simulate_adaptive(model, data, n_sim = 2L, seed = 7L, fit = fit)
    )
  }
  at_fit <- run(model, fit)
  at_init <- run(model)
  literal <- run(.write_fit_adaptive_model(
    c(fit$theta[["TVCL"]], fit$theta[["TVV"]], fit$omega[1, 1], fit$sigma[[1]])
  ))

  for (tbl in c("trajectories", "doses", "decisions", "metrics")) {
    expect_identical(at_fit[[tbl]], literal[[tbl]], label = tbl)
  }
  expect_false(isTRUE(all.equal(at_fit$trajectories$IPRED,
                                at_init$trajectories$IPRED)))
  expect_false(isTRUE(all.equal(at_fit$doses$AMT, at_init$doses$AMT)))
})

# The adaptive engine draws kappa per decision window (occasion = decision
# index, ferx-core `api/adaptive.rs`) and never reads the data's occasion
# labels, so the fit's IOV occasion rule cannot change the output. A fit that
# records no rule (`fit$iov_occasion` absent, as on a .fitrx saved before
# ferx-core #1783), from a model file that sets none and a design with no
# occasion column, must simulate like the same fit with its rule recorded.
# Mutation that reddens this: the from-fit arm calling `derive_fit_occasions`,
# which refuses the unrecorded fit for want of occasion labels.
test_that("fit = ignores the fit's IOV occasion rule (kappa per decision window)", {
  model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
  text <- readLines(model)
  text <- append(text, "  kappa KAPPA_CL ~ 0.04", after = grep("omega ETA_CL", text))
  text <- sub("CL = TVCL * exp(ETA_CL)", "CL = TVCL * exp(ETA_CL + KAPPA_CL)",
              text, fixed = TRUE)
  writeLines(text, model)
  data <- .fit_adaptive_data()
  # One base dose per subject, so `time(6)` (observations at 1 and 4 h, then
  # 12 h) gives the fit the two occasions it needs to identify kappa.
  fit <- ferx_fit(model, data, verbose = FALSE,
                  settings = list(iov_occasion = "time(6)"))
  expect_identical(fit$iov_occasion, "time(6)")
  legacy <- fit
  legacy$iov_occasion <- NULL

  run <- function(fit) {
    .quiet_design_dv(
      ferx_simulate_adaptive(model, data, n_sim = 2L, seed = 5L, fit = fit)
    )
  }
  recorded <- run(fit)
  expect_identical(run(legacy), recorded)
  # The kappa draw reaches the output: a larger fitted kappa variance moves it.
  wider <- fit
  wider$omega_iov <- fit$omega_iov * 4
  expect_false(isTRUE(all.equal(run(wider)$trajectories$IPRED,
                                recorded$trajectories$IPRED)))
})

test_that("fit = with data omitted simulates the fit's own data", {
  model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
  data <- .fit_adaptive_data()
  fit <- ferx_fit(model, data, verbose = FALSE)
  expect_identical(
    .quiet_design_dv(ferx_simulate_adaptive(model, n_sim = 1L, seed = 3L, fit = fit)),
    .quiet_design_dv(ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 3L, fit = fit))
  )
})

# The mixture_override fixture as an ODE model with an [adaptive_dosing] block.
.write_mixture_adaptive_model <- function() {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL1(1.2, 0.01, 100.0)",
    "  theta TVCL2(2.5, 0.01, 100.0)",
    "  theta TVV(10.0, 0.1, 1000.0)",
    "  theta MIXL(0.0, -10.0, 10.0)",
    "  omega ETA_CL ~ 0.09",
    "  sigma EPS ~ 0.04",
    "[mixture]",
    "  nsub = 2",
    "  logit(1) = MIXL",
    "  omega(2) ETA_CL ~ 0.30",
    "[individual_parameters]",
    "  CL = if (MIXNUM == 1) TVCL1 * exp(ETA_CL) else TVCL2 * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  ode(states=[central])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central",
    "[scaling]",
    "  y = central / V",
    "[error_model]",
    "  DV ~ proportional(EPS)",
    "[adaptive_dosing]",
    "  observe     = central / V",
    "  at          = every 12 from 0 to 48",
    "  start_dose  = 100",
    "  route       = bolus(cmt=1)",
    "  dose_bounds = [0, 200]",
    "  confirm     = 1",
    "  when signal < 1 : increase 25%"
  ), path)
  path
}

test_that("ferx_simulate_adaptive(fit = ) refuses a [mixture] fit without override values, coded (#672)", {
  # The fit list carries no `mixture_overrides` (a .fitrx saved before #655),
  # so the class rebuild refuses before anything is simulated. The adaptive
  # prefix on the message must not cost the refusal its code.
  path <- .write_mixture_adaptive_model()
  fit <- list(theta = c(TVCL1 = 1.2, TVCL2 = 2.5, TVV = 10, MIXL = 0),
              omega = matrix(0.02, 1, 1, dimnames = list("ETA_CL", "ETA_CL")),
              sigma = c(EPS = 0.2))
  probe <- engine_error_probe(muffle_mixl_lint(ferx_simulate_adaptive(
    path, test_path("fixtures", "mixture_iv.csv"), n_sim = 1L, seed = 1L, fit = fit
  )))
  expect_coded_refusal(probe, "omega(2) ETA_CL", "E_FIT_INPUT_MISSING")
  expect_true(startsWith(conditionMessage(probe$cond), "ferx_simulate_adaptive: "),
              info = conditionMessage(probe$cond))
  expect_null(probe$value)
})

# ---- supplied etas: `etas =` / `unmatched =` (ferx-core #1594 part 2) ----

# A 1-cpt IV model with two etas, so a column-order mix-up cannot pass. With
# `adaptive = TRUE` it carries an always-hold controller deciding on the dose
# times, so every dose is a base dose from the data and the adaptive run is the
# model's individual prediction; with `covariate_etas = TRUE` (and no
# controller) the etas are read from the data columns ECL / EV, so
# `ferx_predict()`'s PRED (eta = 0) is the IPRED at those etas.
.write_supplied_eta_model <- function(adaptive = TRUE, covariate_etas = FALSE) {
  path <- tempfile(fileext = ".ferx")
  shift <- function(eta, col) if (covariate_etas) paste0(eta, " + ", col) else eta
  writeLines(c(
    "[parameters]",
    "  theta TVCL(4.0, 0.5, 20.0)",
    "  theta TVV(50.0, 5.0, 200.0)",
    "  omega ETA_CL ~ 0.09",
    "  omega ETA_V ~ 0.04",
    "  sigma PROP ~ 0.10 (sd)",
    "[individual_parameters]",
    paste0("  CL = TVCL * exp(", shift("ETA_CL", "ECL"), ")"),
    paste0("  V  = TVV * exp(", shift("ETA_V", "EV"), ")"),
    "[structural_model]",
    "  ode(states=[central])",
    "[odes]",
    "  d/dt(central) = -(CL / V) * central",
    "[scaling]",
    "  y = central / V",
    "[error_model]",
    "  DV ~ proportional(PROP)",
    if (adaptive) c(
      "[adaptive_dosing]",
      "  observe     = central / V",
      "  at          = [0, 24]",
      "  start_dose  = 100",
      "  route       = bolus(cmt=1)",
      "  dose_bounds = [0, 1000]",
      "  when signal >= 0 : hold"
    )
  ), path)
  path
}

# Three subjects with 100 mg at 0 and 24 h. `etas` (rows keyed by ID) adds the
# ECL / EV columns the covariate twin reads.
.supplied_eta_data <- function(etas = NULL) {
  obs <- c(2, 8, 20, 26, 32, 44)
  rows <- do.call(rbind, lapply(c("1", "2", "3"), function(id) {
    data.frame(
      ID = id,
      TIME = c(0, obs[1:3], 24, obs[4:6]),
      DV = 0,
      EVID = c(1, 0, 0, 0, 1, 0, 0, 0),
      AMT = c(100, NA, NA, NA, 100, NA, NA, NA),
      CMT = 1,
      MDV = c(1, 0, 0, 0, 1, 0, 0, 0)
    )
  }))
  if (!is.null(etas)) {
    rows$ECL <- etas$ETA_CL[match(rows$ID, etas$ID)]
    rows$EV <- etas$ETA_V[match(rows$ID, etas$ID)]
  }
  rows
}

.supplied_etas <- data.frame(
  ID = c("1", "2", "3"),
  ETA_CL = c(0.3, -0.4, 0.05),
  ETA_V = c(-0.2, 0.1, 0.25)
)

# The run's IPRED keyed "ID@TIME", sorted, for comparison across engines.
.ipred_by_key <- function(df, col = "IPRED") {
  key <- paste0(as.character(df$ID), "@", format(df$TIME, digits = 15))
  v <- stats::setNames(df[[col]], key)
  v[order(names(v))]
}

# The baseline a zero-eta table must reproduce: the same model with every
# omega declared `0.0 FIX`, which draws eta as exactly 0 (ferx-core #1605).
# The bundled TDM example titrates on the assay-noised trough, so the dose
# ledger, the decisions and the assay draws all have to line up too, not just
# one trajectory. Comparator: `expect_identical()` on every table (bitwise on
# doubles), except `ETA_SUPPLIED`, which by design says TRUE on one run and
# FALSE on the other.
test_that("etas = all zeros equals the run at eta = 0, bit for bit", {
  ex <- ferx_example("adaptive_tdm")
  zero_model <- tempfile(fileext = ".ferx")
  text <- readLines(ex$model)
  expect_length(grep("omega ETA_CL ~ 0.09", text, fixed = TRUE), 1L)
  writeLines(sub("omega ETA_CL ~ 0.09", "omega ETA_CL ~ 0.0 FIX", text,
                 fixed = TRUE), zero_model)
  ids <- unique(as.character(utils::read.csv(ex$data)$ID))
  zeros <- data.frame(ID = ids, ETA_CL = 0)

  supplied <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 11L,
                                     etas = zeros)
  at_zero <- ferx_simulate_adaptive(zero_model, ex$data, n_sim = 2L, seed = 11L)
  for (tbl in c("trajectories", "doses", "decisions")) {
    expect_identical(supplied[[tbl]], at_zero[[tbl]], label = tbl)
  }
  expect_true(all(supplied$metrics$ETA_SUPPLIED))
  expect_false(any(at_zero$metrics$ETA_SUPPLIED))
  keep <- setdiff(names(supplied$metrics), "ETA_SUPPLIED")
  expect_identical(supplied$metrics[keep], at_zero$metrics[keep])
  # The fixture would not notice etas being ignored unless the drawn run differs.
  drawn <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 11L)
  expect_false(isTRUE(all.equal(drawn$trajectories$IPRED,
                                supplied$trajectories$IPRED)))
})

# The oracle is `ferx_predict()` on a twin that reads each subject's etas from
# data columns: its PRED (eta = 0) is exp(0 + eta) = exp(eta), the same
# arithmetic the supplied run does. `ferx_predict()` itself has no per-subject
# eta argument, hence the twin. Mutations that redden this: the glue dropping
# the supplied etas (the drawn run), or flattening the table column-major.
test_that("supplied etas: each subject's IPRED matches ferx_predict at those etas", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  res <- ferx_simulate_adaptive(model, data, n_sim = 2L, seed = 3L,
                                etas = .supplied_etas)
  expect_equal(nrow(res$doses), 0L) # the controller holds throughout
  expect_true(all(res$metrics$ETA_SUPPLIED))
  # A supplied eta stays fixed across replicates; nothing else is drawn here.
  r1 <- res$trajectories[res$trajectories$SIM == 1L, ]
  r2 <- res$trajectories[res$trajectories$SIM == 2L, ]
  expect_identical(r1$IPRED, r2$IPRED)

  pred <- ferx_predict(.write_supplied_eta_model(adaptive = FALSE,
                                                 covariate_etas = TRUE),
                       .supplied_eta_data(.supplied_etas))
  got <- .ipred_by_key(r1)
  want <- .ipred_by_key(pred, "PRED")
  expect_identical(names(got), names(want))
  # Bit-identical, measured (macOS, 2026-11-01): the decisions sit on the dose
  # times, so both engines integrate the same segments at the same etas.
  expect_identical(unname(got), unname(want))
  # The three subjects differ, so a subject mix-up cannot pass.
  expect_equal(length(unique(round(got, 8))), length(got))
})

test_that("etas = fit$ebe_etas reproduces the fit's own IPRED (from-fit path)", {
  model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
  text <- readLines(model)
  text <- sub("at          = every 12 from 12 to 36", "at          = [0]", text,
              fixed = TRUE)
  text <- sub("when signal < 2 : increase 25%", "when signal >= 0 : hold", text,
              fixed = TRUE)
  writeLines(text, model)
  data <- .fit_adaptive_data()
  fit <- ferx_fit(model, data, verbose = FALSE)
  expect_true(all(abs(fit$ebe_etas$ETA_CL) > 1e-6))

  res <- .quiet_design_dv(
    ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 2L, fit = fit,
                           etas = fit$ebe_etas)
  )
  expect_true(all(res$metrics$ETA_SUPPLIED))
  sd <- fit$sdtab
  got <- .ipred_by_key(res$trajectories)
  want <- .ipred_by_key(sd)
  got <- got[names(want)]
  expect_false(anyNA(got))
  # Bit-identical, measured (macOS, 2026-11-01), as ferx-core's own B1 oracle
  # (`tests/adaptive_supplied_etas.rs`) asserts.
  expect_identical(unname(got), unname(want))
  # With data omitted it simulates the fit's own data, as `fit =` alone does.
  expect_identical(
    .quiet_design_dv(ferx_simulate_adaptive(model, n_sim = 1L, seed = 2L,
                                            fit = fit, etas = fit$ebe_etas)),
    res
  )
})

test_that("a data subject without a row: \"refuse\" is coded, \"draw\" draws it alone", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  partial <- .supplied_etas[.supplied_etas$ID != "3", ]
  run <- function(...) {
    ferx_simulate_adaptive(model, data, n_sim = 2L, seed = 9L, ...)
  }

  probe <- engine_error_probe(run(etas = partial))
  expect_coded_refusal(probe, "have no supplied", "E_SUPPLIED_ETAS")
  expect_true(grepl("'3'", conditionMessage(probe$cond), fixed = TRUE),
              info = conditionMessage(probe$cond))
  # The advice is spelled the R way, never as the engine's Rust enum.
  expect_true(grepl("`unmatched = \"draw\"`", conditionMessage(probe$cond),
                    fixed = TRUE),
              info = conditionMessage(probe$cond))
  expect_false(grepl("UnmatchedSubjects", conditionMessage(probe$cond),
                     fixed = TRUE),
               info = conditionMessage(probe$cond))
  expect_false(isTRUE(grepl("UnmatchedSubjects", probe$cond$suggestion,
                            fixed = TRUE)))
  expect_identical(
    engine_error_probe(run(etas = partial, unmatched = "refuse"))$cond$code,
    "E_SUPPLIED_ETAS"
  )
  # The same on the `fit =` route, which goes through the other glue entry.
  fit <- list(theta = c(TVCL = 4, TVV = 50),
              omega = diag(c(0.09, 0.04)), sigma = c(PROP = 0.1))
  from_fit <- engine_error_probe(run(fit = fit, etas = partial))
  expect_coded_refusal(from_fit, "have no supplied", "E_SUPPLIED_ETAS")
  expect_true(grepl("`unmatched = \"draw\"`", conditionMessage(from_fit$cond),
                    fixed = TRUE),
              info = conditionMessage(from_fit$cond))
  expect_false(grepl("UnmatchedSubjects", conditionMessage(from_fit$cond),
                     fixed = TRUE),
               info = conditionMessage(from_fit$cond))

  # Two of three matched: the intended mix, so no no-match warning.
  expect_no_warning(mixed <- run(etas = partial, unmatched = "draw"))
  m <- mixed$metrics
  expect_identical(m$ETA_SUPPLIED, m$ID != "3")
  traj <- function(res, ids) {
    t <- res$trajectories[res$trajectories$ID %in% ids, ]
    rownames(t) <- NULL
    t
  }
  # Subjects 1 and 2 use their rows: as in the run that supplies everyone.
  expect_identical(traj(mixed, c("1", "2")),
                   traj(run(etas = .supplied_etas), c("1", "2")))
  # Subject 3 draws: the supplied subjects still consume their draws, so it
  # gets the eta the same seed gives with no etas at all.
  expect_identical(traj(mixed, "3"), traj(run(), "3"))
})

# Review r1 #1: under "draw" an `etas` keyed differently from the data matches
# nobody, the mistake "refuse" exists to catch. The run goes ahead (every
# subject drawn, as asked) but says so, naming IDs from both sides.
test_that("\"draw\" with etas matching no data subject warns, naming both ID sets", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  padded <- .supplied_etas
  padded$ID <- sprintf("%02d", as.integer(padded$ID))
  expect_warning(
    res <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 4L,
                                  etas = padded, unmatched = "draw"),
    "no row of `etas` matched a subject in `data`.*`etas` has 01, 02, 03; `data` has 1, 2, 3"
  )
  expect_false(any(res$metrics$ETA_SUPPLIED))
  # The run itself is the all-drawn one.
  expect_identical(res$trajectories,
                   ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 4L)$trajectories)
})

# Review r1 #6: an explicit `unmatched` with no `etas` cannot do anything.
test_that("`unmatched` without `etas` warns; the default does not", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  expect_warning(
    ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L, unmatched = "draw"),
    "`unmatched` has no effect without `etas`"
  )
  expect_no_warning(ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L))
  expect_no_warning(ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L,
                                           etas = .supplied_etas,
                                           unmatched = "refuse"))
})

test_that("eta columns that are not the model's are refused with the engine's code", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  refuse <- function(etas) {
    engine_error_probe(ferx_simulate_adaptive(model, data, n_sim = 1L,
                                              seed = 1L, etas = etas))
  }
  e <- .supplied_etas
  wrong <- stats::setNames(e, c("ID", "ETA_CL", "ETA_VC"))
  expect_coded_refusal(refuse(wrong), "[ETA_CL, ETA_VC]", "E_SUPPLIED_ETAS")
  expect_coded_refusal(refuse(e[c("ID", "ETA_CL")]), "[ETA_CL]",
                       "E_SUPPLIED_ETAS")
  # Same names, another order: refused, not silently reordered.
  expect_coded_refusal(refuse(e[c("ID", "ETA_V", "ETA_CL")]), "[ETA_V, ETA_CL]",
                       "E_SUPPLIED_ETAS")
  # The table checks are the engine's as well.
  dup <- rbind(e, e[1, ])
  expect_coded_refusal(refuse(dup), "appears more than once", "E_SUPPLIED_ETAS")
  bad <- e
  bad$ETA_CL[2] <- NA
  expect_coded_refusal(refuse(bad), "non-finite", "E_SUPPLIED_ETAS")
  # The from-fit path refuses the same way.
  fit <- list(theta = c(TVCL = 4, TVV = 50),
              omega = diag(c(0.09, 0.04)), sigma = c(PROP = 0.1))
  probe <- engine_error_probe(ferx_simulate_adaptive(
    model, data, n_sim = 1L, seed = 1L, fit = fit, etas = wrong
  ))
  expect_coded_refusal(probe, "[ETA_CL, ETA_VC]", "E_SUPPLIED_ETAS")
})

test_that("supplied etas on a [mixture] model are refused, coded", {
  ids <- unique(as.character(
    utils::read.csv(test_path("fixtures", "mixture_iv.csv"))$ID
  ))
  probe <- engine_error_probe(muffle_mixl_lint(ferx_simulate_adaptive(
    .write_mixture_adaptive_model(), test_path("fixtures", "mixture_iv.csv"),
    n_sim = 1L, seed = 1L, etas = data.frame(ID = ids, ETA_CL = 0)
  )))
  expect_coded_refusal(probe, "[mixture]", "E_SUPPLIED_ETAS")
})

test_that("a malformed `etas` is refused in R before the engine", {
  model <- .write_supplied_eta_model()
  data <- .supplied_eta_data()
  run <- function(etas, ...) {
    ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L, etas = etas, ...)
  }
  expect_error(run(as.matrix(.supplied_etas[-1])), "must be a data frame")
  expect_error(run(.supplied_etas[-1]), "with an ID")
  na_id <- .supplied_etas
  na_id$ID[1] <- NA
  expect_error(run(na_id), "NA")
  chr <- .supplied_etas
  chr$ETA_V <- as.character(chr$ETA_V)
  expect_error(run(chr), "ETA_V are not numeric")
  expect_error(run(.supplied_etas, unmatched = "skip"), "should be one of")
  # Review r1 #4: a repeated eta column would pass the engine's name check on
  # its first copy and drop the second.
  twice <- cbind(.supplied_etas, data.frame(ETA_V = 0.7))
  expect_error(run(twice), "more than one column named `ETA_V`")
  # Integer IDs match the data's text IDs.
  num <- .supplied_etas
  num$ID <- as.integer(num$ID)
  expect_identical(run(num), run(.supplied_etas))
})

# Review r1 #3 / #5: a double ID is matched as a data.frame `data` spells it
# (`.fitrx_double_text()`). The IDs are chosen where the plausible spellings
# part: `as.character(1e5)` is "1e+05", `format(1e15, digits = 15)` is
# "1000000000000000" and `format(0.1 + 0.2, digits = 15)` "0.3", while the data
# writer gives "100000", "1e+15" and "0.30000000000000004". Under "refuse" a
# mis-spelled ID is a refusal, so the run passing is the check.
test_that("double IDs in etas match a data.frame data's double IDs", {
  model <- .write_supplied_eta_model()
  map <- c("1" = 1e5, "2" = 1e15, "3" = 0.1 + 0.2)
  data <- .supplied_eta_data()
  data$ID <- unname(map[data$ID])
  etas <- .supplied_etas
  etas$ID <- unname(map[etas$ID])
  res <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L, etas = etas)
  expect_true(all(res$metrics$ETA_SUPPLIED))
  expect_setequal(res$metrics$ID, c("100000", "1e+15", "0.30000000000000004"))
  # Each subject on its own row: the IPRED of the character-ID run.
  chr <- ferx_simulate_adaptive(model, .supplied_eta_data(), n_sim = 1L,
                                seed = 1L, etas = .supplied_etas)
  expect_identical(res$trajectories$IPRED, chr$trajectories$IPRED)
})
