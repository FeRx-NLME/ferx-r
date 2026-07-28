
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
