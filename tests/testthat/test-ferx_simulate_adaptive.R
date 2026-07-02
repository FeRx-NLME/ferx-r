
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
