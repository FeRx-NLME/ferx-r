
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




test_that("ferx_simulate_adaptive returns its named data frames", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_type(res, "list")
  expect_named(res, c("trajectories", "doses", "decisions", "metrics", "parameters",
                      "summary_rows", "summary_bands", "prediction_intervals",
                      "variance_components", "variance_skipped"))
  expect_s3_class(res$trajectories, "data.frame")
  expect_s3_class(res$doses, "data.frame")
  expect_s3_class(res$decisions, "data.frame")
  expect_s3_class(res$metrics, "data.frame")
  expect_s3_class(res$parameters, "data.frame")
})
test_that("the result tables carry the documented columns", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "IPRED", "DV_SIM") %in%
    names(res$trajectories)))
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "AMT", "CMT", "RATE", "DECISION",
    "SIGNAL", "LATENT", "RULE") %in% names(res$doses)))
  expect_true(all(c("DRAW", "SIM", "ID", "DECISION", "TIME", "SIGNAL", "LATENT", "OUTCOME",
    "N_DOSED", "INTERVAL") %in% names(res$decisions)))
  expect_true(all(c("DRAW", "SIM", "ID", "CUM_DOSE", "N_DOSES", "N_INCREASES",
    "N_DECREASES", "N_HOLDS", "N_NOT_DUE", "DISCONTINUED", "TIME_TO_DISCONT", "SIGNAL_MIN",
    "SIGNAL_MAX", "SIGNAL_MEAN", "PCT_READINGS_IN_WINDOW", "PCT_TIME_IN_WINDOW") %in%
    names(res$metrics)))
  # The assay replicate (ferx-core #946 PR3) follows SIM; 1 without n_reps.
  for (tab in c("doses", "decisions", "metrics")) {
    nm <- names(res[[tab]])
    expect_identical(nm[match("SIM", nm) + 1L], "REP", info = tab)
    expect_true(all(res[[tab]]$REP == 1L), info = tab)
  }
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
  expect_named(res, c("trajectories", "doses", "decisions", "metrics", "parameters",
                      "summary_rows", "summary_bands", "prediction_intervals",
                      "variance_components", "variance_skipped"))
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
  expect_true(all(res$decisions$OUTCOME %in% c("dosed", "hold", "not_due", "stop")))
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
test_that("metrics: target_window populates PCT_READINGS_IN_WINDOW in [0, 1]", {
  # The bundled adaptive_tdm model declares target_window = [10, 20], so the
  # attainment metric is reported (not NA) and is a valid fraction.
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 5L, seed = 1L)
  expect_false(any(is.na(res$metrics$PCT_READINGS_IN_WINDOW)))
  expect_true(all(res$metrics$PCT_READINGS_IN_WINDOW >= 0 &
    res$metrics$PCT_READINGS_IN_WINDOW <= 1))
  # This example only titrates (no hold/stop rule): nothing discontinues, so the
  # discontinuation time is NA throughout.
  expect_true(all(!res$metrics$DISCONTINUED))
  expect_true(all(is.na(res$metrics$TIME_TO_DISCONT)))
})
# ferx-core #946: `pct_time_in_window` was renamed `pct_readings_in_window` (it
# counts decision readings, never time). The R column follows; the old name stays
# one release as a deprecated duplicate, so it must be the same vector.
test_that("PCT_TIME_IN_WINDOW is a deprecated duplicate of PCT_READINGS_IN_WINDOW", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 3L, seed = 1L)
  expect_identical(res$metrics$PCT_TIME_IN_WINDOW,
                   res$metrics$PCT_READINGS_IN_WINDOW)
})

# ferx-core #946: the outcome metrics are scored on the latent (noise-free) monitor
# value, not on the assay-noised reading the controller saw. adaptive_tdm titrates
# on a `with_assay_error` (Dv) monitor, so LATENT and SIGNAL differ there, and the
# metrics must be recomputable from LATENT and NOT from SIGNAL.
test_that("Dv monitor: metrics are scored on LATENT, not on the noised SIGNAL", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 3L, seed = 1L)
  d   <- res$decisions
  expect_true(all(is.finite(d$LATENT)))
  # The assay draw really moved the reading: otherwise the checks below could not
  # tell LATENT from SIGNAL.
  expect_true(any(d$LATENT != d$SIGNAL))
  # The ledger's LATENT is the decision's LATENT at the issuing decision.
  key_d <- paste(d$DRAW, d$SIM, d$ID, d$DECISION)
  key_l <- paste(res$doses$DRAW, res$doses$SIM, res$doses$ID, res$doses$DECISION)
  expect_identical(res$doses$LATENT, d$LATENT[match(key_l, key_d)])
  expect_identical(res$doses$SIGNAL, d$SIGNAL[match(key_l, key_d)])
  n_signal_scored <- 0L
  for (i in seq_len(nrow(res$metrics))) {
    r    <- res$metrics[i, ]
    rows <- d$DRAW == r$DRAW & d$SIM == r$SIM & d$ID == r$ID
    lat  <- d$LATENT[rows]
    sig  <- d$SIGNAL[rows]
    expect_equal(r$SIGNAL_MIN, min(lat))
    expect_equal(r$SIGNAL_MAX, max(lat))
    expect_equal(r$SIGNAL_MEAN, mean(lat))
    expect_equal(r$PCT_READINGS_IN_WINDOW, mean(lat >= 10 & lat <= 20))
    if (!isTRUE(all.equal(mean(sig), mean(lat)))) {
      expect_false(isTRUE(all.equal(r$SIGNAL_MEAN, mean(sig))))
      n_signal_scored <- n_signal_scored + 1L
    }
  }
  # At least one run where scoring on SIGNAL would have given another mean.
  expect_gt(n_signal_scored, 0L)
})

# An Ipred (latent `observe`) monitor has no assay draw: LATENT is SIGNAL, bit for bit.
test_that("Ipred monitor: LATENT is identical to SIGNAL", {
  ex  <- ferx_example("adaptive_vanco_loading")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_gt(nrow(res$decisions), 0L)
  expect_identical(res$decisions$LATENT, res$decisions$SIGNAL)
  expect_identical(res$doses$LATENT, res$doses$SIGNAL)
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
  expect_true(all(is.na(res$metrics$PCT_READINGS_IN_WINDOW)))
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
test_that("a 24 h interval on a 12 h grid logs every second point not_due (ferx-core #945)", {
  # The review grid stays at 0, 12, 24, 36, 48; a 24 h interval makes only 0, 24
  # and 48 due. The rule never fires, so each due point re-issues the dose.
  model <- .write_adaptive_model(c(
    "  intervals      = [24]",
    "  start_interval = 24",
    "  when signal > 1000000 : hold"
  ))
  data <- ferx_example("adaptive_tdm")$data
  res  <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L)
  d    <- res$decisions
  due  <- d$TIME %% 24 == 0
  expect_true(all(d$OUTCOME[due] == "dosed"))
  expect_true(all(d$OUTCOME[!due] == "not_due"))
  expect_true(all(d$N_DOSED[!due] == 0L))
  expect_true(all(d$INTERVAL == 24))
  # Not-due points still record what the monitor read.
  expect_false(anyNA(d$SIGNAL[!due]))
  # Doses only at the due points.
  expect_gt(nrow(res$doses), 0L)
  expect_true(all(res$doses$TIME %% 24 == 0))
  # N_NOT_DUE counts the not-due points; N_HOLDS does not.
  m <- res$metrics
  n_not_due <- tapply(!due, d$ID, sum)
  expect_equal(m$N_NOT_DUE, as.integer(n_not_due[as.character(m$ID)]))
  expect_true(all(m$N_NOT_DUE == 2L))
  expect_true(all(m$N_HOLDS == 0L))
  expect_true(all(m$N_DOSES == 3L))
})
test_that("decisions$INTERVAL is NA when the block declares no intervals", {
  model <- .write_adaptive_model("  when signal > 1000000 : hold")
  data  <- ferx_example("adaptive_tdm")$data
  res   <- ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L)
  expect_true(all(is.na(res$decisions$INTERVAL)))
  expect_true(all(res$metrics$N_NOT_DUE == 0L))
})
test_that("intervals without start_interval is a parse error naming start_interval", {
  model <- .write_adaptive_model(c(
    "  intervals = [24]",
    "  when signal > 1000000 : hold"
  ))
  data <- ferx_example("adaptive_tdm")$data
  expect_error(ferx_simulate_adaptive(model, data, n_sim = 1L, seed = 1L),
               "start_interval")
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
  expect_named(res, c("trajectories", "doses", "decisions", "metrics", "parameters",
                      "summary_rows", "summary_bands", "prediction_intervals",
                      "variance_components", "variance_skipped"))
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

test_that("auc_target runs on a subject with a system reset (ferx-core #943)", {
  # Refused until ferx-core #2546: the window AUC is now integrated on the
  # run's own trajectory, with a node on each side of the reset.
  auc_model <- .write_vanco_model(dosing_extra = c(
    "  target_window = [10, 15]",
    "  auc_target    = [400, 600]",
    "  when signal < 10 : increase 25%"
  ))
  res <- ferx_simulate_adaptive(auc_model, .write_vanco_data("evid3"),
                                n_sim = 1L, seed = 1L)
  att <- res$metrics$AUC_TARGET_ATTAINMENT
  expect_length(att, 1L)
  expect_true(is.finite(att) && att >= 0 && att <= 1)
})

test_that("AUC_TARGET_ATTAINMENT is NA without an auc_target", {
  model <- .write_vanco_model(dosing_extra = c(
    "  target_window = [10, 15]",
    "  when signal < 10 : increase 25%"
  ))
  res <- ferx_simulate_adaptive(model, .write_vanco_data("none"),
                                n_sim = 1L, seed = 1L)
  expect_identical(res$metrics$AUC_TARGET_ATTAINMENT, NA_real_)
})

test_that("a base regimen with a reset runs under IOV (ferx-core #932)", {
  # Refused with E_ADAPTIVE_UNSUPPORTED before ferx-core #2556.
  iov_model <- .write_vanco_model(
    cl_expr      = "TVCL * exp(ETA_CL + KAPPA_CL)",
    extra_params = "  kappa KAPPA_CL ~ 0.09",
    dosing_extra = c(
      "  target_window = [10, 15]",
      "  when signal < 10 : increase 25%"
    ))
  s_none <- .first_signal(
    ferx_simulate_adaptive(iov_model, .write_vanco_data("none"), n_sim = 1L, seed = 1L))
  res <- ferx_simulate_adaptive(iov_model, .write_vanco_data("evid3"),
                                n_sim = 1L, seed = 1L)
  # The t = 12 reset zeroes CENT under IOV as on the constant path, so the
  # t = 24 trough is ~0. Without the reset the integrated base dose leaves a
  # material trough, which a driver that dropped the base regimen would not.
  expect_lt(.first_signal(res), 1)
  expect_gt(s_none, 3)
  expect_gt(nrow(res$doses), 0L)
})

test_that("an adaptive model the driver cannot simulate raises a typed error", {
  # A closed-form model (no [odes]) is refused before any subject runs.
  model <- .write_adaptive_model("  when signal > 1000000 : hold")
  text  <- readLines(model)
  text  <- sub("  ode(states=[central])",
               "  pk one_cpt_iv(cl=CL, v=V)", text, fixed = TRUE)
  text  <- text[!(text %in% c("[odes]", "  d/dt(central) = -(CL / V) * central",
                              "[scaling]", "  y = central / V"))]
  text  <- sub("observe     = central / V", "observe     = C", text, fixed = TRUE)
  writeLines(text, model)
  expect_error(
    ferx_simulate_adaptive(model, ferx_example("adaptive_tdm")$data,
                           n_sim = 1L, seed = 1L),
    "requires an ODE model.*E_ADAPTIVE_UNSUPPORTED",
    class = "ferx_engine_error"
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
  legacy <- unsealed_fit(legacy)
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
  wider <- unsealed_fit(wider)
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

# The engine builds today's refusals with no suggestion, so the suggestion half
# of the swap is pinned on a condition made here (#690 review r1 #2): for the
# day the engine moves its advice there.
test_that(".ferx_engine_spelling() swaps message and suggestion, for its code only", {
  refusal <- function(code, suggestion) {
    structure(
      class = c("ferx_engine_error", "error", "condition"),
      list(message = "or pass `X::Y` [E_A]", call = NULL, code = code,
           suggestion = suggestion)
    )
  }
  spell <- function(cond) {
    tryCatch(
      ferx:::.ferx_engine_spelling(stop(cond), "E_A", "`X::Y`", "`y = \"z\"`"),
      error = function(e) e
    )
  }
  hit <- spell(refusal("E_A", "pass `X::Y` instead"))
  expect_s3_class(hit, "ferx_engine_error")
  expect_identical(hit$code, "E_A")
  expect_identical(conditionMessage(hit), "or pass `y = \"z\"` [E_A]")
  expect_identical(hit$suggestion, "pass `y = \"z\"` instead")
  expect_identical(spell(refusal("E_A", NA_character_))$suggestion, NA_character_)
  other <- spell(refusal("E_B", "pass `X::Y` instead"))
  expect_identical(conditionMessage(other), "or pass `X::Y` [E_A]")
  expect_identical(other$suggestion, "pass `X::Y` instead")

  # Several pairs in one call, each swapped (replaces nested calls).
  cond <- structure(
    class = c("ferx_engine_error", "error", "condition"),
    list(message = "use `A::B` or `C::D` [E_A]", call = NULL, code = "E_A",
         suggestion = "`C::D` here")
  )
  both <- tryCatch(
    ferx:::.ferx_engine_spelling(stop(cond), "E_A", c("`A::B`", "`C::D`"),
                                 c("`b`", "`d`")),
    error = function(e) e
  )
  expect_identical(conditionMessage(both), "use `b` or `d` [E_A]")
  expect_identical(both$suggestion, "`d` here")
  expect_error(ferx:::.ferx_engine_spelling(1, "E_A", c("a", "b"), "c"))
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

# --- uncertainty = (ferx-core #946) -------------------------------------------

# A fit of `.write_fit_adaptive_model()` with its covariance step, TVV FIX, run
# to convergence so the covariance matrix exists. Fitted once per file.
.unc_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
      txt <- readLines(model)
      txt <- sub("maxiter = 2", "maxiter = 200", txt, fixed = TRUE)
      txt <- sub("covariance = false", "covariance = true", txt, fixed = TRUE)
      txt <- sub("(theta TVV\\(.*\\))", "\\1 FIX", txt)
      writeLines(txt, model)
      data <- .fit_adaptive_data()
      cache <<- list(model = model, data = data,
                     fit = ferx_fit(model, data, verbose = FALSE))
    }
    cache
  }
})

.run_unc <- function(..., n_sim = 2L, seed = 7L) {
  fx <- .unc_fixture()
  .quiet_design_dv(ferx_simulate_adaptive(fx$model, fx$data, n_sim = n_sim,
                                          seed = seed, fit = fx$fit, ...))
}

.adaptive_tables <- c("trajectories", "doses", "decisions", "metrics")

# The contrast (`reference = "fixed"`, ferx-core #946 PR4) on the from-fit path
# with parameter uncertainty: both arms are simulated under the same parameter
# draws, and the adaptive arm is the run without a reference (#785 review,
# finding 3). Mutation that reddens this: the reference arm built without
# `with_uncertainty` (its DRAW is 1 throughout, and core refuses the pair).
test_that("uncertainty: reference = \"fixed\" runs both arms under the same draws", {
  one <- .run_unc(uncertainty = list(n_draws = 3L))
  con <- .run_unc(uncertainty = list(n_draws = 3L), reference = "fixed")
  arm <- function(df, a) {
    out <- df[df$ARM == a, setdiff(names(df), "ARM"), drop = FALSE]
    rownames(out) <- NULL
    out
  }
  for (tbl in c(.adaptive_tables, "parameters", "summary_rows")) {
    expect_equal(arm(con[[tbl]], "adaptive"), one[[tbl]], ignore_attr = TRUE,
                 label = tbl)
  }
  expect_identical(arm(con$parameters, "reference"), arm(con$parameters, "adaptive"))
  for (tbl in .adaptive_tables) {
    expect_identical(sort(unique(arm(con[[tbl]], "reference")$DRAW)), 1:3,
                     label = tbl)
  }
  diff_rows <- arm(con$summary_rows, "difference")
  expect_identical(sort(unique(diff_rows$DRAW)), 1:3)
  bands <- arm(con$summary_bands, "difference")
  expect_true(all(bands$STATUS == "quantiles"))
})

# Mutations that redden this: the glue dropping `with_uncertainty` (every row
# DRAW 1, no parameters), or handing the engine another n_draws than R's.
test_that("uncertainty: n_draws = 3 tags every table DRAW 1..3, one parameters row per draw", {
  res <- .run_unc(uncertainty = list(n_draws = 3L))
  for (tbl in .adaptive_tables) {
    expect_identical(sort(unique(res[[tbl]]$DRAW)), 1:3, label = tbl)
  }
  expect_identical(res$parameters$DRAW, 1:3)
  # Named as `fit$cov_matrix` names its rows: one spelling per parameter.
  expect_identical(names(res$parameters),
                   c("DRAW", rownames(.unc_fixture()$fit$cov_matrix)))
})

# The draws come from their own stream, so adding draws must not reshuffle the
# first. Mutation that reddens this: seeding the parameter draws from the eta
# stream, or from a seed that mixes in n_draws.
test_that("uncertainty: draw 1 does not depend on n_draws", {
  one <- .run_unc(uncertainty = list(n_draws = 1L))
  three <- .run_unc(uncertainty = list(n_draws = 3L))
  for (tbl in c(.adaptive_tables, "parameters")) {
    first <- three[[tbl]][three[[tbl]]$DRAW == 1L, , drop = FALSE]
    rownames(first) <- NULL
    expect_identical(first, one[[tbl]], label = tbl)
  }
})

# Mutations that redden this: a glue that returns the point estimate for every
# draw, or drops the held-at-bound / FIX handling and scatters TVV.
test_that("uncertainty: theta differs across draws, a FIX theta stays at its estimate", {
  fit <- .unc_fixture()$fit
  p <- .run_unc(uncertainty = list(n_draws = 3L))$parameters
  expect_identical(length(unique(p$TVCL)), 3L)
  expect_false(any(p$TVCL == fit$theta[["TVCL"]]))
  expect_identical(p$TVV, rep(fit$theta[["TVV"]], 3L))
})

# The Omega and sigma columns must carry each draw's values, not the model
# file's (ETA_CL 0.09, PROP 0.1). Mutation that reddens this (#763 review,
# finding 1): `packed_parameter_entries` reading Omega / sigma from the
# template instead of the draw.
test_that("uncertainty: the Omega and sigma columns carry each draw's values", {
  p <- .run_unc(uncertainty = list(n_draws = 3L))$parameters
  expect_identical(length(unique(p$ETA_CL)), 3L)
  expect_identical(length(unique(p$PROP)), 3L)
  expect_false(any(p$ETA_CL == 0.09))
  expect_false(any(p$PROP == 0.1))
})

test_that("uncertainty: method = \"sir\" draws from kept SIR samples, refused without them", {
  fx <- .unc_fixture()
  sir_fit <- ferx_sir(fx$fit, sir_samples = 200L, sir_resamples = 50L,
                      sir_keep_samples = TRUE, sir_seed = 1L)
  res <- .quiet_design_dv(ferx_simulate_adaptive(
    fx$model, fx$data, n_sim = 1L, seed = 7L, fit = sir_fit,
    uncertainty = list(method = "sir", n_draws = 2L)
  ))
  expect_identical(res$parameters$DRAW, 1:2)
  expect_identical(sort(unique(res$metrics$DRAW)), 1:2)
  expect_error(.run_unc(uncertainty = list(method = "sir", n_draws = 2L)),
               "sir_keep_samples = TRUE", fixed = TRUE)
})

# A theta named DRAW must not shadow the key: the key stays the one `DRAW`
# column and the parameter becomes `DRAW.1`, so the documented merge holds.
# Mutation that reddens this (#763 review, finding 6): `drawn_parameters_to_df`
# without the rename, which gives two `DRAW` columns.
test_that("uncertainty: a parameter named DRAW is renamed DRAW.1, the key stays unique", {
  model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
  txt <- gsub("TVV", "DRAW", readLines(model), fixed = TRUE)
  txt <- sub("maxiter = 2", "maxiter = 200", txt, fixed = TRUE)
  txt <- sub("covariance = false", "covariance = true", txt, fixed = TRUE)
  writeLines(txt, model)
  data <- .fit_adaptive_data()
  fit <- ferx_fit(model, data, verbose = FALSE)
  res <- .quiet_design_dv(ferx_simulate_adaptive(
    model, data, n_sim = 1L, seed = 1L, fit = fit,
    uncertainty = list(n_draws = 2L)
  ))
  p <- res$parameters
  expect_identical(names(p), c("DRAW", "TVCL", "DRAW.1", "ETA_CL", "PROP"))
  expect_identical(p$DRAW, 1:2)
  expect_true(all(p$DRAW.1 > 1 & p$DRAW.1 < 500))
  joined <- merge(res$metrics, p, by = "DRAW")
  expect_identical(nrow(joined), nrow(res$metrics))
})

# A fit without a covariance step cannot draw asymptotically. The refusal is
# this function's, not the internal helper's (#763 review, finding 3).
test_that("uncertainty: asymptotic on a fit without cov_matrix is refused under this function's name", {
  model <- .write_fit_adaptive_model(c(2, 10, 0.09, 0.1))
  data <- .fit_adaptive_data()
  fit <- ferx_fit(model, data, verbose = FALSE)
  expect_length(fit$cov_matrix, 0L)
  err <- tryCatch(
    ferx_simulate_adaptive(model, data, fit = fit,
                           uncertainty = list(n_draws = 2L)),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err),
               "^ferx_simulate_adaptive: `fit\\$cov_matrix` is empty")
  expect_null(conditionCall(err))
})

test_that("uncertainty without fit, or with etas, is refused", {
  fx <- .unc_fixture()
  expect_error(
    ferx_simulate_adaptive(fx$model, fx$data, uncertainty = list(n_draws = 2L)),
    "`uncertainty` needs `fit`", fixed = TRUE
  )
  # The engine's coded refusal, spelled with this function's arguments.
  probe <- engine_error_probe(
    .run_unc(etas = fx$fit$ebe_etas, uncertainty = list(n_draws = 2L))
  )
  expect_coded_refusal(probe, "cannot be combined", "E_SUPPLIED_ETAS")
  msg <- conditionMessage(probe$cond)
  expect_true(grepl("(`uncertainty`)", msg, fixed = TRUE), info = msg)
  expect_true(grepl("(`etas`)", msg, fixed = TRUE), info = msg)
  expect_false(grepl("with_uncertainty|with_etas", msg), info = msg)
})

test_that("a malformed `uncertainty` is refused in R before the engine", {
  expect_error(.run_unc(uncertainty = list(method = "sir")),
               "`uncertainty$n_draws` is required", fixed = TRUE)
  expect_error(.run_unc(uncertainty = list(n_draws = 0L)), "positive whole")
  expect_error(.run_unc(uncertainty = list(n_draws = 1.5)), "positive whole")
  expect_error(.run_unc(uncertainty = list()), "named list")
  expect_error(.run_unc(uncertainty = list(method = "boot", n_draws = 2L)),
               "\"asymptotic\" or \"sir\"", fixed = TRUE)
  expect_error(.run_unc(uncertainty = list(methd = "sir", n_draws = 2L)),
               "unknown `uncertainty` element(s) `methd`", fixed = TRUE)
})

# Without uncertainty nothing changes: every row is DRAW 1 and `parameters`
# keeps its columns with no rows (the engine's list is empty).
test_that("without uncertainty: DRAW is 1 throughout and parameters has no rows", {
  res <- .run_unc()
  for (tbl in .adaptive_tables) {
    expect_true(all(res[[tbl]]$DRAW == 1L), label = tbl)
  }
  expect_identical(nrow(res$parameters), 0L)
  expect_identical(names(res$parameters),
                   c("DRAW", rownames(.unc_fixture()$fit$cov_matrix)))
  expect_identical(res[.adaptive_tables],
                   .run_unc(uncertainty = NULL)[.adaptive_tables])
})

# --- uncertainty = on every packed block (#763 review, "not checked") ---------
#
# The tests above draw theta, one Omega and one sigma. These run the other
# blocks `packed_parameter_entries` writes - kappa, `[mixture]` overrides,
# `block_sigma` correlations - and a theta-level block, end to end: each draw's
# column is named as `fit$cov_matrix` names the row, scatters around the
# fit's estimate by about its standard error, and a FIX or held-at-bound
# entry stays at its estimate. Every run here also passes the engine's check
# that the run's parameters are the fit's point estimate (E_ADAPTIVE_OPTIONS,
# "the run's parameters differ from it"), which the glue cannot fail: it builds
# the run's parameters and the fit it hands the engine from the same R values.
# So that refusal never reaches R, and needs no R spelling. The check compares
# theta, Omega, Omega_IOV, sigma and rho only - never the `[mixture]` class
# matrices - so for the override values the mixture test's own column check is
# what holds them to the fit (#764 review, finding 5). Mutation that reddens
# every test in this block: `build_fit_result_for_uncertainty` handed a theta
# one ULP off the fit's.

# An ODE 1-cpt adaptive model: `params` the [parameters] lines, `indiv` the
# [individual_parameters] lines, `err` the error model, `opts` extra
# [fit_options] lines.
.unc_shape_model <- function(params, indiv, err = "  DV ~ proportional(PROP)",
                             opts = character()) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]", params,
    "[individual_parameters]", indiv,
    "[structural_model]", "  ode(states=[central])",
    "[odes]", "  d/dt(central) = -(CL / V) * central",
    "[scaling]", "  y = central / V",
    "[error_model]", err,
    "[fit_options]", "  covariance = true", opts,
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

# 24 subjects, two occasions 48 h apart (OCC), each a 100 mg dose observed at
# 1, 4, 12 and 24 h; subjects 1-12 are STUDY 1, 13-24 STUDY 2. The DVs are
# simulated from `truth` (a model file) with a fixed seed.
.unc_shape_data <- function(truth) {
  design <- do.call(rbind, lapply(1:24, function(i) {
    do.call(rbind, lapply(1:2, function(o) {
      t0 <- (o - 1) * 48
      rbind(
        data.frame(ID = i, TIME = t0, DV = ".", EVID = 1, AMT = "100", CMT = 1,
                   MDV = 1, OCC = o, STUDY = if (i <= 12) 1 else 2),
        data.frame(ID = i, TIME = t0 + c(1, 4, 12, 24), DV = "1", EVID = 0,
                   AMT = ".", CMT = 1, MDV = 0, OCC = o,
                   STUDY = if (i <= 12) 1 else 2)
      )
    }))
  }))
  path <- tempfile(fileext = ".csv")
  utils::write.csv(design, path, row.names = FALSE, quote = FALSE)
  sim <- suppressWarnings(ferx_simulate(truth, path, n_sim = 1L, seed = 946L))
  obs <- design$EVID == 0
  key <- paste(design$ID, design$TIME)[obs]
  design$DV[obs] <- sprintf("%.6g", sim$DV_SIM[match(key, paste(sim$ID, sim$TIME))])
  utils::write.csv(design, path, row.names = FALSE, quote = FALSE)
  path
}

# Fit `model` on `data` and run 200 parameter draws through the controller.
.unc_shape_run <- function(model, data, method = "asymptotic", fit = NULL) {
  if (is.null(fit)) fit <- ferx_fit(model, data, verbose = FALSE)
  res <- .quiet_design_dv(ferx_simulate_adaptive(
    model, data, n_sim = 1L, seed = 3L, fit = fit,
    uncertainty = list(method = method, n_draws = 200L)
  ))
  list(fit = fit, p = res$parameters, res = res)
}

# Each `free` column of the draws `p` is centred on `est` and spread by `se`
# (natural scale; 200 draws), each `held` column is exactly its estimate.
.expect_draws_around <- function(p, est, se, free, held = character()) {
  for (nm in free) {
    x <- p[[nm]]
    expect_identical(length(unique(x)), nrow(p), label = nm)
    expect_lt(abs(stats::median(x) - est[[nm]]) / se[[nm]], 0.35)
    expect_gt(stats::sd(x) / se[[nm]], 0.7)
    expect_lt(stats::sd(x) / se[[nm]], 1.45)
  }
  for (nm in held) {
    expect_identical(p[[nm]], rep(est[[nm]], nrow(p)), label = nm)
  }
}

# Mutation that reddens this: `packed_parameter_entries` reading a kappa's
# value from the model template (`p.omega_iov`) instead of the draw, which
# puts the file's 0.05 in every row.
test_that("uncertainty: kappa columns carry each draw's IOV, a FIX kappa its estimate", {
  params <- c("  theta TVCL(2.5, 0.01, 50.0)", "  theta TVV(15, 1.0, 500.0)",
              "  omega ETA_CL ~ 0.09", "  kappa KAPPA_CL ~ 0.05",
              "  kappa KAPPA_V ~ 0.01 FIX", "  sigma PROP ~ 0.1 (sd)")
  indiv <- c("  CL = TVCL * exp(ETA_CL + KAPPA_CL)", "  V  = TVV * exp(KAPPA_V)")
  model <- .unc_shape_model(params, indiv, opts = "  iov_column = OCC")
  r <- .unc_shape_run(model, .unc_shape_data(model))
  fit <- r$fit
  expect_true(fit$converged)
  expect_identical(nrow(fit$held_at_bound), 0L)
  expect_identical(names(r$p), c("DRAW", rownames(fit$cov_matrix)))
  expect_identical(tail(names(r$p), 2L), c("KAPPA_CL", "KAPPA_V"))
  .expect_draws_around(
    r$p,
    est = c(diag(fit$omega_iov), fit$theta),
    se = c(fit$se_kappa, fit$se_theta),
    free = c("KAPPA_CL", "TVCL"), held = "KAPPA_V"
  )
})

# Mutation that reddens this: `packed_parameter_entries` reading an override
# from the base class (`m.omega.get(0)`) instead of its own class, which puts
# the base ETA_CL in the override column.
test_that("uncertainty: a [mixture] override column carries each draw's class value", {
  data <- test_path("fixtures", "mixture_iv.csv")
  model <- .write_mixture_adaptive_model()
  fit <- muffle_mixl_lint(ferx_fit(model, data, verbose = FALSE,
                                   covariance = TRUE))
  r <- muffle_mixl_lint(.unc_shape_run(model, data, fit = fit))
  expect_true(fit$converged)
  expect_identical(names(r$p), c("DRAW", rownames(fit$cov_matrix)))
  expect_identical(tail(names(r$p), 1L), "ETA_CL (class 2)")
  ov <- "ETA_CL (class 2)"
  .expect_draws_around(
    r$p,
    est = c(stats::setNames(fit$mixture_overrides, ov), fit$theta),
    se = c(stats::setNames(fit$se_mixture_overrides, ov), fit$se_theta),
    free = c(ov, "MIXL")
  )
  # Not the base class's draw.
  expect_false(any(r$p[[ov]] == r$p$ETA_CL))
})

# A correlation is drawn on the Fisher-z scale (`atanh(rho)`, the scale of its
# `fit$cov_matrix` row), so it is checked there. The data come from rho = -0.6,
# away from the engine's |z| <= 3 rail: a rho the data push onto it (as a
# PROP/ADD pair with a weak additive part does) has a packed SE of ~10 and
# draws piled on the rail. Mutation that reddens this:
# `packed_parameter_entries` reading a correlation from the model template
# instead of the draw, which puts the file's rho in every row.
test_that("uncertainty: a block_sigma correlation column carries each draw's rho", {
  params <- c("  theta TVCL(2.5, 0.01, 50.0)", "  theta TVV(15, 1.0, 500.0)",
              "  omega ETA_CL ~ 0.09",
              "  block_sigma (PROP, ADD) = [0.04, -0.012, 0.01]")
  indiv <- c("  CL = TVCL * exp(ETA_CL)", "  V  = TVV")
  model <- .unc_shape_model(params, indiv, err = "  DV ~ combined(PROP, ADD)")
  r <- .unc_shape_run(model, .unc_shape_data(model))
  fit <- r$fit
  rc <- fit$residual_correlations
  expect_true(fit$converged)
  expect_identical(names(r$p), c("DRAW", rownames(fit$cov_matrix)))
  expect_identical(tail(names(r$p), 1L), rc$name)
  rho <- r$p[[rc$name]]
  expect_true(all(rho > -1 & rho < 1))
  z <- stats::setNames(data.frame(atanh(rho)), rc$name)
  .expect_draws_around(
    z, est = stats::setNames(atanh(rc$rho), rc$name),
    se = stats::setNames(sqrt(fit$cov_matrix[rc$name, rc$name]), rc$name),
    free = rc$name
  )
})

# A level block's thetas are named per level and drawn per level, on the
# layout bound to the fit. Mutation that reddens this: `packed_parameter_
# entries` reading theta `k` from another position.
test_that("uncertainty: a theta-level block draws each level around its own estimate", {
  params <- c("  theta TVCL[STUDY, contrast = none](2.5, 0.01, 50.0)",
              "  theta TVV(15, 1.0, 500.0)", "  omega ETA_CL ~ 0.09",
              "  sigma PROP ~ 0.1 (sd)")
  indiv <- c("  CL = TVCL * exp(ETA_CL)", "  V  = TVV")
  truth <- .unc_shape_model(sub("TVCL[STUDY, contrast = none]", "TVCL",
                                params, fixed = TRUE), indiv)
  data <- .unc_shape_data(truth)
  # Study 2 clears faster, so the two levels do not share an estimate.
  d <- utils::read.csv(data, colClasses = "character")
  d$DV[d$STUDY == "2" & d$EVID == "0"] <- sprintf(
    "%.6g", as.numeric(d$DV[d$STUDY == "2" & d$EVID == "0"]) *
      exp(-0.8 * as.numeric(d$TIME[d$STUDY == "2" & d$EVID == "0"]) %% 48 / 12))
  utils::write.csv(d, data, row.names = FALSE, quote = FALSE)
  r <- .unc_shape_run(.unc_shape_model(params, indiv), data)
  fit <- r$fit
  lv <- c("TVCL[STUDY=1]", "TVCL[STUDY=2]")
  expect_identical(names(r$p), c("DRAW", rownames(fit$cov_matrix)))
  expect_identical(names(r$p)[2:3], lv)
  expect_gt(abs(diff(fit$theta[lv])) / max(fit$se_theta[lv]), 3)
  .expect_draws_around(r$p, est = fit$theta, se = fit$se_theta,
                       free = c(lv, "TVV"))
})

# `method = "sir"` resamples the fit's kept pool: every draw is one of its
# rows, unpacked (theta and the sd-scale sigma on log, the Omega variance on
# log-sd). Mutation that reddens this: the glue reading the flat pool
# column-major, which mixes coordinates of different resamples.
test_that("uncertainty: method = \"sir\" draws are rows of the fit's kept pool", {
  fx <- .unc_fixture()
  sir_fit <- ferx_sir(fx$fit, sir_samples = 200L, sir_resamples = 30L,
                      sir_keep_samples = TRUE, sir_seed = 1L)
  pool <- matrix(sir_fit$sir_resamples, nrow = sir_fit$sir_resamples_n,
                 byrow = TRUE)
  expect_identical(ncol(pool), 4L)
  # TVCL, TVV (FIX), ETA_CL, PROP.
  nat <- cbind(exp(pool[, 1:2]), exp(2 * pool[, 3]), exp(pool[, 4]))
  res <- .quiet_design_dv(ferx_simulate_adaptive(
    fx$model, fx$data, n_sim = 1L, seed = 7L, fit = sir_fit,
    uncertainty = list(method = "sir", n_draws = 8L)
  ))
  p <- as.matrix(res$parameters[, c("TVCL", "TVV", "ETA_CL", "PROP")])
  for (r in seq_len(nrow(p))) {
    rel <- apply(abs(sweep(nat, 2L, p[r, ])) / rep(p[r, ], each = nrow(nat)),
                 1L, max)
    expect_lt(min(rel), 1e-12)
  }
  # More than one row of the pool is used.
  expect_gt(length(unique(p[, "TVCL"])), 1L)
})

# The population summary under parameter uncertainty (ferx-core #946 PR2): a
# band across draws, read from the per-draw rows. With 3 draws only the median
# is resolvable (the 5% / 95% quantiles need 19 values). Mutation that reddens
# this: the glue dropping the band's quantiles, or mislabelling its status.
test_that("uncertainty: summary bands are the quantiles across draws of the per-draw rows", {
  res <- .run_unc(uncertainty = list(n_draws = 3L))
  rows <- res$summary_rows
  expect_identical(sort(unique(rows$DRAW)), 1:3)
  b <- res$summary_bands
  expect_true(all(b$STATUS == "quantiles"))
  cum <- b[b$METRIC == "CUM_DOSE" & b$STATISTIC == "mean", ]
  per_draw <- rows$VALUE[rows$METRIC == "CUM_DOSE" & rows$STATISTIC == "mean"]
  expect_length(per_draw, 3L)
  expect_identical(cum$N_DRAWS_PRESENT, 3L)
  expect_equal(cum$MED, stats::median(per_draw))
  expect_true(is.na(cum$LO) && is.na(cum$HI))
  # Each draw's mean is its own runs' mean.
  m <- res$metrics
  expect_equal(per_draw, as.vector(tapply(m$CUM_DOSE, m$DRAW, mean)))
})

# n_reps reaches the from-fit / uncertainty glue path too, and with two draws
# the decomposition estimates a parameter layer. Mutation that reddens this:
# the from-fit entry point not passing `run_options` on.
test_that("uncertainty: n_reps and decompose on the from-fit path estimate a parameter layer", {
  res <- .run_unc(uncertainty = list(n_draws = 2L), n_reps = 2L, decompose = TRUE)
  for (tbl in c("doses", "decisions", "metrics")) {
    expect_identical(sort(unique(res[[tbl]]$REP)), 1:2, label = tbl)
  }
  vc <- res$variance_components
  expect_identical(attr(vc, "n_draws"), 2L)
  expect_identical(attr(vc, "n_reps"), 2L)
  par <- vc[vc$LAYER == "parameter" & vc$METRIC == "CUM_DOSE", ]
  expect_true(all(par$STATUS == "estimated"))
  # One draw under uncertainty cannot estimate it: refused before the run (a
  # million replicates would take minutes), in R's spelling, not core's
  # after-the-run message (#779 review, finding 5).
  t0 <- proc.time()[["elapsed"]]
  e <- tryCatch(.run_unc(uncertainty = list(n_draws = 1L), n_reps = 1000000L,
                         decompose = TRUE), error = function(e) e)
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
  expect_identical(e$code, "E_ADAPTIVE_OPTIONS")
  expect_match(conditionMessage(e),
               "`decompose = TRUE` with `uncertainty` needs `uncertainty$n_draws` >= 2 to estimate the parameter layer; got 1",
               fixed = TRUE)
  expect_false(grepl("decompose_adaptive", conditionMessage(e), fixed = TRUE))
})
