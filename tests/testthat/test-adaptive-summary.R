# Tests for the adaptive population summary, the variance decomposition and
# assay replicates of ferx_simulate_adaptive() (ferx-core #946 PR2-4).

.tdm_run <- function(...) {
  ex <- ferx_example("adaptive_tdm")
  ferx_simulate_adaptive(ex$model, ex$data, seed = 1L, ...)
}

# The run's own tables: what n_reps must leave untouched.
.run_tables <- c("trajectories", "doses", "decisions", "metrics", "parameters")

# A 1-cpt IV adaptive model titrating on the noise-free (Ipred) signal and with
# no IOV: nothing is redrawn per assay replicate, so every replicate of a
# (SIM, ID) must repeat the first exactly if its eta is held fixed.
.write_ipred_model <- function(rules = "  when signal < 10 : increase 25%") {
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
    rules
  ), path)
  path
}

test_that("the summary tables are attached, with the documented columns", {
  res <- .tdm_run(n_sim = 2L)
  expect_named(res, c(.run_tables, "summary_rows", "summary_bands",
                      "prediction_intervals", "variance_components",
                      "variance_skipped"))
  expect_identical(names(res$summary_rows),
    c("DRAW", "METRIC", "ANALYTE", "DECISION", "STATISTIC", "THRESHOLD",
      "VALUE", "N_RUNS", "N_MISSING"))
  expect_identical(names(res$summary_bands),
    c("METRIC", "ANALYTE", "DECISION", "STATISTIC", "THRESHOLD",
      "N_DRAWS_PRESENT", "N_DRAWS_MISSING", "STATUS", "LO", "MED", "HI"))
  expect_identical(names(res$prediction_intervals),
    c("METRIC", "ANALYTE", "DECISION", "LO", "MED", "HI", "N_RUNS",
      "N_MISSING"))
  expect_identical(attr(res$summary_bands, "probs"), c(0.05, 0.5, 0.95))
  expect_identical(attr(res$prediction_intervals, "probs"), c(0.05, 0.5, 0.95))
  # Without decompose = TRUE there is no decomposition.
  expect_null(res$variance_components)
  expect_null(res$variance_skipped)
  # Every scalar metric is a metrics column; SIGNAL rows name the monitor and
  # a decision index of the decisions table.
  scalar <- setdiff(unique(res$summary_rows$METRIC), "SIGNAL")
  expect_true(all(scalar %in% names(res$metrics)))
  sig <- res$summary_rows[res$summary_rows$METRIC == "SIGNAL", ]
  expect_true(nrow(sig) > 0L)
  expect_false(anyNA(sig$ANALYTE))
  expect_setequal(unique(sig$DECISION), unique(res$decisions$DECISION))
  expect_true(all(is.na(res$summary_rows$ANALYTE[res$summary_rows$METRIC != "SIGNAL"])))
})

test_that("n_reps = 1 is the default run, bit for bit", {
  a <- .tdm_run(n_sim = 2L)
  b <- .tdm_run(n_sim = 2L, n_reps = 1L)
  expect_identical(a, b)
})

test_that("n_reps = 3 tags REP 1..3; replicate 1 is the run without replicates", {
  one <- .tdm_run(n_sim = 2L)
  three <- .tdm_run(n_sim = 2L, n_reps = 3L)
  for (tab in c("doses", "decisions", "metrics")) {
    expect_identical(sort(unique(three[[tab]]$REP)), 1:3, label = tab)
  }
  m <- three$metrics
  expect_identical(nrow(m), 3L * nrow(one$metrics))
  expect_false(anyDuplicated(m[c("SIM", "ID", "REP")]) > 0L)
  # trajectories hold replicate 1 only, and it is the n_reps = 1 run's.
  expect_identical(three$trajectories, one$trajectories)
  for (tab in c("doses", "decisions", "metrics")) {
    first <- three[[tab]][three[[tab]]$REP == 1L, , drop = FALSE]
    rownames(first) <- NULL
    expect_identical(first, one[[tab]], label = tab)
  }
  # The assay noise is redrawn per replicate, so the measured signals differ.
  d <- three$decisions
  expect_false(identical(d$SIGNAL[d$REP == 1L], d$SIGNAL[d$REP == 2L]))
})

test_that("n_reps holds each subject's eta fixed across replicates", {
  ex <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(.write_ipred_model(), ex$data, n_sim = 2L,
                                seed = 3L, n_reps = 3L)
  # No assay error and no IOV: only eta could make replicates differ.
  d <- res$decisions
  key <- c("SIM", "ID", "DECISION")
  r1 <- d[d$REP == 1L, c(key, "LATENT")]
  for (r in 2:3) {
    rr <- d[d$REP == r, c(key, "LATENT")]
    expect_identical(rr$LATENT[do.call(order, rr[key])],
                     r1$LATENT[do.call(order, r1[key])], label = paste("rep", r))
  }
  # ...while the SIMs (new etas) do differ.
  expect_false(identical(d$LATENT[d$SIM == 1L], d$LATENT[d$SIM == 2L]))
})

test_that("a malformed n_reps is refused in R before the engine", {
  for (bad in list(0L, -1L, 1.5, NA, c(1L, 2L), "2")) {
    expect_error(.tdm_run(n_reps = bad), "`n_reps` must be a single whole number",
                 label = format(bad))
  }
})

test_that("the per-draw summary matches the mean / median / fraction recomputed from metrics", {
  res <- .tdm_run(n_sim = 5L, n_reps = 2L,
                  summary_thresholds = list(PCT_READINGS_IN_WINDOW = 0.5,
                                            CUM_DOSE = c(8000, 12000)))
  m <- res$metrics
  rows <- res$summary_rows
  get <- function(metric, stat, thr = NA) {
    r <- rows[rows$METRIC == metric & rows$STATISTIC == stat &
                (is.na(thr) | (!is.na(rows$THRESHOLD) & rows$THRESHOLD == thr)), ]
    expect_identical(nrow(r), 1L, label = paste(metric, stat, thr))
    r
  }
  for (col in c("CUM_DOSE", "N_DOSES", "N_INCREASES", "N_DECREASES", "N_HOLDS",
                "N_NOT_DUE", "SIGNAL_MIN", "SIGNAL_MAX", "SIGNAL_MEAN",
                "PCT_READINGS_IN_WINDOW")) {
    v <- as.numeric(m[[col]])
    present <- v[!is.na(v)]
    mean_row <- get(col, "mean")
    expect_equal(mean_row$VALUE, mean(present), tolerance = 1e-12, label = col)
    expect_identical(mean_row$N_RUNS, nrow(m))
    expect_identical(mean_row$N_MISSING, sum(is.na(v)))
    expect_equal(get(col, "median")$VALUE, stats::median(present),
                 tolerance = 1e-12, label = col)
  }
  # DISCONTINUED as 1/0: its mean is the fraction discontinued.
  expect_equal(get("DISCONTINUED", "mean")$VALUE, mean(m$DISCONTINUED))
  # No auc_target: every run is missing, the value is NA, never imputed.
  auc <- get("AUC_TARGET_ATTAINMENT", "mean")
  expect_true(is.na(auc$VALUE))
  expect_identical(auc$N_MISSING, nrow(m))
  # FractionAtLeast: >= the threshold, over the runs carrying the metric.
  p <- m$PCT_READINGS_IN_WINDOW
  expect_equal(get("PCT_READINGS_IN_WINDOW", "fraction_at_least", 0.5)$VALUE,
               mean(p[!is.na(p)] >= 0.5))
  for (a in c(8000, 12000)) {
    expect_equal(get("CUM_DOSE", "fraction_at_least", a)$VALUE,
                 mean(m$CUM_DOSE >= a), label = a)
  }
  # A SIGNAL row is the mean of that decision's LATENT over the runs.
  d <- res$decisions
  sig <- rows[rows$METRIC == "SIGNAL" & rows$STATISTIC == "mean", ]
  for (k in seq_len(nrow(sig))) {
    lat <- d$LATENT[d$DECISION == sig$DECISION[k]]
    expect_equal(sig$VALUE[k], mean(lat), tolerance = 1e-12)
  }
  # The pooled prediction interval: type-7 quantiles over every run.
  pi <- res$prediction_intervals[res$prediction_intervals$METRIC == "CUM_DOSE", ]
  expect_equal(unname(unlist(pi[c("LO", "MED", "HI")])),
               unname(stats::quantile(m$CUM_DOSE, c(0.05, 0.5, 0.95), type = 7)))
})

test_that("a SIGNAL threshold is given by data frame, with its analyte and decision", {
  base <- .tdm_run(n_sim = 3L)
  analyte <- unique(stats::na.omit(base$summary_rows$ANALYTE))
  expect_length(analyte, 1L)
  thr <- data.frame(METRIC = c("SIGNAL", "PCT_READINGS_IN_WINDOW"),
                    THRESHOLD = c(10, 0.5), ANALYTE = c(analyte, NA),
                    DECISION = c(2L, NA))
  res <- .tdm_run(n_sim = 3L, summary_thresholds = thr)
  r <- res$summary_rows
  hit <- r[r$STATISTIC == "fraction_at_least" & r$METRIC == "SIGNAL", ]
  expect_identical(nrow(hit), 1L)
  expect_identical(hit$DECISION, 2L)
  lat <- res$decisions$LATENT[res$decisions$DECISION == 2L]
  expect_equal(hit$VALUE, mean(lat >= 10))
  expect_identical(sum(r$STATISTIC == "fraction_at_least"), 2L)
})

test_that("summary options are refused: bad probs and non-finite thresholds coded, unknown metrics named", {
  # Refused in R's spelling - the argument and the metrics column, never the
  # Rust enum or array (#779 review, finding 1) - and still coded.
  coded <- function(expr, pattern) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "ferx_engine_error")
    expect_identical(e$code, "E_ADAPTIVE_OPTIONS")
    expect_match(conditionMessage(e), pattern, fixed = TRUE)
    expect_false(grepl("summarize_adaptive|CumulativeDose|Signal \\{|\\[0",
                       conditionMessage(e)))
  }
  coded(.tdm_run(summary_probs = c(0.5, 0.1, 0.9)),
        "`summary_probs` c(0.5, 0.1, 0.9) must be three probabilities")
  coded(.tdm_run(summary_probs = c(0, 0.5, 1)), "`summary_probs` c(0, 0.5, 1)")
  coded(.tdm_run(summary_thresholds = list(CUM_DOSE = Inf)),
        "the threshold for CUM_DOSE is Inf; it must be finite")
  coded(.tdm_run(summary_thresholds = list(PCT_READINGS_IN_WINDOW = NA_real_)),
        "the threshold for PCT_READINGS_IN_WINDOW is NA")
  coded(.tdm_run(summary_thresholds = data.frame(
          METRIC = "SIGNAL", THRESHOLD = -Inf, ANALYTE = "signal", DECISION = 1L)),
        "the threshold for SIGNAL (ANALYTE \"signal\", DECISION 1) is -Inf")
  # A negative decision is named as such, not as a missing one (finding 7).
  expect_error(.tdm_run(summary_thresholds = data.frame(
                 METRIC = "SIGNAL", THRESHOLD = 1, ANALYTE = "signal",
                 DECISION = -1L)),
               "`DECISION` is -1; decisions are numbered from 0")
  expect_error(.tdm_run(summary_thresholds = list(FOO = 1)),
               "`FOO` is not a summary metric")
  expect_error(.tdm_run(summary_thresholds = data.frame(METRIC = "SIGNAL",
                                                        THRESHOLD = 1)),
               "needs an `ANALYTE`")
  expect_error(.tdm_run(summary_probs = c(0.1, 0.9)), "three probabilities")
  expect_error(.tdm_run(summary_thresholds = list(1)), "must be NULL, a named list")
  expect_error(.tdm_run(summary_thresholds = list(CUM_DOSE = "1")),
               "must be NULL, a named list")
  expect_error(.tdm_run(decompose = NA), "`decompose` must be a single")
})

test_that("a run with one draw has no band: STATUS not_drawn, never a zero-width band", {
  res <- .tdm_run(n_sim = 5L)
  b <- res$summary_bands
  expect_true(all(b$STATUS == "not_drawn"))
  expect_true(all(is.na(b$LO) & is.na(b$MED) & is.na(b$HI)))
  expect_true(all(b$N_DRAWS_PRESENT + b$N_DRAWS_MISSING == 1L))
  expect_true(all(res$summary_rows$DRAW == 1L))
})

test_that("a read.csv threshold table (\"\" for no analyte) is accepted", {
  thr <- data.frame(METRIC = c("CUM_DOSE", "SIGNAL"), THRESHOLD = c(9000, 10),
                    ANALYTE = c("", "signal"), DECISION = c(NA, 3L))
  res <- expect_silent(.tdm_run(n_sim = 3L, summary_thresholds = thr))
  fr <- res$summary_rows[res$summary_rows$STATISTIC == "fraction_at_least", ]
  expect_setequal(fr$METRIC, c("CUM_DOSE", "SIGNAL"))
  expect_equal(fr$VALUE[fr$METRIC == "CUM_DOSE"], mean(res$metrics$CUM_DOSE >= 9000))
})

test_that("a SIGNAL threshold no run recorded warns, naming what was recorded", {
  # A mistyped analyte, or a decision past the schedule: NA rows, and a
  # warning rather than silence (#779 review, finding 2).
  for (thr in list(
    data.frame(METRIC = "SIGNAL", THRESHOLD = 10, ANALYTE = "signal_typo", DECISION = 1L),
    data.frame(METRIC = "SIGNAL", THRESHOLD = 10, ANALYTE = "signal", DECISION = 999L)
  )) {
    expect_warning(
      res <- .tdm_run(n_sim = 3L, summary_thresholds = thr),
      paste0("no run recorded a signal for the `summary_thresholds` row\\(s\\) ",
             "ANALYTE \"", thr$ANALYTE, "\" at DECISION ", thr$DECISION,
             ".*recorded ANALYTE \"signal\" at DECISION 0-8")
    )
    r <- res$summary_rows
    hit <- r[r$METRIC == "SIGNAL" & r$ANALYTE == thr$ANALYTE &
               r$DECISION == thr$DECISION, ]
    expect_true(all(is.na(hit$VALUE) & hit$N_MISSING == hit$N_RUNS))
  }
  # A recorded one does not warn.
  expect_silent(.tdm_run(n_sim = 3L, summary_thresholds = data.frame(
    METRIC = "SIGNAL", THRESHOLD = 10, ANALYTE = "signal", DECISION = 8L)))
})

test_that("custom summary_probs reach the prediction intervals and the bands", {
  res <- .tdm_run(n_sim = 5L, summary_probs = c(0.25, 0.5, 0.75))
  expect_identical(attr(res$prediction_intervals, "probs"), c(0.25, 0.5, 0.75))
  expect_identical(attr(res$summary_bands, "probs"), c(0.25, 0.5, 0.75))
  pi <- res$prediction_intervals[res$prediction_intervals$METRIC == "CUM_DOSE", ]
  expect_equal(unname(unlist(pi[c("LO", "MED", "HI")])),
               unname(stats::quantile(res$metrics$CUM_DOSE, c(0.25, 0.5, 0.75))))
})

test_that("decompose = TRUE: the layers sum to the total, per subject and for the population", {
  res <- .tdm_run(n_sim = 4L, n_reps = 3L, decompose = TRUE)
  vc <- res$variance_components
  expect_identical(names(vc),
    c("METRIC", "ANALYTE", "DECISION", "SCOPE", "ID", "LAYER", "STATUS",
      "VALUE", "MC_SE", "N_SUBJECTS"))
  expect_identical(names(res$variance_skipped),
                   c("METRIC", "ANALYTE", "DECISION", "ID", "N_MISSING"))
  expect_identical(attr(vc, "n_draws"), 1L)
  expect_identical(attr(vc, "n_sims"), 4L)
  expect_identical(attr(vc, "n_reps"), 3L)
  # No uncertainty: the parameter layer is not drawn, never a zero.
  par <- vc[vc$LAYER == "parameter", ]
  expect_true(all(par$STATUS == "not_drawn" & is.na(par$VALUE)))
  expect_setequal(unique(vc$SCOPE), c("subject", "population"))
  expect_true(all(is.na(vc$ID[vc$SCOPE == "population"])))
  expect_setequal(unique(vc$ID[vc$SCOPE == "subject"]), unique(res$metrics$ID))
  # A subject's total is the sum of its estimated layers; the population total
  # is the sum of its lines (no subject's eta was supplied).
  groups <- split(vc, list(vc$METRIC, vc$ANALYTE, vc$DECISION, vc$SCOPE, vc$ID),
                  drop = TRUE)
  pop <- vc[vc$SCOPE == "population", ]
  pop_groups <- split(pop, list(pop$METRIC, pop$ANALYTE, pop$DECISION),
                      drop = TRUE)
  pop_scalar <- split(pop[pop$METRIC != "SIGNAL", ],
                      pop$METRIC[pop$METRIC != "SIGNAL"])
  checked <- 0L
  for (g in c(groups, pop_groups, pop_scalar)) {
    est <- g[g$STATUS == "estimated", ]
    if (!"total" %in% est$LAYER) next
    parts <- sum(est$VALUE[est$LAYER != "total"])
    expect_equal(est$VALUE[est$LAYER == "total"], parts, tolerance = 1e-9)
    checked <- checked + 1L
  }
  expect_true(checked > 0L)
  # The population lines: parameter, bsv, within, between_subject, total.
  cum <- pop[pop$METRIC == "CUM_DOSE", ]
  expect_identical(cum$LAYER,
                   c("parameter", "bsv", "within", "between_subject", "total"))
})

# A design decompose_adaptive() cannot separate is refused before the run
# (#779 review, finding 5), so each call below is one that would take minutes
# if it ran (a million runs or more): it must come back at once, coded, in R's
# spelling and not core's after-the-run message (`decompose_adaptive: ...`).
# Mutation that reddens this: dropping a condition from the glue's
# check_decomposition_design(), which leaves the refusal to core after the run.
.expect_refused_up_front <- function(expr, pattern) {
  t0 <- proc.time()[["elapsed"]]
  e <- tryCatch(expr, error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - t0
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_ADAPTIVE_OPTIONS")
  expect_match(conditionMessage(e), pattern, fixed = TRUE)
  expect_false(grepl("decompose_adaptive", conditionMessage(e), fixed = TRUE))
  expect_lt(elapsed, 5)
}

test_that("decompose = TRUE refuses an unseparable design before the run, in R's spelling", {
  .expect_refused_up_front(
    .tdm_run(n_sim = 1000000L, decompose = TRUE),
    "`decompose = TRUE` needs `n_reps` >= 2 to separate the within-subject layer from BSV; got 1"
  )
  .expect_refused_up_front(
    .tdm_run(n_sim = 1L, n_reps = 1000000L, decompose = TRUE),
    "`decompose = TRUE` needs `n_sim` >= 2 to separate BSV from the parameter layer; got 1"
  )
  # Core checks n_reps first; so does the glue.
  .expect_refused_up_front(.tdm_run(n_sim = 1L, decompose = TRUE),
                           "needs `n_reps` >= 2")
})

test_that("decompose = TRUE accepts the smallest separable design, n_sim = 2 and n_reps = 2", {
  res <- .tdm_run(n_sim = 2L, n_reps = 2L, decompose = TRUE)
  expect_identical(attr(res$variance_components, "n_sims"), 2L)
  expect_identical(attr(res$variance_components, "n_reps"), 2L)
  expect_gt(nrow(res$variance_components), 0L)
})

test_that("a supplied-eta subject's BSV is `supplied`", {
  res <- .tdm_run(n_sim = 3L, n_reps = 2L, decompose = TRUE,
                  etas = data.frame(ID = "1", ETA_CL = 0.2), unmatched = "draw")
  vc <- res$variance_components
  bsv <- vc[vc$LAYER == "bsv" & vc$SCOPE == "subject", ]
  expect_true(all(bsv$STATUS[bsv$ID == "1"] == "supplied"))
  expect_true(all(bsv$STATUS[bsv$ID != "1"] == "estimated"))
})

# Runs that stop early record no signal at the later decisions, so those
# SIGNAL metrics are partly missing: the summary is over the present runs and
# counts the rest, and the decomposition skips the (metric, subject) pairs
# with a gap. Mutations that redden this: MC_SE always NA, the skipped table's
# N_MISSING always 0, N_MISSING not reaching R (#779 review, finding 3).
test_that("partly missing metrics: summary over present runs, skipped pairs counted, MC_SE finite", {
  ex <- ferx_example("adaptive_tdm")
  model <- .write_ipred_model(c("  when signal > 11 : stop",
                                "  when signal < 10 : increase 25%"))
  res <- ferx_simulate_adaptive(model, ex$data, n_sim = 4L, seed = 3L,
                                n_reps = 2L, decompose = TRUE)
  m <- res$metrics
  expect_true(any(m$DISCONTINUED) && nrow(m) == 5L * 4L * 2L)
  d <- res$decisions
  rows <- res$summary_rows
  sig <- rows[rows$METRIC == "SIGNAL" & rows$STATISTIC == "mean", ]
  partial <- sig[sig$N_MISSING > 0L & sig$N_MISSING < sig$N_RUNS, ]
  expect_gt(nrow(partial), 0L)
  for (k in seq_len(nrow(sig))) {
    lat <- d$LATENT[d$DECISION == sig$DECISION[k]]
    expect_identical(sig$N_RUNS[k], nrow(m))
    expect_identical(sig$N_MISSING[k], nrow(m) - length(lat))
    expect_equal(sig$VALUE[k], mean(lat), tolerance = 1e-12)
  }
  # The decomposition skips a subject's SIGNAL at a decision some of its runs
  # never reached, counting the runs without it.
  sk <- res$variance_skipped
  expect_gt(nrow(sk), 0L)
  expect_true(all(sk$METRIC == "SIGNAL"))
  for (k in seq_len(nrow(sk))) {
    present <- sum(d$ID == sk$ID[k] & d$DECISION == sk$DECISION[k])
    expect_identical(sk$N_MISSING[k], 4L * 2L - present)
  }
  vc <- res$variance_components
  # A skipped pair has no subject rows for that metric.
  for (k in seq_len(nrow(sk))) {
    expect_false(any(vc$SCOPE == "subject" & vc$ID == sk$ID[k] &
                       vc$METRIC == "SIGNAL" & vc$DECISION %in% sk$DECISION[k]))
  }
  # A jackknife MC SE for every estimated BSV line of CUM_DOSE. The model has
  # no assay error and no IOV, so the within layer is exactly 0.
  cum <- vc[vc$METRIC == "CUM_DOSE", ]
  bsv <- cum[cum$LAYER == "bsv", ]
  expect_true(all(bsv$STATUS == "estimated"))
  expect_true(all(is.finite(bsv$MC_SE) & bsv$MC_SE > 0))
  expect_true(all(cum$VALUE[cum$LAYER == "within"] == 0))
})
