# Tests for the controller-effect contrast of ferx_simulate_adaptive()
# (`reference = "fixed"`, ferx-core #946 PR4) and the summary rows' Monte Carlo
# SE (MC_SE, ferx-core #3215).

.tdm <- function(...) {
  ex <- ferx_example("adaptive_tdm")
  ferx_simulate_adaptive(ex$model, ex$data, seed = 1L, ...)
}

.run_tabs <- c("trajectories", "doses", "decisions", "metrics", "parameters")
.outcome_tabs <- c("summary_rows", "summary_bands", "prediction_intervals",
                   "variance_components", "variance_skipped")

# One arm's rows of a stacked table, without ARM: comparable to a single run.
.arm <- function(df, arm) {
  out <- df[df$ARM == arm, setdiff(names(df), "ARM"), drop = FALSE]
  rownames(out) <- NULL
  out
}

# A 1-cpt IV model whose one rule can never fire (the signal is never
# negative): its controller re-issues the start dose at every due point, which
# is exactly the fixed reference. So the adaptive arm IS the reference arm and
# every difference must be zero: the contrast of a regimen with itself.
.write_inert_model <- function() {
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
    "  when signal < -1 : increase 25%"
  ), path)
  path
}

test_that("summary_rows carries MC_SE: NA on median rows, set on mean rows", {
  res <- .tdm(n_sim = 5L)
  rows <- res$summary_rows
  expect_true(all(is.na(rows$MC_SE[rows$STATISTIC == "median"])))
  mean_cd <- rows[rows$STATISTIC == "mean" & rows$METRIC == "CUM_DOSE", ]
  expect_false(anyNA(mean_cd$MC_SE))
  expect_true(all(mean_cd$MC_SE >= 0))
  # It is the jackknife over sims: n_sim = 5 sims of one subject, so the SE of
  # the mean is sd(per-sim value) / sqrt(5).
  m <- res$metrics
  per_sim <- tapply(m$CUM_DOSE, m$SIM, mean)
  expect_equal(mean_cd$MC_SE, stats::sd(per_sim) / sqrt(length(per_sim)),
               tolerance = 1e-10)
})

test_that("MC_SE is NA throughout with a single sim", {
  rows <- .tdm(n_sim = 1L)$summary_rows
  expect_true(all(is.na(rows$MC_SE)))
  expect_false(all(is.na(rows$VALUE)))
})

test_that("reference = \"fixed\" returns a single run's tables, stacked by ARM", {
  con <- .tdm(n_sim = 3L, reference = "fixed")
  one <- .tdm(n_sim = 3L)
  expect_named(con, c(names(one), "arm_labels"))
  expect_identical(con$arm_labels,
                   c(adaptive = "[adaptive_dosing] rules",
                     reference = "start dose every due point"))
  for (tab in c(.run_tabs, .outcome_tabs[1:3])) {
    expect_identical(names(con[[tab]]), c("ARM", names(one[[tab]])), label = tab)
  }
  # Run tables: the two arms and nothing else; outcome tables add the
  # difference, and no stray arm.
  for (tab in c("doses", "decisions", "metrics", "trajectories")) {
    expect_setequal(unique(con[[tab]]$ARM), c("adaptive", "reference"))
    expect_false("difference" %in% con[[tab]]$ARM, label = tab)
  }
  for (tab in .outcome_tabs[1:3]) {
    expect_setequal(unique(con[[tab]]$ARM),
                    c("adaptive", "reference", "difference"))
  }
  expect_null(con$variance_components)
  # Attributes survive the stacking.
  expect_identical(attr(con$summary_bands, "probs"), c(0.05, 0.5, 0.95))
  # The adaptive arm is exactly the run without a reference.
  for (tab in c(.run_tabs, .outcome_tabs[1:3])) {
    expect_equal(.arm(con[[tab]], "adaptive"), one[[tab]],
                 ignore_attr = TRUE, label = tab)
  }
})

test_that("the controller changes the reference arm's doses, and the difference is paired", {
  con <- .tdm(n_sim = 5L, reference = "fixed")
  ref_doses <- .arm(con$doses, "reference")
  # The fixed reference gives the start dose at every due point.
  expect_identical(length(unique(ref_doses$AMT)), 1L)
  expect_gt(length(unique(.arm(con$doses, "adaptive")$AMT)), 1L)
  # With no metric missing, the mean of the per-run differences is the
  # difference of the arms' means.
  rows <- con$summary_rows
  key <- rows$STATISTIC == "mean" & rows$METRIC == "CUM_DOSE"
  v <- stats::setNames(rows$VALUE[key], rows$ARM[key])
  expect_equal(v[["difference"]], v[["adaptive"]] - v[["reference"]],
               tolerance = 1e-10)
  expect_false(is.na(rows$MC_SE[key & rows$ARM == "difference"]))
})

test_that("a regimen contrasted with itself has zero difference, exactly", {
  ex <- ferx_example("adaptive_tdm")
  con <- ferx_simulate_adaptive(.write_inert_model(), ex$data, n_sim = 4L,
                                seed = 2L, n_reps = 2L, decompose = TRUE,
                                reference = "fixed")
  for (tab in c(.run_tabs, .outcome_tabs)) {
    if (is.null(con[[tab]])) next
    expect_identical(.arm(con[[tab]], "adaptive"), .arm(con[[tab]], "reference"),
                     label = tab)
  }
  rows <- .arm(con$summary_rows, "difference")
  mean_rows <- rows[rows$STATISTIC == "mean" & !is.na(rows$VALUE), ]
  expect_gt(nrow(mean_rows), 0L)
  expect_identical(mean_rows$VALUE, rep(0, nrow(mean_rows)))
  vc <- .arm(con$variance_components, "difference")
  est <- vc[vc$STATUS == "estimated", ]
  expect_gt(nrow(est), 0L)
  expect_identical(est$VALUE, rep(0, nrow(est)))
  expect_setequal(unique(con$variance_components$ARM),
                  c("adaptive", "reference", "difference"))
})

test_that("arm_labels are passed through; an empty one is refused before the run", {
  labs <- c(reference = "1000 mg q12h", adaptive = "TDM rule")
  con <- .tdm(n_sim = 2L, reference = "fixed", arm_labels = labs)
  expect_identical(con$arm_labels,
                   c(adaptive = "TDM rule", reference = "1000 mg q12h"))

  e <- tryCatch(.tdm(n_sim = 2L, reference = "fixed",
                     arm_labels = c(adaptive = "TDM", reference = " ")),
                error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_ADAPTIVE_OPTIONS")
  expect_match(conditionMessage(e),
               "`arm_labels`: each arm needs a label that says what it is",
               fixed = TRUE)
  expect_false(grepl("summarize_contrast", conditionMessage(e), fixed = TRUE))
})

test_that("malformed reference / arm_labels are refused in R", {
  expect_error(.tdm(reference = "loading"), "`reference` must be NULL or \"fixed\"")
  expect_error(.tdm(reference = TRUE), "`reference` must be NULL or \"fixed\"")
  expect_error(.tdm(reference = "fixed", arm_labels = c("a", "b")),
               "`arm_labels` must be a character vector")
  expect_error(.tdm(reference = "fixed", arm_labels = c(adaptive = "a")),
               "`arm_labels` must be a character vector")
  expect_warning(.tdm(arm_labels = c(adaptive = "a", reference = "b")),
                 "`arm_labels` has no effect without `reference`")
})
