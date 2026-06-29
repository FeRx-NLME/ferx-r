# Tests for ferx_simulate_adaptive() - the declarative [adaptive_dosing] /
# feedback-dosing simulation path (ferx-core #585, epic #391).

test_that("ferx_simulate_adaptive returns three named data frames", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  expect_type(res, "list")
  expect_named(res, c("trajectories", "doses", "decisions"))
  expect_s3_class(res$trajectories, "data.frame")
  expect_s3_class(res$doses, "data.frame")
  expect_s3_class(res$decisions, "data.frame")
})

test_that("the result tables carry the documented columns", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_true(all(c("DRAW", "SIM", "ID", "TIME", "IPRED", "DV_SIM") %in%
    names(res$trajectories)))
  expect_true(all(c("SIM", "ID", "TIME", "AMT", "CMT", "RATE", "DECISION",
    "SIGNAL", "RULE") %in% names(res$doses)))
  expect_true(all(c("SIM", "ID", "DECISION", "TIME", "SIGNAL", "OUTCOME",
    "N_DOSED") %in% names(res$decisions)))
})

test_that("n_sim = 3 produces three replicates", {
  ex  <- ferx_example("adaptive_tdm")
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 3L, seed = 1L)
  expect_equal(length(unique(res$trajectories$SIM)), 3L)
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
  # The engine prints an error and returns NULL when the model has no
  # [adaptive_dosing] block (this entry point requires one).
  res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 1L, seed = 1L)
  expect_null(res)
})
