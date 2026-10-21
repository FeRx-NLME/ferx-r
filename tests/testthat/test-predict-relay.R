# The live side of the findings `ferx_predict()` promises to relay (ferx-r
# #434, a #431 follow-up). Core covers each emitter; these tests catch the glue
# dropping one, or core's reader filter withholding one everywhere.
#
# The third promise from #431, skipped flip-flop draws on
# `ferx_simulate_with_uncertainty()`, is pinned in
# test-ferx_simulate_with_uncertainty.R ("skipped flip-flop draws reach R as
# one coded warning").

predict_codes <- function(x) {
  w <- attr(x, "simulation_warnings", exact = TRUE)
  expect_type(w, "character")
  sub(":.*$", "", w)
}

# A copy of the bundled Michaelis-Menten model with a solver budget too small
# for any dosing interval. The system is nonlinear, so the exact linear-ODE
# path (ferx-core #2021) cannot step around the budget.
mm_starved <- function(max_steps, env = parent.frame()) {
  lines <- readLines(ferx_example("mm_oral")$model)
  at <- grep("^\\[fit_options\\]", lines)
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(append(lines, sprintf("  ode_max_steps = %d", max_steps), after = at),
             path)
  normalizePath(path)
}

test_that("ferx_predict relays the engine's ODE-solver diagnostic", {
  data <- ferx_example("mm_oral")$data
  pred <- NULL
  w <- testthat::capture_warnings(pred <- ferx_predict(mm_starved(3L), data))
  expect_s3_class(pred, "data.frame")
  expect_true("W_ODE_SOLVER_DIAGNOSTICS" %in% predict_codes(pred))
  expect_true(any(grepl("W_ODE_SOLVER_DIAGNOSTICS", w, fixed = TRUE)))
})

test_that("an ample solver budget raises no ODE-solver diagnostic on predict", {
  # Control: the finding comes from the starved budget, not the model.
  data <- ferx_example("mm_oral")$data
  pred <- suppressWarnings(ferx_predict(mm_starved(100000L), data))
  expect_false("W_ODE_SOLVER_DIAGNOSTICS" %in% predict_codes(pred))
})

# The bundled warfarin data without its CMT column, so every row's compartment
# is defaulted.
warfarin_no_cmt <- function(env = parent.frame()) {
  data <- utils::read.csv(ferx_example("warfarin")$data, na.strings = ".")
  data$CMT <- NULL
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
  utils::write.csv(data, path, row.names = FALSE, quote = FALSE, na = ".")
  normalizePath(path)
}

test_that("ferx_predict on a compartment model still raises W_CMT_DEFAULTED", {
  # test-diag-warnings.R checks the finding is withheld where it is moot; this
  # is the other side, so a filter withholding it everywhere fails here.
  ex <- ferx_example("warfarin")
  pred <- NULL
  w <- testthat::capture_warnings(pred <- ferx_predict(ex$model, warfarin_no_cmt()))
  expect_gt(nrow(pred), 0L)
  expect_true("W_CMT_DEFAULTED" %in% predict_codes(pred))
  expect_true(any(grepl("W_CMT_DEFAULTED", w, fixed = TRUE)))
})

test_that("ferx_predict(fit = ) relays W_CMT_DEFAULTED too", {
  # The separate from-fit glue entry point (`ferx_rust_predict_from_fit`).
  ex <- ferx_example("warfarin")
  fit <- warfarin_fit()
  pred <- suppressWarnings(ferx_predict(ex$model, warfarin_no_cmt(), fit = fit))
  expect_true("W_CMT_DEFAULTED" %in% predict_codes(pred))
})
