# A user's own `data.frame` in the global environment must not reach the glue
# (#780). extendr's `data_frame!` evaluated `data.frame(..)` from there, so a
# shadowing function made each of these entry points panic and print the whole
# table. The glue now builds every table with `r_data_frame()`, which never
# calls R; `tools/check-glue-r-lookups.sh` keeps the macro out. Each test
# compares the shadowed result with the unshadowed one: identical, not merely
# "no error", since the glue's frame is what the user gets back.

# Evaluate `expr` with a `data.frame` in the global environment that refuses
# to run. Lazy: `expr` is forced only after the shadow is in place.
with_shadowed_data_frame <- function(expr) {
  stopifnot(!exists("data.frame", envir = globalenv(), inherits = FALSE))
  assign("data.frame", function(...) stop("shadowed"), envir = globalenv())
  on.exit(rm("data.frame", envir = globalenv()), add = TRUE)
  force(expr)
}

test_that("the shadow is what data_frame! would have found, and is removed after", {
  # Without this, a helper that assigned somewhere else would pass every test
  # below without testing anything.
  expect_error(
    with_shadowed_data_frame(eval(quote(data.frame(a = 1)), globalenv())),
    "shadowed"
  )
  expect_false(exists("data.frame", envir = globalenv(), inherits = FALSE))
})

test_that("ferx_predict() ignores a data.frame in the global environment (#780)", {
  ex <- ferx_example("warfarin")
  ref <- ferx_predict(ex$model, ex$data)
  expect_gt(nrow(ref), 0L)
  expect_identical(with_shadowed_data_frame(ferx_predict(ex$model, ex$data)), ref)
})

test_that("ferx_predict(fit = ) ignores a data.frame in the global environment (#780)", {
  ex <- ferx_example("warfarin")
  fit <- warfarin_fit()
  ref <- ferx_predict(ex$model, ex$data, fit = fit)
  expect_gt(nrow(ref), 0L)
  expect_identical(
    with_shadowed_data_frame(ferx_predict(ex$model, ex$data, fit = fit)),
    ref
  )
})

test_that("ferx_simulate_adaptive() ignores a data.frame in the global environment (#780)", {
  ex <- ferx_example("adaptive_tdm")
  run <- function() ferx_simulate_adaptive(ex$model, ex$data, n_sim = 2L, seed = 1L)
  ref <- run()
  for (tbl in c("doses", "decisions", "metrics")) expect_gt(nrow(ref[[tbl]]), 0L)
  expect_identical(with_shadowed_data_frame(run()), ref)
})

test_that("ferx_calc_npde() ignores a data.frame in the global environment (#780)", {
  fit <- warfarin_fit()
  run <- function() ferx_calc_npde(fit, nsim = 100L, seed = 3L)
  ref <- run()
  expect_true(all(is.finite(ref$sdtab$NPDE)))
  expect_identical(with_shadowed_data_frame(run())$sdtab, ref$sdtab)
})

pktte_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("pktte_joint")
      fit <<- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                        verbose = FALSE))
    }
    fit
  }
})

test_that("ferx_predict_survival() ignores a data.frame in the global environment (#780)", {
  ex <- ferx_example("pktte_joint")
  fit <- pktte_fit()
  times <- c(0, 6, 12, 24)
  pop <- ferx_predict_survival(ex$model, ex$data, times = times)
  at_fit <- ferx_predict_survival(ex$model, ex$data, times = times, fit = fit)
  expect_gt(nrow(pop), 0L)
  expect_gt(nrow(at_fit), 0L)
  expect_identical(
    with_shadowed_data_frame(ferx_predict_survival(ex$model, ex$data, times = times)),
    pop
  )
  expect_identical(
    with_shadowed_data_frame(
      ferx_predict_survival(ex$model, ex$data, times = times, fit = fit)
    ),
    at_fit
  )
})

# -- The layout `r_data_frame()` builds, checked against what `data_frame!`
# (that is, `data.frame()`) gave on `main` (#781 review finding 6). The
# shadow tests above compare one build with itself, so a wrong `nrow` or a
# column of the wrong type passes them; these pin the frame as R would have
# built it. The tables are taken straight from the glue, before any R code
# reshapes them, by wrapping each `ferx_rust_*` binding.

# Run `expr`, collecting every value the namespace function `fn` returns.
capture_glue <- function(fn, expr) {
  orig <- get(fn, envir = asNamespace("ferx"))
  # The wrapper appends to this environment rather than with `<<-`.
  seen <- new.env(parent = emptyenv())
  seen$values <- list()
  wrapper <- function(...) {
    out <- orig(...)
    seen$values[[length(seen$values) + 1L]] <- out
    out
  }
  do.call(local_mocked_bindings, c(setNames(list(wrapper), fn), .package = "ferx"))
  force(expr)
  seen$values
}

# A frame laid out as `data.frame()` lays it out: these columns in this order
# with these classes, all of one length, compact automatic row names counted
# from the columns themselves (not from `nrow()`, which reads the row names).
expect_data_frame_layout <- function(x, classes) {
  expect_identical(class(x), "data.frame")
  expect_identical(vapply(x, function(col) class(col)[1], ""), classes)
  n <- length(x[[1]])
  expect_true(all(lengths(x) == n))
  expect_identical(.row_names_info(x, 0L),
                   if (n == 0L) integer(0) else c(NA_integer_, -n))
}

test_that("every glue table has data.frame()'s layout (#780, #781 review)", {
  ex <- ferx_example("warfarin")
  pred_cols <- c(ID = "character", TIME = "numeric", PRED = "numeric")
  pop <- capture_glue("ferx_rust_predict", ferx_predict(ex$model, ex$data))
  expect_length(pop, 1L)
  expect_gt(length(pop[[1]][[1]]), 0L)
  expect_data_frame_layout(pop[[1]], pred_cols)
  at_fit <- capture_glue("ferx_rust_predict_from_fit",
                         ferx_predict(ex$model, ex$data, fit = warfarin_fit()))
  expect_data_frame_layout(at_fit[[1]], pred_cols)

  surv_cols <- c(ID = "character", CMT = "integer", TIME = "numeric",
                 survival = "numeric", cum_hazard = "numeric", hazard = "numeric",
                 cif = "numeric", survival_all = "numeric",
                 median_survival = "numeric", mean_survival = "numeric")
  pk <- ferx_example("pktte_joint")
  surv <- capture_glue("ferx_rust_predict_survival",
                       ferx_predict_survival(pk$model, pk$data, times = c(0, 12)))
  expect_gt(length(surv[[1]][[1]]), 0L)
  expect_data_frame_layout(surv[[1]], surv_cols)
  surv_fit <- capture_glue("ferx_rust_predict_survival_from_fit",
                           ferx_predict_survival(pk$model, pk$data, times = c(0, 12),
                                                 fit = pktte_fit()))
  expect_data_frame_layout(surv_fit[[1]], surv_cols)
  # No TTE endpoint: 0 rows, which data.frame() gives `integer(0)` row names.
  iv <- ferx_example("one_cpt_iv")
  empty <- capture_glue("ferx_rust_predict_survival",
                        ferx_predict_survival(iv$model, iv$data, times = c(0, 12)))
  expect_length(empty[[1]][[1]], 0L)
  expect_data_frame_layout(empty[[1]], surv_cols)

  ad <- ferx_example("adaptive_tdm")
  run <- capture_glue("ferx_rust_simulate_adaptive",
                      ferx_simulate_adaptive(ad$model, ad$data, n_sim = 2L, seed = 1L))[[1]]
  for (tbl in c("doses", "decisions", "metrics")) expect_gt(length(run[[tbl]][[1]]), 0L)
  expect_data_frame_layout(run$doses, c(
    DRAW = "integer", SIM = "integer", REP = "integer", ID = "character",
    TIME = "numeric", AMT = "numeric", CMT = "integer", RATE = "numeric",
    DECISION = "integer", SIGNAL = "numeric", LATENT = "numeric", RULE = "character"))
  expect_data_frame_layout(run$decisions, c(
    DRAW = "integer", SIM = "integer", REP = "integer", ID = "character",
    DECISION = "integer", TIME = "numeric", SIGNAL = "numeric", LATENT = "numeric",
    OUTCOME = "character", N_DOSED = "integer", INTERVAL = "numeric"))
  expect_data_frame_layout(run$metrics, c(
    DRAW = "integer", SIM = "integer", REP = "integer", ID = "character",
    CUM_DOSE = "numeric", N_DOSES = "integer", N_INCREASES = "integer",
    N_DECREASES = "integer", N_HOLDS = "integer", N_NOT_DUE = "integer",
    DISCONTINUED = "logical", TIME_TO_DISCONT = "numeric", SIGNAL_MIN = "numeric",
    SIGNAL_MAX = "numeric", SIGNAL_MEAN = "numeric",
    PCT_READINGS_IN_WINDOW = "numeric", PCT_TIME_IN_WINDOW = "numeric",
    AUC_TARGET_ATTAINMENT = "numeric", ETA_SUPPLIED = "logical",
    DOSE_ORACLE = "logical"))

  npde <- capture_glue("ferx_rust_npde_from_fit",
                       suppressWarnings(ferx_calc_npde(pktte_fit(), nsim = 50L, seed = 5L)))[[1]]
  expect_data_frame_layout(npde, c(
    ID = "character", TIME = "numeric", NPDE = "numeric", NPD = "numeric",
    EPRED = "numeric", EWRES = "numeric", PWRES = "numeric"))
  tte <- attr(npde, "tte_npde", exact = TRUE)
  expect_gt(length(tte[[1]]), 0L)
  expect_data_frame_layout(tte, c(ID = "character", TIME = "numeric",
                                  CMT = "integer", NPDE = "numeric"))
})

test_that("ferx_calc_npde() on a TTE fit ignores a data.frame in the global environment (#780)", {
  # The time-to-event records are scored in a table of their own
  # (`tte_npde_table` in the glue) and placed on the sdtab's TTE rows.
  fit <- pktte_fit()
  run <- function() suppressWarnings(ferx_calc_npde(fit, nsim = 100L, seed = 5L))
  ref <- run()
  tte <- ferx:::.fitrx_sdtab_blocks(ref)$tte
  expect_gt(length(tte), 0L)
  expect_true(all(is.finite(ref$sdtab$NPDE[tte])))
  expect_identical(with_shadowed_data_frame(run())$sdtab, ref$sdtab)
})
