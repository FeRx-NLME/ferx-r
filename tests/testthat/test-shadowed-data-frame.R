# A user's own `data.frame` in the global environment must not reach the glue
# (#780). extendr's `data_frame!` evaluated `data.frame(..)` from there, so a
# shadowing function made each of these entry points panic and print the whole
# table. The glue now builds every table with `r_data_frame()`, which never
# calls R; `tools/check-glue-data-frame.sh` keeps the macro out. Each test
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
