# The simulate/predict warnings channel, and the missing-value marker in the
# returned frame (ferx-r #283).
#
# Two defects, both about something the caller cannot see:
#   1. Every data-reader diagnostic the engine raised for the dataset was
#      dropped on the simulate path, and the predict path had no warnings
#      channel at all. Only a count re-derived in the R glue, for the one
#      empty-DV case it was written for, ever reached the caller.
#   2. Columns with no value carried a bare `NaN` instead of `NA`, so an
#      `OBSERVED` that simply has no event flag (every non-TTE row) printed as
#      an arithmetic failure.

# A dataset with scored observations and no dose records at all: the engine
# answers with `W_NO_DOSES`, which is a data-reader diagnostic with nothing to do
# with a missing DV -- so it shows the whole channel is open, not just the one
# case the old glue counted.
write_no_dose_data <- function() {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      ID   = rep(1:2, each = 3L),
      TIME = rep(c(0.5, 4, 24), times = 2L),
      DV   = c(5.3, 3.1, 1.2, 4.8, 2.9, 1.0),
      EVID = 0L,
      AMT  = ".",
      CMT  = 1L,
      MDV  = 0L
    ),
    path,
    row.names = FALSE,
    quote = FALSE
  )
  path
}

test_that("ferx_simulate surfaces a data-reader diagnostic unrelated to the DV", {
  ex <- ferx_example("warfarin")
  expect_warning(
    sim <- ferx_simulate(ex$model, write_no_dose_data(), n_sim = 1L, seed = 1L),
    "W_NO_DOSES"
  )
  w <- attr(sim, "simulation_warnings", exact = TRUE)
  expect_true(any(grepl("^W_NO_DOSES", w)))
})

test_that("ferx_predict surfaces data-reader diagnostics too", {
  # `ferx_predict()` returned a bare frame with no warnings channel whatsoever,
  # so the same dataset problem reached the caller nowhere at all.
  ex <- ferx_example("warfarin")
  expect_warning(
    pred <- ferx_predict(ex$model, write_no_dose_data()),
    "W_NO_DOSES"
  )
  expect_true(any(grepl("^W_NO_DOSES", attr(pred, "simulation_warnings", exact = TRUE))))
})

test_that("ferx_predict(fit = ) surfaces data-reader diagnostics too", {
  # The `fit =` path is a separate glue entry point, and the common one after a
  # fit; it must not be the quiet half of the pair.
  ex <- ferx_example("warfarin")
  fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
  expect_warning(
    pred <- ferx_predict(ex$model, write_no_dose_data(), fit = fit),
    "W_NO_DOSES"
  )
  expect_true(any(grepl("^W_NO_DOSES", attr(pred, "simulation_warnings", exact = TRUE))))
})

test_that("a clean dataset raises no simulation warning", {
  # The channel must stay quiet on data the engine has nothing to say about,
  # or every VPC run trains its user to ignore it.
  ex <- ferx_example("warfarin")
  expect_warning(
    sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L),
    NA
  )
  expect_length(attr(sim, "simulation_warnings", exact = TRUE), 0L)
  expect_warning(pred <- ferx_predict(ex$model, ex$data), NA)
  expect_length(attr(pred, "simulation_warnings", exact = TRUE), 0L)
})

# Every empty cell in the returned frame is asserted twice, and both halves are
# load-bearing (Codex review on ferx-r#381):
#   is.na(x)  alone passes for a NaN, which is the bug.
#   !is.nan(x) alone passes for any ordinary number, e.g. a wrong 0.
# Only the pair pins the cell to NA_real_.
expect_na_real <- function(x) {
  expect_true(all(is.na(x)))
  expect_false(any(is.nan(x)))
}

test_that("OBSERVED is NA, not NaN, on a continuous row", {
  # Specifically NA_real_: `is.na()` was already TRUE for a NaN, so the
  # documented `is.na(OBSERVED)` idiom never broke -- but the column printed as
  # `NaN` for every row of an ordinary PK simulation, which reads as a failed
  # computation rather than "this row has no event flag".
  ex <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L)

  expect_na_real(sim$OBSERVED)
  # DV_SIM and IPRED are real numbers here; nothing was turned into NA.
  expect_true(all(is.finite(sim$DV_SIM)))
  expect_true(all(is.finite(sim$IPRED)))
})

test_that("a TTE row's empty DV_SIM/IPRED are NA, not NaN", {
  ex <- ferx_example("pktte_joint")
  sim <- ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, horizon = 24)

  ev <- sim[!is.na(sim$OBSERVED), ]
  expect_gt(nrow(ev), 0)
  # The event row has no Gaussian prediction; that absence is NA, not a NaN and
  # not a stand-in number.
  expect_na_real(ev$DV_SIM)
  expect_na_real(ev$IPRED)
  # OBSERVED itself is a real 0/1 on these rows...
  expect_true(all(ev$OBSERVED %in% c(0, 1)))
  # ...and NA_real_ on the continuous rows of the same frame. This is the case a
  # warfarin-only check cannot reach: OBSERVED is all-empty there, while here the
  # column mixes real flags with empty cells.
  expect_na_real(sim$OBSERVED[is.na(sim$OBSERVED)])
  expect_true(all(is.finite(sim$DV_SIM[is.na(sim$OBSERVED)])))
})
