# Simulating from a *design* template: dosing plus sampling times, with `DV`
# left missing (`.` / `NA`) because the DV is the column the run is about to
# produce (ferx-core #957, ferx#286).
#
# The engine's fitting reader skips an EVID=0 / MDV=0 row whose DV is missing —
# right when the DV is an input, since a forgotten MDV must not be scored as
# DV=0 (ferx-core #258). Simulation now reads such a row through
# `read_population_for_simulation()` instead, which keeps it as a sampling time.
# Before that swap `ferx_simulate()` returned zero rows for the most natural
# template there is, and only a placeholder number in the DV column made it work.

# Two subjects, one bolus dose each, three sampling times each, no DV values.
write_design_template <- function(dv = ".") {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      ID   = rep(1:2, each = 4L),
      TIME = rep(c(0, 0.5, 4, 24), times = 2L),
      DV   = dv,
      EVID = rep(c(1L, 0L, 0L, 0L), times = 2L),
      AMT  = rep(c(100, NA, NA, NA), times = 2L),
      CMT  = 1L,
      MDV  = rep(c(1L, 0L, 0L, 0L), times = 2L)
    ),
    path,
    row.names = FALSE,
    na = "."
  )
  path
}

test_that("ferx_simulate simulates a design template with no DV values", {
  ex <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, write_design_template(), n_sim = 1L, seed = 1L)

  expect_s3_class(sim, "data.frame")
  # Two subjects x three sampling times: the rows that used to vanish.
  expect_equal(nrow(sim), 6L)
  expect_equal(sort(unique(sim$TIME)), c(0.5, 4, 24))
  expect_setequal(as.character(sim$ID), c("1", "2"))
  expect_true(all(is.finite(sim$DV_SIM)))
  expect_true(all(is.finite(sim$IPRED)))
})

test_that("an NA DV template simulates the same as a `.` one", {
  # `.` and `NA` are the same missing cell to the reader; an R user writing the
  # frame with `NA` should not have to know which sentinel gets written.
  ex <- ferx_example("warfarin")
  dot <- ferx_simulate(ex$model, write_design_template("."), n_sim = 1L, seed = 7L)
  na <- ferx_simulate(ex$model, write_design_template(NA), n_sim = 1L, seed = 7L)

  expect_equal(nrow(na), nrow(dot))
  expect_equal(na$DV_SIM, dot$DV_SIM)
})

test_that("MDV = 1 still excludes a sampling row from the simulation", {
  # MDV=1 is the user explicitly saying the record is not an observation, which
  # is unambiguous on either path. Only the *forgotten* MDV is reread.
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      ID   = 1L,
      TIME = c(0, 0.5, 4, 24),
      DV   = ".",
      EVID = c(1L, 0L, 0L, 0L),
      AMT  = c(100, NA, NA, NA),
      CMT  = 1L,
      MDV  = c(1L, 1L, 0L, 0L)
    ),
    path,
    row.names = FALSE,
    na = "."
  )
  ex <- ferx_example("warfarin")
  sim <- ferx_simulate(ex$model, path, n_sim = 1L, seed = 1L)

  expect_equal(nrow(sim), 2L)
  expect_equal(sort(sim$TIME), c(4, 24))
})

test_that("a template with DV values is unaffected", {
  # The workaround users resorted to — a placeholder number in the column about
  # to be overwritten — must keep working, and give the same simulated values as
  # the DV-less template it replaces.
  ex <- ferx_example("warfarin")
  placeholder <- ferx_simulate(ex$model, write_design_template("0"),
                               n_sim = 1L, seed = 3L)
  missing <- ferx_simulate(ex$model, write_design_template("."),
                           n_sim = 1L, seed = 3L)

  expect_equal(nrow(placeholder), 6L)
  expect_equal(placeholder$DV_SIM, missing$DV_SIM)
})
