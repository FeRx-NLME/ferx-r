# Shared by test-engine-errors-raise.R (#385).

# Run `expr`, recording everything a caller could observe about a failure: the
# condition (if one was raised), the value (if one came back), and whatever
# reached the console through R's output stream. The third is the point - the
# old glue `rprintln!`-ed the engine's message and returned NULL, so "nothing
# was raised" and "something was printed" have to be seen separately. Warnings
# are muffled: the data-reader diagnostics are not what these tests are about.
engine_error_probe <- function(expr) {
  cond <- NULL
  value <- NULL
  printed <- utils::capture.output(
    value <- tryCatch(
      suppressWarnings(expr),
      error = function(e) {
        cond <<- e
        NULL
      }
    )
  )
  list(cond = cond, value = value, printed = printed)
}

# Every R entry point whose glue used to print-and-return-NULL, in both the
# model + data form and the `fit =` form, as functions of (model, data).
engine_entry_points <- function(fit) {
  list(
    "ferx_predict()" = function(m, d) ferx_predict(m, d),
    "ferx_predict(fit = )" = function(m, d) ferx_predict(m, d, fit = fit),
    "ferx_simulate()" = function(m, d) ferx_simulate(m, d, n_sim = 1L, seed = 1L),
    "ferx_simulate(fit = )" = function(m, d) {
      ferx_simulate(m, d, n_sim = 1L, seed = 1L, fit = fit)
    },
    "ferx_simulate_with_uncertainty()" = function(m, d) {
      ferx_simulate_with_uncertainty(m, d, fit = fit, n_uncertainty_draws = 2L,
                                     n_sim_per_draw = 1L, seed = 1L)
    },
    "ferx_predict_survival()" = function(m, d) {
      ferx_predict_survival(m, d, times = c(1, 2))
    },
    "ferx_predict_survival(fit = )" = function(m, d) {
      ferx_predict_survival(m, d, times = c(1, 2), fit = fit)
    },
    "ferx_calc_npde()" = function(m, d) {
      ferx_calc_npde(fit, nsim = 20L, seed = 1L, model = m, data = d)
    },
    "ferx_inits_from_nca()" = function(m, d) ferx_inits_from_nca(m, d)
  )
}

# -- Refused inputs, all built on the bundled warfarin example ----------------

# #385's report: a block the engine does not know.
unknown_block_model <- function(env = parent.frame()) {
  ex   <- ferx_example("warfarin")
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(readLines(ex$model), "", "[not_a_block]", "  x = 1"), path)
  path
}

# A `[data_selection]` clause with no right-hand side.
bad_selection_model <- function(env = parent.frame()) {
  ex   <- ferx_example("warfarin")
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(readLines(ex$model), "", "[data_selection]", "  ignore = DV <"), path)
  path
}

warfarin_rows <- function() {
  utils::read.csv(ferx_example("warfarin")$data, stringsAsFactors = FALSE,
                  na.strings = c(".", "NA"))
}

write_nonmem_csv <- function(d, env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
  utils::write.csv(d, path, row.names = FALSE, quote = FALSE, na = ".")
  path
}

# A dataset the reader refuses: no TIME column.
no_time_data <- function(env = parent.frame()) {
  d <- warfarin_rows()
  write_nonmem_csv(d[, setdiff(names(d), "TIME")], env)
}

# A dataset that reads, but breaks a precondition of predict / simulate: a
# zero-order infusion into `CMT = 0`. At the pinned engine this is a panic that
# extendr turns into an R error; FeRx-NLME/ferx-core#898 turns it into an `Err`.
# The contract below has to hold on either channel.
infusion_into_cmt0_data <- function(env = parent.frame()) {
  d <- warfarin_rows()
  d$RATE[d$EVID == 1] <- 10
  d$CMT[d$EVID == 1]  <- 0
  write_nonmem_csv(d, env)
}

# What holds for every refusal, whichever channel the engine used to fail:
# a condition and no value, nothing on the console in its place, and a message
# that is the engine's own text rather than a wrapper's.
expect_refusal <- function(probe, phrase) {
  expect_s3_class(probe$cond, "error")
  expect_null(probe$value)
  expect_identical(probe$printed, character(0))

  msg <- conditionMessage(probe$cond)
  expect_true(grepl(phrase, msg, fixed = TRUE), info = msg)
  # extendr's stand-in text for a panic it could not describe
  expect_false(grepl("worker thread panicked", msg, fixed = TRUE), info = msg)
  # a prefix stacked on a prefix
  expect_false(grepl("Error: Error", msg, fixed = TRUE), info = msg)
  expect_false(grepl("^(Error [a-z ]+: ){2}", msg), info = msg)
}

# ... and, where the engine's validation pass names the finding, the same
# classed condition a refused `ferx_fit()` raises, code in the message too.
expect_coded_refusal <- function(probe, phrase, code) {
  expect_refusal(probe, phrase)
  expect_s3_class(probe$cond, "ferx_engine_error")
  expect_identical(probe$cond$code, code)
  expect_true(
    grepl(sprintf("[%s]", code), conditionMessage(probe$cond), fixed = TRUE),
    info = conditionMessage(probe$cond)
  )
}
