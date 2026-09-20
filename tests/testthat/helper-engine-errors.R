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

# The same names, for generating one test per entry point at file level without
# running the fit there: the closures touch `fit` only when called.
engine_entry_point_names <- function() names(engine_entry_points(NULL))

# One entry point, with the cached covariance fit behind its `fit =` form.
# Called inside `test_that()`, so a fit that fails is that test's failure.
engine_entry_point <- function(name) engine_entry_points(warfarin_fit_cov())[[name]]

# -- Refused inputs, built on the bundled warfarin example unless said ---------

model_with_lines <- function(lines, env) {
  ex   <- ferx_example("warfarin")
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(c(readLines(ex$model), "", lines), path)
  path
}

# #385's report: a block the engine does not know.
unknown_block_model <- function(env = parent.frame()) {
  model_with_lines(c("[not_a_block]", "  x = 1"), env)
}

# A `[data_selection]` clause with no right-hand side.
bad_selection_model <- function(env = parent.frame()) {
  model_with_lines(c("[data_selection]", "  ignore = DV <"), env)
}

# A `[data_selection]` clause whose refusal quotes a `%` back at the user. The
# glue used to hand its message to `Rf_error()` as the printf format (#388), so
# `5%': r` was read as a conversion; it now goes through a `"%s"` the glue owns.
percent_selection_model <- function(env = parent.frame()) {
  model_with_lines(c("[data_selection]", "  ignore = DV < 5%"), env)
}

# The same, with a literal `%%` in the user's own text: a glue that escaped `%`
# as `%%` to suit a printf format would print this back as `5%`.
percent_escape_model <- function(env = parent.frame()) {
  model_with_lines(c("[data_selection]", "  ignore = DV < 5%%"), env)
}

# The same, with a `%s` chain: printf would read a pointer off the stack and
# dereference it, ending the session. Only ever run in a `callr` child.
percent_chain_model <- function(env = parent.frame()) {
  model_with_lines(c("[data_selection]", "  ignore = DV < %s%s%s%s"), env)
}

# A `.ferxsearch` file whose MFL carries a `%`: the third input space, text the
# user typed rather than a model the engine parsed.
percent_search_config <- function(env = parent.frame()) {
  ex   <- ferx_example("warfarin")
  path <- withr::local_tempfile(fileext = ".ferxsearch", .local_envir = env)
  writeLines(c("[run]",
               sprintf('model = "%s"', ex$model),
               sprintf('data = "%s"', ex$data),
               "mfl = X5%dY(1)"), path)
  path
}

example_rows <- function(example = "warfarin") {
  utils::read.csv(ferx_example(example)$data, stringsAsFactors = FALSE,
                  na.strings = c(".", "NA"))
}

write_nonmem_csv <- function(d, env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".csv", .local_envir = env)
  utils::write.csv(d, path, row.names = FALSE, quote = FALSE, na = ".")
  path
}

# A dataset the reader refuses: no TIME column.
no_time_data <- function(env = parent.frame()) {
  d <- example_rows()
  write_nonmem_csv(d[, setdiff(names(d), "TIME")], env)
}

# A dataset that reads, but breaks a precondition of predict / simulate: a
# zero-order infusion into `CMT = 0`. At the pinned engine this is a panic that
# extendr turns into an R error; FeRx-NLME/ferx-core#898 turns it into an `Err`.
# It is also the one validation error (`E_DOSE_CMT_NOT_INFUSABLE`) the
# mislabelling tests put next to an unrelated failure.
infusion_into_cmt0_data <- function(example = "warfarin", env = parent.frame()) {
  d <- example_rows(example)
  if (!"RATE" %in% names(d)) d$RATE <- 0
  d$RATE[d$EVID == 1] <- 10
  d$CMT[d$EVID == 1]  <- 0
  write_nonmem_csv(d, env)
}

# The same refusal, on a dataset whose subject id carries the `%`: the second
# input space of #388, and the one that arrives through the panic path rather
# than through a glue `Err` (the diagnostic names the subject).
percent_id_data <- function(id = "s1%d", example = "warfarin",
                            env = parent.frame()) {
  d <- example_rows(example)
  if (!"RATE" %in% names(d)) d$RATE <- 0
  d$RATE[d$EVID == 1] <- 10
  d$CMT[d$EVID == 1]  <- 0
  d$ID <- id
  write_nonmem_csv(d, env)
}

# What holds for every refusal, whichever channel the engine used to fail: a
# condition, nothing on the console in its place, and a message that is the
# engine's own text rather than a wrapper's.
expect_refusal <- function(probe, phrase) {
  expect_s3_class(probe$cond, "error")
  expect_identical(probe$printed, character(0))

  msg <- conditionMessage(probe$cond)
  expect_true(grepl(phrase, msg, fixed = TRUE), info = msg)
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

# A refusal the validation pass cannot name: the engine's prose and no code -
# in particular not the code of some other finding in the same model or data.
expect_uncoded_refusal <- function(probe, phrase) {
  expect_refusal(probe, phrase)
  expect_false(inherits(probe$cond, "ferx_engine_error"))
  expect_false(grepl("[E_", conditionMessage(probe$cond), fixed = TRUE),
               info = conditionMessage(probe$cond))
}
