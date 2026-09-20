# A `%` in an engine error message is ordinary text in this domain - `CV%`,
# `5%`, a column name, a path like `my%20data.csv` - and until #388 it was read
# as a printf conversion: `Rf_error()` takes its first argument as its *format*
# string, and the glue handed it the engine's message. `%d` printed whatever
# happened to be on the stack in place of the text; a `%s` chain ended the R
# session with "An irrecoverable exception occurred".
#
# Every `#[extendr]` body now returns a `Result` and `entry()` raises once,
# after that body has returned, through a `"%s"` format the glue owns. What is
# tested here is the four places a `%` can enter from: the model file, the
# dataset, MFL / `.ferxsearch` text, and an R argument.
#
# Two rules for the fixtures:
#
#   * in-process cases carry `%d` only. A `%s` chain against a build that
#     regressed segfaults, which would take the whole test run with it, so the
#     `%s` cases run in a `callr` child where only the child dies.
#   * nothing asserts what the bug printed. Those digits are whatever was on
#     that platform's stack; only the verbatim text is portable.
#
# The assertion itself is `expect_refusal()` from helper-engine-errors.R, the
# same one the twin file uses: a `%` must reach R verbatim *and* the refusal
# must still leave the console clean (#385), and these 15 entry points appear
# in no other test that would notice if it stopped doing so.

# The engine's own words for `[data_selection] ignore = DV < 5%`, which quote
# the `%` back twice.
PERCENT_PARSE <- "malformed filter expression 'DV < 5%': right-hand side '5%' is not a number"

# -- The model file -----------------------------------------------------------
#
# One test per glue function that quotes a refused model back. The nine entry
# points of #386 are covered by test-engine-errors-raise.R; these are the ones
# that were raising through `throw_r_error` on `main`.

percent_model_callers <- function() {
  ex <- ferx_example("warfarin")
  repointed <- function(m) {
    fit <- warfarin_fit_cov()
    fit$model_path <- m
    fit$model_hash <- NULL
    fit$data_hash <- NULL
    fit
  }
  list(
    "ferx_fit()" = function(m) ferx_fit(m, ex$data, verbose = FALSE),
    "ferx_simulate_adaptive()" = function(m) {
      ferx_simulate_adaptive(m, ex$data, n_sim = 1L)
    },
    "ferx_bootstrap()" = function(m) ferx_bootstrap(m, ex$data, samples = 1),
    "ferx_model_to_frem()" = function(m) {
      ferx_model_to_frem(m, ex$data, covariates = "WT")
    },
    "ferx_allometry()" = function(m) ferx_allometry(m, ex$data, covariate = "WT"),
    "ferx_ruvsearch()" = function(m) ferx_ruvsearch(m, ex$data),
    "ferx_iovsearch()" = function(m) ferx_iovsearch(m, ex$data),
    "ferx_covariance()" = function(m) ferx_covariance(repointed(m)),
    "ferx_sir()" = function(m) {
      ferx_sir(repointed(m), sir_samples = 20L, sir_resamples = 5L)
    }
  )
}

for (nm in names(percent_model_callers())) {
  local({
    nm <- nm

    test_that(paste(nm, "quotes a '%' from the model back as written"), {
      probe <- engine_error_probe(percent_model_callers()[[nm]](percent_selection_model()))
      expect_refusal(probe, PERCENT_PARSE)
    })
  })
}

# -- MFL and .ferxsearch text -------------------------------------------------
#
# The second input space: a search space the user typed. The MFL lexer names
# the character it stopped on, and that character is the `%` - which is exactly
# what used to be eaten, leaving `unexpected character `` `` at offset 2`. The
# `.ferxsearch` case below is the one that quotes the whole line back.

percent_mfl_callers <- function() {
  ex <- ferx_example("warfarin")
  list(
    "ferx_covsearch()" = function(s) ferx_covsearch(ex$model, ex$data, search_space = s),
    "ferx_modelsearch()" = function(s) ferx_modelsearch(ex$model, ex$data, search_space = s),
    "ferx_iivsearch()" = function(s) ferx_iivsearch(ex$model, ex$data, search_space = s),
    "ferx_amd()" = function(s) ferx_amd(ex$model, ex$data, search_space = s),
    "ferx_globalsearch()" = function(s) ferx_globalsearch(ex$model, ex$data, search_space = s),
    "ferx_search_space()" = function(s) ferx_search_space(s)
  )
}

for (nm in names(percent_mfl_callers())) {
  local({
    nm <- nm

    test_that(paste(nm, "quotes a '%' from the search space back as written"), {
      probe <- engine_error_probe(percent_mfl_callers()[[nm]]("X5%dY(1)"))
      expect_refusal(probe, "unexpected character `%` at offset 2")
    })
  })
}

test_that("ferx_search_config() quotes a '%' from the .ferxsearch file back as written", {
  probe <- engine_error_probe(ferx_search_config(percent_search_config()))
  expect_refusal(probe, "X5%dY(1)")
})

# -- An R argument ------------------------------------------------------------
#
# The input space the issue did not have: the text is not the engine quoting a
# file, it is the name the caller passed.

test_that("ferx_fit() quotes a '%' in a settings key back as written", {
  ex    <- ferx_example("warfarin")
  probe <- engine_error_probe(
    ferx_fit(ex$model, ex$data, verbose = FALSE, settings = list("bad%dkey" = 1))
  )
  expect_refusal(probe, "`bad%dkey`")
})

# -- The dataset --------------------------------------------------------------
#
# The input space that does not arrive through a glue `Err`: the engine's dose
# diagnostic names the subject, and at the pinned engine it is a panic, which
# extendr raised through `throw_r_error` in turn.
# FeRx-NLME/ferx-core#1487 turns this particular one into an `Err`; the
# assertion is the same either way, which is the point of making it here.

test_that("a '%' in a subject id comes back as written", {
  probe <- engine_error_probe(
    ferx_predict(ferx_example("warfarin")$model, percent_id_data())
  )
  expect_refusal(probe, "subject s1%d, time 0")
})

# -- A literal `%%` in the user's text ----------------------------------------
#
# The other half of "verbatim": escaping `%` as `%%` before handing the text to
# a printf format prints `5%%` as `5%`. extendr main (b0cb8a81, unreleased)
# switches `throw_r_error` to `"%s"`, at which point a glue that escaped would
# print `5%%%%` - so this is the test that makes that upgrade safe either way.

test_that("a literal '%%' in the model survives as '%%'", {
  probe <- engine_error_probe(
    ferx_fit(percent_escape_model(), ferx_example("warfarin")$data, verbose = FALSE)
  )
  expect_refusal(probe, "'DV < 5%%'")
  msg <- conditionMessage(probe$cond)
  expect_false(grepl("5%%%%", msg, fixed = TRUE), info = msg)
})

# -- The cases that kill the process ------------------------------------------
#
# `%s` reads a pointer off the stack and dereferences it, so a regressed build
# aborts rather than printing nonsense. The child loads the build under test:
# under `devtools::test()` a bare `library(ferx)` would load whatever is
# installed, which on a development machine is some older build - it would pass
# or crash for the wrong reason.

# Not `skip_on_cran()`: these three are the regression tests for the abort
# itself, and a run that skips them says nothing. `callr` is in Imports, so the
# child is always available, and each one is a single failed parse.
child_message <- function(expr_text) {
  dev  <- requireNamespace("pkgload", quietly = TRUE) &&
    pkgload::is_dev_package("ferx")
  path <- if (dev) pkgload::pkg_path() else NULL
  callr::r(
    function(dev, path, expr_text) {
      if (dev) {
        pkgload::load_all(path, compile = FALSE, quiet = TRUE)
      } else {
        library(ferx)
      }
      tryCatch(
        {
          eval(parse(text = expr_text), envir = globalenv())
          "no error"
        },
        error = function(e) conditionMessage(e)
      )
    },
    args = list(dev = dev, path = path, expr_text = expr_text)
  )
}

test_that("a '%s' chain in the model raises instead of ending the session", {
  model <- percent_chain_model()
  data  <- ferx_example("warfarin")$data
  msg <- child_message(sprintf(
    'ferx::ferx_fit("%s", "%s", verbose = FALSE)', model, data
  ))
  expect_true(grepl("%s%s%s%s", msg, fixed = TRUE), info = msg)
})

test_that("a '%s' chain in a settings key raises instead of ending the session", {
  ex <- ferx_example("warfarin")
  msg <- child_message(sprintf(
    'ferx::ferx_fit("%s", "%s", verbose = FALSE, settings = list("%s" = 1))',
    ex$model, ex$data, "%s%s%s%s"
  ))
  expect_true(grepl("%s%s%s%s", msg, fixed = TRUE), info = msg)
})

# -- The panic path -----------------------------------------------------------
#
# extendr's own wrapper raises a panic's payload through `throw_r_error`, so a
# panic whose text carries a `%` had the same two failures. `entry()` catches
# the unwind before extendr sees it. The fixture is a glue function that panics
# on demand: the panics that are reachable through a dataset today
# (FeRx-NLME/ferx-core#1487 turns one of them into an `Err`) are not a durable
# way to reach this path.

test_that("a panic message with a '%' arrives verbatim", {
  probe <- engine_error_probe(ferx:::ferx_rust_test_panic("50% of %d subjects"))
  expect_refusal(probe, "50% of %d subjects")
})

test_that("a panic message with a '%s' chain raises instead of ending the session", {
  msg <- child_message('ferx:::ferx_rust_test_panic("%s%s%s%s")')
  expect_true(grepl("%s%s%s%s", msg, fixed = TRUE), info = msg)
})

# A payload that is neither `&str` nor `String` carries no text, and this is the
# one raise-position message ferx writes itself rather than passing through from
# the engine. extendr names the function in that case; `entry()` is one function
# for all 44, so it is `#[track_caller]` and the caller's line stands in for the
# name.
#
# The line is what makes the assertion worth anything, so read the line the
# fixture actually calls `entry(` on out of lib.rs and compare. The package
# source is not there for an installed package or under `R CMD check`, and the
# fallback is a bound: `ferx_rust_test_panic` calls `entry(` around lib.rs:3370
# while `entry()` itself is at lib.rs:166, so a dropped `#[track_caller]` makes
# `Location::caller()` report its own line, ~170. The bound is the weaker test -
# it would go false-green if `entry()` ever moved past line 1000 - which is why
# it is only the fallback.
glue_entry_line <- function(fn) {
  p <- file.path("..", "..", "src", "rust", "src", "lib.rs")
  if (!file.exists(p)) return(NA_integer_)
  src   <- readLines(p, warn = FALSE)
  start <- grep(sprintf("^fn %s\\(", fn), src)
  if (length(start) != 1L) return(NA_integer_)
  window <- seq(start, min(start + 20L, length(src)))
  hit    <- grep("^\\s*entry\\(", src[window])
  if (length(hit) < 1L) return(NA_integer_)
  window[hit[1]]
}

test_that("a panic carrying no text names the entry point it came out of", {
  probe <- engine_error_probe(ferx:::ferx_rust_test_panic("<non-string payload>"))
  expect_refusal(
    probe, "ferx: the engine panicked without a message, in the entry point at "
  )
  msg  <- conditionMessage(probe$cond)
  expect_true(grepl("src/lib.rs:", msg, fixed = TRUE), info = msg)
  line <- as.integer(sub("^.*src/lib\\.rs:([0-9]+).*$", "\\1", msg))
  expect_false(is.na(line), info = msg)

  want <- glue_entry_line("ferx_rust_test_panic")
  if (is.na(want)) expect_gt(line, 1000L) else expect_identical(line, want)
})

# The other half of the same restructure - that a refused call no longer leaks
# what the body was holding (#389) - is a memory measurement, not an assertion:
# 3000 calls take ~13 minutes and RSS is too noisy to threshold in CI. It lives
# in tools/measure-refusal-rss.R, to be run by hand when the raise path changes.
