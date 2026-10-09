# A `ferx_engine_error` shows its suggestion in its message, once, the same way
# whichever constructor built it (#504): the engine's text, ` [CODE]` on the
# same line, then `hint: <suggestion>` on the next - unless the text already
# gives that advice. All of it is `.ferx_engine_message()`.

# How often `needle` occurs in `haystack`, literally.
count_fixed <- function(needle, haystack) {
  sum(gregexpr(needle, haystack, fixed = TRUE)[[1L]] > 0L)
}

# A copy of warfarin whose CL line is `edit`ed, and what fitting it raised.
ses_model <- function(edit, env = parent.frame()) {
  ex <- ferx_example("warfarin")
  txt <- readLines(ex$model)
  i <- grep("^\\s*CL\\s*=", txt)[1L]
  txt[i] <- edit(txt[i])
  path <- withr::local_tempfile(fileext = ".ferx", .local_envir = env)
  writeLines(txt, path)
  path
}
ses_err <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
ses_fit <- function(model) {
  ses_err(ferx_fit(model, ferx_example("warfarin")$data, covariance = FALSE,
                   verbose = FALSE))
}

# -- Re-validation path: ferx_fit(), and a non-fit entry point -----------------

# E_ETA_NOT_DECLARED's suggestion (`declare \`omega ETA_FOO ~ <variance>\``) is
# not in the engine's text, so it is the case that must gain a hint line.
test_that("ferx_fit() and ferx_predict() show a suggestion the text lacks, once, on its own line", {
  model <- ses_model(function(l) paste0(l, " * exp(ETA_FOO)"))
  got <- list(
    fit = ses_fit(model),
    predict = ses_err(ferx_predict(model, ferx_example("warfarin")$data))
  )
  for (k in names(got)) {
    e <- got[[k]]
    expect_s3_class(e, "ferx_engine_error")
    expect_identical(e$code, "E_ETA_NOT_DECLARED", info = k)
    expect_identical(e$suggestion, "declare `omega ETA_FOO ~ <variance>`",
                     info = k)
    lines <- strsplit(conditionMessage(e), "\n", fixed = TRUE)[[1L]]
    expect_length(lines, 2L)
    # The engine's prose and the code stay on the first line ...
    expect_match(lines[1L], "Model references ETA_FOO as a random effect",
                 fixed = TRUE, info = k)
    expect_true(endsWith(lines[1L], "[E_ETA_NOT_DECLARED]"), info = k)
    # ... the hint follows, and the suggestion appears exactly once.
    expect_identical(lines[2L], "hint: declare `omega ETA_FOO ~ <variance>`",
                     info = k)
    expect_identical(count_fixed(e$suggestion, conditionMessage(e)), 1L,
                     info = k)
  }
  # Both entry points render the same message after their stage prefix, byte
  # for byte (one arrives marked UTF-8, the other as unmarked bytes).
  tail_of <- function(e) {
    msg <- conditionMessage(e)
    at <- regexpr("Model references", msg, fixed = TRUE, useBytes = TRUE)
    charToRaw(substring(msg, at))
  }
  expect_identical(tail_of(got$fit), tail_of(got$predict))
})

# The parser writes "did you mean `[fit_options]`" into its own text, and the
# data check "Available covariate columns: ..." (capitalised): appending the
# suggestion would say each twice.
test_that("a suggestion the engine's text already gives is not repeated", {
  ex <- ferx_example("warfarin")
  block <- withr::local_tempfile(fileext = ".ferx")
  writeLines(c(readLines(ex$model), "", "[fit_option]", "  method = focei"),
             block)
  got <- list(
    unknown_block = ses_fit(block),
    missing_cov = ses_fit(ses_model(function(l) paste0(l, " * (WTX / 70)")))
  )
  expect_identical(got$unknown_block$code, "E_UNKNOWN_BLOCK")
  expect_identical(got$missing_cov$code, "E_MISSING_COVARIATE")
  for (k in names(got)) {
    e <- got[[k]]
    expect_false(is.na(e$suggestion), info = k)
    msg <- conditionMessage(e)
    expect_false(grepl("\n", msg, fixed = TRUE), info = k)
    expect_true(endsWith(msg, sprintf("[%s]", e$code)), info = k)
  }
})

# -- Coded path: a refusal carrying ferx-core's own record ---------------------

# No coded refusal R can reach carries a suggestion yet (E_COV_LEVEL_UNKNOWN and
# E_COVARIATE_STATS_BINDING write their advice into the message), so the record
# is stubbed. `text` is what the glue raised; `message` is it without the
# folded suggestion (ferx-core `EngineError::with_suggestion_in_display`).
ses_coded <- function(text, message, suggestion) {
  testthat::local_mocked_bindings(
    ferx_rust_take_engine_diagnostic = function() {
      list(text = text, message = message, code = "E_TEST", block = "",
           line = 0L, suggestion = suggestion)
    },
    .package = "ferx"
  )
  ferx:::.ferx_engine_coded_error(simpleError(text), text)
}

test_that("coded path: a suggestion is shown once, built from the record's message", {
  plain <- ses_coded("Error predicting: no theta for level 4",
                     "Error predicting: no theta for level 4",
                     "bind the level first")
  expect_s3_class(plain, "ferx_engine_error")
  expect_identical(
    conditionMessage(plain),
    "Error predicting: no theta for level 4 [E_TEST]\nhint: bind the level first"
  )

  # The engine folded the advice into the raised text: the message is built
  # from `message`, so the hint is the only place the advice appears.
  folded <- ses_coded("Error predicting: no theta for level 4 Bind the level first.",
                      "Error predicting: no theta for level 4",
                      "bind the level first")
  expect_identical(conditionMessage(folded), conditionMessage(plain))
})

test_that("a refusal without a suggestion is unchanged: text, then the code", {
  none <- ses_coded("Error predicting: no theta", "Error predicting: no theta", "")
  expect_identical(conditionMessage(none), "Error predicting: no theta [E_TEST]")
  expect_true(is.na(none$suggestion))
  # The re-validation path hands over NA for a diagnostic without one.
  expect_identical(ferx:::.ferx_engine_message("x", "E_X", NA_character_),
                   "x [E_X]")
})
