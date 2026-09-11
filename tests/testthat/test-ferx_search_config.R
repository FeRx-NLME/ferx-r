# `.ferxsearch` loading and validation. The loader is the engine's, so what is
# asserted here is that the R surface reports what it validated - and that a
# space the engine refuses is an R error naming the offender, before any fit.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

minimal_cfg <- function(mfl = "COVARIATE?(CL, WT, [pow, lin])", ...) {
  write_cfg(
    'base = "model.ferx"',
    "[space]",
    paste0('mfl = "', mfl, '"'),
    ...
  )
}

test_that("a minimal configuration loads with the engine's defaults", {
  cfg <- ferx_search_config(minimal_cfg())

  expect_s3_class(cfg, "ferx_search_config")
  expect_equal(basename(cfg$base), "model.ferx")
  expect_null(cfg$data)
  expect_equal(cfg$mfl, "COVARIATE?(CL, WT, [pow, lin])")

  # One statement, exploratory, reported under its keyword.
  expect_equal(nrow(cfg$space), 1L)
  expect_equal(cfg$space$keyword, "COVARIATE")
  expect_true(cfg$space$optional)

  # Defaults come from the engine, not from R: an unstated rank is left to
  # whichever tool runs the file, and there is no cutoff.
  expect_true(is.na(cfg$rank$type))
  expect_true(is.na(cfg$rank$cutoff))
  expect_true(cfg$strictness$require_converged)
  expect_equal(cfg$strictness$max_condition_number, 1000)
  expect_length(cfg$strictness_set, 0L)
  expect_null(cfg$run$threads)
  expect_false(cfg$run$resume)
  expect_length(cfg$tools, 0L)
})

test_that("an unstated [rank] type is NA, not a guess at the tool's default", {
  # `[rank] type` became optional in the engine (covsearch tests rather than
  # ranks, and refuses a BIC ranking outright), so "the file did not say" has to
  # survive the trip to R as a missing value rather than as a plausible-looking
  # "bic".
  cfg <- ferx_search_config(minimal_cfg())
  expect_true(is.na(cfg$rank$type))
  expect_true(is.na(cfg$rank$cutoff))
  expect_match(paste(capture.output(print(cfg)), collapse = "\n"),
               "Rank: (tool default)", fixed = TRUE)

  stated <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, [pow, lin])",
    "[rank]", 'type = "bic"', "cutoff = 3.84"
  ))
  expect_equal(stated$rank$type, "bic")
  expect_equal(stated$rank$cutoff, 3.84)
})

test_that("base, data and cache_dir resolve against the file's own directory", {
  dir <- tempfile()
  dir.create(dir)
  # macOS puts the session temp directory behind a symlink, and the engine
  # resolves paths against the config file's own (normalized) directory.
  dir <- normalizePath(dir)
  path <- file.path(dir, "space.ferxsearch")
  writeLines(c(
    'base = "models/m.ferx"',
    'data = "data/d.csv"',
    "[space]",
    'mfl = "COVARIATE?(CL, WT, pow)"',
    "[run]",
    'cache_dir = "cache"'
  ), path)

  cfg <- ferx_search_config(path)
  expect_equal(cfg$base, file.path(dir, "models", "m.ferx"))
  expect_equal(cfg$data, file.path(dir, "data", "d.csv"))
  expect_equal(cfg$run$cache_dir, file.path(dir, "cache"))
})

test_that("the strictness gate is the file's keys over the engine defaults", {
  cfg <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[strictness]",
    "require_covariance = true",
    "max_correlation = 0.99"
  ))

  expect_true(cfg$strictness$require_covariance)
  expect_equal(cfg$strictness$max_correlation, 0.99)
  # Untouched keys keep the engine's defaults, and are not reported as set.
  expect_true(cfg$strictness$require_converged)
  expect_equal(cfg$strictness$max_condition_number, 1000)
  expect_setequal(cfg$strictness_set,
                  c("require_covariance", "max_correlation"))
})

test_that("[run] settings come back typed, with auto threads as NULL", {
  cfg <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[run]",
    "threads = 4",
    "retries = 0",
    "resume = true"
  ))
  expect_identical(cfg$run$threads, 4L)
  expect_identical(cfg$run$retries, 0L)
  expect_true(cfg$run$resume)
})

test_that("tool sections are kept for their tool and reported by name", {
  cfg <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[covsearch]",
    "p_forward = 0.01"
  ))
  expect_equal(cfg$tools, "covsearch")
})

test_that("a coverage gap is an error naming the feature, before any fit", {
  # Needs a feature the engine does not implement yet - see the note in
  # test-ferx_search_space.R for how to repoint this when SEQ-ZO-FO lands.
  err <- expect_error(
    ferx_search_config(minimal_cfg("ABSORPTION([INST, SEQ-ZO-FO])")),
    "SEQ-ZO-FO"
  )
  # The engine explains why, rather than silently narrowing the space.
  expect_match(conditionMessage(err), "ABSORPTION\\(SEQ-ZO-FO\\)")
})

test_that("a misspelt section is an error rather than a silently ignored one", {
  expect_error(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[strictnes]",
      "require_converged = false"
    )),
    "strictnes"
  )
})

test_that("an unparseable or empty space is an error", {
  expect_error(ferx_search_config(minimal_cfg("COVARIATE?(CL, WT")), "\\[space\\]")
  expect_error(ferx_search_config(minimal_cfg("LET(X, [a, b])")), "empty")
})

test_that("the penalized rank type loads and reports itself", {
  # This test used to assert the opposite - that `penalized` was refused as an
  # unimplemented rank type. ferx-core #1185 implemented it for *every* search
  # tool, so the refusal is gone and the assertion had to flip.
  #
  # It cannot be repointed at some other unimplemented type the way the
  # SEQ-ZO-FO coverage-gap test above can: all eight `RankType` variants (ofv,
  # aic, bic, bic_mixed, bic_iiv, bic_random, bic_fixed, penalized) are now
  # implemented, so the "declared but unimplemented" category is empty. What
  # survives is the pair below - the type loads, and an unrecognised one is
  # still refused.
  cfg <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[rank]",
    'type = "penalized"'
  ))
  expect_s3_class(cfg, "ferx_search_config")
  expect_equal(cfg$rank$type, "penalized")
})

test_that("[rank.penalties] comes back as the effective schedule", {
  # The file validated the schedule and then dropped it (#348): `cfg$rank` held
  # `type` and `cutoff` only, so a config with an overlaid penalty printed
  # byte-identically to one with the defaults. What comes back now is the
  # *effective* schedule - the file's keys over pyDarwin's defaults - which is
  # what the run would charge.
  defaults <- ferx_search_config(minimal_cfg("COVARIATE?(CL, WT, pow)"))
  expect_type(defaults$rank$penalties, "double")
  expect_true(all(c("theta", "omega", "sigma", "convergence", "covariance",
                    "correlation", "max_correlation", "condition_number",
                    "max_condition_number", "non_influential", "crash",
                    "gate") %in% names(defaults$rank$penalties)))
  # pyDarwin's defaults, from the engine - never a copy kept in R.
  expect_equal(defaults$rank$penalties[["theta"]], 10)
  expect_equal(defaults$rank$penalties[["crash"]], 99999999)
  expect_length(defaults$rank$penalties_set, 0L)

  overlaid <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[rank]", 'type = "penalized"',
    "[rank.penalties]", "theta = 5.0", "gate = 0.0"
  ))
  expect_equal(overlaid$rank$penalties[["theta"]], 5)
  expect_equal(overlaid$rank$penalties[["gate"]], 0)
  # Untouched charges keep the engine's default rather than going missing.
  expect_equal(overlaid$rank$penalties[["omega"]], 10)
  expect_setequal(overlaid$rank$penalties_set, c("theta", "gate"))
})

test_that("print() shows the penalty schedule only when it matters", {
  quiet <- ferx_search_config(minimal_cfg("COVARIATE?(CL, WT, pow)"))
  expect_false(any(grepl("Penalties", capture.output(print(quiet)))))

  loud <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[rank]", 'type = "penalized"',
    "[rank.penalties]", "theta = 5.0"
  ))
  out <- paste(capture.output(print(loud)), collapse = "\n")
  expect_match(out, "Penalties")
  # The overlaid charge is starred, an untouched one is not - the same idiom
  # the strictness block uses.
  expect_match(out, "\\* theta +5")
  expect_match(out, "  omega +10")
  # The crash charge is the number the user would write, not a rounded 1e+08.
  expect_match(out, "crash +99999999")

  # A schedule overlaid under a non-penalized type still prints: a global
  # search charges the search-level penalties whatever the criterion.
  bic <- ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[rank]", 'type = "bic"',
    "[rank.penalties]", "crash = 1000.0"
  ))
  expect_match(paste(capture.output(print(bic)), collapse = "\n"), "Penalties")
})

test_that("an invalid penalty charge is still an error at load", {
  expect_error(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[rank.penalties]", "theta = -1.0"
    )),
    "non-negative"
  )
  expect_error(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[rank.penalties]", "not_a_penalty = 5.0"
    )),
    "not_a_penalty"
  )
})

test_that("a section no R tool runs warns at load rather than loading silently", {
  # #347: the engine's `TOOL_SECTIONS` admits `[globalsearch]` and
  # `[structsearch]`, neither of which has an R binding. The file loads, the
  # section is ignored by whichever tool gets the file, and the user gets a
  # stepwise search having asked for a global one. The load says so.
  path <- minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[globalsearch]", 'algorithm = "exhaustive"'
  )
  w <- expect_warning(ferx_search_config(path), "globalsearch")
  expect_s3_class(w, "ferx_search_unconsumed_section")
  expect_match(conditionMessage(w), "ferx_search_config")
  # `[globalsearch]` is the one section with somewhere else to run it.
  expect_match(conditionMessage(w), "ferx globalsearch", fixed = TRUE)

  cfg <- suppressWarnings(ferx_search_config(path))
  expect_equal(cfg$tools, "globalsearch")
  expect_match(paste(capture.output(print(cfg)), collapse = "\n"),
               "no R tool runs: globalsearch")

  # Both orphan sections at once, in one warning.
  both <- expect_warning(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[globalsearch]", 'algorithm = "exhaustive"',
      "[structsearch]", "dummy = 1"
    )),
    "structsearch"
  )
  expect_match(conditionMessage(both), "globalsearch")
  expect_match(conditionMessage(both), "have")
})

test_that("the remediation is per section, not a blanket pointer at the CLI", {
  # `[structsearch]` is accepted vocabulary and nothing more: the engine has no
  # structsearch module and the CLI has no structsearch command, so telling the
  # user to run it there would name a command that does not exist.
  w <- expect_warning(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[structsearch]", "dummy = 1"
    )),
    "structsearch"
  )
  expect_no_match(conditionMessage(w), "command-line tool")
  expect_no_match(conditionMessage(w), "ferx globalsearch", fixed = TRUE)
})

test_that("a section the engine adds later warns without a binding for it", {
  # The check is the complement of the sections this package has a tool for,
  # not a list of the two known-bad names - so a `TOOL_SECTIONS` entry added to
  # a future ferx-core is reported here from the day a file can carry it. The
  # tools vector is synthetic because the engine of the day refuses the name at
  # load; that is the point, and the behaviour under test is what happens once
  # it does not.
  w <- tryCatch(
    ferx:::.ferx_search_warn_unconsumed(c("covsearch", "sometoolfrom2027"),
                                        "ferx_search_config"),
    warning = function(w) w
  )
  expect_s3_class(w, "ferx_search_unconsumed_section")
  expect_match(conditionMessage(w), "sometoolfrom2027")
  # The section that does have a tool is not swept up with it.
  expect_no_match(conditionMessage(w), "covsearch")
  # And nothing is invented about where to run it instead.
  expect_no_match(conditionMessage(w), "command-line tool")

  # Every section this package binds stays silent.
  expect_silent(ferx:::.ferx_search_warn_unconsumed(
    c("allometry", "amd", "covsearch", "iivsearch", "iovsearch",
      "modelsearch", "ruvsearch"),
    "ferx_search_config"
  ))
  expect_silent(ferx:::.ferx_search_warn_unconsumed(character(0),
                                                    "ferx_search_config"))
})

test_that("a section every R tool does run is not warned about", {
  expect_silent(ferx_search_config(minimal_cfg(
    "COVARIATE?(CL, WT, pow)",
    "[covsearch]", "p_forward = 0.01"
  )))
})

test_that("an unrecognised rank type is refused at load, naming the offender", {
  expect_error(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[rank]",
      'type = "nonsense"'
    )),
    "nonsense"
  )
})

test_that("a missing file is a caller error", {
  expect_error(ferx_search_config(tempfile(fileext = ".ferxsearch")),
               "not found")
  expect_error(ferx_search_config(c("a", "b")), "single file path")
})

test_that("the bundled example configuration loads and prints", {
  ex <- ferx_example("two_cpt_oral_cov")
  expect_false(is.null(ex$search))
  expect_true(file.exists(ex$search))

  cfg <- ferx_search_config(ex$search)
  expect_equal(basename(cfg$base), "two_cpt_oral_cov.ferx")
  expect_true(file.exists(cfg$base))
  expect_true(file.exists(cfg$data))
  expect_equal(cfg$tools, "covsearch")
  expect_output(print(cfg), "ferx search configuration")
  expect_output(print(cfg), "COVARIATE")
})

test_that("an example without a search space has no `search` element", {
  # `$search` is per example and optional; two_cpt_iv ships no .ferxsearch.
  expect_null(ferx_example("two_cpt_iv")$search)
})

test_that("the warfarin example ships a structural search space", {
  path <- ferx_example("warfarin")$search
  expect_true(file.exists(path))
  cfg <- ferx_search_config(path)
  expect_equal(basename(cfg$base), "warfarin.ferx")
  expect_equal(cfg$tools, "modelsearch")
  expect_output(print(cfg), "PERIPHERALS")
})
