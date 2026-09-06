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

  # Defaults come from the engine, not from R: mixed BIC, no cutoff.
  expect_equal(cfg$rank$type, "bic")
  expect_true(is.na(cfg$rank$cutoff))
  expect_true(cfg$strictness$require_converged)
  expect_equal(cfg$strictness$max_condition_number, 1000)
  expect_length(cfg$strictness_set, 0L)
  expect_null(cfg$run$threads)
  expect_false(cfg$run$resume)
  expect_length(cfg$tools, 0L)
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
  err <- expect_error(
    ferx_search_config(minimal_cfg("ELIMINATION([FO, MM])")),
    "MM"
  )
  # The engine explains why, rather than silently narrowing the space.
  expect_match(conditionMessage(err), "ELIMINATION\\(MM\\)")
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

test_that("an unimplemented rank type is refused at load", {
  expect_error(
    ferx_search_config(minimal_cfg(
      "COVARIATE?(CL, WT, pow)",
      "[rank]",
      'type = "penalized"'
    )),
    "penalized"
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
  expect_null(ferx_example("warfarin")$search)
})
