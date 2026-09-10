# Stepwise covariate modelling from R.
#
# The trajectory itself is the engine's (ferx-core #1180 anchors it against PsN
# `scm`); what is asserted here is the R surface's own contract: the two entry
# forms cannot disagree, the object reports what the engine ran rather than
# anything recomputed in R, and a space covsearch cannot honour is an error
# before the first fit.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[space]", 'mfl = "COVARIATE?(CL, WT, pow)"')

  expect_error(
    ferx_covsearch(config = cfg, search_space = "COVARIATE?(CL, WT, pow)"),
    "search_space"
  )
  expect_error(ferx_covsearch(config = cfg, model = "model.ferx"), "model")
  # `data` names the dataset the search runs on, so the file states it too and
  # a second one here would be silently dropped.
  expect_error(ferx_covsearch(config = cfg, data = "other.csv"), "`data`")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that does
  # not exist - but not on the entry form.
  msg <- tryCatch(ferx_covsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_covsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_covsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(
    ferx_covsearch(model = tempfile(fileext = ".ferx"),
                   search_space = "COVARIATE?(CL, WT, pow)"),
    "model file not found"
  )
})

test_that("the inline form requires a search space", {
  ex <- ferx_example("two_cpt_oral_base")
  expect_error(ferx_covsearch(model = ex$model, data = ex$data),
               "`search_space` is required")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("two_cpt_oral_base")
  args <- list(model = ex$model, data = ex$data,
               search_space = "COVARIATE?(CL, WT, pow)")

  expect_error(do.call(ferx_covsearch, c(args, list(p_forward = -1))), "positive")
  expect_error(do.call(ferx_covsearch, c(args, list(p_forward = c(0.01, 0.05)))),
               "single number")
  # An infinite scalar cannot reach the config file - the binding emits a key
  # only when its value is finite - so it is refused rather than vanishing
  # into "keep the engine default".
  expect_error(do.call(ferx_covsearch, c(args, list(p_forward = Inf))), "finite")
  expect_error(do.call(ferx_covsearch, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_covsearch, c(args, list(max_steps = 2.5))), "whole number")
  expect_error(do.call(ferx_covsearch, c(args, list(retries = -1))), "at least 0")
  expect_error(do.call(ferx_covsearch, c(args, list(resume = NA))), "TRUE or FALSE")
  expect_error(do.call(ferx_covsearch, c(args, list(algorithm = "scm-backward"))),
               "should be one of")
  expect_error(do.call(ferx_covsearch, c(args, list(rank = 1))), "single string")
})

test_that("a structural feature in a covariate search is refused by name", {
  ex <- ferx_example("two_cpt_oral_base")
  # covsearch takes COVARIATE / COVARIATE? only; a structural statement is a
  # file meant for modelsearch, and running its covariate half silently would
  # answer a structural question with a covariate model.
  expect_error(
    ferx_covsearch(model = ex$model, data = ex$data,
                   search_space = "PERIPHERALS(0..1)"),
    "PERIPHERALS|not a covariate statement"
  )
})

test_that("a coverage gap is an error naming the feature, before any fit", {
  ex <- ferx_example("two_cpt_oral_base")
  expect_error(
    ferx_covsearch(model = ex$model, data = ex$data,
                   search_space = "ELIMINATION(MM)"),
    "ELIMINATION"
  )
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: slow, so it runs once and every assertion about the result
# object reads that one run.

covsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("two_cpt_oral_base")
    dir <- file.path(tempdir(), "covsearch-run")
    cached <<- ferx_covsearch(
      model        = ex$model,
      data         = ex$data,
      search_space = "COVARIATE?(CL, WT, pow)",
      max_steps    = 1,
      retries      = 0,
      directory    = dir,
      progress     = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- covsearch_run()

  expect_s3_class(res, "ferx_covsearch")
  expect_s3_class(res, "ferx_search_result")

  # The step table is the engine's own, column for column.
  expect_named(res$steps,
               c("step", "phase", "candidate", "parameter", "covariate", "form",
                 "parent_ofv", "ofv", "dofv", "df", "p_value", "alpha",
                 "significant", "selected", "converged", "passed", "failures"))
  expect_gt(nrow(res$steps), 0L)
  expect_true(is.logical(res$steps$passed))
  # Termination status is a column, never omitted - an init-stalled candidate
  # losing a step is a selection error, not a result.
  expect_true(is.logical(res$steps$converged))

  expect_true(is.finite(res$base_ofv))
  expect_true(is.finite(res$final_ofv))
  expect_true(file.exists(res$final_model_path))
})

test_that("the R step table equals the run's own steps.csv", {
  skip_on_cran()
  res <- covsearch_run()
  csv <- read.csv(file.path(res$directory, "steps.csv"), stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$steps))
  expect_equal(names(csv), names(res$steps))
  expect_equal(csv$candidate, res$steps$candidate)
  expect_equal(csv$ofv, res$steps$ofv, tolerance = 1e-6)
  # `read.csv` types the engine's booleans as text; the R object types them.
  expect_equal(as.logical(csv$selected), res$steps$selected)
  expect_equal(as.logical(csv$passed), res$steps$passed)
})

test_that("the winning fit is a ferx_fit, not a refit", {
  skip_on_cran()
  res <- covsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  expect_equal(unname(res$fit$ofv), res$final_ofv, tolerance = 1e-8)
  expect_true(!is.null(res$fit$theta))
  expect_true(is.data.frame(res$fit$sdtab))
})

test_that("the candidate table comes back with the run", {
  skip_on_cran()
  res <- covsearch_run()
  expect_true(is.data.frame(res$candidates))
  expect_true(all(c("run", "id", "criterion", "passed") %in% names(res$candidates)))
  # Each phase of the search runs its candidates in its own directory, and the
  # stacked table says which phase a row came from.
  expect_true("base" %in% res$candidates$run)
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  inline <- covsearch_run()

  ex <- ferx_example("two_cpt_oral_base")
  cfg <- write_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[space]",
    'mfl = "COVARIATE?(CL, WT, pow)"',
    "[run]",
    "retries = 0",
    "[covsearch]",
    "max_steps = 1"
  )
  from_file <- ferx_covsearch(config = cfg,
                              directory = file.path(tempdir(), "covsearch-file"),
                              progress = FALSE)

  expect_equal(from_file$steps$candidate, inline$steps$candidate)
  expect_equal(from_file$steps$ofv, inline$steps$ofv, tolerance = 1e-6)
  expect_equal(from_file$final_ofv, inline$final_ofv, tolerance = 1e-6)
  expect_equal(from_file$included, inline$included)
})

test_that("print() shows the step table and the final relations", {
  skip_on_cran()
  res <- covsearch_run()
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx covariate search")
  expect_match(out, "Final covariate relations")
})
