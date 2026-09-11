# The AMD pipeline from R.
#
# The pipeline itself is the engine's (ferx-core #1184 anchors its steps
# against Pharmpy's `amd`); what is asserted here is the R surface's own
# contract: the two entry forms cannot disagree, the plan is the engine's own
# rather than a second derivation of it, the step order is the strategy's, the
# tables are the engine's column for column, and the strictness verdict and the
# termination status are columns at both levels.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[space]",
                   'mfl = "IIV?(@PK, exp)"')

  expect_error(ferx_amd(config = cfg, search_space = "IIV?(CL, exp)"),
               "search_space")
  expect_error(ferx_amd(config = cfg, strategy = "SIR"), "strategy")
  expect_error(ferx_amd(config = cfg, retries_on = "final"), "retries_on")
  expect_error(ferx_amd(config = cfg, skip = "residual"), "skip")
  expect_error(ferx_amd(config = cfg, model = "model.ferx"), "model")
  # `data` states *which dataset* the pipeline runs on, so it belongs with
  # `model` rather than with the run knobs.
  expect_error(ferx_amd(config = cfg, data = "other.csv"), "`data`")
  expect_error(ferx_amd_plan(config = cfg, strategy = "RSI"), "strategy")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that does
  # not exist - but not on the entry form.
  msg <- tryCatch(ferx_amd(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_amd(), "`config`.*or `model`")
  expect_error(ferx_amd_plan(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_amd(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(ferx_amd(model = tempfile(fileext = ".ferx"),
                        search_space = "IIV?(@PK, exp)"),
               "model file not found")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("warfarin_amd")
  args <- list(model = ex$model, data = ex$data,
               search_space = "IIV?(@PK, exp)")

  expect_error(do.call(ferx_amd, c(args, list(strategy = "greedy"))),
               "should be one of")
  expect_error(do.call(ferx_amd, c(args, list(retries_on = "sometimes"))),
               "should be one of")
  expect_error(do.call(ferx_amd, c(args, list(skip = "structrual"))),
               "does not name a pipeline step")
  expect_error(do.call(ferx_amd, c(args, list(skip = NA_character_))),
               "character vector")
  expect_error(do.call(ferx_amd, c(args, list(rank = c("bic", "ofv")))),
               "single string")
  expect_error(do.call(ferx_amd, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_amd, c(args, list(retries = -1))), "at least 0")
  expect_error(do.call(ferx_amd, c(args, list(resume = NA))), "TRUE or FALSE")
  # Resuming needs the directory the fits were journalled in; without one the
  # pipeline would run in a temporary directory and reuse nothing.
  expect_error(do.call(ferx_amd, c(args, list(resume = TRUE))),
               "needs the `directory`")
  # A pipeline with no space is a residual-error search with five skipped
  # steps, which is ferx_ruvsearch() under another name.
  expect_error(ferx_amd(model = ex$model, data = ex$data),
               "`search_space` is required")
})

test_that("the engine's own validation runs before the dataset is read", {
  # Both name a file key rather than an R argument, and both are refused while
  # the config is being read - the base model named below does not exist, so
  # anything that got as far as loading it would fail differently.
  #
  # A `skip` naming a step the chosen strategy does not run.
  expect_error(ferx_amd(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "IIV?(@PK, exp)"',
    "[amd]", 'strategy = "SIR"', 'skip = ["covariates"]'
  )), "not a step of the SIR strategy")
  # A space carrying a statement no step of this pipeline can route. AMD's own
  # refusal ("a PD / metabolite statement, which the AMD pipeline has no step
  # for") sits behind the loader's coverage check today, since every PD
  # statement is a feature the engine builds no candidate for at all - so what
  # this asserts is that such a space is refused by name, not which layer says
  # so.
  expect_error(ferx_amd(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "DIRECTEFFECT(LINEAR)"'
  )), "DIRECTEFFECT")
  # The same file with nothing wrong in it gets as far as the base model and
  # fails on that instead, which is what says the refusals above were not this
  # one.
  expect_error(ferx_amd(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "IIV?(@PK, exp)"'
  )), "read model file")
})

# -- The plan ----------------------------------------------------------------

test_that("the plan is the pipeline the engine would run", {
  ex <- ferx_example("warfarin_amd")
  plan <- ferx_amd_plan(config = ex$search)

  expect_equal(names(plan),
               c("index", "step", "tool", "rerun", "directory", "skipped"))
  expect_equal(plan$index, seq_len(nrow(plan)))
  # The default order, and the tools that run it - which is not the step's own
  # name for three of the six.
  expect_equal(plan$step,
               c("structural", "iivsearch", "residual", "iovsearch",
                 "allometry", "covariates"))
  expect_equal(plan$tool,
               c("modelsearch", "iivsearch", "ruvsearch", "iovsearch",
                 "allometry", "covsearch"))
  expect_false(any(plan$rerun))
  expect_equal(plan$directory[1L], "01-modelsearch")

  # A skipped step says why, in the engine's own words. The bundled example
  # has no occasion column and no covariates, so three of the six are skipped.
  skipped <- plan[!is.na(plan$skipped), , drop = FALSE]
  expect_equal(skipped$step, c("iovsearch", "allometry", "covariates"))
  expect_match(skipped$skipped[1L], "iov_column")
  expect_match(skipped$skipped[2L], "ALLOMETRY")
  expect_match(skipped$skipped[3L], "COVARIATE")

  opts <- attr(plan, "options")
  expect_equal(opts$strategy, "default")
  expect_equal(opts$retries_on, "all_final")
  expect_true(is.na(opts$iov_column))
})

test_that("each strategy plans its own order, reruns marked", {
  ex <- ferx_example("warfarin_amd")
  space <- "ABSORPTION(FO); PERIPHERALS(0..1); IIV?(@PK, exp)"
  order_of <- function(strategy) {
    ferx_amd_plan(model = ex$model, data = ex$data, search_space = space,
                  strategy = strategy)
  }

  expect_equal(order_of("SIR")$step, c("structural", "iivsearch", "residual"))
  expect_equal(order_of("SRI")$step, c("structural", "residual", "iivsearch"))
  expect_equal(order_of("RSI")$step, c("residual", "structural", "iivsearch"))

  # `reevaluation` is the default order with IIV and residual run a second
  # time; the second occurrence is a rerun and has its own directory.
  re <- order_of("reevaluation")
  expect_equal(re$step,
               c("structural", "iivsearch", "residual", "iovsearch",
                 "allometry", "covariates", "iivsearch", "residual"))
  expect_equal(re$rerun, c(rep(FALSE, 6L), TRUE, TRUE))
  expect_equal(re$directory[7L], "07-rerun-iivsearch")
})

test_that("skip leaves a step out, and says that is why", {
  ex <- ferx_example("warfarin_amd")
  plan <- ferx_amd_plan(model = ex$model, data = ex$data,
                        search_space = "ABSORPTION(FO); IIV?(@PK, exp)",
                        strategy = "SIR", skip = "residual")
  expect_equal(plan$step, c("structural", "iivsearch", "residual"))
  expect_true(is.na(plan$skipped[2L]))
  expect_match(plan$skipped[3L], "skip")
})

test_that("the config file and the inline arguments plan the same pipeline", {
  ex <- ferx_example("warfarin_amd")
  # The space the bundled file states, spelled inline.
  space <- c("ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
             "IIV?(@PK, exp); COVARIANCE?(IIV, *)")

  from_file <- ferx_amd_plan(config = ex$search)
  inline <- ferx_amd_plan(model = ex$model, data = ex$data,
                          search_space = space, strategy = "default",
                          retries_on = "all_final")

  expect_equal(inline[, c("index", "step", "tool", "rerun", "directory",
                          "skipped")],
               from_file[, c("index", "step", "tool", "rerun", "directory",
                             "skipped")])
  expect_equal(attr(inline, "options")$strategy,
               attr(from_file, "options")$strategy)
  expect_equal(attr(inline, "options")$retries_on,
               attr(from_file, "options")$retries_on)
})

# -- The step verdict --------------------------------------------------------

test_that("a step's verdict is the verdict on the model it selected", {
  candidates <- data.frame(
    step      = c(1L, 1L, 1L, 2L, 2L, 3L),
    tool      = c("modelsearch", "modelsearch", "retries",
                  "iivsearch", "iivsearch", "ruvsearch"),
    id        = c("base", "run1", "retry1", "base", "run1", "input"),
    converged = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
    passed    = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    selected  = c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  verdict <- ferx:::.ferx_amd_step_verdict(1:3, candidates)

  # Step 1 selected a candidate the gate rejected and then improved on it in
  # the retries pass, which is filed under the same step and selected in its
  # turn: the later row is the model the step handed on.
  expect_equal(verdict$converged, c(TRUE, TRUE, NA))
  expect_equal(verdict$passed, c(TRUE, FALSE, NA))

  # A step whose winner fails the gate is visible as such rather than reported
  # by its rank alone.
  one_step <- ferx:::.ferx_amd_step_verdict(
    2L, candidates[candidates$step == 2L, , drop = FALSE])
  expect_false(one_step$passed)
})

# -- The end-to-end run ------------------------------------------------------
#
# One real pipeline: slow, so it runs once and every assertion about the result
# object reads that one run. `SIR` is the three-step strategy - structural,
# IIV, residual - which is every kind of step the pipeline has (a ranked
# search, a two-stage ranked search, and one that selects on a likelihood-ratio
# test) without the three the bundled example skips.

amd_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("warfarin_amd")
    cached <<- ferx_amd(
      model        = ex$model,
      data         = ex$data,
      search_space = "ABSORPTION(FO); PERIPHERALS(0..1); IIV?(@PK, exp)",
      strategy     = "SIR",
      retries_on   = "final",
      retries      = 0,
      directory    = file.path(tempdir(), "amd-run"),
      progress     = FALSE
    )
    cached
  }
})

test_that("the result reports the pipeline the engine ran", {
  skip_on_cran()
  res <- amd_run()

  expect_s3_class(res, "ferx_amd")
  expect_s3_class(res, "ferx_search_result")

  # Both tables are the engine's own, column for column: the step table with
  # the pipeline position in front and the joined verdict behind, the candidate
  # table exactly as `CANDIDATE_COLUMNS` orders it.
  expect_equal(names(res$steps),
               c("index", ferx_rust_amd_step_columns(),
                 "converged", "passed", "notes"))
  expect_equal(names(res$candidates), ferx_rust_amd_candidate_columns())

  # The order is the strategy's, asserted on the run rather than assumed.
  expect_equal(res$steps$step, c("structural", "iivsearch", "residual"))
  expect_equal(res$steps$tool, c("modelsearch", "iivsearch", "ruvsearch"))
  expect_equal(res$steps$index, 1:3)
  expect_true(all(res$steps$status == "ran"))

  # The options come back as the engine read them, not as R passed them.
  expect_equal(res$options$strategy, "SIR")
  expect_equal(res$options$retries_on, "final")
  expect_equal(res$options$skip, character(0))

  expect_true(is.finite(res$input_ofv))
  expect_true(is.finite(res$final_ofv))
  expect_equal(res$d_ofv, res$final_ofv - res$input_ofv)
  expect_s3_class(res$fit, "ferx_fit")
  expect_true(file.exists(res$final_model_path))
  expect_false(res$cancelled)
})

test_that("the run reproduces the engine's own report, rather than recomputing it", {
  skip_on_cran()
  res <- amd_run()

  # The two tables the engine wrote beside the run are the two tables on the
  # object: R reads the same numbers out of the same result, it does not build
  # a second version of them.
  steps_csv <- utils::read.csv(res$steps_csv, stringsAsFactors = FALSE)
  expect_equal(steps_csv$step, res$steps$step)
  expect_equal(steps_csv$status, res$steps$status)
  expect_equal(steps_csv$selected, res$steps$selected)
  expect_equal(steps_csv$value_after, res$steps$value_after, tolerance = 1e-6)
  # `d_value` is the engine's own subtraction; R does the same one rather than
  # a different one.
  expect_equal(steps_csv$d_value, res$steps$d_value, tolerance = 1e-6)

  cand_csv <- utils::read.csv(res$candidates_csv, stringsAsFactors = FALSE)
  expect_equal(nrow(cand_csv), nrow(res$candidates))
  expect_equal(cand_csv$id, res$candidates$id)
  expect_equal(cand_csv$value, res$candidates$value, tolerance = 1e-6)
  # The CSV spells a logical the engine's way ("true" / "false"), which is the
  # same verdict the object carries as a logical.
  expect_equal(ferx:::.ferx_search_lgl(cand_csv$passed), res$candidates$passed)
})

test_that("the strictness verdict and the termination status are columns at both levels", {
  skip_on_cran()
  res <- amd_run()

  expect_true(is.logical(res$candidates$converged))
  expect_true(is.logical(res$candidates$passed))
  expect_true(is.logical(res$steps$converged))
  expect_true(is.logical(res$steps$passed))

  # A candidate the gate rejected keeps the gate's reasons, in the gate's own
  # words - a winner-only table is what this pipeline must not be.
  rejected <- res$candidates[!res$candidates$passed, , drop = FALSE]
  expect_gt(nrow(rejected), 0L)
  expect_true(all(nzchar(rejected$failures) | !is.na(rejected$error)))

  # Each step's verdict is the verdict on the model that step selected.
  for (i in which(!is.na(res$steps$passed))) {
    rows <- res$candidates[res$candidates$step == res$steps$index[i] &
                             res$candidates$selected, , drop = FALSE]
    expect_gt(nrow(rows), 0L)
    expect_equal(res$steps$passed[i], rows$passed[nrow(rows)])
    expect_equal(res$steps$converged[i], rows$converged[nrow(rows)])
  }
})

test_that("every step's candidates are under it, the pipeline's own fits beside them", {
  skip_on_cran()
  res <- amd_run()

  expect_equal(names(res$tools),
               c("01-structural", "02-iivsearch", "03-residual"))
  for (tool in res$tools) {
    expect_equal(unique(tool$candidates$step), tool$index)
    expect_gt(nrow(tool$candidates), 0L)
    expect_true(dir.exists(file.path(res$directory, tool$directory)))
  }

  # The pipeline's own start fit is a row of the candidate table, filed at step
  # 0 - it belongs to no step, and it is the model every d_ofv is measured from.
  start <- res$candidates[res$candidates$tool == "start", , drop = FALSE]
  expect_equal(nrow(start), 1L)
  expect_equal(start$step, 0L)
  expect_equal(start$ofv, res$input_ofv, tolerance = 1e-6)

  # Each step's candidates were fitted by that step's tool.
  for (i in seq_len(nrow(res$steps))) {
    rows <- res$candidates[res$candidates$step == res$steps$index[i], ,
                           drop = FALSE]
    expect_true(all(rows$tool %in% c(res$steps$tool[i], "retries")))
  }
})

test_that("a run with no directory keeps nothing but its tables", {
  skip_on_cran()
  ex <- ferx_example("warfarin_amd")
  res <- ferx_amd(
    model        = ex$model,
    data         = ex$data,
    search_space = "IIV?(@PK, exp)",
    strategy     = "SIR",
    skip         = c("structural", "residual"),
    retries_on   = "skip",
    retries      = 0,
    progress     = FALSE
  )

  expect_true(is.na(res$directory))
  expect_true(is.na(res$steps_csv))
  expect_true(is.na(res$candidates_csv))
  # The tables are on the object either way, and the final model still has a
  # file to be read from.
  expect_equal(nrow(res$steps), 3L)
  expect_gt(nrow(res$candidates), 0L)
  expect_true(file.exists(res$final_model_path))
  # The two skipped steps say they were skipped by name, not by the space.
  skipped <- res$steps[res$steps$status == "skipped", , drop = FALSE]
  expect_equal(skipped$step, c("structural", "residual"))
  expect_true(all(grepl("skip", skipped$reason)))
})

test_that("print and summary show the step table and the candidates under it", {
  skip_on_cran()
  res <- amd_run()

  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(out, "SIR strategy")
  expect_match(out, "Pipeline:")
  expect_match(out, "structural")
  expect_match(out, "passed")
  expect_match(out, "strictness gate")

  out2 <- paste(utils::capture.output(summary(res)), collapse = "\n")
  expect_match(out2, "Step 1 - structural")
  expect_match(out2, "Step 2 - iivsearch")
  expect_match(out2, "excluded, and why")
})
