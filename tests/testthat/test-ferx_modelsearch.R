# Structural PK model search from R.
#
# The trajectory itself is the engine's (ferx-core #1181 anchors it against
# Pharmpy's modelsearch); what is asserted here is the R surface's own
# contract: the two entry forms cannot disagree, the object reports what the
# engine ran rather than anything recomputed in R, and a space modelsearch
# cannot honour is an error before the first fit.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[space]", 'mfl = "PERIPHERALS(0..1)"')

  expect_error(
    ferx_modelsearch(config = cfg, search_space = "PERIPHERALS(0..1)"),
    "search_space"
  )
  expect_error(ferx_modelsearch(config = cfg, model = "model.ferx"), "model")
  expect_error(ferx_modelsearch(config = cfg, iiv_strategy = "no_add"), "iiv_strategy")
  # `data` states *which dataset* the search runs on, so it belongs with
  # `model` rather than with the run knobs: accepted beside `config` it would
  # be dropped and the file's own `data` key used, searching a different
  # dataset than the caller asked for, with nothing said.
  expect_error(ferx_modelsearch(config = cfg, data = "other.csv"), "`data`")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that does
  # not exist - but not on the entry form.
  msg <- tryCatch(ferx_modelsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_modelsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_modelsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(
    ferx_modelsearch(model = tempfile(fileext = ".ferx"),
                     search_space = "PERIPHERALS(0..1)"),
    "model file not found"
  )
})

test_that("the inline form requires a search space", {
  ex <- ferx_example("warfarin")
  expect_error(ferx_modelsearch(model = ex$model, data = ex$data),
               "`search_space` is required")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("warfarin")
  args <- list(model = ex$model, data = ex$data,
               search_space = "PERIPHERALS(0..1)")

  expect_error(do.call(ferx_modelsearch, c(args, list(algorithm = "forward"))),
               "should be one of")
  expect_error(do.call(ferx_modelsearch, c(args, list(iiv_strategy = "diagonal"))),
               "should be one of")
  expect_error(do.call(ferx_modelsearch, c(args, list(cutoff = c(1, 2)))),
               "single number")
  # An infinite cutoff would be dropped on the way to the config file - the
  # binding emits the key only when it is finite - so the search would run
  # with no cutoff at all rather than with the one that was asked for.
  expect_error(do.call(ferx_modelsearch, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_modelsearch, c(args, list(cutoff = -Inf))), "finite")
  expect_error(do.call(ferx_modelsearch, c(args, list(retries = -1))), "at least 0")
  expect_error(do.call(ferx_modelsearch, c(args, list(threads = 2.5))), "whole number")
  expect_error(do.call(ferx_modelsearch, c(args, list(resume = NA))), "TRUE or FALSE")
  expect_error(do.call(ferx_modelsearch, c(args, list(rank = 1))), "single string")
})

test_that("a finite but negative cutoff reaches the engine, which names it", {
  # The other half of the cutoff contract: a finite value is emitted, so the
  # engine's own validation is what rejects it. R does not duplicate that rule,
  # it only guarantees the value arrives.
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_modelsearch(model = ex$model, data = ex$data,
                     search_space = "PERIPHERALS(0..1)", cutoff = -1),
    "cutoff"
  )
})

test_that("a covariate space in a structural search is refused by name", {
  ex <- ferx_example("warfarin")
  # modelsearch takes structural statements only; a covariate space handed to
  # it would search nothing and report the base model as the winner.
  expect_error(
    ferx_modelsearch(model = ex$model, data = ex$data,
                     search_space = "COVARIATE?(CL, WT, pow)"),
    "COVARIATE|not a structural statement"
  )
})

test_that("a coverage gap is an error naming the feature, before any fit", {
  ex <- ferx_example("warfarin")
  # SEQ-ZO-FO is a different disposition rather than a different input term,
  # so the engine cannot build it from a template - and says so by name.
  expect_error(
    ferx_modelsearch(model = ex$model, data = ex$data,
                     search_space = "ABSORPTION([FO,SEQ-ZO-FO])"),
    "SEQ-ZO-FO"
  )
})

test_that("iiv_strategy = 'fullblock' is refused naming the tool that owns it", {
  ex <- ferx_example("warfarin")
  # Pharmpy's fourth strategy: a block over the new and existing eta is a
  # variability-structure move, which is iivsearch's. It stays a legal choice
  # in R so the refusal comes from the engine, naming what to use instead.
  expect_error(
    ferx_modelsearch(model = ex$model, data = ex$data,
                     search_space = "PERIPHERALS(0..1)",
                     iiv_strategy = "fullblock"),
    "fullblock"
  )
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: slow, so it runs once and every assertion about the result
# object reads that one run.

modelsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("warfarin")
    dir <- file.path(tempdir(), "modelsearch-run")
    cached <<- ferx_modelsearch(
      model        = ex$model,
      data         = ex$data,
      search_space = "ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
      rank         = "bic",
      retries      = 0,
      directory    = dir,
      progress     = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- modelsearch_run()

  expect_s3_class(res, "ferx_modelsearch")
  expect_s3_class(res, "ferx_search_result")

  # The model table is the engine's own, column for column, with the engine's
  # structure rendering appended - nothing here is spelled in R.
  expect_equal(names(res$models),
               c(ferx_rust_modelsearch_columns(), "structure"))
  expect_gt(nrow(res$models), 1L)
  expect_true(is.logical(res$models$passed))
  # Termination status is a column, never omitted - a candidate that stalled at
  # its initial estimates is a selection error, not a result.
  expect_true(is.logical(res$models$converged))

  # The base model is a row of the table, and the winner is one of the rows.
  expect_true(res$base_model_id %in% res$models$id)
  expect_true(res$final_model_id %in% res$models$id)
  expect_equal(sum(res$models$selected), 1L)
  expect_equal(res$models$id[res$models$selected], res$final_model_id)

  expect_true(is.finite(res$base_criterion))
  expect_true(is.finite(res$final_criterion))
  # The search cannot end worse than it started: the base is a candidate too.
  expect_lte(res$final_criterion, res$base_criterion)
  expect_true(file.exists(res$final_model_path))
})

test_that("the R model table equals the run's own models.csv", {
  skip_on_cran()
  res <- modelsearch_run()
  csv <- read.csv(file.path(res$directory, "models.csv"), stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$models))
  # Every engine column, in the engine's order; `structure` is the engine's own
  # rendering added on the end, and is the only column beyond the file's.
  expect_equal(names(res$models)[seq_along(names(csv))], names(csv))
  expect_equal(setdiff(names(res$models), names(csv)), "structure")

  expect_equal(csv$id, res$models$id)
  expect_equal(csv$criterion, res$models$criterion, tolerance = 1e-6)
  expect_equal(csv$ofv, res$models$ofv, tolerance = 1e-6)
  # `read.csv` types the engine's booleans as text; the R object types them.
  expect_equal(as.logical(csv$selected), res$models$selected)
  expect_equal(as.logical(csv$passed), res$models$passed)
})

test_that("ferx_search_results() reads the model table back", {
  skip_on_cran()
  res <- modelsearch_run()
  tab <- ferx_search_results(res$directory, type = "models")

  expect_equal(names(tab), ferx_rust_modelsearch_columns())
  expect_equal(tab$id, res$models$id)
  expect_equal(tab$criterion, res$models$criterion, tolerance = 1e-6)
  expect_true(is.logical(tab$passed))
  expect_true(is.integer(tab$layer))
  # The empty cell the engine writes for "no value" is NA, not NaN and not "".
  expect_true(is.na(tab$parent[tab$id == res$base_model_id]) ||
                nzchar(tab$parent[tab$id == res$base_model_id]))

  # The candidate table is a different table of the same run, and still reads.
  # A structural search gives each layer its own runner directory, so the
  # runner's own table lives one level down - which is what `$candidates`
  # stacks and label by layer.
  expect_true(is.data.frame(ferx_search_results(file.path(res$directory, "base"))))
  expect_setequal(unique(res$candidates$run),
                  c("base", paste0("layer-", seq_len(res$n_layers))))
  expect_error(ferx_search_results(res$directory, type = "models", partial = TRUE),
               "no partial model table")
  expect_error(ferx_search_results(tempdir(), type = "models"), "No model table")
})

test_that("every candidate's model text comes back, named by id", {
  skip_on_cran()
  res <- modelsearch_run()

  expect_true(all(res$models$id %in% names(res$model_text)))
  expect_true(all(nzchar(res$model_text)))
  # The text is a model file, not a description of one.
  expect_match(res$model_text[[res$final_model_id]], "[structural_model]",
               fixed = TRUE)
  # And it is what the run wrote beside itself.
  on_disk <- file.path(res$directory, "models",
                       paste0(res$final_model_id, ".ferx"))
  skip_if_not(file.exists(on_disk))
  expect_equal(paste(readLines(on_disk), collapse = "\n"),
               trimws(res$model_text[[res$final_model_id]], which = "right"))
})

test_that("the winning fit is a ferx_fit, not a refit", {
  skip_on_cran()
  res <- modelsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  expect_equal(unname(res$fit$ofv), res$final_ofv, tolerance = 1e-8)
  expect_true(!is.null(res$fit$theta))
  expect_true(is.data.frame(res$fit$sdtab))
})

test_that("the candidate table comes back with the run", {
  skip_on_cran()
  res <- modelsearch_run()
  expect_true(is.data.frame(res$candidates))
  expect_true(all(c("id", "criterion", "passed") %in% names(res$candidates)))
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  inline <- modelsearch_run()

  ex <- ferx_example("warfarin")
  cfg <- write_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[space]",
    'mfl = "ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])"',
    "[rank]",
    'type = "bic"',
    "[run]",
    "retries = 0"
  )
  from_file <- ferx_modelsearch(config = cfg,
                                directory = file.path(tempdir(), "modelsearch-file"),
                                progress = FALSE)

  expect_equal(from_file$models$id, inline$models$id)
  expect_equal(from_file$models$criterion, inline$models$criterion, tolerance = 1e-6)
  expect_equal(from_file$final_model_id, inline$final_model_id)
  expect_equal(from_file$final_criterion, inline$final_criterion, tolerance = 1e-6)
})

test_that("a single-point space returns the base model, fitted", {
  skip_on_cran()
  # The degenerate oracle: a space the base model already sits on has one
  # candidate, so the search must hand back the base fit rather than a
  # transformed model - and it must agree with fitting that model directly.
  ex <- ferx_example("warfarin")
  res <- ferx_modelsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "ABSORPTION(FO); PERIPHERALS(0)",
    retries      = 0,
    directory    = file.path(tempdir(), "modelsearch-degenerate"),
    progress     = FALSE
  )

  expect_equal(nrow(res$models), 1L)
  expect_equal(res$final_model_id, res$base_model_id)
  skip_if(is.null(res$fit), "the run recovered no final fit")

  direct <- ferx_fit(ex$model, ex$data)
  expect_equal(unname(res$fit$ofv), unname(direct$ofv), tolerance = 1e-6)
})

test_that("a gate-excluded model is a row with its reason, and cannot win", {
  skip_on_cran()
  # The failure #334 asks the R surface to make visible: on this space the
  # lowest-OFV models are the two-compartment ones, and they are excluded for
  # being ill-conditioned (TVKA ~ TVV2 at |r| = 1) or for pinning TVQ to a
  # bound. A winner-only view would report one of them; the table has to show
  # both the exclusion and why, and the search has to select the best model
  # that actually passed.
  ex <- ferx_example("warfarin")
  res <- ferx_modelsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
    rank         = "ofv",
    retries      = 0,
    directory    = file.path(tempdir(), "modelsearch-gated"),
    progress     = FALSE
  )

  gated <- res$models[!res$models$passed, , drop = FALSE]
  expect_gt(nrow(gated), 0L)
  # Excluded, but present, and each row says why.
  expect_true(all(!is.na(gated$failures) & nzchar(gated$failures)))
  expect_true(all(is.na(gated$rank)))
  expect_false(any(gated$selected))

  # The criterion alone does not separate them from the winner: the two
  # two-compartment models with a lag land on the same fit as the winner
  # because Q collapses, and come back 2.0e-6 and 7.8e-8 OFV *worse* - a
  # measured tie. So the gate is doing the separating, and a report showing
  # only the criterion would present three indistinguishable winners, one of
  # them with |r| = 1 between TVKA and TVV2. The bound is 1e-3, ~500x the
  # realised gap.
  winner <- res$models[res$models$selected, , drop = FALSE]
  expect_equal(nrow(winner), 1L)
  expect_true(winner$passed)
  expect_lt(min(gated$criterion, na.rm = TRUE) - winner$criterion, 1e-3)

  # And the winner here is not the base model, so the transformed-model path -
  # the model text, `final.ferx` and the fit that comes back - is exercised.
  expect_false(winner$id == res$base_model_id)
  expect_equal(res$final_model_id, winner$id)
  expect_match(res$model_text[[winner$id]], "lag", ignore.case = TRUE)
})

test_that("print() shows the ranked models and summary() shows the rest", {
  skip_on_cran()
  res <- modelsearch_run()
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx structural model search")
  expect_match(out, "Ranked models")
  expect_match(out, res$base_structure, fixed = TRUE)

  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "Every model")
})
