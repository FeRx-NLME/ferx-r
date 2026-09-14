# Global model search from R.
#
# The trajectory itself is the engine's (ferx-core #1185 anchors the genetic
# algorithm against exhaustive enumeration on a known landscape); what is
# asserted here is the R surface's own contract: the two entry forms cannot
# disagree, the object reports what the engine ran rather than anything
# recomputed in R, the per-candidate charges add up to the fitness the search
# ranked on, and a space globalsearch cannot lay out as a grid is an error
# before the first fit.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[space]", 'mfl = "PERIPHERALS(0..1)"')

  expect_error(
    ferx_globalsearch(config = cfg, search_space = "PERIPHERALS(0..1)"),
    "search_space"
  )
  expect_error(ferx_globalsearch(config = cfg, model = "model.ferx"), "model")
  expect_error(ferx_globalsearch(config = cfg, algorithm = "ga"), "algorithm")
  expect_error(ferx_globalsearch(config = cfg, ga = list(seed = 1)), "ga")
  expect_error(ferx_globalsearch(config = cfg, penalties = list(gate = 1)),
               "penalties")
  # `data` states *which dataset* the search runs on, so it belongs with
  # `model` rather than with the run knobs.
  expect_error(ferx_globalsearch(config = cfg, data = "other.csv"), "`data`")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that
  # does not exist - but not on the entry form.
  msg <- tryCatch(ferx_globalsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_globalsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_globalsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(
    ferx_globalsearch(model = tempfile(fileext = ".ferx"),
                      search_space = "PERIPHERALS(0..1)"),
    "model file not found"
  )
})

test_that("the inline form requires a search space", {
  ex <- ferx_example("two_cpt_oral_global")
  expect_error(ferx_globalsearch(model = ex$model, data = ex$data),
               "`search_space` is required")
})

test_that("an unknown algorithm is refused by name, not partially matched", {
  # #364 asks for this explicitly: `match.arg()` would accept `"e"` and pick
  # `exhaustive`, which on a grid nobody sized is not a typo worth guessing
  # at. Both the unknown value and the valid ones are in the message.
  ex <- ferx_example("two_cpt_oral_global")
  args <- list(model = ex$model, data = ex$data,
               search_space = "PERIPHERALS(0..1)")

  err <- tryCatch(do.call(ferx_globalsearch, c(args, list(algorithm = "genetic"))),
                  error = conditionMessage)
  expect_match(err, "genetic", fixed = TRUE)
  expect_match(err, '"ga"', fixed = TRUE)
  expect_match(err, '"exhaustive"', fixed = TRUE)

  expect_error(do.call(ferx_globalsearch, c(args, list(algorithm = "e"))), "\"e\"")
  expect_error(do.call(ferx_globalsearch, c(args, list(algorithm = "ex"))), "\"ex\"")
  expect_error(do.call(ferx_globalsearch, c(args, list(algorithm = 1))),
               "single string")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("two_cpt_oral_global")
  args <- list(model = ex$model, data = ex$data,
               search_space = "PERIPHERALS(0..1)")

  expect_error(do.call(ferx_globalsearch, c(args, list(iiv_strategy = "diagonal"))),
               "iiv_strategy")
  expect_error(do.call(ferx_globalsearch, c(args, list(cutoff = c(1, 2)))),
               "single number")
  # An infinite cutoff would be dropped on the way to the config file - the
  # binding emits the key only when it is finite - so the search would run
  # with no cutoff at all rather than with the one that was asked for.
  expect_error(do.call(ferx_globalsearch, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_globalsearch, c(args, list(max_models = 0))),
               "at least 1")
  expect_error(do.call(ferx_globalsearch, c(args, list(max_models = 2.5))),
               "whole number")
  expect_error(do.call(ferx_globalsearch, c(args, list(retries = -1))), "at least 0")
  expect_error(do.call(ferx_globalsearch, c(args, list(threads = 2.5))), "whole number")
  expect_error(do.call(ferx_globalsearch, c(args, list(resume = NA))), "TRUE or FALSE")
  expect_error(do.call(ferx_globalsearch, c(args, list(rank = 1))), "single string")
})

test_that("resume with nothing to resume from is refused before any fit", {
  # A `resume = TRUE` the engine cannot honour reaches it as an ordinary run:
  # the flag is set, there is no journal, and the search silently refits
  # everything. On a tool whose whole cost is fits that is an expensive way to
  # learn the argument did nothing, so both halves are refused up front
  # (PR #375 review).
  ex <- ferx_example("two_cpt_oral_global")
  args <- list(model = ex$model, data = ex$data,
               search_space = "PERIPHERALS(0..1)")

  # No directory at all: nowhere a journal could be.
  expect_error(do.call(ferx_globalsearch, c(args, list(resume = TRUE))),
               "needs the `directory`", fixed = TRUE)
  # A directory no earlier run wrote.
  expect_error(
    do.call(ferx_globalsearch,
            c(args, list(resume = TRUE,
                         directory = file.path(tempdir(), "no-such-run")))),
    "there is nothing at"
  )
  # `resume = NA` is still answered by the type check, not by this one.
  expect_error(do.call(ferx_globalsearch, c(args, list(resume = NA))),
               "TRUE or FALSE")
})

test_that("the resume guard covers the whole search family, not just this tool", {
  # The rule lives in `.ferx_search_directory()`, which every tool routes its
  # `directory` through, so a sibling cannot quietly keep the old behaviour.
  # `ferx_allometry()` is absent on purpose: it fits a base and a scaled arm
  # and journals nothing, so it is the one search tool with no `resume`.
  for (what in c("ferx_modelsearch", "ferx_covsearch", "ferx_iivsearch",
                 "ferx_iovsearch", "ferx_ruvsearch", "ferx_amd",
                 "ferx_globalsearch")) {
    expect_error(ferx:::.ferx_search_directory(NULL, what, resume = TRUE),
                 "needs the `directory`", fixed = TRUE)
  }
  expect_false("resume" %in% names(formals(ferx_allometry)))
  # Without `resume` an absent directory is still the in-memory run.
  expect_equal(ferx:::.ferx_search_directory(NULL, "ferx_globalsearch"), "")
  expect_equal(
    ferx:::.ferx_search_directory(NULL, "ferx_globalsearch", resume = FALSE), ""
  )
  # And an existing directory resumes fine.
  d <- file.path(tempdir(), "resume-guard-ok")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  expect_equal(ferx:::.ferx_search_directory(d, "ferx_globalsearch", resume = TRUE),
               normalizePath(d))
})

test_that("an unknown ga or penalty knob is refused by name, from the engine's list", {
  ex <- ferx_example("two_cpt_oral_global")
  args <- list(model = ex$model, data = ex$data,
               search_space = "PERIPHERALS(0..1)")

  err <- tryCatch(
    do.call(ferx_globalsearch, c(args, list(ga = list(populaton_size = 8)))),
    error = conditionMessage
  )
  expect_match(err, "populaton_size", fixed = TRUE)
  # The valid names come from the engine, not from a list restated in R.
  expect_match(err, "population_size", fixed = TRUE)

  expect_error(
    do.call(ferx_globalsearch, c(args, list(penalties = list(gates = 1)))),
    "gates"
  )
  expect_error(
    do.call(ferx_globalsearch, c(args, list(ga = list(8)))),
    "must be named"
  )
  expect_error(
    do.call(ferx_globalsearch, c(args, list(ga = list(seed = 1, seed = 2)))),
    "twice"
  )
  # `population_size` is a count and `final_downhill` a flag; the kinds come
  # from the engine beside the names.
  expect_error(
    do.call(ferx_globalsearch, c(args, list(ga = list(population_size = 2.5)))),
    "whole number"
  )
  expect_error(
    do.call(ferx_globalsearch, c(args, list(ga = list(final_downhill = 1)))),
    "TRUE or FALSE"
  )
  expect_error(
    do.call(ferx_globalsearch, c(args, list(penalties = list(gate = Inf)))),
    "finite"
  )
})

test_that("the engine's option keys are the ones R validates against", {
  keys <- ferx_rust_globalsearch_option_keys()
  expect_true(all(c("population_size", "generations", "seed", "final_downhill")
                  %in% keys$ga_name))
  expect_true(all(c("theta", "crash", "gate", "non_influential")
                  %in% keys$penalty_name))
  expect_equal(length(keys$ga_name), length(keys$ga_kind))
  expect_setequal(unique(keys$ga_kind), c("count", "number", "flag"))
})

test_that("a variability or allometry space is refused by name, before any fit", {
  ex <- ferx_example("two_cpt_oral_global")
  # The grid takes the structural statements and COVARIATE?; IIV? is
  # iivsearch's axis and ALLOMETRY is ferx_allometry()'s.
  expect_error(
    ferx_globalsearch(model = ex$model, data = ex$data,
                      search_space = "IIV?(@PK, exp)"),
    "IIV|iivsearch"
  )
})

test_that("a coverage gap is an error naming the feature, before any fit", {
  ex <- ferx_example("two_cpt_oral_global")
  # SEQ-ZO-FO is a different disposition rather than a different input term,
  # so the engine cannot build it from a template - and says so by name.
  expect_error(
    ferx_globalsearch(model = ex$model, data = ex$data,
                      search_space = "ABSORPTION([FO,SEQ-ZO-FO])"),
    "SEQ-ZO-FO"
  )
})

test_that("an exhaustive grid above max_models is refused, not truncated", {
  ex <- ferx_example("two_cpt_oral_global")
  expect_error(
    ferx_globalsearch(model = ex$model, data = ex$data,
                      search_space = "PERIPHERALS(0..1); LAGTIME([OFF,ON])",
                      algorithm = "exhaustive", max_models = 1),
    "max_models"
  )
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: it runs once and every assertion about the result object
# reads that one run. Exhaustive rather than the GA default, so the table is
# deterministic and every grid point is present; four axes rather than two so
# the run reaches models the strictness gate refuses, which is what the gate
# charge has to be visible on.

GLOBALSEARCH_SPACE <- paste(
  "PERIPHERALS(0..1); LAGTIME([OFF,ON])",
  "COVARIATE?(CL, WT, pow); COVARIATE?(CL, CRCL, pow)",
  sep = "\n"
)

globalsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("two_cpt_oral_global")
    dir <- file.path(tempdir(), "globalsearch-run")
    cached <<- ferx_globalsearch(
      model        = ex$model,
      data         = ex$data,
      search_space = GLOBALSEARCH_SPACE,
      algorithm    = "exhaustive",
      retries      = 0,
      directory    = dir,
      progress     = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- globalsearch_run()

  expect_s3_class(res, "ferx_globalsearch")
  expect_s3_class(res, "ferx_search_result")

  # The model table is the engine's own, column for column, with the three
  # charges appended - nothing here is spelled in R.
  expect_equal(names(res$models),
               c(ferx_rust_globalsearch_columns(),
                 "charge_non_influential", "charge_gate", "charge_crash"))
  expect_gt(nrow(res$models), 1L)
  expect_true(is.logical(res$models$passed))
  expect_true(is.logical(res$models$converged))

  # Four axes, sixteen points, and the input is a row of the table beside them.
  expect_equal(length(res$axes), 4L)
  expect_equal(res$space_size, 16)
  # A double, not an integer: 31 binary axes already overflow R's integer
  # range, and `as.integer()` would answer a searchable grid with NA
  # (PR #375 review).
  expect_type(res$space_size, "double")
  expect_setequal(names(res$axes),
                  c("PERIPHERALS", "LAGTIME", "CL-WT", "CL-CRCL"))
  expect_true("input" %in% res$models$id)
  expect_true(res$final_model_id %in% res$models$id)
  expect_equal(sum(res$models$selected), 1L)
  expect_equal(res$models$id[res$models$selected], res$final_model_id)

  expect_true(is.finite(res$input_fitness))
  expect_true(is.finite(res$final_fitness))
  # The search cannot end worse than it started: the input is a candidate too.
  expect_lte(res$final_fitness, res$input_fitness)
  expect_true(file.exists(res$final_model_path))

  # Exhaustive enumeration writes no generations.
  expect_equal(res$algorithm, "exhaustive")
  expect_null(res$generations)
  # The tool's own default criterion, which is not the family's.
  expect_match(res$criterion, "penal", ignore.case = TRUE)
})

test_that("the charges decompose the fitness the search ranked on", {
  skip_on_cran()
  res <- globalsearch_run()

  charges <- res$models$charge_non_influential + res$models$charge_gate +
    res$models$charge_crash
  # A row with a usable criterion: fitness is the criterion plus the charges.
  ok <- is.finite(res$models$criterion)
  expect_true(any(ok))
  expect_equal(res$models$fitness[ok],
               res$models$criterion[ok] + charges[ok], tolerance = 1e-9)
  # A row without one is charged the crash value outright.
  if (any(!ok)) {
    expect_equal(res$models$fitness[!ok], charges[!ok], tolerance = 1e-9)
    expect_equal(unname(res$models$charge_crash[!ok]),
                 rep(unname(res$penalties[["crash"]]), sum(!ok)))
  }

  # The charges are the schedule the run reports, not numbers invented here.
  expect_true(all(c("crash", "gate", "non_influential") %in% names(res$penalties)))
  gated <- !res$models$passed & is.finite(res$models$criterion)
  if (any(gated)) {
    expect_equal(unname(res$models$charge_gate[gated]),
                 rep(unname(res$penalties[["gate"]]), sum(gated)))
  }
  expect_equal(res$models$charge_non_influential,
               res$models$non_influential * unname(res$penalties[["non_influential"]]),
               tolerance = 1e-12)
  expect_true(all(res$models$charge_gate[res$models$passed] == 0))
})

test_that("the R model table equals the run's own models.csv", {
  skip_on_cran()
  res <- globalsearch_run()
  csv <- read.csv(file.path(res$directory, "models.csv"), stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$models))
  # Every engine column, in the engine's order; the three charges are the only
  # columns beyond the file's.
  expect_equal(names(res$models)[seq_along(names(csv))], names(csv))
  expect_equal(setdiff(names(res$models), names(csv)),
               c("charge_non_influential", "charge_gate", "charge_crash"))

  expect_equal(csv$id, res$models$id)
  expect_equal(csv$criterion, res$models$criterion, tolerance = 1e-6)
  expect_equal(csv$fitness, res$models$fitness, tolerance = 1e-6)
  # `read.csv` types the engine's booleans as text; the R object types them.
  expect_equal(as.logical(csv$selected), res$models$selected)
  expect_equal(as.logical(csv$passed), res$models$passed)
})

test_that("ferx_search_results() reads the global search's model table back", {
  skip_on_cran()
  res <- globalsearch_run()
  tab <- ferx_search_results(res$directory, type = "models")

  expect_equal(names(tab), ferx_rust_globalsearch_columns())
  expect_equal(attr(tab, "tool"), "globalsearch")
  expect_equal(tab$id, res$models$id)
  expect_equal(tab$fitness, res$models$fitness, tolerance = 1e-6)
  expect_true(is.logical(tab$passed))
  expect_true(is.integer(tab$non_influential))
  # `step` is a batch label here, not the layer index it is in the
  # variability tables - it must survive as text rather than become NA.
  expect_true(is.character(tab$step))
  expect_true(all(nzchar(tab$step)))

  # The candidate table is a different table of the same run: the runner gives
  # each batch its own directory, which is what `$candidates` stacks.
  expect_true(is.data.frame(res$candidates))
  expect_true(all(c("id", "criterion", "passed") %in% names(res$candidates)))
  expect_true(all(unique(res$candidates$run) %in% unique(res$models$step)))
})

test_that("every evaluated model's text comes back, named by id", {
  skip_on_cran()
  res <- globalsearch_run()

  expect_true(all(res$models$id %in% names(res$model_text)))
  expect_true(all(nzchar(res$model_text)))
  expect_match(res$model_text[[res$final_model_id]], "[structural_model]",
               fixed = TRUE)
  on_disk <- file.path(res$directory, "models",
                       paste0(res$final_model_id, ".ferx"))
  skip_if_not(file.exists(on_disk))
  expect_equal(paste(readLines(on_disk), collapse = "\n"),
               trimws(res$model_text[[res$final_model_id]], which = "right"))
})

test_that("the winning fit is a ferx_fit, not a refit", {
  skip_on_cran()
  res <- globalsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  final_row <- res$models[res$models$id == res$final_model_id, , drop = FALSE]
  expect_equal(unname(res$fit$ofv), final_row$ofv, tolerance = 1e-8)
  expect_true(!is.null(res$fit$theta))
  expect_true(is.data.frame(res$fit$sdtab))
})

test_that("a gate-failed genome is a row with its verdict and its gate charge", {
  skip_on_cran()
  # #364's bar: never silently dropped. Whatever the run's trajectory, a model
  # the gate refused must be present, say why, carry the gate charge, hold no
  # rank and not be selected.
  res <- globalsearch_run()
  gated <- res$models[!res$models$passed & is.finite(res$models$criterion), ,
                      drop = FALSE]
  expect_gt(nrow(gated), 0L)

  expect_true(all(!is.na(gated$failures) & nzchar(gated$failures)))
  expect_true(all(is.na(gated$rank)))
  expect_false(any(gated$selected))
  expect_true(all(gated$charge_gate > 0))
  # And the charge is what put its fitness above its criterion, so the table
  # does not read as though it lost on OFV.
  expect_true(all(gated$fitness > gated$criterion))
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  inline <- globalsearch_run()

  ex <- ferx_example("two_cpt_oral_global")
  cfg <- write_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[space]",
    "mfl = '''",
    GLOBALSEARCH_SPACE,
    "'''",
    "[run]",
    "retries = 0",
    "[globalsearch]",
    'algorithm = "exhaustive"'
  )
  from_file <- ferx_globalsearch(
    config    = cfg,
    directory = file.path(tempdir(), "globalsearch-file"),
    progress  = FALSE
  )

  expect_equal(from_file$models$id, inline$models$id)
  expect_equal(from_file$models$genome, inline$models$genome)
  expect_equal(from_file$models$fitness, inline$models$fitness, tolerance = 1e-6)
  expect_equal(from_file$final_model_id, inline$final_model_id)
  expect_equal(from_file$axes, inline$axes)
})

test_that("a single-point space returns the input model, fitted", {
  skip_on_cran()
  # The degenerate oracle: a grid the input model already sits on has one
  # point, so the search must hand back that model's fit rather than a
  # transformed one, and must land on the same optimum as fitting it
  # directly.
  #
  # Not bit-identical, and cannot be: the grid point is the input model
  # refitted, not a transformed one, and the fit handed back is that row's.
  #
  # What is deliberately NOT asserted is numeric agreement with a standalone
  # `ferx_fit()` of the same model. That comparison is not a property of this
  # package: the cold optimizer trajectory on this model is platform
  # dependent, and on Linux / R 4.6 the standalone fit lands in a far worse
  # basin than the one the search reaches (OFV -62.0 against -172.0), while
  # on macOS the two agree to 9e-4 (-172.049 against -172.196). A tolerance
  # calibrated on either platform is false on the other - CI on PR #375
  # failed exactly there. So the oracle is the structural claim below plus
  # the one numeric direction that means something: the search, which fits
  # the same model, must not come back *worse* than fitting it directly.
  # Same lesson as the status-aware AMD assertions in #356.
  ex <- ferx_example("two_cpt_oral_global")
  res <- ferx_globalsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "PERIPHERALS(0)",
    algorithm    = "exhaustive",
    retries      = 0,
    directory    = file.path(tempdir(), "globalsearch-degenerate"),
    progress     = FALSE
  )

  expect_equal(res$space_size, 1)
  # The input and the one grid point, nothing else.
  expect_equal(nrow(res$models), 2L)
  expect_equal(sum(res$models$selected), 1L)

  # The winner is the grid's one point, and that point is the input model
  # unchanged - one compartment, no covariate relation, nothing transformed.
  winner <- res$models[res$models$selected, , drop = FALSE]
  expect_equal(winner$id, res$final_model_id)
  expect_equal(winner$genome, "PERIPHERALS=0")
  expect_equal(winner$peripherals, 0L)
  expect_true(is.na(winner$covariates))

  skip_if(is.null(res$fit), "the run recovered no final fit")
  # The fit belongs to the row the table selected - that is the whole of
  # "returns the base fit" that this package is responsible for.
  expect_s3_class(res$fit, "ferx_fit")
  expect_equal(unname(res$fit$ofv), winner$ofv, tolerance = 1e-8)

  direct <- ferx_fit(ex$model, ex$data)
  expect_lte(unname(res$fit$ofv), unname(direct$ofv) + 1e-6)
})

test_that("the genetic algorithm reports its trajectory and finds the same winner", {
  skip_on_cran()
  exhaustive <- globalsearch_run()

  ex <- ferx_example("two_cpt_oral_global")
  res <- ferx_globalsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = GLOBALSEARCH_SPACE,
    algorithm    = "ga",
    ga           = list(population_size = 4, generations = 2, seed = 20250914),
    retries      = 0,
    directory    = file.path(tempdir(), "globalsearch-ga"),
    progress     = FALSE
  )

  expect_equal(res$algorithm, "ga")
  expect_true(is.data.frame(res$generations))
  expect_true(all(c("index", "best", "best_fitness", "mean_fitness", "polished")
                  %in% names(res$generations)))
  expect_gt(nrow(res$generations), 0L)
  # The knobs the call passed are the ones the engine read back.
  expect_equal(unname(res$ga[["population_size"]]), 4)
  expect_equal(unname(res$ga[["generations"]]), 2)
  expect_equal(unname(res$ga[["seed"]]), 20250914)

  # The enumeration is the grid's optimum by construction, so the algorithm
  # cannot beat it - that, not equality, is the invariant: a population of
  # four over sixteen points may or may not reach the best cell, and a test
  # that demanded it would be asserting the GA's luck.
  expect_gte(res$final_fitness, exhaustive$final_fitness - 1e-6)
  # It does have to beat doing nothing.
  expect_lte(res$final_fitness, res$input_fitness)
  # And every generation's best is a model in the table.
  expect_true(all(res$generations$best %in% res$models$id))
})

test_that("a penalties override reaches the schedule the run charged", {
  skip_on_cran()
  ex <- ferx_example("two_cpt_oral_global")
  res <- ferx_globalsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "PERIPHERALS(0)",
    algorithm    = "exhaustive",
    rank         = "bic",
    penalties    = list(gate = 250, crash = 12345678),
    retries      = 0,
    directory    = file.path(tempdir(), "globalsearch-penalties"),
    progress     = FALSE
  )

  expect_equal(unname(res$penalties[["gate"]]), 250)
  # Eight significant digits: `format()`'s 7-digit default would have written
  # `1.234568e+07` into the config file.
  expect_equal(unname(res$penalties[["crash"]]), 12345678)
  # The three search-level charges apply under any criterion, so naming a BIC
  # does not turn them off.
  expect_match(res$criterion, "BIC", ignore.case = TRUE)
})

test_that("a [globalsearch] file no longer warns that no R tool runs it", {
  # #364 removes the #347 warning for this one section: it now has a tool.
  cfg <- write_cfg(
    'base = "model.ferx"',
    "[space]",
    'mfl = "COVARIATE?(CL, WT, pow)"',
    "[globalsearch]",
    'algorithm = "exhaustive"'
  )
  expect_silent(ferx_search_config(cfg))
  expect_equal(ferx_search_config(cfg)$tools, "globalsearch")
})

test_that("print() shows the grid and the charges, and summary() shows the rest", {
  skip_on_cran()
  res <- globalsearch_run()
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx global model search")
  expect_match(out, "Grid:")
  expect_match(out, "Ranked models")

  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "Every model")
  expect_match(full, "Penalty schedule")
})

test_that("print() renders a grid larger than R's integer range", {
  skip_on_cran()
  # 31 binary axes is 2^31 points - a grid the GA searches perfectly well
  # (it evaluates only its population), and one `%d` would refuse outright.
  # The object is doctored rather than run: the point under test is the
  # formatting, not a 2-billion-fit search.
  res <- globalsearch_run()
  res$space_size <- 2^31

  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "2,147,483,648 points", fixed = TRUE)
  expect_no_match(out, "NA points", fixed = TRUE)
})
