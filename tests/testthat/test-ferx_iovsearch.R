# Inter-occasion variability search from R.
#
# The trajectory itself is the engine's (ferx-core #1183 anchors it against
# Pharmpy's iovsearch); what is asserted here is the R surface's own contract:
# the two entry forms cannot disagree, the object reports what the engine ran
# rather than anything recomputed in R, both removal steps stay visible, and
# every structure row is labelled with the model's own declared kappa names.

iov_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- iov_cfg('base = "model.ferx"', "[iovsearch]", 'distribution = "joint"')

  expect_error(ferx_iovsearch(config = cfg, distribution = "disjoint"),
               "distribution")
  expect_error(ferx_iovsearch(config = cfg, column = "OCC"), "column")
  expect_error(ferx_iovsearch(config = cfg, groups = list("CL")), "groups")
  expect_error(ferx_iovsearch(config = cfg, model = "model.ferx"), "model")
  expect_error(ferx_iovsearch(config = cfg, data = "other.csv"), "`data`")

  msg <- tryCatch(ferx_iovsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_iovsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_iovsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(ferx_iovsearch(model = tempfile(fileext = ".ferx")),
               "model file not found")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("warfarin_iov")
  args <- list(model = ex$model, data = ex$data)

  expect_error(do.call(ferx_iovsearch, c(args, list(distribution = "same_as_iiv"))),
               "should be one of")
  expect_error(do.call(ferx_iovsearch, c(args, list(column = c("OCC", "VISIT")))),
               "single column name")
  expect_error(do.call(ferx_iovsearch, c(args, list(groups = list()))),
               "list of character vectors")
  expect_error(do.call(ferx_iovsearch, c(args, list(groups = list(1)))),
               "list of character vectors")
  expect_error(do.call(ferx_iovsearch, c(args, list(block_retries = 1.5))),
               "whole number")
  expect_error(do.call(ferx_iovsearch, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_iovsearch, c(args, list(rank = 1))), "single string")
  expect_error(do.call(ferx_iovsearch, c(args, list(resume = NA))),
               "TRUE or FALSE")
})

test_that("the engine's own validation runs before the dataset is read", {
  expect_error(ferx_iovsearch(config = iov_cfg(
    'base = "model.ferx"', "[iovsearch]", 'distribution = "explicit"'
  )), "groups")
  expect_error(ferx_iovsearch(config = iov_cfg(
    'base = "model.ferx"', "[iovsearch]", 'distribution = "joint"',
    'groups = [["CL", "V"]]'
  )), "explicit")
  expect_error(ferx_iovsearch(config = iov_cfg(
    'base = "model.ferx"', "[iovsearch]", 'distribution = "explicit"',
    'groups = [["CL", "V"], ["V", "KA"]]'
  )), "two groups")
})

test_that("a base model without an occasion column is refused by name", {
  skip_on_cran()
  # The occasions are read by the engine from the base's own `iov_column`; a
  # base without one has a population with no occasions, so no candidate could
  # be fitted. Saying so beats an empty search reporting the input as a winner.
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_iovsearch(model = ex$model, data = ex$data,
                   directory = file.path(tempdir(), "iovsearch-no-column"),
                   progress = FALSE),
    "iov_column"
  )
})

test_that("a column that disagrees with the base model is refused", {
  skip_on_cran()
  ex <- ferx_example("warfarin_iov")
  expect_error(
    ferx_iovsearch(model = ex$model, data = ex$data, column = "VISIT",
                   directory = file.path(tempdir(), "iovsearch-wrong-column"),
                   progress = FALSE),
    "VISIT"
  )
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: slow, so it runs once and every assertion about the result
# object reads that one run.
#
# warfarin_iov reads its occasions from OCC and already carries KAPPA_CL, so
# the full-IOV model adds a kappa to V and KA and the two removal steps have
# something to decide.

iovsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("warfarin_iov")
    cached <<- ferx_iovsearch(
      model     = ex$model,
      data      = ex$data,
      retries   = 0,
      directory = file.path(tempdir(), "iovsearch-run"),
      progress  = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- iovsearch_run()

  expect_s3_class(res, "ferx_iovsearch")
  expect_s3_class(res, "ferx_search_result")

  expect_equal(names(res$models),
               c(ferx_rust_iovsearch_columns(),
                 "step_kind", "eta_labels", "kappa_labels",
                 "kappa_block_labels", "structure"))
  expect_gt(nrow(res$models), 2L)
  expect_true(is.logical(res$models$converged))
  expect_true(is.logical(res$models$passed))
  expect_true(is.integer(res$models$starts))

  expect_equal(res$models$id[1L], "input")
  expect_equal(res$models$step[1L], 0L)

  # The occasion column comes back as the engine resolved it, not as R guessed.
  expect_equal(res$column, "OCC")
  expect_equal(res$distribution, "same-as-iiv")
  expect_equal(res$criterion, "bic_random")
  expect_equal(res$options$block_retries, 2L)
  expect_equal(res$options$starts, 1L)

  expect_true(is.finite(res$input_criterion))
  expect_true(is.finite(res$final_criterion))
  expect_true(file.exists(res$final_model_path))
})

test_that("both removal steps are visible, each with its parent", {
  skip_on_cran()
  res <- iovsearch_run()

  expect_equal(names(res$steps),
               c("step", "kind", "parent", "id", "criterion", "d_criterion",
                 "rank", "best"))
  # Step 1 removes kappas from the full-IOV model; step 2 removes the etas of
  # the parameters that kept one, taking step 1's winner as its parent.
  expect_true(all(c(1L, 2L) %in% res$steps$step))
  first <- res$steps[res$steps$step == 1L, , drop = FALSE]
  second <- res$steps[res$steps$step == 2L, , drop = FALSE]
  expect_equal(unique(second$parent), first$id[first$best])
  for (s in unique(res$steps$step)) {
    expect_equal(sum(res$steps$best[res$steps$step == s]), 1L)
  }
  # The stage is on the model rows too.
  expect_true(all(res$models$step %in% c(0L, res$steps$step)))
})

test_that("every structure row is labelled with the declared random effects", {
  skip_on_cran()
  res <- iovsearch_run()

  # The label convention (CLAUDE.md): the bare declared names, never
  # OMEGA(i,i) / KAPPA<i> on a model that names its random effects.
  etas <- unlist(strsplit(stats::na.omit(res$models$eta_labels), ";",
                          fixed = TRUE))
  kappas <- unlist(strsplit(stats::na.omit(res$models$kappa_labels), ";",
                            fixed = TRUE))
  expect_true(all(etas %in% c("ETA_CL", "ETA_V", "ETA_KA")))
  expect_true(all(kappas %in% c("KAPPA_CL", "KAPPA_V", "KAPPA_KA")))
  expect_false(any(grepl("OMEGA(", res$models$structure, fixed = TRUE)))
  expect_false(any(grepl("KAPPA1", res$models$structure, fixed = TRUE)))

  # The labelled structure is the engine's `description` with the parameters
  # replaced by their random effects, so the two say the same in two spellings.
  stripped <- gsub("KAPPA_", "", gsub("ETA_", "", res$models$structure,
                                      fixed = TRUE), fixed = TRUE)
  expect_equal(stripped, res$models$description)

  expect_match(res$input_structure, "IOV\\(\\[KAPPA_CL\\]\\)")
  expect_equal(res$input_description, "IIV([CL]+[KA]+[V]);IOV([CL])")
})

test_that("the R model table equals the run's own models.csv", {
  skip_on_cran()
  res <- iovsearch_run()
  csv <- read.csv(file.path(res$directory, "models.csv"),
                  stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$models))
  expect_equal(names(res$models)[seq_along(names(csv))], names(csv))
  expect_equal(setdiff(names(res$models), names(csv)),
               c("step_kind", "eta_labels", "kappa_labels",
                 "kappa_block_labels", "structure"))

  expect_equal(csv$id, res$models$id)
  expect_equal(csv$description, res$models$description)
  expect_equal(csv$kappas, ifelse(is.na(res$models$kappas), "",
                                  res$models$kappas))
  expect_equal(csv$criterion, res$models$criterion, tolerance = 1e-6)
  expect_equal(as.logical(csv$selected), res$models$selected)
})

test_that("ferx_search_results() tells an iovsearch table from an iivsearch one", {
  skip_on_cran()
  res <- iovsearch_run()
  tab <- ferx_search_results(res$directory, type = "models")

  expect_equal(names(tab), ferx_rust_iovsearch_columns())
  expect_equal(attr(tab, "tool"), "iovsearch")
  expect_equal(tab$id, res$models$id)
  expect_true(is.integer(tab$step))
  expect_true(is.logical(tab$passed))

  # The runner directories this tool writes are its own.
  expect_true(all(c("input", "iov-all") %in% unique(res$candidates$run)))
})

test_that("every fitted model's text comes back, named by model id", {
  skip_on_cran()
  res <- iovsearch_run()

  expect_true(all(res$models$id %in% names(res$model_text)))
  expect_true(all(nzchar(res$model_text)))
  # The full-IOV model is the one with a kappa on every candidate parameter,
  # and it stays readable whether or not the search kept it.
  expect_true(any(grepl("kappa", res$model_text, ignore.case = TRUE)))
})

test_that("the winning fit is a ferx_fit", {
  skip_on_cran()
  res <- iovsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  expect_true(all(nzchar(res$final_kappas)))
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  ex <- ferx_example("warfarin_iov")
  inline <- ferx_iovsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "IOV?(*, exp)",
    distribution = "disjoint",
    retries      = 0,
    directory    = file.path(tempdir(), "iovsearch-inline"),
    progress     = FALSE
  )
  cfg <- iov_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[space]",
    'mfl = "IOV?(*, exp)"',
    "[iovsearch]",
    'distribution = "disjoint"',
    "[run]",
    "retries = 0"
  )
  from_file <- ferx_iovsearch(config = cfg,
                              directory = file.path(tempdir(), "iovsearch-file"),
                              progress = FALSE)

  expect_equal(inline$models$id, from_file$models$id)
  expect_equal(inline$models$description, from_file$models$description)
  expect_equal(inline$models$structure, from_file$models$structure)
  expect_equal(inline$final_model_id, from_file$final_model_id)
  expect_equal(inline$final_criterion, from_file$final_criterion,
               tolerance = 1e-6)
  expect_equal(inline$options, from_file$options)
  expect_equal(from_file$config, normalizePath(cfg))
  expect_true(is.na(inline$config))
})

test_that("a space with nothing to add is refused by name, not searched", {
  skip_on_cran()
  # The degenerate oracle. CL already carries the only kappa this space allows,
  # so the space holds exactly one structure - the model's own. The engine says
  # so and stops, rather than fitting the input twice and reporting it as a
  # decision; a search that "selected" the model it was given, with no
  # alternative ever built, is the failure this refusal exists to prevent.
  ex <- ferx_example("warfarin_iov")
  expect_error(
    ferx_iovsearch(
      model        = ex$model,
      data         = ex$data,
      search_space = "IOV(CL, exp)",
      retries      = 0,
      directory    = file.path(tempdir(), "iovsearch-degenerate"),
      progress     = FALSE
    ),
    "nothing to add"
  )
})

test_that("print() shows the per-step tables and summary() the rest", {
  skip_on_cran()
  res <- iovsearch_run()

  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx inter-occasion variability search")
  expect_match(out, "Occasions: OCC")
  expect_match(out, "KAPPA_CL")
  expect_match(out, "Selected variability structure")

  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "Structures not selected, and why")
  expect_match(full, "Every model:")
  expect_true(nchar(full) > nchar(out))
})
