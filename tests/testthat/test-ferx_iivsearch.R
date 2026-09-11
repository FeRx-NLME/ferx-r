# Variability-structure search from R.
#
# The trajectory itself is the engine's (ferx-core #1183 anchors it against
# Pharmpy's iivsearch); what is asserted here is the R surface's own contract:
# the two entry forms cannot disagree, the object reports what the engine ran
# rather than anything recomputed in R, the two stages stay two stages, and
# every structure row is labelled with the model's own declared eta names.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[space]",
                   'mfl = "IIV?(@PK, exp)"')

  expect_error(ferx_iivsearch(config = cfg, search_space = "IIV?(CL, exp)"),
               "search_space")
  expect_error(ferx_iivsearch(config = cfg, algorithm = "bottom_up_stepwise"),
               "algorithm")
  expect_error(ferx_iivsearch(config = cfg, as_fullblock = TRUE),
               "as_fullblock")
  expect_error(ferx_iivsearch(config = cfg, model = "model.ferx"), "model")
  # `data` states *which dataset* the search runs on, so it belongs with
  # `model` rather than with the run knobs.
  expect_error(ferx_iivsearch(config = cfg, data = "other.csv"), "`data`")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that does
  # not exist - but not on the entry form.
  msg <- tryCatch(ferx_iivsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_iivsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_iivsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(ferx_iivsearch(model = tempfile(fileext = ".ferx"),
                              search_space = "IIV?(@PK, exp)"),
               "model file not found")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("warfarin_block_omega")
  args <- list(model = ex$model, data = ex$data,
               search_space = "IIV?(@PK, exp)")

  expect_error(do.call(ferx_iivsearch, c(args, list(algorithm = "stepwise"))),
               "should be one of")
  expect_error(do.call(ferx_iivsearch,
                       c(args, list(correlation_algorithm = "bottom_up"))),
               "should be one of")
  expect_error(do.call(ferx_iivsearch, c(args, list(as_fullblock = NA))),
               "TRUE, FALSE or NULL")
  expect_error(do.call(ferx_iivsearch, c(args, list(block_retries = 1.5))),
               "whole number")
  expect_error(do.call(ferx_iivsearch, c(args, list(cutoff = c(1, 2)))),
               "single number")
  # An infinite cutoff would be dropped on the way to the config file - the
  # binding emits the key only when it is finite - so the search would run at
  # the engine default rather than at the level that was asked for.
  expect_error(do.call(ferx_iivsearch, c(args, list(cutoff = Inf))), "finite")
  expect_error(do.call(ferx_iivsearch, c(args, list(rank = c("bic", "aic")))),
               "single string")
  expect_error(do.call(ferx_iivsearch, c(args, list(resume = NA))),
               "TRUE or FALSE")
  expect_error(do.call(ferx_iivsearch, c(args, list(retries = -1))),
               "at least 0")
  # The inline form has no default space: iivsearch without one would search
  # nothing.
  expect_error(ferx_iivsearch(model = ex$model, data = ex$data),
               "`search_space` is required")
})

test_that("the engine's own validation runs before the dataset is read", {
  # These all name a file key rather than an R argument, and all of them are
  # refused while the config is being read - the base model named below does
  # not exist, so anything that got as far as loading it would fail differently.
  expect_error(ferx_iivsearch(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "IIV?(@PK, exp)"',
    "[iivsearch]", 'algorithm = "skip"'
  )), "correlation_algorithm")
  expect_error(ferx_iivsearch(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "IIV?(@PK, exp)"',
    "[iivsearch]", 'algorithm = "simultaneous_stepwise"',
    'correlation_algorithm = "top_down_exhaustive"'
  )), "simultaneous_stepwise")
  # A structural space is another tool's.
  expect_error(ferx_iivsearch(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "PERIPHERALS(0..1)"'
  )), "iivsearch")
  # Every one of those failed on the file's own keys. The same file with
  # nothing wrong in it gets as far as the base model and fails on that
  # instead, which is what says the refusals above were not this one.
  expect_error(ferx_iivsearch(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "IIV?(@PK, exp)"'
  )), "read model file")
})

test_that("an inline cutoff reaches the engine, which names it", {
  # The other half of the sentinel contract: a value R accepts is emitted into
  # the config file, so the engine's own validation is what rejects it. R does
  # not duplicate that rule, it only guarantees the value arrives.
  ex <- ferx_example("warfarin_block_omega")
  expect_error(ferx_iivsearch(model = ex$model, data = ex$data,
                              search_space = "IIV?(@PK, exp)", cutoff = -1),
               "cutoff")
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: slow, so it runs once and every assertion about the result
# object reads that one run.
#
# warfarin_block_omega carries a block over ETA_CL and ETA_V plus a diagonal
# ETA_KA on 10 subjects, so both stages have something to decide: which of the
# three etas survive, and which of the survivors end up correlated.

iivsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("warfarin_block_omega")
    cached <<- ferx_iivsearch(
      model        = ex$model,
      data         = ex$data,
      search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
      retries      = 0,
      directory    = file.path(tempdir(), "iivsearch-run"),
      progress     = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- iivsearch_run()

  expect_s3_class(res, "ferx_iivsearch")
  expect_s3_class(res, "ferx_search_result")

  # The model table is the engine's own, column for column, with the stage and
  # the label columns appended - nothing here is spelled in R.
  expect_equal(names(res$models),
               c(ferx_rust_iivsearch_columns(),
                 "step_kind", "eta_labels", "block_labels", "structure"))
  expect_gt(nrow(res$models), 2L)
  # Termination status and the gate verdict are columns, never omitted.
  expect_true(is.logical(res$models$converged))
  expect_true(is.logical(res$models$passed))
  # The starts a candidate was fitted with are on its row: a block candidate is
  # fitted from more starting points than a diagonal one.
  expect_true(is.integer(res$models$starts))
  expect_true(all(res$models$starts >= 1L))

  # The input is row zero and belongs to no stage.
  expect_equal(res$models$id[1L], "input")
  expect_equal(res$models$step[1L], 0L)
  expect_true(is.na(res$models$step_kind[1L]))

  # The options come back as the engine read them, not as R passed them.
  expect_equal(res$options$algorithm, "top_down_exhaustive")
  expect_equal(res$options$block_retries, 2L)
  expect_equal(res$options$starts, 1L)
  expect_false(res$options$as_fullblock)
  expect_true(res$block_stage)
  expect_equal(res$criterion, "bic_iiv")

  expect_true(is.finite(res$final_criterion))
  expect_true(file.exists(res$final_model_path))
})

test_that("both stages are visible, each with its StepKind", {
  skip_on_cran()
  res <- iivsearch_run()

  # The number of etas and the block structure are two decisions, taken against
  # two parents. A table that collapsed them would report a search that was
  # never run.
  expect_true(all(c("no_of_etas", "block_structure") %in% res$steps$kind))
  expect_true("compare_to_input" %in% res$steps$kind)
  expect_equal(names(res$steps),
               c("step", "kind", "parent", "id", "criterion", "d_criterion",
                 "rank", "best"))

  eta_step <- res$steps[res$steps$kind == "no_of_etas", , drop = FALSE]
  blk_step <- res$steps[res$steps$kind == "block_structure", , drop = FALSE]
  expect_gt(nrow(eta_step), 1L)
  expect_gt(nrow(blk_step), 0L)
  # The block stage starts from the eta stage's winner, not from the input.
  expect_equal(unique(blk_step$parent),
               eta_step$id[eta_step$best & eta_step$step == min(eta_step$step)])
  # Exactly one winner per step, and it is a row of that step.
  for (s in unique(res$steps$step)) {
    one <- res$steps[res$steps$step == s, , drop = FALSE]
    expect_equal(sum(one$best), 1L)
  }

  # The stage is on the model rows too, so a table filtered to one stage is
  # still the engine's table.
  expect_setequal(
    unique(res$models$step_kind[res$models$step > 0L]),
    unique(res$steps$kind[res$steps$step %in% res$models$step])
  )
})

test_that("every structure row is labelled with the model's declared etas", {
  skip_on_cran()
  res <- iivsearch_run()

  # The label convention (CLAUDE.md): the bare declared name, never OMEGA(i,i)
  # on a model that names its etas. warfarin_block_omega declares ETA_CL,
  # ETA_V and ETA_KA, so every label is one of those.
  labels <- unlist(strsplit(stats::na.omit(res$models$eta_labels), ";",
                            fixed = TRUE))
  expect_gt(length(labels), 0L)
  expect_true(all(labels %in% c("ETA_CL", "ETA_V", "ETA_KA")))
  expect_false(any(grepl("OMEGA(", res$models$structure, fixed = TRUE)))

  # The labelled structure is the engine's `description` with the parameters
  # replaced by their etas, so the two say the same thing in two spellings.
  for (i in seq_len(nrow(res$models))) {
    expect_equal(
      gsub("ETA_", "", res$models$structure[i], fixed = TRUE),
      if (is.na(res$models$description[i])) "" else res$models$description[i]
    )
  }

  # A block is a block in both columns.
  blocked <- !is.na(res$models$block_labels)
  expect_equal(
    gsub("ETA_", "", res$models$block_labels[blocked], fixed = TRUE),
    res$models$blocks[blocked]
  )

  # The input carries the block the model file declares, under its own names.
  expect_equal(res$input_structure, "[ETA_CL,ETA_V]+[ETA_KA]")
  expect_equal(res$input_description, "[CL,V]+[KA]")
})

test_that("a model whose etas are not named ETA_<P> is labelled as written", {
  skip_on_cran()
  # The trap the label convention exists to catch: a table that spelled
  # `ETA_<parameter>` for itself would be right on every bundled model and
  # wrong on this one.
  ex <- ferx_example("warfarin_block_omega")
  text <- readLines(ex$model)
  text <- gsub("ETA_CL", "E1", text, fixed = TRUE)
  text <- gsub("ETA_V", "E2", text, fixed = TRUE)
  text <- gsub("ETA_KA", "E3", text, fixed = TRUE)
  model <- tempfile(fileext = ".ferx")
  writeLines(text, model)

  res <- ferx_iivsearch(
    model        = model,
    data         = ex$data,
    search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
    algorithm             = "skip",
    correlation_algorithm = "top_down_exhaustive",
    retries      = 0,
    directory    = file.path(tempdir(), "iivsearch-renamed-etas"),
    progress     = FALSE
  )
  labels <- unlist(strsplit(stats::na.omit(res$models$eta_labels), ";",
                            fixed = TRUE))
  expect_true(all(labels %in% c("E1", "E2", "E3")))
  expect_false(any(grepl("ETA_", res$models$structure, fixed = TRUE)))
  expect_equal(res$input_structure, "[E1,E2]+[E3]")
})

test_that("the R model table equals the run's own models.csv", {
  skip_on_cran()
  res <- iivsearch_run()
  csv <- read.csv(file.path(res$directory, "models.csv"),
                  stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$models))
  # Every engine column, in the engine's order; the four added here are the
  # only columns beyond the file's.
  expect_equal(names(res$models)[seq_along(names(csv))], names(csv))
  expect_equal(setdiff(names(res$models), names(csv)),
               c("step_kind", "eta_labels", "block_labels", "structure"))

  expect_equal(csv$id, res$models$id)
  expect_equal(csv$description, res$models$description)
  expect_equal(csv$etas, ifelse(is.na(res$models$etas), "", res$models$etas))
  expect_equal(csv$criterion, res$models$criterion, tolerance = 1e-6)
  expect_equal(csv$starts, res$models$starts)
  # `read.csv` types the engine's booleans as text; the R object types them.
  expect_equal(as.logical(csv$selected), res$models$selected)
  expect_equal(as.logical(csv$passed), res$models$passed)
})

test_that("ferx_search_results() reads the model table back", {
  skip_on_cran()
  res <- iivsearch_run()
  tab <- ferx_search_results(res$directory, type = "models")

  expect_equal(names(tab), ferx_rust_iivsearch_columns())
  # Three tools write a `models.csv`; the header says which one did.
  expect_equal(attr(tab, "tool"), "iivsearch")
  expect_equal(tab$id, res$models$id)
  expect_equal(tab$criterion, res$models$criterion, tolerance = 1e-6)
  expect_true(is.logical(tab$passed))
  expect_true(is.integer(tab$step))
  expect_true(is.integer(tab$starts))

  # The candidate table is a different table of the same run, and still reads.
  # Each step gets its own runner directory, which is what `$candidates`
  # stacks and labels by step.
  expect_true(is.data.frame(ferx_search_results(file.path(res$directory, "input"))))
  expect_true(all(c("input", "step-1") %in% unique(res$candidates$run)))
  expect_error(ferx_search_results(res$directory, type = "models", partial = TRUE),
               "no partial model table")
})

test_that("a model table that is no tool's is refused rather than mistyped", {
  dir <- file.path(tempdir(), "iivsearch-not-a-model-table")
  dir.create(dir, showWarnings = FALSE)
  write.csv(data.frame(a = 1, b = 2), file.path(dir, "models.csv"),
            row.names = FALSE)
  expect_error(ferx_search_results(dir, type = "models"),
               "not a search model table")
})

test_that("every fitted model's text comes back, named by model id", {
  skip_on_cran()
  res <- iivsearch_run()

  expect_true(all(res$models$id %in% names(res$model_text)))
  expect_true(all(nzchar(res$model_text)))
  # A structure the search rejected is still readable, which is the point of
  # keeping the texts: it can be refitted without re-running the search.
  rejected <- setdiff(res$models$id, res$final_model_id)
  expect_gt(length(rejected), 0L)
  expect_true(all(grepl("[individual_parameters]",
                        res$model_text[rejected], fixed = TRUE)))
})

test_that("the winning fit is a ferx_fit whose omega carries the eta names", {
  skip_on_cran()
  res <- iivsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  # The label convention again, on the fit this time: omega dimnames are the
  # declared eta names, and they are the labels the table reported.
  expect_equal(rownames(res$fit$omega), res$final_etas)
  expect_equal(colnames(res$fit$omega), res$final_etas)
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  ex <- ferx_example("warfarin_block_omega")
  inline <- ferx_iivsearch(
    model        = ex$model,
    data         = ex$data,
    search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
    algorithm             = "skip",
    correlation_algorithm = "top_down_exhaustive",
    retries      = 0,
    directory    = file.path(tempdir(), "iivsearch-inline"),
    progress     = FALSE
  )
  cfg <- write_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[space]",
    'mfl = "IIV?(@PK, exp); COVARIANCE?(IIV, *)"',
    "[iivsearch]",
    'algorithm = "skip"',
    'correlation_algorithm = "top_down_exhaustive"',
    "[run]",
    "retries = 0"
  )
  from_file <- ferx_iivsearch(config = cfg,
                              directory = file.path(tempdir(), "iivsearch-file"),
                              progress = FALSE)

  expect_equal(inline$models$id, from_file$models$id)
  expect_equal(inline$models$description, from_file$models$description)
  expect_equal(inline$models$structure, from_file$models$structure)
  expect_equal(inline$final_model_id, from_file$final_model_id)
  expect_equal(inline$final_criterion, from_file$final_criterion,
               tolerance = 1e-6)
  expect_equal(inline$options, from_file$options)
  # The file names itself on the object; the inline form has nothing to name.
  expect_equal(from_file$config, normalizePath(cfg))
  expect_true(is.na(inline$config))
})

test_that("a single-structure space returns the base model, fitted", {
  skip_on_cran()
  # The degenerate oracle: `IIV(...)` keeps every eta and `correlation_algorithm
  # = "skip"` leaves no block to decide, so the only structure in the space is
  # the one the model already has. A search with nothing to choose returns what
  # it started from, and says so with a real fit rather than a promise of one.
  ex <- ferx_example("warfarin_block_omega")
  res <- ferx_iivsearch(
    model                 = ex$model,
    data                  = ex$data,
    search_space          = "IIV(@PK, exp)",
    correlation_algorithm = "skip",
    retries               = 0,
    directory             = file.path(tempdir(), "iivsearch-degenerate"),
    progress              = FALSE
  )

  expect_equal(res$final_model_id, "input")
  expect_equal(res$final_structure, res$input_structure)
  expect_equal(res$final_description, res$input_description)
  expect_false(res$block_stage)
  expect_false("block_structure" %in% res$steps$kind)
  expect_s3_class(res$fit, "ferx_fit")
})

test_that("a shorter run does not inherit the last run's candidate tables", {
  skip_on_cran()
  # #336 review: the engine rewrites the steps it executes but removes nothing,
  # so a two-stage run leaves a `step-2/` directory that a later one-stage run
  # in the same directory would otherwise fold into its `$candidates` - a
  # result object contradicting its own model table.
  ex <- ferx_example("warfarin_block_omega")
  dir <- file.path(tempdir(), "iivsearch-reused-directory")
  args <- list(model = ex$model, data = ex$data,
               search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
               retries = 0, directory = dir, progress = FALSE)

  long <- do.call(ferx_iivsearch, args)
  skip_if(max(long$models$step) < 2L, "the long run did not reach a second stage")
  expect_true("step-2" %in% unique(long$candidates$run))

  short <- do.call(ferx_iivsearch,
                   c(args, list(correlation_algorithm = "skip")))
  expect_equal(max(short$models$step), 1L)
  expect_false("step-2" %in% unique(short$candidates$run))
  # The stale directory is still on disk - this is the object refusing to read
  # it, not the run deleting a previous run's journal.
  expect_true(dir.exists(file.path(dir, "step-2")))
})

test_that("print() shows the per-stage tables and summary() the rest", {
  skip_on_cran()
  res <- iivsearch_run()

  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx variability-structure search")
  expect_match(out, "no_of_etas")
  expect_match(out, "block_structure")
  # The structures printed are the labelled ones, per the output convention.
  expect_match(out, "ETA_CL")
  expect_match(out, "Selected variability structure")

  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "Structures not selected, and why")
  expect_match(full, "Every model:")
  expect_true(nchar(full) > nchar(out))
})

test_that("an off-diagonal is printed with the convention's tilde", {
  skip_on_cran()
  res <- iivsearch_run()
  skip_if(length(res$final_blocks) == 0L,
          "the search selected a diagonal structure")

  out <- paste(capture.output(print(res)), collapse = "\n")
  # `ETA_V ~ ETA_CL`, not `OMEGA(2,1)`: the whole point of this search's output
  # is saying which parameters ended up correlated.
  expect_match(out, "Correlated:")
  expect_match(out, "ETA_[A-Z0-9_]+ ~ ETA_[A-Z0-9_]+")
})
