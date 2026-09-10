# Residual-error model search from R.
#
# The trajectory itself is the engine's (ferx-core #1182 anchors it against
# Pharmpy's ruvsearch); what is asserted here is the R surface's own contract:
# the two entry forms cannot disagree, the object reports what the engine ran
# rather than anything recomputed in R, and a file ruvsearch cannot honour is
# an error before the first fit.

write_cfg <- function(...) {
  path <- tempfile(fileext = ".ferxsearch")
  writeLines(c(...), path)
  path
}

test_that("the entry forms are mutually exclusive", {
  cfg <- write_cfg('base = "model.ferx"', "[ruvsearch]", "max_iter = 1")

  expect_error(ferx_ruvsearch(config = cfg, p_value = 0.05), "p_value")
  expect_error(ferx_ruvsearch(config = cfg, skip = "power"), "skip")
  expect_error(ferx_ruvsearch(config = cfg, groups = 2), "groups")
  expect_error(ferx_ruvsearch(config = cfg, model = "model.ferx"), "model")
  # `data` states *which dataset* the search runs on, so it belongs with
  # `model` rather than with the run knobs.
  expect_error(ferx_ruvsearch(config = cfg, data = "other.csv"), "`data`")

  # The run knobs say *how* to run, not *what* to search, so they stay legal
  # beside a file. This one still fails - the file names a base model that does
  # not exist - but not on the entry form.
  msg <- tryCatch(ferx_ruvsearch(config = cfg, threads = 2),
                  error = conditionMessage)
  expect_false(grepl("cannot be given beside it", msg, fixed = TRUE))
})

test_that("neither entry form is an error naming both", {
  expect_error(ferx_ruvsearch(), "`config`.*or `model`")
})

test_that("a missing file is refused before anything runs", {
  expect_error(ferx_ruvsearch(config = tempfile(fileext = ".ferxsearch")),
               "config file not found")
  expect_error(ferx_ruvsearch(model = tempfile(fileext = ".ferx")),
               "model file not found")
})

test_that("argument validation happens in R, before the engine is called", {
  ex <- ferx_example("one_cpt_transit")
  args <- list(model = ex$model, data = ex$data)

  expect_error(do.call(ferx_ruvsearch, c(args, list(skip = "powerr"))),
               "does not name a residual-error family")
  expect_error(do.call(ferx_ruvsearch, c(args, list(skip = 1))),
               "character vector")
  expect_error(do.call(ferx_ruvsearch, c(args, list(p_value = c(0.01, 0.05)))),
               "single number")
  # An infinite p_value would be dropped on the way to the config file - the
  # binding emits the key only when it is finite - so the search would run at
  # the engine default rather than at the level that was asked for.
  expect_error(do.call(ferx_ruvsearch, c(args, list(p_value = Inf))), "finite")
  expect_error(do.call(ferx_ruvsearch, c(args, list(p_value = 0))), "positive")
  expect_error(do.call(ferx_ruvsearch, c(args, list(max_iter = 1.5))),
               "whole number")
  expect_error(do.call(ferx_ruvsearch, c(args, list(cwres_prescreen = NA))),
               "TRUE, FALSE or NULL")
  expect_error(do.call(ferx_ruvsearch, c(args, list(resume = NA))),
               "TRUE or FALSE")
  expect_error(do.call(ferx_ruvsearch, c(args, list(retries = -1))),
               "at least 0")
})

test_that("the engine's own validation runs before the dataset is read", {
  # These all name a file key rather than an R argument, and all of them are
  # refused while the config is being read - the base model named below does
  # not exist, so anything that got as far as loading it would fail differently.
  expect_error(ferx_ruvsearch(config = write_cfg(
    'base = "model.ferx"', "[ruvsearch]", "groups = 1"
  )), "groups")
  expect_error(ferx_ruvsearch(config = write_cfg(
    'base = "model.ferx"', "[ruvsearch]", "max_iter = 4"
  )), "max_iter")
  # A search space is a different tool's: ruvsearch's candidates are the four
  # residual-error families, not anything MFL can name.
  expect_error(ferx_ruvsearch(config = write_cfg(
    'base = "model.ferx"', "[space]", 'mfl = "PERIPHERALS(0..1)"'
  )), "space")
  # And it selects on the likelihood-ratio test, so a BIC ranking or a dOFV
  # cutoff is refused by name rather than ignored.
  expect_error(ferx_ruvsearch(config = write_cfg(
    'base = "model.ferx"', "[rank]", 'type = "bic"'
  )), "ruvsearch selects by the likelihood-ratio test")
  expect_error(ferx_ruvsearch(config = write_cfg(
    'base = "model.ferx"', "[rank]", "cutoff = 3.84"
  )), "p_value")
})

test_that("an inline groups/max_iter reaches the engine, which names it", {
  # The other half of the sentinel contract: a value R accepts is emitted into
  # the config file, so the engine's own validation is what rejects it. R does
  # not duplicate that rule, it only guarantees the value arrives.
  ex <- ferx_example("one_cpt_transit")
  expect_error(ferx_ruvsearch(model = ex$model, data = ex$data, groups = 1),
               "groups")
})

# -- The end-to-end run ------------------------------------------------------
#
# One real search: slow, so it runs once and every assertion about the result
# object reads that one run.
#
# transit_oral.csv is simulated with a plain proportional error, but the model
# fitted here is `ka`-free while the data was generated with a transit chain
# *and* a first-order `ka` step - so the absorption-phase residuals carry that
# misspecification, and a time-varying magnitude is accepted at p = 0.05.

ruvsearch_run <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    ex <- ferx_example("one_cpt_transit")
    cached <<- ferx_ruvsearch(
      model     = ex$model,
      data      = ex$data,
      p_value   = 0.05,
      max_iter  = 1,
      retries   = 0,
      directory = file.path(tempdir(), "ruvsearch-run"),
      progress  = FALSE
    )
    cached
  }
})

test_that("the result reports the search the engine ran", {
  skip_on_cran()
  res <- ruvsearch_run()

  expect_s3_class(res, "ferx_ruvsearch")
  expect_s3_class(res, "ferx_search_result")

  # The step table is the engine's own, column for column, with the feature's
  # family and its note appended - nothing here is spelled in R.
  expect_equal(names(res$steps),
               c(ferx_rust_ruvsearch_columns(), "family", "note"))
  expect_gt(nrow(res$steps), 1L)
  # Termination status and the gate verdict are columns, never omitted.
  expect_true(is.logical(res$steps$converged))
  expect_true(is.logical(res$steps$passed))
  # The p-value is what this search decides on, so it is a column of numbers.
  tested <- res$steps[res$steps$iteration > 0L, , drop = FALSE]
  expect_true(all(is.finite(tested$p_value)))
  expect_true(all(tested$family %in%
                    c("IIV_on_RUV", "power", "combined", "time_varying")))

  # The input is row zero and is never a candidate.
  expect_equal(res$steps$candidate[1L], "input")
  expect_equal(res$steps$iteration[1L], 0L)
  expect_true(is.na(res$steps$feature[1L]))

  # The options come back as the engine read them, not as R passed them.
  expect_equal(res$options$p_value, 0.05)
  expect_equal(res$options$max_iter, 1L)
  expect_equal(res$options$groups, 4L)
  expect_true(is.finite(res$options$cutoff))

  expect_true(is.finite(res$input_ofv))
  expect_true(is.finite(res$final_ofv))
  # The search cannot end worse than it started: the input is a candidate too.
  expect_lte(res$final_ofv, res$input_ofv)
  expect_true(file.exists(res$final_model_path))
})

test_that("the accepted form is the one the table selected", {
  skip_on_cran()
  res <- ruvsearch_run()

  selected <- res$steps[!is.na(res$steps$selected) & res$steps$selected, ,
                        drop = FALSE]
  expect_equal(nrow(selected), length(res$final_features))
  skip_if(length(res$final_features) == 0L,
          "no form was accepted on this run")

  expect_equal(selected$feature, res$final_features)
  expect_equal(selected$family, res$final_families)
  # Accepted means significant at the level asked for, and past the gate.
  expect_true(all(selected$significant))
  expect_true(all(selected$p_value < res$options$p_value))
  expect_true(all(selected$passed))
  # The winning model carries the accepted form in its error model.
  expect_match(res$final_model, "TAD", fixed = TRUE)
})

test_that("the R step table equals the run's own steps.csv", {
  skip_on_cran()
  res <- ruvsearch_run()
  csv <- read.csv(file.path(res$directory, "steps.csv"), stringsAsFactors = FALSE)

  expect_equal(nrow(csv), nrow(res$steps))
  # Every engine column, in the engine's order; `family` and `note` are added
  # on the end and are the only columns beyond the file's.
  expect_equal(names(res$steps)[seq_along(names(csv))], names(csv))
  expect_equal(setdiff(names(res$steps), names(csv)), c("family", "note"))

  expect_equal(csv$candidate, res$steps$candidate)
  expect_equal(csv$ofv, res$steps$ofv, tolerance = 1e-6)
  expect_equal(csv$p_value, res$steps$p_value, tolerance = 1e-4)
  # `read.csv` types the engine's booleans as text; the R object types them.
  expect_equal(as.logical(csv$selected), res$steps$selected)
  expect_equal(as.logical(csv$passed), res$steps$passed)
})

test_that("ferx_search_results() reads the step table back", {
  skip_on_cran()
  res <- ruvsearch_run()
  tab <- ferx_search_results(res$directory, type = "steps")

  expect_equal(names(tab), ferx_rust_ruvsearch_columns())
  # Both stepwise tools write a `steps.csv`; the header says which one did.
  expect_equal(attr(tab, "tool"), "ruvsearch")
  expect_equal(tab$candidate, res$steps$candidate)
  expect_equal(tab$ofv, res$steps$ofv, tolerance = 1e-6)
  expect_true(is.logical(tab$passed))
  expect_true(is.integer(tab$iteration))
  # The empty cell the engine writes for "no value" is NA, not NaN and not "".
  expect_true(is.na(tab$feature[1L]))
  expect_true(is.na(tab$p_value[1L]))

  # The candidate table is a different table of the same run, and still reads.
  # Each step gets its own runner directory, which is what `$candidates`
  # stacks and labels by step.
  expect_true(is.data.frame(ferx_search_results(file.path(res$directory, "input"))))
  expect_true(all(c("input", "iteration-1") %in% unique(res$candidates$run)))
  expect_error(ferx_search_results(res$directory, type = "steps", partial = TRUE),
               "no partial step table")
  expect_error(ferx_search_results(tempdir(), type = "steps"), "No step table")
})

test_that("a table that is neither tool's is refused rather than mistyped", {
  dir <- file.path(tempdir(), "ruvsearch-not-a-step-table")
  dir.create(dir, showWarnings = FALSE)
  write.csv(data.frame(a = 1, b = 2), file.path(dir, "steps.csv"),
            row.names = FALSE)
  expect_error(ferx_search_results(dir, type = "steps"),
               "not a search step table")
})

test_that("every fitted model's text comes back, named by candidate id", {
  skip_on_cran()
  res <- ruvsearch_run()

  expect_true(all(res$steps$candidate %in% names(res$model_text)))
  expect_true(all(nzchar(res$model_text)))
  # The text is a model file, not a description of one.
  expect_match(res$model_text[["input"]], "[error_model]", fixed = TRUE)
  # And it is what the run wrote beside itself.
  on_disk <- file.path(res$directory, "models", "input.ferx")
  skip_if_not(file.exists(on_disk))
  expect_equal(paste(readLines(on_disk), collapse = "\n"),
               trimws(res$model_text[["input"]], which = "right"))
})

test_that("the winning fit is a ferx_fit, not a refit", {
  skip_on_cran()
  res <- ruvsearch_run()
  skip_if(is.null(res$fit), "the run recovered no final fit")

  expect_s3_class(res$fit, "ferx_fit")
  expect_equal(unname(res$fit$ofv), res$final_ofv, tolerance = 1e-8)
  expect_true(!is.null(res$fit$theta))
  expect_true(is.data.frame(res$fit$sdtab))
})

test_that("both entry forms run the same search", {
  skip_on_cran()
  inline <- ruvsearch_run()

  ex <- ferx_example("one_cpt_transit")
  cfg <- write_cfg(
    sprintf('base = "%s"', ex$model),
    sprintf('data = "%s"', ex$data),
    "[ruvsearch]",
    "p_value = 0.05",
    "max_iter = 1",
    "[run]",
    "retries = 0"
  )
  from_file <- ferx_ruvsearch(config = cfg,
                              directory = file.path(tempdir(), "ruvsearch-file"),
                              progress = FALSE)

  expect_equal(from_file$steps$candidate, inline$steps$candidate)
  expect_equal(from_file$steps$ofv, inline$steps$ofv, tolerance = 1e-6)
  expect_equal(from_file$final_model_id, inline$final_model_id)
  expect_equal(from_file$final_features, inline$final_features)
  expect_equal(from_file$final_ofv, inline$final_ofv, tolerance = 1e-6)
})

test_that("a space with every family skipped returns the input model, fitted", {
  skip_on_cran()
  # The degenerate oracle: with nothing left to test the search has one model
  # to report - the one it was given - and its fit must agree with fitting that
  # model directly.
  ex <- ferx_example("one_cpt_transit")
  res <- ferx_ruvsearch(
    model     = ex$model,
    data      = ex$data,
    skip      = c("IIV_on_RUV", "power", "combined", "time_varying"),
    max_iter  = 1,
    retries   = 0,
    directory = file.path(tempdir(), "ruvsearch-degenerate"),
    progress  = FALSE
  )

  expect_equal(nrow(res$steps), 1L)
  expect_equal(res$final_model_id, "input")
  expect_length(res$final_features, 0L)
  # Nothing was tested, and the run says so rather than reporting an empty
  # iteration as a decision.
  expect_true(any(grepl("no candidate", res$notes)))
  skip_if(is.null(res$fit), "the run recovered no final fit")

  direct <- ferx_fit(ex$model, ex$data)
  expect_equal(unname(res$fit$ofv), unname(direct$ofv), tolerance = 1e-6)
})

test_that("a form the test rejects is a row carrying its verdict", {
  skip_on_cran()
  res <- ruvsearch_run()

  rejected <- res$steps[res$steps$iteration > 0L & !res$steps$selected, ,
                        drop = FALSE]
  expect_gt(nrow(rejected), 0L)
  # Excluded, but present: every rejected form carries its own comparison and
  # its own verdict, which is what a winner-only report would hide.
  expect_true(all(is.finite(rejected$p_value)))
  expect_true(all(!is.na(rejected$significant)))
  expect_true(all(!rejected$significant | rejected$iteration > 0L))
  # A row the strictness gate excluded says why in the same column.
  gated <- rejected[!rejected$passed, , drop = FALSE]
  if (nrow(gated) > 0L) {
    expect_true(all(!is.na(gated$failures) & nzchar(gated$failures)))
    expect_false(any(gated$selected))
  }
})

test_that("a family the method cannot take is a note, not a missing row", {
  skip_on_cran()
  # `iiv_on_ruv` needs eta-epsilon interaction. On a `foce` model it is not
  # tested at all - and the run says so, rather than leaving a candidate
  # unaccounted for. Everything else is skipped here so this costs one fit.
  ex <- ferx_example("warfarin")
  res <- ferx_ruvsearch(
    model     = ex$model,
    data      = ex$data,
    skip      = c("power", "combined", "time_varying"),
    max_iter  = 1,
    retries   = 0,
    directory = file.path(tempdir(), "ruvsearch-nointeraction"),
    progress  = FALSE
  )

  expect_true(any(grepl("IIV_on_RUV not tested", res$notes)))
  expect_false(any(res$steps$family %in% "IIV_on_RUV", na.rm = TRUE))
  expect_length(res$final_features, 0L)
})

test_that("the CWRES pre-screen's rows are marked as what they are", {
  skip_on_cran()
  # Pharmpy's cheap path: every candidate is fitted to the parent's CWRES and
  # only the winner is refitted to the data. Those rows carry an OFV on the
  # CWRES scale, so they have to be distinguishable from the fits that decide
  # the search - otherwise the table holds two objective functions in one
  # column with nothing to tell them apart.
  ex <- ferx_example("one_cpt_transit")
  res <- ferx_ruvsearch(
    model           = ex$model,
    data            = ex$data,
    p_value         = 0.05,
    max_iter        = 1,
    retries         = 0,
    cwres_prescreen = TRUE,
    directory       = file.path(tempdir(), "ruvsearch-prescreen"),
    progress        = FALSE
  )

  expect_true(res$options$cwres_prescreen)
  screened <- res$steps[res$steps$screened, , drop = FALSE]
  expect_gt(nrow(screened), 0L)
  # A screened row is compared on the CWRES scale, so it carries `cwres_dofv`
  # and no likelihood-ratio test at all.
  expect_true(all(is.na(screened$p_value)))
  expect_true(all(is.na(screened$dofv)))
  expect_true(any(is.finite(screened$cwres_dofv)))
  expect_true(all(grepl("^cwres-", screened$candidate)))

  # Only the feature the screen picked was refitted to the data, and that fit
  # is what the likelihood-ratio test decided on.
  refit <- res$steps[res$steps$iteration > 0L & !res$steps$screened, ,
                     drop = FALSE]
  expect_lte(nrow(refit), 1L)
  skip_if(nrow(refit) == 0L, "the pre-screen found nothing to refit")
  expect_true(is.finite(refit$p_value))
  expect_equal(res$final_features, refit$feature[refit$selected])

  # print() says which rows are the screen's, and does not count them as
  # candidates the strictness gate excluded.
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "CWRES pre-screen fit")
  expect_match(out, "cwres_dofv", fixed = TRUE)
  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "CWRES pre-screen only")
})

test_that("print() shows the iteration table and summary() the rest", {
  skip_on_cran()
  res <- ruvsearch_run()
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "ferx residual-error search")
  expect_match(out, "Candidates by iteration")
  expect_match(out, "Selected residual-error model")
  # The selected form is named, and so is the error model it produced.
  if (length(res$final_features) > 0L) {
    expect_match(out, res$final_features[1L], fixed = TRUE)
    expect_match(out, "DV ~", fixed = TRUE)
  }

  full <- paste(capture.output(summary(res)), collapse = "\n")
  expect_match(full, "Forms not selected")
  expect_match(full, "Every model fitted")
})
