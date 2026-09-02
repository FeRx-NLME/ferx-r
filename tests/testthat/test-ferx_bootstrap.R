# Tests for ferx_bootstrap() / ferx_bootstrap_summarize().
#
# The statistics themselves are the engine's (ferx-tools), and are tested there.
# What is checked here is the wrapper: argument validation, the shape of the
# returned object, that the run is reproducible from its seed, that the disk
# artefacts round-trip through the summarize path, and the print/plot methods.
#
# Runs use a handful of samples on the bundled warfarin example, which fits in
# well under a second - enough to exercise every path without a slow test file.

boot_small <- local({
  bs <- NULL
  function() {
    if (is.null(bs)) {
      ex <- ferx_example("warfarin")
      bs <<- ferx_bootstrap(ex$model, ex$data, samples = 4L, seed = 42,
                            threads = 2L)
    }
    bs
  }
})

test_that("ferx_bootstrap returns the documented structure", {
  bs <- boot_small()

  expect_s3_class(bs, "ferx_bootstrap")
  expect_true(all(c("parameters", "raw", "diagnostics", "parameter_names",
                    "n_completed", "n_included", "chi_square_df",
                    "confidence_level") %in% names(bs)))

  expect_s3_class(bs$parameters, "data.frame")
  expect_identical(
    names(bs$parameters),
    c("parameter", "original", "mean", "bias", "standard_error", "median",
      "ci_lower", "ci_upper", "ci_lower_normal", "ci_upper_normal")
  )
  expect_gt(nrow(bs$parameters), 0L)
  expect_identical(bs$parameters$parameter, bs$parameter_names)

  # One raw row per fit, the original dataset first.
  expect_s3_class(bs$raw, "data.frame")
  expect_identical(nrow(bs$raw), 5L)
  expect_identical(bs$raw$sample, 0:4)
  expect_true(all(c("minimization_successful", "estimate_near_boundary",
                    "covariance_step_successful", "covariance_step_warnings",
                    "ofv", "seconds", "error") %in% names(bs$raw)))
  expect_true(all(bs$parameter_names %in% names(bs$raw)))
  expect_true(is.logical(bs$raw$minimization_successful))
  # No covariance step for the replicates by default, so no se_ columns.
  expect_false(any(startsWith(names(bs$raw), "se_")))
  # "" from the CSV convention is NA on the R side.
  expect_true(all(is.na(bs$raw$error)))

  expect_s3_class(bs$diagnostics, "data.frame")
  expect_identical(names(bs$diagnostics), c("statistic", "value"))
  expect_true(all(c("samples_requested", "samples_completed",
                    "samples_included", "chi_square_df") %in%
                    bs$diagnostics$statistic))

  # dofv is off by default, and nothing is written without a `directory`.
  expect_null(bs$delta_ofv)
  expect_null(bs$directory)
  expect_identical(bs$confidence_level, 95)
  expect_identical(bs$n_completed, 4L)
})

test_that("the original dataset's row carries the base fit", {
  bs <- boot_small()
  original <- bs$raw[bs$raw$sample == 0L, ]
  expect_identical(nrow(original), 1L)
  # `original` in the parameter table is that row's estimate, by construction.
  for (p in bs$parameter_names) {
    expect_equal(bs$parameters$original[bs$parameters$parameter == p],
                 original[[p]], tolerance = 1e-8)
  }
})

test_that("a run is reproducible from its seed and moves with it", {
  ex <- ferx_example("warfarin")
  a <- ferx_bootstrap(ex$model, ex$data, samples = 3L, seed = 99, threads = 2L)
  b <- ferx_bootstrap(ex$model, ex$data, samples = 3L, seed = 99, threads = 1L)
  c <- ferx_bootstrap(ex$model, ex$data, samples = 3L, seed = 100, threads = 2L)

  # Same seed, different thread count: the draw is derived from the seed and
  # the replicate index, so it must not depend on the schedule.
  expect_equal(a$parameters$mean, b$parameters$mean, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(a$parameters$mean, c$parameters$mean)))
})

test_that("directory writes the artefacts and summarize reads them back", {
  ex <- ferx_example("warfarin")
  d <- file.path(tempdir(), "ferx-bs-dir")
  unlink(d, recursive = TRUE)

  bs <- ferx_bootstrap(ex$model, ex$data, samples = 4L, seed = 7, threads = 2L,
                       dofv = TRUE, directory = d, ci = 90)
  expect_identical(bs$directory, normalizePath(d))
  expect_true(all(c("raw_results.csv", "bootstrap_results.csv",
                    "bootstrap_diagnostics.csv", "delta_ofv.csv") %in%
                    list.files(d)))
  expect_s3_class(bs$delta_ofv, "data.frame")
  expect_identical(names(bs$delta_ofv), c("sample", "delta_ofv"))
  expect_identical(nrow(bs$delta_ofv), 4L)
  expect_identical(bs$confidence_level, 90)

  again <- ferx_bootstrap_summarize(d, ci = 90)
  expect_s3_class(again, "ferx_bootstrap")
  expect_equal(again$parameters$mean, bs$parameters$mean, tolerance = 1e-8)
  expect_identical(again$n_included, bs$n_included)
  # `$raw` is read back off disk, delta_ofv included.
  expect_identical(nrow(again$raw), nrow(bs$raw))
  expect_identical(nrow(again$delta_ofv), nrow(bs$delta_ofv))

  # Relaxing a filter cannot drop samples, and the counts move with it.
  relaxed <- ferx_bootstrap_summarize(
    d, skip_minimization_terminated = FALSE,
    skip_estimate_near_boundary = FALSE, ci = 90
  )
  expect_gte(relaxed$n_included, again$n_included)
})

test_that("ferx_bootstrap_summarize rejects a directory that is not a run", {
  d <- file.path(tempdir(), "ferx-bs-empty")
  dir.create(d, showWarnings = FALSE)
  expect_error(ferx_bootstrap_summarize(d), "no raw_results.csv")
  expect_error(ferx_bootstrap_summarize(file.path(tempdir(), "nope-1144")),
               "not found")
})

test_that("ferx_bootstrap_summarize validates its flags before touching disk", {
  d <- file.path(tempdir(), "ferx-bs-empty")
  dir.create(d, showWarnings = FALSE)
  for (nm in c("skip_minimization_terminated", "skip_estimate_near_boundary",
               "skip_covariance_step_terminated", "skip_with_covstep_warnings")) {
    args <- list(d); args[[nm]] <- NA
    expect_error(do.call(ferx_bootstrap_summarize, args), "TRUE or FALSE")
    args[[nm]] <- "TRUE"
    expect_error(do.call(ferx_bootstrap_summarize, args), "TRUE or FALSE")
  }
  # A bad flag is reported even when the directory itself is unusable: the
  # argument error is the one the caller can act on.
  expect_error(ferx_bootstrap_summarize(d, ci = 100), "confidence level")
})

test_that("stratified resampling and sample_size reach the engine", {
  ex <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  ids <- unique(dat$ID)
  dat$STUDY <- ifelse(dat$ID %in% ids[seq_len(length(ids) %/% 2L)], 1001, 1002)
  path <- file.path(tempdir(), "warfarin-strata.csv")
  write.csv(dat, path, row.names = FALSE)

  expect_s3_class(
    ferx_bootstrap(ex$model, path, samples = 2L, seed = 3, threads = 2L,
                   stratify_on = "STUDY"),
    "ferx_bootstrap"
  )
  # PsN's `-sample_size=1001=>3,1002=>4`, spelled as a named vector.
  expect_s3_class(
    ferx_bootstrap(ex$model, path, samples = 2L, seed = 3, threads = 2L,
                   stratify_on = "STUDY",
                   sample_size = c("1001" = 3, "1002" = 4)),
    "ferx_bootstrap"
  )
  expect_s3_class(
    ferx_bootstrap(ex$model, path, samples = 2L, seed = 3, threads = 2L,
                   sample_size = 8),
    "ferx_bootstrap"
  )
  expect_error(
    ferx_bootstrap(ex$model, path, samples = 2L, stratify_on = "NOPE"),
    "NOPE"
  )
})

test_that("argument validation happens before anything is fitted", {
  ex <- ferx_example("warfarin")

  expect_error(ferx_bootstrap("no-such-model.ferx"), "Model file not found")
  expect_error(ferx_bootstrap(ex$model, "no-such-data.csv"),
               "Data file not found")
  expect_error(ferx_bootstrap(ex$model, ex$data, samples = 0),
               "whole number >= 1")
  expect_error(ferx_bootstrap(ex$model, ex$data, ci = 100),
               "confidence level")
  expect_error(ferx_bootstrap(ex$model, ex$data, threads = 0),
               "whole number >= 1")
  expect_error(ferx_bootstrap(ex$model, ex$data, update_inits = TRUE,
                              run_base_model = FALSE),
               "run_base_model")
  expect_error(ferx_bootstrap(ex$model, ex$data, dofv = NA), "TRUE or FALSE")
  expect_error(ferx_bootstrap(ex$model, ex$data, progress = NA), "TRUE or FALSE")
  expect_error(ferx_bootstrap(ex$model, ex$data, stratify_on = NA),
               "single column name")
  expect_error(ferx_bootstrap(ex$model, ex$data, sample_size = c(3, 4)),
               "single number")
  expect_error(ferx_bootstrap(ex$model, ex$data, sample_size = c("1001" = 12)),
               "needs `stratify_on`")
  # The covariance-step filters read a diagnostic that only exists when the
  # step ran; the engine says so rather than dropping the filter silently.
  expect_error(
    ferx_bootstrap(ex$model, ex$data, samples = 2L,
                   skip_covariance_step_terminated = TRUE),
    "keep-covariance"
  )
})

test_that("sample_size shapes map onto the engine's spec", {
  expect_identical(.ferx_bootstrap_sample_size(NULL),
                   list(keys = character(0), values = numeric(0)))
  expect_identical(.ferx_bootstrap_sample_size(12),
                   list(keys = character(0), values = 12))
  expect_identical(
    .ferx_bootstrap_sample_size(c("1001" = 12, "1002" = 24)),
    list(keys = c("1001", "1002"), values = c(12, 24))
  )
  expect_error(.ferx_bootstrap_sample_size("12"), "must be a number")
  expect_error(.ferx_bootstrap_sample_size(c(a = 1, 2)), "every element")
  # 0 reaches the engine as SampleSize::Total(0) - replicates drawing no
  # subjects at all, which surfaces as a fit error far from its cause.
  expect_error(.ferx_bootstrap_sample_size(0), "whole numbers >= 1")
  expect_error(.ferx_bootstrap_sample_size(c("1001" = 12, "1002" = 0)),
               "whole numbers >= 1")
  expect_error(.ferx_bootstrap_sample_size(2.5), "whole numbers >= 1")
})

test_that("keep_covariance adds the per-replicate standard errors", {
  ex <- ferx_example("warfarin")
  bs <- ferx_bootstrap(ex$model, ex$data, samples = 2L, seed = 5, threads = 2L,
                       keep_covariance = TRUE)
  expect_true(all(paste0("se_", bs$parameter_names) %in% names(bs$raw)))
})

test_that("`$raw` has the same column types from either entry point", {
  ex <- ferx_example("warfarin")
  d <- file.path(tempdir(), "ferx-bs-types")
  unlink(d, recursive = TRUE)
  bs <- ferx_bootstrap(ex$model, ex$data, samples = 3L, seed = 11, threads = 2L,
                       directory = d)
  again <- ferx_bootstrap_summarize(d)

  # raw_results.csv is untyped text: the flags are written 0/1 and `error` is
  # blank for a replicate that succeeded, so read.csv would hand back integer
  # flags and - with no replicate errored - a logical NA `error` column.
  flags <- c("minimization_successful", "estimate_near_boundary",
             "covariance_step_successful", "covariance_step_warnings")
  for (nm in flags) {
    expect_type(bs$raw[[nm]], "logical")
    expect_type(again$raw[[nm]], "logical")
    expect_identical(again$raw[[nm]], bs$raw[[nm]])
  }
  expect_type(bs$raw$error, "character")
  expect_type(again$raw$error, "character")
  expect_identical(is.na(again$raw$error), is.na(bs$raw$error))
})

test_that("plot survives a parameter with a single finite estimate", {
  bs <- boot_small()
  one <- bs
  # Everything but the first replicate lost its estimate: `breaks = "FD"` needs
  # a spread to derive a bin width, and hist() aborts without one.
  p1 <- bs$parameter_names[1]
  keep <- one$raw$sample > 0L
  one$raw[[p1]][keep][-1] <- NA_real_
  one$delta_ofv <- data.frame(sample = 1L, delta_ofv = 3.1)
  one$chi_square_df <- 1L

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(one, parameters = p1), p1)

  # And with no finite estimate at all, the panel is a placeholder.
  none <- one
  none$raw[[p1]] <- NA_real_
  none$delta_ofv <- NULL
  expect_identical(plot(none, parameters = p1), p1)
})

test_that("print and plot work on a bootstrap result", {
  bs <- boot_small()
  out <- capture.output(print(bs))
  expect_true(any(grepl("ferx bootstrap", out, fixed = TRUE)))
  expect_true(any(grepl("included", out, fixed = TRUE)))
  expect_true(any(grepl(bs$parameter_names[1], out, fixed = TRUE)))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(bs), bs$parameter_names)
  expect_identical(plot(bs, parameters = bs$parameter_names[1]),
                   bs$parameter_names[1])
  expect_error(plot(bs, parameters = "NOPE"), "Unknown parameter")
  expect_error(plot(bs, y = 1), "does not use `y`")
})

# -- progress ---------------------------------------------------------------

test_that("a watched run reports its fits and returns the same numbers", {
  skip_if_not_installed("mockery")

  ex <- ferx_example("warfarin")

  seen <- list()
  handler <- function(stage, completed, total, base_fit) {
    seen[[length(seen) + 1L]] <<- list(stage = stage, completed = completed,
                                       total = total, base_fit = base_fit)
    invisible(NULL)
  }
  mockery::stub(ferx_bootstrap, ".ferx_bootstrap_progress_handler",
                function(...) handler)

  watched <- ferx_bootstrap(ex$model, ex$data, samples = 4L, seed = 42,
                            threads = 2L, progress = TRUE)
  quiet <- ferx_bootstrap(ex$model, ex$data, samples = 4L, seed = 42,
                          threads = 2L, progress = FALSE)

  # The events are drawn on a timer, so which ones arrive is a race - but the
  # last one always is, and it is the one that clears the bar.
  expect_gt(length(seen), 0L)
  expect_identical(seen[[length(seen)]]$stage, "finished")
  stages <- vapply(seen, function(e) e$stage, character(1))
  expect_true(all(stages %in% c("started", "base_done", "replicate", "dofv",
                                "finished")))

  # Watching a run must not change it: same seed, same estimates.
  expect_equal(watched$parameters, quiet$parameters)
})

test_that("the progress handler tolerates a coalesced event stream", {
  # Events are coalesced on the way over, so a stage may arrive without the one
  # that would normally have opened its bar. Each has to stand on its own.
  h <- .ferx_bootstrap_progress_handler(environment())
  expect_no_error(h("replicate", 3L, 10L, FALSE))
  expect_no_error(h("dofv", 1L, 10L, FALSE))
  expect_no_error(h("finished", 0L, 0L, FALSE))
  # And a second "finished" - closing a bar that is already closed is what a
  # duplicated final event does.
  expect_no_error(h("finished", 0L, 0L, FALSE))
})

test_that("the progress handler falls back to a text bar without cli", {
  skip_if_not_installed("mockery")

  mockery::stub(.ferx_bootstrap_progress_handler, "requireNamespace", FALSE)
  h <- .ferx_bootstrap_progress_handler(environment())

  out <- capture.output(
    suppressMessages({
      h("started", 0L, 4L, FALSE)
      h("replicate", 2L, 4L, FALSE)
      h("finished", 0L, 0L, FALSE)
    })
  )
  expect_true(any(grepl("50%", out, fixed = TRUE)))
})

test_that("the handler holds the frame it was made in", {
  # The bars are drawn long after the factory returns, so the `envir` default
  # has to be forced while its caller is still on the stack - cli otherwise
  # ties the bar to the global environment and nothing ever closes it.
  captured <- NULL
  make <- function() {
    captured <<- environment()
    .ferx_bootstrap_progress_handler()
  }
  h <- make()
  expect_identical(environment(h)$envir, captured)
})

test_that("a bar opened without a total is reopened once the total arrives", {
  skip_if_not_installed("mockery")

  # "started" and "base_done" can coalesce into one poll window, leaving the
  # replicate bar opened against an unknown total. The first counted event has
  # to correct it, or the run draws no progress at all.
  mockery::stub(.ferx_bootstrap_progress_handler, "requireNamespace", FALSE)
  h <- .ferx_bootstrap_progress_handler(environment())

  out <- capture.output(
    suppressMessages({
      h("base_done", 0L, 0L, FALSE)
      h("replicate", 2L, 4L, FALSE)
      h("finished", 0L, 0L, FALSE)
    })
  )
  expect_true(any(grepl("50%", out, fixed = TRUE)))
})

test_that("a fully resumed run does not open a zero-width text bar", {
  skip_if_not_installed("mockery")

  # A `directory` resume holding every replicate reports 0 left to fit, and
  # txtProgressBar() stops on `max <= min`.
  mockery::stub(.ferx_bootstrap_progress_handler, "requireNamespace", FALSE)
  h <- .ferx_bootstrap_progress_handler(environment())

  msgs <- capture_messages({
    h("started", 4L, 0L, FALSE)
    h("finished", 0L, 0L, FALSE)
  })
  expect_true(any(grepl("Reusing 4 replicates", msgs, fixed = TRUE)))
  expect_true(any(grepl("Bootstrap replicates", msgs, fixed = TRUE)))
})

test_that("a handler error cannot take the run down with it", {
  # The engine drops what the handler throws, but the handler is ours: it must
  # not throw in the first place, whatever it is handed.
  h <- .ferx_bootstrap_progress_handler(environment())
  expect_no_error(h("nonsense", NA_integer_, NA_integer_, NA))
})
