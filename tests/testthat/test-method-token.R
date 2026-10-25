# `ferx_fit(method = ...)` takes the engine's method tokens (#506). The R side
# and the glue used to keep their own method lists, neither of which had
# `"vi"`, so `method = "vi"` was refused twice before it reached the engine.

vi_settings <- list(vi_iters = 30, vi_seed = 42)

vi_call_fit <- function(method, settings = list()) {
  ex <- ferx_example("warfarin")
  suppressWarnings(ferx_fit(
    ex$model, ex$data,
    method = method, covariance = FALSE, verbose = FALSE,
    settings = c(vi_settings, settings)
  ))
}

test_that("the engine's method parser reads `vi`", {
  labels <- getFromNamespace(".ferx_method_labels", "ferx")
  expect_identical(labels("vi"), "VI")
  expect_identical(labels("VI"), "VI")
  expect_identical(labels(c("vi", "laplace")), c("VI", "LAPLACE"))
  expect_identical(labels(c("focei", "vi")), c("FOCEI", "VI"))
  expect_identical(labels(character(0L)), character(0L))
})

test_that("an unknown method token is refused with the engine's message", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, method = "advi"),
    "unknown estimation method.*advi"
  )
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "nope")),
    "unknown estimation method.*nope"
  )
})

test_that("ferx_fit(method = \"vi\") runs on warfarin", {
  fit <- vi_call_fit("vi")
  expect_s3_class(fit, "ferx_fit")
  expect_identical(fit$method, "VI")
  expect_false(is.null(fit$vi))
  # A bound, not a likelihood: no OFV under the default `vi_final_ofv`.
  expect_true(is.nan(fit$ofv))
})

test_that("c(\"vi\", \"laplace\") with agq_eval_only finishes on a -2 log L", {
  fit <- vi_call_fit(c("vi", "laplace"), list(agq_eval_only = TRUE))
  expect_s3_class(fit, "ferx_fit")
  expect_identical(toupper(fit$method_chain), c("VI", "LAPLACE"))
  expect_true(is.finite(fit$ofv))
})

test_that("c(\"vi\", \"imp\") with imp_eval_only runs", {
  fit <- vi_call_fit(c("vi", "imp"), list(imp_eval_only = TRUE))
  expect_s3_class(fit, "ferx_fit")
  expect_identical(toupper(fit$method_chain), c("VI", "IMP"))
  # Not asserted: a finite `fit$ofv`. At ferx-core 5bd862b8 this chain reports
  # NaN whether it is selected here or by a model file's `method = [vi, imp]`,
  # so that is the engine's result, not this argument's.
})

test_that("an evaluate-only stage that is not last is refused", {
  ex <- ferx_example("warfarin")
  # R-side for `imp_eval_only` passed through `settings`.
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("imp", "vi"),
             settings = list(imp_eval_only = TRUE)),
    "final stage"
  )
  # The engine's own chain check for `agq_eval_only`, as for a model file.
  expect_error(
    suppressWarnings(ferx_fit(
      ex$model, ex$data, method = c("laplace", "vi"),
      covariance = FALSE, verbose = FALSE,
      settings = c(vi_settings, list(agq_eval_only = TRUE))
    )),
    "agq_eval_only"
  )
})

# The warfarin model with its own `method` line replaced by `method_line`
# (dropped when NULL) and the VI settings `vi_call_fit()` passes.
vi_model_file <- function(method_line, extra = character(0L)) {
  ex <- ferx_example("warfarin")
  src <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", src)
  opts <- src[-seq_len(at)]
  opts <- opts[!grepl("^\\s*(method|covariance)\\s*=", opts)]
  model <- tempfile(fileext = ".ferx")
  writeLines(c(
    src[seq_len(at - 1L)], "[fit_options]", opts,
    if (!is.null(method_line)) paste0("  method     = ", method_line),
    "  vi_iters   = 30",
    "  vi_seed    = 42",
    extra,
    "  covariance = false"
  ), model)
  model
}

test_that("a model-file `method = vi` runs with `method = NULL`", {
  ex <- ferx_example("warfarin")
  model <- vi_model_file("vi")
  on.exit(unlink(model), add = TRUE)
  fit <- suppressWarnings(ferx_fit(model, ex$data, method = NULL, verbose = FALSE))
  expect_identical(fit$method, "VI")
  expect_false(is.null(fit$vi))
})

test_that("a call-time method fits as the same model-file method line does", {
  # The call's `method` takes the place of the file's `method =` line, so the
  # two routes must leave the engine the same options - `interaction` included,
  # which core sets from the final stage (on for focei, off for foce, the
  # default otherwise) and which scores CWRES, SIR and the reported objective.
  # The warfarin model's own `method = foce` must not leak into the call route.
  ex <- ferx_example("warfarin")
  cases <- list(
    list(r = "vi", file = "vi", extra = character(0L)),
    list(r = c("vi", "laplace"), file = "[vi, laplace]",
         extra = "  agq_eval_only = true"),
    list(r = c("vi", "foce"), file = "[vi, foce]", extra = character(0L))
  )
  for (cs in cases) {
    model <- vi_model_file(cs$file, cs$extra)
    by_file <- suppressWarnings(ferx_fit(model, ex$data, verbose = FALSE))
    unlink(model)
    settings <- if (length(cs$extra)) list(agq_eval_only = TRUE) else list()
    by_call <- vi_call_fit(cs$r, settings)
    label <- paste(cs$r, collapse = ", ")
    expect_identical(by_call$interaction, by_file$interaction, info = label)
    expect_identical(by_call$method_chain, by_file$method_chain, info = label)
    expect_identical(by_call$ofv, by_file$ofv, info = label)
  }
  expect_true(vi_call_fit("vi")$interaction)
  expect_false(vi_call_fit(c("vi", "foce"))$interaction)
})

test_that("one `method` element cannot smuggle in a second stage", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, method = "gn_hybrid, saem"),
    "not one method token|unknown estimation method"
  )
})

test_that("the glue keeps no method parser of its own", {
  # A copy of the engine's token list in the glue is what let `vi` go missing;
  # every token has to go through ferx-core's `parse_method_token`. The package
  # source is not there for an installed package or under `R CMD check`.
  p <- file.path("..", "..", "src", "rust", "src", "lib.rs")
  skip_if_not(file.exists(p), "package source not available")
  src <- readLines(p, warn = FALSE)
  code <- src[!grepl("^\\s*//", src)]
  expect_false(any(grepl("\\bfn\\s+parse_method\\b", code)))
  # Aliases only a hand-written method parser would spell.
  expect_false(any(grepl(
    "\"(importance_sampling|importance-sampling|laplacian|bayesian|mcmc|foce-i)\"",
    code
  )))
  expect_true(any(grepl("parse_method_token", code, fixed = TRUE)))
})
