# `ferx_fit(method = "vi")` (#506). The method token used to be refused twice
# before it reached the engine - by an R-side `match.arg()` list and by a copy
# of the engine's token parser in the glue - neither of which had "vi". Both
# now defer to ferx-core's own reader, so VI runs from the call as it does from
# a model file, and so do its documented finishing chains.

vi_settings <- function(...) list(vi_iters = 30L, vi_seed = 42L, ...)

test_that("ferx_fit(method = \"vi\") runs on warfarin", {
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(
    ex$model, ex$data, method = "vi", covariance = FALSE, verbose = FALSE,
    settings = vi_settings()
  ))
  expect_identical(fit$method, "VI")
  # VI reports a bound, not an objective: NaN OFV under `vi_final_ofv = none`.
  expect_true(is.nan(fit$ofv))
  expect_true(is.finite(fit$vi$neg_two_elbo))
  expect_true(all(is.finite(fit$theta)))
})

test_that("an explicit method = \"vi\" gets the interaction flag a model-file vi does", {
  # VI has no R-side history (R refused it until #506), so it follows the
  # model file's rule, not the legacy R one (off for any final stage but
  # FOCEI) - and it overrides the file's own method line: warfarin's
  # `method = foce` clears the flag, which must not carry over to VI.
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(
    ex$model, ex$data, method = "vi", covariance = FALSE, verbose = FALSE,
    settings = vi_settings()
  ))
  expect_true(fit$interaction)
})

test_that("a method alias in the call does not warn as overriding the same method", {
  warn_conflicts <- getFromNamespace(".ferx_warn_fit_option_conflicts", "ferx")
  file_opts <- list(method = "[focei, imp]")
  expect_no_warning(warn_conflicts(
    file_opts, list(method = "focei, importance_sampling"), list()
  ))
  expect_no_warning(warn_conflicts(
    list(method = "laplace"), list(method = "laplacian"), list()
  ))
  expect_warning(
    warn_conflicts(file_opts, list(method = "vi"), list()),
    "overrides it with `vi`"
  )
  # A model-file value the engine does not read still compares as text.
  expect_warning(
    warn_conflicts(list(method = "nonsense"), list(method = "vi"), list()),
    "sets `method = nonsense`"
  )
})

test_that("c(\"vi\", \"laplace\") with agq_eval_only = TRUE scores the VI estimate", {
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(
    ex$model, ex$data, method = c("vi", "laplace"), covariance = FALSE,
    verbose = FALSE, settings = vi_settings(agq_eval_only = TRUE, n_agq = 3L)
  ))
  # An evaluate-only readout does not estimate: the fit stays a VI fit, and the
  # quadrature -2 log L at the VI estimate is its OFV.
  expect_identical(fit$method_chain, c("VI", "LAPLACE"))
  expect_identical(fit$method, "VI")
  expect_true(is.finite(fit$ofv))
  expect_false(is.null(fit$vi))
})

test_that("c(\"vi\", \"imp\") with imp_eval_only = TRUE scores the VI estimate", {
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(
    ex$model, ex$data, method = c("vi", "imp"), covariance = FALSE,
    verbose = FALSE, settings = vi_settings(imp_eval_only = TRUE)
  ))
  # The IMP readout's -2 log L at the VI estimate is on
  # fit$importance_sampling; the fit itself stays a VI fit.
  expect_identical(fit$method_chain, c("VI", "IMP"))
  expect_identical(fit$method, "VI")
  expect_true(is.finite(fit$importance_sampling$minus2_log_likelihood))
  expect_false(is.null(fit$vi))
})

test_that("an evaluate-only stage that is not last is refused, as in a model file", {
  ex <- ferx_example("warfarin")
  expect_error(
    suppressWarnings(ferx_fit(
      ex$model, ex$data, method = c("vi", "laplace", "focei"),
      covariance = FALSE, verbose = FALSE,
      settings = vi_settings(agq_eval_only = TRUE)
    )),
    "agq_eval_only.*final stage|final stage.*agq_eval_only"
  )
})

test_that("a model-file `method = vi` runs with method = NULL", {
  ex <- ferx_example("warfarin")
  src <- readLines(ex$model)
  src <- src[!grepl("^\\s*(method|covariance)\\s*=", src)]
  at <- grep("^\\[fit_options\\]", src)
  src <- append(
    src,
    c("  method     = vi", "  vi_iters   = 30", "  vi_seed    = 42",
      "  covariance = false"),
    after = at
  )
  model <- tempfile(fileext = ".ferx")
  writeLines(src, model)
  fit <- suppressWarnings(ferx_fit(model, ex$data, verbose = FALSE))
  expect_identical(fit$method, "VI")
  expect_true(fit$interaction)
  expect_true(is.nan(fit$ofv))
  expect_false(is.null(fit$vi))
})
