# fit$status / fit$status_reasons (ferx-core #2799): the usability verdict beside
# `converged`, read from core's `FitResult::status()` and never recomputed in R.

status_tokens <- c("converged", "converged_with_issues", "not_converged")

status_line <- function(fit) {
  grep("STATUS:", capture.output(print(fit)), value = TRUE)
}

# `key` removed wherever it occurs in a parsed fit.json, so a bundle reads as one
# written before the key existed.
drop_key <- function(x, key) {
  if (!is.list(x)) return(x)
  x[[key]] <- NULL
  lapply(x, drop_key, key = key)
}

eval_only_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                       covariance = FALSE, settings = list(maxiter = 0L))
    }
    fit
  }
})

test_that("a converged fit reports status `converged` with no reasons", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  expect_true(fit$converged)
  expect_identical(fit$status, "converged")
  expect_identical(fit$status_label, "CONVERGED")
  expect_identical(fit$status_reasons, character(0))
  expect_false(fit$eval_budget_spent)
  expect_false(fit$evaluation_only)
  line <- status_line(fit)
  expect_length(line, 1L)
  expect_match(line, "STATUS: CONVERGED", fixed = TRUE)
  expect_false(grepl("(", line, fixed = TRUE))
})

test_that("a maxiter = 0 fit is `not_converged` with reason evaluation_only", {
  skip_on_cran()
  fit <- eval_only_fit()
  expect_false(fit$converged)
  expect_true(fit$evaluation_only)
  expect_identical(fit$status, "not_converged")
  expect_identical(fit$status_label, "NOT CONVERGED")
  # Core's fixed order puts evaluation_only first.
  expect_identical(fit$status_reasons[1], "evaluation_only")
  line <- status_line(fit)
  expect_match(line, "STATUS: NOT CONVERGED   (evaluation_only", fixed = TRUE)
  # The two halves of the verdict: the converged fit above prints no reason, so
  # the reason text here comes from `status_reasons`, not from the print method.
  expect_false(grepl("evaluation_only", status_line(warfarin_fit_cov()), fixed = TRUE))
})

test_that("status and its inputs survive a save/load round-trip", {
  skip_on_cran()
  for (fit in list(warfarin_fit_cov(), eval_only_fit())) {
    path <- tempfile(fileext = ".fitrx")
    on.exit(unlink(path), add = TRUE)
    ferx_save_fit(fit, path)
    back <- ferx_load_fit(path)
    expect_identical(back$status, fit$status)
    expect_identical(back$status_label, fit$status_label)
    expect_identical(back$status_reasons, fit$status_reasons)
    expect_identical(back$eval_budget_spent, fit$eval_budget_spent)
    expect_identical(back$evaluation_only, fit$evaluation_only)
  }
})

test_that("a bundle written before #2799 loads with both inputs FALSE", {
  skip_on_cran()
  fit <- eval_only_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  writeBin(fit$bundle, path)
  # Precondition: the unedited bundle carries the flag, so the edit below is
  # what drops `evaluation_only` from the reasons.
  expect_true(ferx_load_fit(path)$evaluation_only)
  old <- fitrx_edit_bundle(path, function(w) {
    drop_key(drop_key(w, "evaluation_only"), "eval_budget_spent")
  })
  on.exit(unlink(old), add = TRUE)
  back <- ferx_load_fit(old)
  expect_identical(back$evaluation_only, FALSE)
  expect_identical(back$eval_budget_spent, FALSE)
  # Without the flag core no longer knows the run was evaluation-only; it is
  # still not converged, so the verdict stays `not_converged`.
  expect_identical(back$status, "not_converged")
  expect_false("evaluation_only" %in% back$status_reasons)
})

test_that("a list without a status falls back to `converged` in the STATUS line", {
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1), converged = TRUE)
  fit[["status"]] <- NULL
  expect_match(status_line(fit), "STATUS: CONVERGED", fixed = TRUE)
})

test_that("the STATUS line prints core's label and reasons as given", {
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1), converged = TRUE)
  fit$status <- "converged_with_issues"
  fit$status_label <- "CONVERGED WITH ISSUES"
  fit$status_reasons <- c("condition_number", "covariance_failed")
  expect_match(status_line(fit),
               "STATUS: CONVERGED WITH ISSUES   (condition_number, covariance_failed)",
               fixed = TRUE)
  expect_true(all(fit$status %in% status_tokens))
})
