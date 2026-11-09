# fit$status / fit$status_reasons (ferx-core #2799): the usability verdict beside
# `converged`, read from core's `FitResult::status()` and never recomputed in R.

quiet <- function(x) suppressWarnings(suppressMessages(x))

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
  expect_true(fit$status %in% status_tokens)
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
  expect_true(fit$status %in% status_tokens)
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

test_that("a status without its label reads `converged` alone", {
  # A list that carries the token but not the label was not built by the glue;
  # R does not spell the label itself, so it falls back as a status-less list.
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1), converged = FALSE,
                       status = "converged_with_issues",
                       status_reasons = "covariance_failed")
  line <- status_line(fit)
  expect_match(line, "STATUS: NOT CONVERGED", fixed = TRUE)
  expect_false(grepl("(", line, fixed = TRUE))
})

test_that("summary() and ferx_runlog() show the verdict", {
  skip_on_cran()
  fit <- eval_only_fit()
  s <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl(paste0("^Status:    NOT CONVERGED \\(",
                               fit$status_reasons[1]), s)))
  rl <- utils::capture.output(quiet(ferx_runlog(fit)))
  expect_true(any(startsWith(rl, "  Status:      NOT CONVERGED (evaluation_only")))
})

test_that("ferx_covariance() gives a fit the verdict its own covariance step would", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  plain <- quiet(ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                          covariance = FALSE, settings = list(maxiter = 30L)))
  expect_identical(plain$covariance_status, "not_requested")
  cv <- quiet(ferx_covariance(plain))
  ref <- warfarin_fit_cov()
  expect_identical(cv$covariance_status, ref$covariance_status)
  expect_identical(cv$status, ref$status)
  expect_identical(cv$status_reasons, ref$status_reasons)
})

test_that("the STATUS line prints core's label and reasons as given", {
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1), converged = TRUE)
  fit$status <- "converged_with_issues"
  fit$status_label <- "CONVERGED WITH ISSUES"
  fit$status_reasons <- c("condition_number", "covariance_failed")
  expect_match(status_line(fit),
               "STATUS: CONVERGED WITH ISSUES   (condition_number, covariance_failed)",
               fixed = TRUE)
})

# R's own critical condition-number flag (#737 review, finding 1): ferx-core
# computes the number but raises no warning for it, so its `status()` cannot see
# the critical row R adds. The verdict counts it.
cn_case <- function(status, reasons, cn) {
  list(status = status, status_label = toupper(gsub("_", " ", status)),
       status_reasons = reasons, condition_number = cn)
}

test_that("R's condition-number flag joins the verdict in core's order", {
  f <- ferx:::.ferx_status_with_r_flags
  r <- f(cn_case("converged", character(0), 1397.7))
  expect_identical(r$status, "converged_with_issues")
  expect_identical(r$status_label, "CONVERGED WITH ISSUES")
  expect_identical(r$status_reasons, "condition_number")
  # After evaluation_only / budget_spent, before covariance_failed / collapsed.
  r <- f(cn_case("not_converged",
                 c("budget_spent", "covariance_failed", "collapsed_at_guard: ETA_KA"), 2e4))
  expect_identical(r$status, "not_converged")
  expect_identical(r$status_reasons, c("budget_spent", "condition_number",
                                       "covariance_failed", "collapsed_at_guard: ETA_KA"))
  # At or below 1000, non-finite, missing, or already listed by core: unchanged.
  for (case in list(cn_case("converged", character(0), 1000),
                    cn_case("converged", character(0), NaN),
                    cn_case("converged", character(0), NULL),
                    cn_case("converged_with_issues", "condition_number", 5e3))) {
    expect_identical(f(case), case)
  }
})

test_that("a fit with a condition number above 1000 is CONVERGED WITH ISSUES", {
  skip_on_cran()
  ex <- ferx_example("bioavailability")
  fit <- quiet(ferx_fit(ex$model, ex$data, verbose = FALSE, settings = list(maxiter = 50L)))
  skip_if_not(isTRUE(fit$condition_number > 1000), "condition number not above 1000 here")
  w <- ferx_get_warnings(fit, as_df = TRUE)
  expect_true("condition_number" %in% w$category[w$severity == "critical"])
  expect_true("condition_number" %in% fit$status_reasons)
  expect_false(identical(fit$status, "converged"))
  expect_match(status_line(fit), "condition_number", fixed = TRUE)
  # The same verdict after a save / load round trip.
  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$status, fit$status)
  expect_identical(back$status_reasons, fit$status_reasons)
})
