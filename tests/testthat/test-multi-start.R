# fit$starts: what each start of a multi-start fit did (ferx-core #2867). Every
# number is the engine's `FitResult::starts`; these tests check that the glue
# carries it into R unchanged, against the same rows in the bundle's fit.json
# (core's own serialisation of the struct), and that it survives a save/load.

warfarin_multistart_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- suppressWarnings(ferx_fit(
        ex$model, ex$data, method = "focei", verbose = FALSE,
        covariance = FALSE,
        settings = list(maxiter = 30L, n_starts = 3L, multi_start_seed = 7L)
      ))
    }
    fit
  }
})

starts_columns <- c("start", "ofv", "penalized", "objective", "converged",
                    "status", "status_reasons", "n_iterations",
                    "wall_time_secs", "theta", "selected", "error",
                    "error_code")

test_that("a single-start fit has zero start rows and no multi-start line", {
  fit <- warfarin_fit()
  expect_s3_class(fit$starts, "data.frame")
  expect_identical(names(fit$starts), starts_columns)
  expect_identical(nrow(fit$starts), 0L)
  expect_false(any(grepl("Multi-start:", capture.output(print(fit)))))
})

test_that("fit$starts has one row per start, equal to the engine's records", {
  fit <- warfarin_multistart_fit()
  s <- fit$starts
  expect_identical(names(s), starts_columns)
  expect_identical(nrow(s), 3L)
  expect_identical(s$start, 0:2)
  expect_identical(sum(s$selected), 1L)

  # The engine's own rows, as it wrote them into the bundle.
  wire <- fitrx_fit_json(fit$bundle)$starts
  expect_length(wire, 3L)
  num <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
  chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)
  for (k in seq_along(wire)) {
    w <- wire[[k]]
    expect_identical(s$start[k], as.integer(w$start))
    expect_equal(s$ofv[k], num(w$ofv), tolerance = 1e-12)
    expect_identical(s$penalized[k], isTRUE(w$penalized))
    expect_equal(s$objective[k], num(w$objective))
    expect_identical(s$converged[k], isTRUE(w$converged))
    expect_identical(s$status[k], chr(w$status))
    expect_identical(s$status_reasons[[k]],
                     as.character(unlist(w$status_reasons)))
    expect_identical(s$n_iterations[k], as.integer(w$n_iterations))
    expect_equal(s$wall_time_secs[k], as.numeric(w$wall_time_secs))
    expect_equal(unname(s$theta[[k]]), as.numeric(unlist(w$theta)),
                 tolerance = 1e-12)
    expect_identical(s$selected[k], isTRUE(w$selected))
    expect_identical(s$error[k], chr(w$error))
    expect_identical(s$error_code[k], chr(w$error_code))
  }

  # Each start's theta is named like fit$theta; the selected one is the fit's.
  sel <- which(s$selected)
  expect_identical(names(s$theta[[sel]]), names(fit$theta))
  expect_equal(s$theta[[sel]], fit$theta, tolerance = 1e-10)
  expect_equal(s$ofv[sel], fit$ofv, tolerance = 1e-8)
  expect_true(all(s$status %in% c("converged", "converged_with_issues",
                                  "not_converged")))
})

test_that("print() shows the engine's multi-start line", {
  fit <- warfarin_multistart_fit()
  line <- grep("^Multi-start:", fit$warnings, value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "^Multi-start: 3 starts, ")
  out <- capture.output(print(fit))
  expect_true(any(out == paste0(" ", line, "  (fit$starts)")))
})

test_that("fit$starts survives ferx_save_fit() / ferx_load_fit()", {
  fit <- warfarin_multistart_fit()
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path))
  ferx_save_fit(fit, path)
  back <- ferx_load_fit(path)
  expect_identical(back$starts, fit$starts)
})

test_that(".ferx_starts_frame keeps a failed start's empty fields NA", {
  frame <- getFromNamespace(".ferx_starts_frame", "ferx")
  s <- frame(list(
    start = c(0L, 1L), ofv = c(-10.5, NA), penalized = c(FALSE, FALSE),
    objective = c(NA, NA), converged = c(TRUE, FALSE),
    status = c("converged", NA), status_reasons = list(character(), character()),
    n_iterations = c(12L, 0L), wall_time_secs = c(1.5, 0.1),
    theta = list(c(1, 2), numeric()), selected = c(TRUE, FALSE),
    error = c(NA, "boom"), error_code = c(NA, "E_X")
  ), theta_names = c("TVCL", "TVV"))
  expect_identical(names(s$theta[[1]]), c("TVCL", "TVV"))
  expect_identical(s$theta[[2]], numeric())
  expect_true(is.na(s$status[2]))
  expect_identical(s$error_code, c(NA, "E_X"))
  # Zero rows from an empty record list.
  expect_identical(nrow(frame(list())), 0L)
})
