# fit$held_at_bound (ferx-core #2797): every coordinate the covariance step
# reports no standard error for. A real guard collapse is hard to provoke on a
# bundled example, so these tests put one in the fit's bytes with
# fitrx_edit_fit() and let the engine's loader and verdict read it back.

quiet <- function(x) suppressWarnings(suppressMessages(x))

# `fit` with ETA `eta` of warfarin_fit_cov() listed as collapsed at the guard
# (`reason`, `kind` as the wire spells them), held by the step when `held`.
with_collapse <- function(fit, eta, held, reason = "collapsed_at_guard",
                          kind = "Omega") {
  i <- match(eta, fit$eta_names)
  packed <- match(eta, rownames(fit$cov_matrix))
  stopifnot(!is.na(i), !is.na(packed))
  fitrx_edit_fit(fit, function(w) {
    w$held_at_bound <- list(list(
      name = eta, kind = kind, index = i - 1L, packed_index = packed - 1L,
      reason = reason
    ))
    if (held) {
      mask <- as.logical(unlist(w$cov_fixed))
      if (length(mask) == 0L) mask <- rep(FALSE, nrow(fit$cov_matrix))
      mask[packed] <- TRUE
      w$cov_fixed <- as.list(mask)
    }
    w
  })
}

se_cell <- function(lines, label) {
  row <- grep(paste0("^  ", label, " "), lines, value = TRUE)
  expect_length(row, 1L)
  sub(".*SE = ([^ ]+).*", "\\1", row)
}

test_that("a real guard collapse prints N/A and is held quietly (ferx-core #2797)", {
  skip_on_cran()
  # ferx-core's tests/guard_hold.rs "stationary" case: ETA_KA sits on the floor
  # of its log packing at the start, `maxiter = 0` keeps it there, and the
  # covariance step holds it.
  dir <- test_path("fixtures", "guard_hold")
  fit <- expect_no_warning(ferx_fit(file.path(dir, "guard_hold.ferx"),
                                    file.path(dir, "guard_hold.csv"),
                                    verbose = FALSE))
  expect_identical(
    fit$held_at_bound,
    data.frame(parameter = "ETA_KA", kind = "omega", index = 3L,
               reason = "collapsed_at_guard", stringsAsFactors = FALSE)
  )
  # Held by the step: a FIX-like row, so the engine's 0 stays in se_omega,
  # cov_fixed marks it, and cor_matrix does not warn about a zero variance
  # (expect_no_warning above).
  expect_identical(unname(fit$se_omega[3]), 0)
  expect_true(fit$cov_fixed[["ETA_KA"]])
  expect_false(any(startsWith(fit$status_reasons, "collapsed_at_guard")))
  # No surface reports that 0 as a standard error.
  out <- utils::capture.output(print(fit))
  expect_identical(se_cell(out, "ETA_KA"), "N/A")
  expect_true(is.na(fit$estimates$se[fit$estimates$param == "ETA_KA"]))
  rl <- utils::capture.output(ferx_runlog(fit))
  row <- strsplit(trimws(grep("^  ETA_KA ", rl, value = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(tail(row, 2L), c("N/A", "N/A"))
})

test_that("a fit that holds nothing has a zero-row held_at_bound", {
  skip_on_cran()
  fit <- warfarin_fit_cov()
  h <- fit$held_at_bound
  expect_s3_class(h, "data.frame")
  expect_identical(names(h), c("parameter", "kind", "index", "reason"))
  expect_identical(nrow(h), 0L)
})

test_that("a collapse the step held prints N/A and leaves the status clean", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  fit <- with_collapse(base, eta, held = TRUE)
  expect_identical(
    fit$held_at_bound,
    data.frame(parameter = eta, kind = "omega", index = 2L,
               reason = "collapsed_at_guard", stringsAsFactors = FALSE)
  )
  # Held by the step: core lists no reason for it.
  expect_identical(fit$status, base$status)
  expect_false(any(startsWith(fit$status_reasons, "collapsed_at_guard")))
  # The step's mask reaches cov_fixed, so cor_matrix reads the row as held.
  expect_true(fit$cov_fixed[[eta]])
  # Every surface shows no SE for it; the other ETA keeps its own.
  est <- fit$estimates
  expect_true(all(is.na(est[est$param == eta, c("se", "rse_pct", "lower_95", "upper_95")])))
  expect_true(is.finite(est$se[est$param == base$eta_names[1]]))
  out <- utils::capture.output(quiet(print(fit)))
  expect_identical(se_cell(out, eta), "N/A")
  expect_false(identical(se_cell(out, base$eta_names[1]), "N/A"))
  rl <- utils::capture.output(quiet(ferx_runlog(fit)))
  row <- strsplit(trimws(grep(paste0("^  ", eta, " "), rl, value = TRUE)), "[[:space:]]+")[[1]]
  expect_identical(tail(row, 2L), c("N/A", "N/A"))
  # The fit's se_omega keeps the engine's number, as for a FIX parameter.
  expect_identical(fit$se_omega, base$se_omega)
})

test_that("a collapse the step did not hold is a status reason", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  fit <- with_collapse(base, eta, held = FALSE)
  expect_true(fit$converged)
  expect_identical(fit$status, "converged_with_issues")
  expect_identical(fit$status_label, "CONVERGED WITH ISSUES")
  expect_identical(tail(fit$status_reasons, 1L), paste0("collapsed_at_guard: ", eta))
  line <- grep("STATUS:", utils::capture.output(quiet(print(fit))), value = TRUE)
  expect_match(line, paste0("CONVERGED WITH ISSUES   (", paste(fit$status_reasons, collapse = ", "), ")"),
               fixed = TRUE)
  # Not held, still no SE.
  expect_true(is.na(fit$estimates$se[fit$estimates$param == eta]))
})

test_that("held_at_bound survives ferx_save_fit() / ferx_load_fit()", {
  skip_on_cran()
  fit <- with_collapse(warfarin_fit_cov(), warfarin_fit_cov()$eta_names[2], held = FALSE)
  path <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, path))
  back <- quiet(ferx_load_fit(path))
  expect_identical(back$held_at_bound, fit$held_at_bound)
  expect_identical(back$status_reasons, fit$status_reasons)
  expect_identical(back$estimates, fit$estimates)
})

test_that("a held reason or family a later ferx added still loads (ferx-core #2817)", {
  skip_on_cran()
  base <- warfarin_fit_cov()
  eta <- base$eta_names[2]
  later <- with_collapse(base, eta, held = FALSE, reason = "some_future_reason")
  expect_identical(later$held_at_bound$reason, "unknown")
  # Still held: no SE, and no guard-collapse reason it is not.
  expect_true(is.na(later$estimates$se[later$estimates$param == eta]))
  expect_false(any(startsWith(later$status_reasons, "collapsed_at_guard")))
  # A family this build does not know is dropped, not a load failure.
  alien <- with_collapse(base, eta, held = FALSE, kind = "Hyperprior")
  expect_identical(nrow(alien$held_at_bound), 0L)
})
