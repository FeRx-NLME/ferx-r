# print.ferx_fit()'s OMEGA rows are ferx-core's own console rows (#543):
# format_omega_rows_from over the R fit's fields, through
# ferx_rust_omega_rows(). The oracle never touches that glue: the fit is saved,
# ferx-core's .fitrx loader rebuilds a FitResult, and the rows come from
# OmegaRowsInput::from_result with the same options
# (ferx_rust_fitrx_omega_rows). Any field the glue drops, misreads or
# reorders shows up as a difference between the two. The KAPPA twin is
# test-kappa-rows.R.

omega_rows <- getFromNamespace(".ferx_omega_rows", "ferx")
kappa_rows <- getFromNamespace(".ferx_kappa_rows", "ferx")
fitrx_omega_rows <- getFromNamespace("ferx_rust_fitrx_omega_rows", "ferx")

# The engine's rows for `fit`, via a saved bundle.
core_omega_rows <- function(fit) {
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  strsplit(fitrx_omega_rows(f), "\n", fixed = TRUE)[[1L]]
}

# The diagonal rows of print(fit)'s OMEGA block: the lines between its rule
# and the next blank line or `--- Correlations ---`.
printed_omega_rows <- function(fit) {
  out <- utils::capture.output(print(fit))
  start <- grep("OMEGA  (between-subject variability)", out, fixed = TRUE)
  expect_length(start, 1L)
  body <- out[-seq_len(start + 1L)]
  body[seq_len(match(TRUE, body == "" | startsWith(body, "  ---"), nomatch = length(body) + 1L) - 1L)]
}

# Both sides, plus print() agreeing with the helper it calls.
expect_core_omega_rows <- function(fit, info) {
  core <- core_omega_rows(fit)
  expect_length(core, nrow(fit$omega))
  expect_identical(omega_rows(fit), core, info = info)
  expect_identical(printed_omega_rows(fit), core, info = info)
  invisible(core)
}

# One fit per model variant per file.
omega_fit <- local({
  cache <- list()
  function(name, example, edit = identity, ...) {
    if (is.null(cache[[name]])) {
      ex <- ferx_example(example)
      txt <- edit(readLines(ex$model))
      mod <- tempfile(fileext = ".ferx")
      writeLines(txt, mod)
      cache[[name]] <<- suppressWarnings(ferx_fit(mod, ex$data, verbose = FALSE, ...))
    }
    cache[[name]]
  }
})

# Replace exactly one line of a model file.
edit_line <- function(old, new) function(txt) {
  stopifnot(sum(txt == old) == 1L)
  txt[txt == old] <- new
  txt
}

# warfarin_logit_f with an additive ETA_V (beside log-normal ETA_CL and
# logit-probability ETA_F), and ETA_CL's initial value written as an SD.
scales_fit <- function() {
  omega_fit("scales", "warfarin_logit_f", function(txt) {
    txt <- edit_line("  V  = TVV  * exp(ETA_V)", "  V  = TVV + ETA_V")(txt)
    edit_line("  omega ETA_CL ~ 0.07                 # BSV on CL (log-normal)",
              "  omega ETA_CL ~ 0.26 (sd)")(txt)
  })
}

lognormal_cv <- function(v) sprintf("CV%% = %.1f", sqrt(exp(v) - 1) * 100)

test_that("diagonal omega (warfarin): print shows ferx-core's rows, exact CV%", {
  skip_on_cran()
  fit <- omega_fit("diag", "warfarin")
  skip_if(!identical(fit$covariance_status, "computed"), "covariance step did not run")
  core <- expect_core_omega_rows(fit, "diag")
  for (i in seq_len(nrow(fit$omega))) {
    expect_match(core[i], sprintf("^  %s +\\[log-normal\\] +=", fit$eta_names[i]))
    expect_match(core[i], lognormal_cv(fit$omega[i, i]), fixed = TRUE)
    expect_match(core[i], sprintf("SE = %.6f$", fit$se_omega[[i]]))
  }
})

test_that("block omega: packed SEs land on the right rows, correlations stay", {
  skip_on_cran()
  fit <- omega_fit("block", "warfarin_block_omega")
  skip_if(is.null(fit$se_omega), "covariance step did not run")
  n <- nrow(fit$omega)
  expect_length(fit$se_omega, n * (n + 1L) / 2L)
  core <- expect_core_omega_rows(fit, "block")
  se_at <- getFromNamespace(".omega_se_at", "ferx")
  for (i in seq_len(n)) {
    expect_match(core[i], sprintf("SE = %.6f$", se_at(fit$se_omega, n, i, i)))
  }
  out <- utils::capture.output(print(fit))
  expect_true(any(startsWith(out, "  ETA_V ~ ETA_CL : cov = ")))
})

test_that("a FIX omega prints `NAME [FIX]` and `SE = ---`, like a FIX kappa", {
  skip_on_cran()
  fit <- omega_fit("fix", "warfarin_block_omega",
                   edit_line("  omega ETA_KA ~ 0.40", "  omega ETA_KA ~ 0.40 FIX"))
  expect_identical(unname(as.logical(fit$omega_fixed)), c(FALSE, FALSE, TRUE))
  core <- expect_core_omega_rows(fit, "fix")
  expect_match(core[3L], "^  ETA_KA \\[FIX\\] +\\[log-normal\\] += 0\\.400000")
  expect_match(core[3L], "SE = ---$")
  expect_false(any(grepl("[FIX]", core[1:2], fixed = TRUE)))
})

test_that("additive and logit-probability ETAs: SD notes with the linked theta", {
  skip_on_cran()
  fit <- scales_fit()
  skip_if(!identical(fit$covariance_status, "computed"), "covariance step did not run")
  expect_identical(unname(fit$eta_param_types), c("log_normal", "additive", "logit_probability"))
  core <- expect_core_omega_rows(fit, "scales")
  sd_v <- sqrt(fit$omega[2L, 2L])
  expect_match(core[2L], "^  ETA_V +\\[additive\\] +=")
  expect_match(core[2L], sprintf("(SD = %.4f, CV%% = %.1f at TVV = %.4f)", sd_v,
                                 sd_v / fit$theta[["TVV"]] * 100, fit$theta[["TVV"]]),
               fixed = TRUE)
  # `[initial specified as SD]` on ETA_CL alone.
  expect_identical(unname(as.logical(fit$omega_init_as_sd)), c(TRUE, FALSE, FALSE))
  expect_true(endsWith(core[1L], "  [initial specified as SD]"), info = core[1L])
  expect_false(any(grepl("initial specified", core[2:3], fixed = TRUE)))
  sd_f <- sqrt(fit$omega[3L, 3L])
  centre <- stats::qlogis(fit$theta[["THETA_F"]])
  expect_match(core[3L], "^  ETA_F +\\[logit\\] +=")
  expect_match(core[3L], sprintf("(SD = %.4f, logit scale, +/-1 SD = [%.3f, %.3f])", sd_f,
                                 stats::plogis(centre - sd_f), stats::plogis(centre + sd_f)),
               fixed = TRUE)
})

test_that("a logit ETA on a logit-scale theta centres its +/-1 SD on the theta", {
  skip_on_cran()
  fit <- omega_fit("logit", "warfarin_logit_f", function(txt) {
    txt <- edit_line("  theta THETA_F(0.80, 0.001, 0.999)   # typical bioavailability, (0,1) scale",
                     "  theta THETA_F(1.4, -5, 5)")(txt)
    edit_line("  F  = inv_logit(logit(THETA_F) + ETA_F)", "  F  = inv_logit(THETA_F + ETA_F)")(txt)
  })
  expect_identical(fit$eta_param_types[[3L]], "logit")
  core <- expect_core_omega_rows(fit, "logit")
  sd_f <- sqrt(fit$omega[3L, 3L])
  centre <- fit$theta[["THETA_F"]]
  expect_match(core[3L], sprintf("+/-1 SD = [%.3f, %.3f])", stats::plogis(centre - sd_f),
                                 stats::plogis(centre + sd_f)), fixed = TRUE)
})

test_that("the theta's transform, not the ETA's label, decides the +/-1 SD centre", {
  # A drifted fit (#371): ETA labelled `logit`, theta a probability. The glue
  # must hand over the transform, or core falls back to the label.
  fit <- make_fake_fit(
    theta = c(THETA_F = 0.7), theta_transforms = c(THETA_F = "logit_probability"),
    omega = matrix(0.10, 1, 1), eta_names = "ETA_F", eta_param_types = "logit",
    eta_linked_theta = "THETA_F"
  )
  sd <- sqrt(0.10)
  expect_match(printed_omega_rows(fit), sprintf("+/-1 SD = [%.3f, %.3f])",
               stats::plogis(stats::qlogis(0.7) - sd), stats::plogis(stats::qlogis(0.7) + sd)),
               fixed = TRUE)
})

test_that("a failed or SIR-fallback covariance step hides every CV%/SD", {
  skip_on_cran()
  fit <- scales_fit()
  for (status in c("failed", "sir_fallback", "computed", "not_requested")) {
    f <- fit
    f$covariance_status <- status
    core <- expect_core_omega_rows(f, status)
    shown <- !status %in% c("failed", "sir_fallback")
    # Both sides of the gate, on every row.
    expect_identical(grepl("(", core, fixed = TRUE), rep(shown, 3L), info = status)
  }
})

test_that("a reloaded .fitrx prints the same ferx-core rows", {
  skip_on_cran()
  for (name in c("fix", "scales")) {
    fit <- switch(name,
      fix = omega_fit("fix", "warfarin_block_omega",
                      edit_line("  omega ETA_KA ~ 0.40", "  omega ETA_KA ~ 0.40 FIX")),
      scales = scales_fit()
    )
    f <- withr::local_tempfile(fileext = ".fitrx")
    suppressWarnings(ferx_save_fit(fit, f))
    loaded <- suppressWarnings(ferx_load_fit(f))
    expect_identical(omega_rows(loaded),
                     strsplit(fitrx_omega_rows(f), "\n", fixed = TRUE)[[1L]], info = name)
    expect_identical(printed_omega_rows(loaded), printed_omega_rows(fit), info = name)
  }
})

# -- OMEGA and KAPPA agree (#543) ----------------------------------------------

test_that("OMEGA and KAPPA rows agree on CV%, FIX and the covariance gate", {
  fit <- make_fake_fit(
    omega = matrix(0.2, 1, 1), eta_names = "ETA_CL", eta_param_types = "log_normal",
    se_omega = 0.05, omega_fixed = FALSE,
    omega_iov = matrix(0.2, 1, 1), kappa_names = "KAPPA_CL", kappa_param_types = "log_normal",
    se_kappa = 0.05, kappa_fixed = FALSE, covariance_status = "computed"
  )
  om <- printed_omega_rows(fit)
  kp <- kappa_rows(fit)
  # The exact CV, sqrt(exp(0.2) - 1) * 100 = 47.1, on both.
  expect_match(om, "(CV% = 47.1)  SE = 0.050000", fixed = TRUE)
  expect_match(kp, "(CV% = 47.1)  SE = 0.050000", fixed = TRUE)
  fit$omega_fixed <- TRUE
  fit$kappa_fixed <- TRUE
  expect_match(printed_omega_rows(fit), "^  ETA_CL \\[FIX\\] .*SE = ---$")
  expect_match(kappa_rows(fit), "^  KAPPA_CL \\[FIX\\] .*SE = ---$")
  fit$covariance_status <- "failed"
  expect_false(grepl("CV%", printed_omega_rows(fit), fixed = TRUE))
  expect_false(grepl("CV%", kappa_rows(fit), fixed = TRUE))
})

test_that("a custom ETA prints no note; unnamed ETAs print OMEGA(i,i)", {
  fit <- make_fake_fit(
    omega = diag(c(0.1, 0.2)), eta_param_types = c("custom", "log_normal"),
    se_omega = c(0.01, 0.02)
  )
  rows <- printed_omega_rows(fit)
  expect_identical(rows, c(
    "  OMEGA(1,1)           [custom]     = 0.100000  SE = 0.010000",
    sprintf("  OMEGA(2,2)           [log-normal] = 0.200000  (%s)  SE = 0.020000",
            lognormal_cv(0.2))
  ))
  # A fit without scales (from before ferx-core #1643): log-normal CV%, no tag.
  fit$eta_param_types <- NULL
  expect_match(printed_omega_rows(fit)[1L], "^  OMEGA\\(1,1\\) {24}= 0\\.100000  \\(CV% = ")
})

test_that("print refuses an unreadable eta scale by name, accepts NA as unknown", {
  base <- make_fake_fit(omega = diag(c(0.1, 0.2)), eta_names = c("ETA_CL", "ETA_V"),
                        eta_param_types = c("additive", "log_normal"))
  f <- base; f$eta_param_types <- c("lognormal", "additive")
  err <- expect_error(print(f), class = "ferx_fit_field_error")
  expect_identical(err$field, "eta_param_types")
  expect_match(conditionMessage(err), "has \"lognormal\" at position 1", fixed = TRUE)
  f <- base; f$eta_param_types <- c(NA, "additive")
  expect_match(printed_omega_rows(f)[1L], "(CV% = ", fixed = TRUE)
})

test_that("the glue refuses what it cannot read rather than guessing", {
  glue <- getFromNamespace("ferx_rust_omega_rows", "ferx")
  call <- function(...) {
    args <- utils::modifyList(list(
      omega = 0.04, n_eta = 1L, eta_names = "E", omega_fixed = NULL, se_omega = NULL,
      eta_param_types = "additive", omega_init_as_sd = TRUE, linked_theta_names = "TV",
      linked_theta_values = 2, linked_theta_transforms = "identity",
      covariance_status = "computed"
    ), list(...))
    do.call(glue, args)
  }
  expect_identical(call(), paste0(
    "  E                    [additive]   = 0.040000  (SD = 0.2000, CV% = 10.0 at TV = 2.0000)",
    "  SE = N/A  [initial specified as SD]\n"
  ))
  expect_error(call(covariance_status = "done"), "covariance_status")
  expect_error(call(eta_param_types = "lognormal"), "random-effect scale")
  expect_error(call(linked_theta_transforms = "exp"), "theta transform")
  expect_error(call(omega = c(0.04, 0)), "1 x 1")
  # No linked theta: the SD alone.
  expect_match(call(linked_theta_names = NA_character_, linked_theta_values = NA_real_),
               "(SD = 0.2000)  SE", fixed = TRUE)
})
