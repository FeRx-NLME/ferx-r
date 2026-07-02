
# ---- header from test-model-format-structural.R ----
# Unit tests for .ferx_format_structural(): builds the one-line structural
# summary from a model_structure list. Pure string formatting, no fit needed.



test_that(".ferx_format_structural combines type and theta names", {
  ms <- list(model_type = "1-cpt oral", theta_names = c("CL", "V", "KA"))
  expect_identical(ferx:::.ferx_format_structural(ms), "1-cpt oral  (CL, V, KA)")
})
test_that(".ferx_format_structural uses type alone when no thetas", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = "ODE", theta_names = character(0))),
    "ODE"
  )
})
test_that(".ferx_format_structural uses theta names alone when no type", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = NULL, theta_names = c("CL", "V"))),
    "CL, V"
  )
})
test_that(".ferx_format_structural falls back to 'unknown' when nothing is known", {
  expect_identical(
    ferx:::.ferx_format_structural(list(model_type = NULL, theta_names = character(0))),
    "unknown"
  )
})

# ---- header from test-fit-options-conflicts.R ----
# Tests for the model-file vs. call-time [fit_options] conflict detection
# (#31 / PR #66). Covers the parser ferx:::.ferx_parse_model_fit_options() (now
# tested alongside its caller in test-ferx_fit.R) and the warning helper
# ferx:::.ferx_warn_fit_option_conflicts() in isolation, plus a small end-to-end
# check via ferx_fit() against the bundled warfarin example.





# ---------------------------------------------------------------------------
# .ferx_warn_fit_option_conflicts
# ---------------------------------------------------------------------------

empty_settings <- function() list(keys = character(0), values = character(0))










# ---------------------------------------------------------------------------
# End-to-end: ferx_fit() integration with match.call()-based detection
# ---------------------------------------------------------------------------

# Helper: copy warfarin example to a temp dir and rewrite [fit_options] so the
# fit runs in just a few iterations. The model file's `method = foce` and
# `covariance = true` lines are preserved so we can test conflict detection
# against them without paying for a full estimation.
make_fast_warfarin <- function() {
  ex <- ferx_example("warfarin")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  model_path <- file.path(tmpdir, "warfarin.ferx")
  data_path  <- file.path(tmpdir, "warfarin.csv")
  file.copy(ex$data, data_path)
  src <- readLines(ex$model)
  src <- sub("^(\\s*maxiter\\s*=\\s*)\\d+", "\\15", src)
  writeLines(src, model_path)
  list(model = model_path, data = data_path, dir = tmpdir)
}






# Helper: like make_fast_warfarin(), but force `covariance = false` in the model
# file so we can test that the R-side `covariance` default no longer overrides it.
make_fast_warfarin_no_cov <- function() {
  fx <- make_fast_warfarin()
  src <- readLines(fx$model)
  src <- src[!grepl("^\\s*covariance\\s*=", src)]
  src <- sub("(\\[fit_options\\])", "\\1\n  covariance = false", src)
  writeLines(src, fx$model)
  fx
}








# ---- header from test-omega-se.R ----
# Unit tests for internal omega-SE lookup and IMPMAP trace reconstruction.
# Local aliases avoid the ::: operator (undesirable_operator_linter).
omega_se_at <- getFromNamespace(".omega_se_at", "ferx")
reconstruct_impmap_trace <- getFromNamespace(".reconstruct_impmap_trace", "ferx")










test_that(".omega_se_at returns NA when se_omega is NULL", {
  expect_identical(omega_se_at(NULL, 2L, 1L, 1L), NA_real_)
})
test_that(".omega_se_at reads diagonal-only se_omega by position", {
  se <- c(0.11, 0.22)  # length n_eta -> diagonal-only
  expect_equal(omega_se_at(se, 2L, 1L, 1L), 0.11)
  expect_equal(omega_se_at(se, 2L, 2L, 2L), 0.22)
})
test_that(".omega_se_at returns NA for off-diagonal of a diagonal-only se_omega", {
  se <- c(0.11, 0.22)
  expect_identical(omega_se_at(se, 2L, 2L, 1L), NA_real_)
})
test_that(".omega_se_at returns NA for diagonal index beyond se_omega length", {
  # n_eta = 2 but only one SE present -> (2,2) is out of range.
  expect_identical(omega_se_at(c(0.11), 2L, 2L, 2L), NA_real_)
})
test_that(".omega_se_at reads a 2x2 full lower-triangle (column-major)", {
  # Layout for n_eta=2: [(1,1), (2,1), (2,2)]
  se <- c(0.10, 0.20, 0.30)
  expect_equal(omega_se_at(se, 2L, 1L, 1L), 0.10)  # (1,1)
  expect_equal(omega_se_at(se, 2L, 2L, 1L), 0.20)  # (2,1) off-diagonal
  expect_equal(omega_se_at(se, 2L, 2L, 2L), 0.30)  # (2,2)
})
test_that(".omega_se_at is symmetric in (i, j)", {
  se <- c(0.10, 0.20, 0.30)
  expect_equal(omega_se_at(se, 2L, 1L, 2L), omega_se_at(se, 2L, 2L, 1L))
})
test_that(".omega_se_at reads a 3x3 full lower-triangle (column-major)", {
  # Layout for n_eta=3: [(1,1),(2,1),(3,1),(2,2),(3,2),(3,3)]
  se <- c(1, 2, 3, 4, 5, 6)
  expect_equal(omega_se_at(se, 3L, 1L, 1L), 1)
  expect_equal(omega_se_at(se, 3L, 3L, 1L), 3)  # (3,1)
  expect_equal(omega_se_at(se, 3L, 2L, 2L), 4)  # (2,2)
  expect_equal(omega_se_at(se, 3L, 3L, 2L), 5)  # (3,2)
  expect_equal(omega_se_at(se, 3L, 3L, 3L), 6)  # (3,3)
})

test_that("fit carries a structured warnings data frame", {
  fit <- warfarin_fit_cov()
  expect_true(is.data.frame(fit$warnings_structured))
  expect_named(fit$warnings_structured, c("severity", "category", "message", "source_method"))
  # Severity values are always from the fixed vocabulary
  if (nrow(fit$warnings_structured) > 0L) {
    expect_true(all(fit$warnings_structured$severity %in%
      c("critical", "warning", "info")))
  }
})
test_that("structured warnings include R-side condition number flag", {
  fake <- structure(
    list(
      model_name = "m",
      condition_number = 5000,
      eta_normality = NULL
    ),
    class = "ferx_fit"
  )
  raw <- list(
    warnings_severity = character(0),
    warnings_category = character(0),
    warnings_message  = character(0)
  )
  df <- ferx:::.ferx_assemble_structured_warnings(raw, fake)
  expect_true(any(df$category == "condition_number"))
  expect_equal(df$severity[df$category == "condition_number"], "critical")
})
test_that("multiple non-normal ETAs fold into one structured warning (#163)", {
  eta_norm <- data.frame(
    eta   = c("ETA_CL", "ETA_V", "ETA_KA"),
    W     = c(0.90, 0.80, 0.99),
    p_val = c(0.0001, 0.0000, 0.6000),
    flag  = c("[!]", "[!]", ""),
    stringsAsFactors = FALSE
  )
  fake <- structure(
    list(model_name = "m", condition_number = NULL, eta_normality = eta_norm),
    class = "ferx_fit"
  )
  raw <- list(
    warnings_severity = character(0),
    warnings_category = character(0),
    warnings_message  = character(0)
  )
  df <- ferx:::.ferx_assemble_structured_warnings(raw, fake)
  norm_rows <- df[df$category == "eta_normality", , drop = FALSE]
  # Exactly one eta_normality warning, regardless of how many ETAs were flagged
  expect_equal(nrow(norm_rows), 1L)
  # The single message lists every flagged ETA but not the normal one
  expect_match(norm_rows$message, "ETA_CL")
  expect_match(norm_rows$message, "ETA_V")
  expect_false(grepl("ETA_KA", norm_rows$message))
  # Plural noun for >1 flagged ETA
  expect_match(norm_rows$message, "2 ETAs:")
})
test_that(".ferx_eta_normality_warning uses singular noun for 1 ETA", {
  one <- data.frame(
    eta = "ETA_CL", W = 0.9, p_val = 0.0001, flag = "[!]",
    stringsAsFactors = FALSE
  )
  msg <- ferx:::.ferx_eta_normality_warning(one)
  expect_match(msg, "1 ETA:")
  expect_false(grepl("ETA(s)", msg, fixed = TRUE))
})
test_that(".ferx_eta_normality_warning returns NULL when nothing is flagged", {
  expect_null(ferx:::.ferx_eta_normality_warning(NULL))
  all_ok <- data.frame(
    eta = c("ETA_CL", "ETA_V"), W = c(0.99, 0.98),
    p_val = c(0.7, 0.6), flag = c("", ""), stringsAsFactors = FALSE
  )
  expect_null(ferx:::.ferx_eta_normality_warning(all_ok))
})

# ---- header from test-diagnostics-helpers.R ----
# Unit tests for the pure-R diagnostic/formatting helpers in diagnostics.R:
# the warning-guidance lookup and the cli styling shims. No model fit required.

# ---------------------------------------------------------------------------
# .ferx_warning_guidance
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# .ferx_use_cli / .ferx_style
# ---------------------------------------------------------------------------




test_that(".ferx_use_cli follows cli's colour detection", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 1)
  expect_false(ferx:::.ferx_use_cli())
  withr::local_options(cli.num_colors = 256)
  expect_true(ferx:::.ferx_use_cli())
})
test_that(".ferx_style returns text unchanged when cli is off", {
  for (s in c("bold", "green", "red", "yellow", "dim", "unknownstyle")) {
    expect_identical(ferx:::.ferx_style("hi", s, use_cli = FALSE), "hi")
  }
})
test_that(".ferx_style applies ANSI for each known style when cli is on", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 256)
  for (s in c("bold", "green", "red", "yellow", "dim")) {
    styled <- ferx:::.ferx_style("hi", s, use_cli = TRUE)
    expect_true(cli::ansi_has_any(styled), info = s)   # colour actually applied
    expect_identical(cli::ansi_strip(styled), "hi", info = s) # text intact
  }
  # Unknown style falls through to the raw text even with cli on.
  expect_identical(ferx:::.ferx_style("hi", "nope", use_cli = TRUE), "hi")
})

# ---- header from test-transform-output.R ----
# Tests for parameter transform display: CI helpers, ferx_estimates(), print.ferx_fit
# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# .ferx_inv_logit helper ------------------------------------------------


# ferx_estimates() — identity transform ---------------------------------


# ferx_estimates() — log transform --------------------------------------


# ferx_estimates() — logit transform ------------------------------------


# ferx_estimates() — sigma types ----------------------------------------


# ferx_estimates() — omega has variance transform -----------------------


# ferx_estimates() — bare eta/sigma names when declared ------------------



# ferx_estimates() — no SE → NA columns ---------------------------------


# print.ferx_fit — THETA section transform tags -------------------------



# print.ferx_fit — OMEGA section ETA type labels ------------------------




# print.ferx_fit — SIGMA section type labels ----------------------------



# backward compatibility — absent transform fields ----------------------



test_that(".ferx_inv_logit gives expected values", {
  expect_equal(ferx:::.ferx_inv_logit(0),    0.5,     tolerance = 1e-10)
  expect_equal(ferx:::.ferx_inv_logit(Inf),  1,       tolerance = 1e-10)
  expect_equal(ferx:::.ferx_inv_logit(-Inf), 0,       tolerance = 1e-10)
  # logit(0.7) ≈ 0.8473; inv_logit(0.8473) ≈ 0.7
  expect_equal(ferx:::.ferx_inv_logit(log(0.7 / 0.3)), 0.7, tolerance = 1e-6)
})

# ---- header from test-runlog.R ----
# test-runlog.R
# Tests for ferx_runlog() -- shape, content, and graceful degradation.
# Uses warfarin_fit() (no covariance) and warfarin_fit_cov() (with covariance)
# from helper-warfarin-fit.R, which cap iterations at 30 for CI speed.

# -- Tier 1: basic contract ---------------------------------------------------




# -- Tier 2: header fields ----------------------------------------------------



# -- Tier 2: data summary -----------------------------------------------------







# -- Tier 2: parameter table --------------------------------------------------







# -- Tier 2: covariance step --------------------------------------------------




# -- Tier 2: estimation settings ----------------------------------------------






# -- Tier 2: diagnostics ------------------------------------------------------






# -- Tier 2: gradient ---------------------------------------------------------







# -- Tier 2: runtime ----------------------------------------------------------


# -- Tier 3: graceful degradation ---------------------------------------------



# -- Tier 3: NA-guard regression tests ----------------------------------------
# extendr maps Option<String>::None to NA_character_, not NULL.
# Verify that ferx_runlog() never emits NA tokens for these optional fields,
# even when they arrive as NA and another field (covariate_names) triggers
# has_settings = TRUE.





# -- Tier 2: iteration history -----------------------------------------------





# Local alias avoids the ::: operator which lintr flags.
.iter_table <- getFromNamespace(".runlog_iter_table", "ferx")













test_that(".runlog_iter_table: FOCE/FOCEI header has GRAD_NORM and STEP_NORM", {
  path <- write_fake_trace(n = 5, method = "focei")
  on.exit(unlink(path))
  tr  <- read.csv(path, stringsAsFactors = FALSE)
  tbl <- .iter_table(tr, truncate = FALSE)
  expect_match(tbl[1], "GRAD_NORM", fixed = TRUE)
  expect_match(tbl[1], "STEP_NORM", fixed = TRUE)
  expect_no_match(tbl[1], "LM_LAMBDA",  fixed = TRUE)
  expect_no_match(tbl[1], "MH_ACCEPT",  fixed = TRUE)
  expect_match(tbl[1], "dOFV",      fixed = TRUE)
})
test_that(".runlog_iter_table: GN header has LM_LAMBDA and ACC, not GRAD_NORM", {
  path <- write_fake_trace(n = 5, method = "gn")
  on.exit(unlink(path))
  tr <- read.csv(path, stringsAsFactors = FALSE)
  tr$lm_lambda     <- c(1, 2, 4, 8, 16)
  tr$ofv_delta     <- c(NA, -5, -3, -1, -0.5)
  tr$step_accepted <- c(NA_integer_, 1L, 1L, 1L, 0L)
  tr$grad_norm     <- NA_real_
  tr$step_norm     <- NA_real_
  tbl <- .iter_table(tr, truncate = FALSE)
  expect_match(tbl[1], "LM_LAMBDA", fixed = TRUE)
  expect_match(tbl[1], "ACC",       fixed = TRUE)
  expect_no_match(tbl[1], "GRAD_NORM", fixed = TRUE)
  # YES/NO in data rows
  rows <- tbl[-(1:2)]
  expect_true(any(grepl("YES", rows, fixed = TRUE)))
  expect_true(any(grepl("NO",  rows, fixed = TRUE)))
})
test_that(".runlog_iter_table: SAEM header uses COND_NLL and dCOND_NLL", {
  path <- write_fake_trace(n = 6, method = "saem")
  on.exit(unlink(path))
  tr <- read.csv(path, stringsAsFactors = FALSE)
  tr$phase          <- c(rep("explore", 3), rep("accumulate", 3))
  tr$gamma          <- seq(1.0, by = -0.1, length.out = 6)
  tr$mh_accept_rate <- rep(0.35, 6)
  tr$cond_nll       <- tr$ofv
  tbl <- .iter_table(tr, truncate = FALSE)
  expect_match(tbl[1], "COND_NLL",  fixed = TRUE)
  expect_match(tbl[1], "dCOND_NLL", fixed = TRUE)
  expect_match(tbl[1], "MH_ACCEPT", fixed = TRUE)
  expect_no_match(tbl[1], "dOFV",   fixed = TRUE)
})
test_that(".runlog_iter_table truncates long traces correctly", {
  path <- write_fake_trace(n = 40)
  on.exit(unlink(path))
  tr  <- read.csv(path, stringsAsFactors = FALSE)
  tbl <- .iter_table(tr, truncate = TRUE,
                                   trunc_total = 30L, trunc_head = 10L,
                                   trunc_tail  = 10L)
  expect_true(any(grepl("rows not shown", tbl, fixed = TRUE)))
  # head rows (iter 1-10) and tail rows (iter 31-40) present
  expect_true(any(grepl("^  +1 ", tbl)))
  expect_true(any(grepl("^  +40 ", tbl)))
})
test_that(".runlog_iter_table shows all rows when truncate = FALSE", {
  path <- write_fake_trace(n = 40)
  on.exit(unlink(path))
  tr  <- read.csv(path, stringsAsFactors = FALSE)
  tbl <- .iter_table(tr, truncate = FALSE)
  expect_false(any(grepl("rows not shown", tbl, fixed = TRUE)))
  # 40 data rows + 1 header + 1 bar = 42 elements
  expect_equal(length(tbl), 42L)
})
test_that(".runlog_iter_table output is pure ASCII", {
  path <- write_fake_trace(n = 5)
  on.exit(unlink(path))
  tr  <- read.csv(path, stringsAsFactors = FALSE)
  tbl <- .iter_table(tr, truncate = FALSE)
  out <- paste(tbl, collapse = "\n")
  expect_false(nchar(out, type = "bytes") != nchar(out, type = "chars"))
})
