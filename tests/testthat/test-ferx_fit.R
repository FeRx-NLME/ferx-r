
# write_pipe_test_model()/modifying_editor() come from helper-model-pipe.R

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------



test_that("ferx_fit() on ferx_model returns a ferx_fit (not a ferx_model)", {
  ex     <- ferx_example("warfarin")
  result <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_fit(verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
  expect_false(inherits(result, "ferx_model"))
})
test_that("ferx_fit() on ferx_model via named argument works (ferx_fit(model = m))", {
  ex     <- ferx_example("warfarin")
  m      <- ferx_model(ex$model, data = ex$data)
  result <- suppressWarnings(
    ferx_fit(model = m, verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
})
test_that("ferx_fit() on ferx_model produces same theta names as path-style call", {
  ex      <- ferx_example("warfarin")
  by_path <- suppressWarnings(
    ferx_fit(ex$model, data = ex$data, verbose = FALSE, settings = list(maxiter = 30L))
  )
  by_pipe <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_fit(verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_equal(names(by_path$theta), names(by_pipe$theta))
})
test_that("ferx_fit() data argument overrides data stored in ferx_model", {
  ex <- ferx_example("warfarin")
  result <- suppressWarnings(
    ferx_model(ex$model) |>
      ferx_fit(data = ex$data, verbose = FALSE, settings = list(maxiter = 30L))
  )
  expect_s3_class(result, "ferx_fit")
})
test_that("ferx_fit() errors with clear message when ferx_model has no data and none supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_error(
    ferx_model(model = path) |> ferx_fit(verbose = FALSE),
    regexp = "data"
  )
})
test_that("ferx_fit() exposes all estimation arguments in its formals (#52)", {
  # The whole point of dropping the S3 generic: IDE argument completion and
  # inline help introspect formals(ferx_fit), which previously showed only
  # model/data/... and hid every real option from the user.
  fmls <- names(formals(ferx_fit))
  expect_true(all(c(
    "model", "data", "method", "covariance", "verbose", "bloq_method",
    "threads", "mu_referencing", "sir", "gradient", "optimizer_trace",
    "scale_params", "inits_from_nca", "settings", "output", "include_data"
  ) %in% fmls))
  # Convergence-tolerance knobs flow through `settings`, not the signature (#51).
  expect_false("max_unconverged_frac" %in% fmls)
  expect_false("min_obs_for_convergence_check" %in% fmls)
  # No longer an S3 generic: the .default / .ferx_model methods are gone.
  expect_false(exists("ferx_fit.default"))
  expect_false(exists("ferx_fit.ferx_model"))
})
test_that("ferx_fit() rejects unrecognised arguments instead of silently ignoring them", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, methdo = "focei"),
    regexp = "unused argument"
  )
})
test_that("ferx_fit() reports an unrecognised arg without evaluating its value", {
  # The `...` guard must inspect names via ...names()/...length(), never force
  # the promises with list(...): forcing a side-effecting / erroring value would
  # mask the clear "unused argument(s)" message. Regression for the #67 review.
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, bogus = stop("must not be evaluated")),
    regexp = "unused argument\\(s\\) passed to ferx_fit\\(\\): bogus"
  )
})

# ---- header from test-fit.R ----
# Local aliases avoid the ::: operator (undesirable_operator_linter).
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# Return structure — Tier 1
















# These tests require the covariance step to have succeeded. With maxiter = 30L
# the outer optimisation may not converge on all machines, so we skip rather
# than fail — a skip here means "covariance step needs more iterations", not
# a bug. A full-convergence run (no maxiter cap) is tested manually / locally.









# fd_hessian_step argument — Tier 1 (R-side validation, no Rust call)








# ferx_check_init — Tier 1




# ferx_cor_matrix — Tier 1, requires covariance = TRUE





# Gradient correctness — [ENZYME ONLY]
#
# Two blockers before this test can be fully implemented:
#
# 1. Build tier: requires the Enzyme Rust toolchain (a custom Rust fork with
#    automatic differentiation). The Enzyme build takes ~1.5 h and is not
#    available on normal dev machines. The skip_if() guard below makes this
#    test inert on all Tier 1 (FERX_NO_AUTODIFF=1) machines — CI with the
#    Enzyme toolchain is required to exercise it.
#
# 2. Missing API: comparing autodiff vs finite-difference gradients requires
#    a gradient inspection entry point to be exposed from the Rust side. No
#    such API exists yet in ferx-core. Once it lands, replace the inner
#    skip() with a real comparison (e.g. relative error < 1e-4 per element).


# -- Diagnostic fields (Step 3 feature-parity) --------------------------------





# -- IOV: shrinkage_kappa_by_occ (Step 1 feature-parity) ---------------------




# SAEM HMC proposals — [ENZYME ONLY]
#
# HMC proposals in the SAEM E-step are gated on the Enzyme autodiff build
# (`hmc_step` is `#[cfg(feature = "autodiff")]` in ferx-core). On a stable /
# FERX_NO_AUTODIFF=1 build, `n_leapfrog > 0` is silently ignored and the
# sampler uses Metropolis-Hastings, so `saem_n_subjects_hmc` is NA. The
# guard below makes this test inert on Tier 1 machines; CI with the Enzyme
# toolchain exercises the behavioural assertions.




test_that("returns class ferx_fit", {
  expect_s3_class(warfarin_fit(), "ferx_fit")
})
test_that("$theta is a named numeric vector", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$theta))
  expect_false(is.null(names(fit$theta)))
  expect_true(length(names(fit$theta)) > 0L)
})
test_that("$theta names match model theta declarations", {
  fit <- warfarin_fit()
  ex  <- ferx_example("warfarin")
  s   <- ferx_model_inspect(ex$model)
  expect_setequal(names(fit$theta), s$theta_names)
})
test_that("$omega is a square numeric matrix", {
  fit <- warfarin_fit()
  expect_true(is.matrix(fit$omega))
  expect_true(is.numeric(fit$omega))
  expect_equal(nrow(fit$omega), ncol(fit$omega))
})
test_that("$sigma is a finite positive numeric", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$sigma))
  expect_true(is.finite(fit$sigma))
  expect_true(fit$sigma > 0)
})
test_that("$ofv is finite numeric", {
  fit <- warfarin_fit()
  expect_true(is.numeric(fit$ofv))
  expect_true(is.finite(fit$ofv))
})
test_that("$aic > $ofv and $bic > $ofv", {
  fit <- warfarin_fit()
  expect_gt(fit$aic, fit$ofv)
  expect_gt(fit$bic, fit$ofv)
})
test_that("$method equals requested method", {
  fit <- warfarin_fit()
  expect_equal(tolower(fit$method), "focei")
})
test_that("$model_name falls back to basename when engine returns Unnamed", {
  fit <- warfarin_fit()
  expect_equal(fit$model_name, "warfarin")
})
test_that("$sdtab is a data frame", {
  fit <- warfarin_fit()
  expect_s3_class(fit$sdtab, "data.frame")
})
test_that("$sdtab has required columns", {
  fit <- warfarin_fit()
  required <- c("ID", "TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES")
  expect_true(all(required %in% names(fit$sdtab)))
})
test_that("$sdtab always includes TAFD and TAD", {
  fit <- warfarin_fit()
  expect_true("TAFD" %in% names(fit$sdtab),
              label = "TAFD must be present even without [derived] block")
  expect_true("TAD"  %in% names(fit$sdtab),
              label = "TAD must be present even without [derived] block")
  expect_true(all(is.finite(fit$sdtab$TAFD)),
              label = "TAFD must be finite for all rows")
  expect_true(all(is.finite(fit$sdtab$TAD)),
              label = "TAD must be finite for all rows")
})
test_that("$sdtab does not contain ETA columns (ETAs live in ebe_etas)", {
  # Regression guard against ferx-core#188: ETAs were briefly written into
  # sdtab, creating duplicate information with ebe_etas and breaking the
  # per-observation-row contract of sdtab.
  fit <- warfarin_fit()
  eta_pattern <- "^ETA_|^ETA[0-9]"
  eta_cols <- grep(eta_pattern, names(fit$sdtab), value = TRUE)
  expect_length(eta_cols, 0L)
})
test_that("$sdtab has one row per observation", {
  fit <- warfarin_fit()
  ex  <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(fit$sdtab), n_obs)
})
test_that("$se_theta absent when covariance = FALSE", {
  expect_null(warfarin_fit()$se_theta)
})
test_that("$se_theta present and named when covariance = TRUE", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$se_theta), "covariance step did not converge — skipping")
  expect_true(is.numeric(fit$se_theta))
  expect_false(is.null(names(fit$se_theta)))
})
test_that("$se_theta length matches $theta length when covariance = TRUE", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$se_theta), "covariance step did not converge — skipping")
  expect_equal(length(fit$se_theta), length(fit$theta))
})
test_that("errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit("no_such_model.ferx", ex$data, method = "focei"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})
test_that("method = 'imp' passes the ferx_fit validation step", {
  # Tier 1: validate that the method-normalisation block in `ferx_fit()`
  # accepts the new `imp` token (single or chained, plus documented aliases)
  # without erroring. The actual IMP run is exercised in the integration tests
  # once the engine supports it; here we only cover the R-side allowlist + the
  # IMP alias fold. Mirrors the logic in `R/fit.R::ferx_fit()`.
  # Exercises the real `.normalize_method_token` used by `ferx_fit()` so the
  # allowlist and alias folds register coverage instead of a test-local copy.
  normalize_token <- getFromNamespace(".normalize_method_token", "ferx")
  normalize <- function(m) vapply(m, normalize_token, character(1L), USE.NAMES = FALSE)
  expect_equal(normalize("imp"), "imp")
  expect_equal(normalize("IMP"), "imp")
  # The two documented aliases — partial-prefix matching in `match.arg` can't
  # resolve these to `"imp"`, so the wrapper has to fold them explicitly.
  expect_equal(normalize("importance_sampling"), "imp")
  expect_equal(normalize("importance-sampling"), "imp")
  expect_equal(normalize("Importance_Sampling"), "imp")
  expect_equal(normalize(c("focei", "imp")), c("focei", "imp"))
  expect_equal(normalize(c("saem", "imp")), c("saem", "imp"))
  expect_equal(
    normalize(c("focei", "importance_sampling")),
    c("focei", "imp")
  )
  # IMPMAP: exact token, the long alias, and case-insensitivity. Exact-match
  # semantics keep "imp" and "impmap" unambiguous despite the shared prefix.
  expect_equal(normalize("impmap"), "impmap")
  expect_equal(normalize("IMPMAP"), "impmap")
  expect_equal(normalize("importance_sampling_map"), "impmap")
  expect_equal(normalize("importance-sampling-map"), "impmap")
  expect_equal(normalize(c("focei", "impmap")), c("focei", "impmap"))
})
test_that("method = 'bayes' passes the ferx_fit validation step", {
  # Tier 1: R-side allowlist + alias fold for the Bayesian estimator. Mirrors
  # the normalisation block in `R/fit.R::ferx_fit()`. The actual MCMC run is
  # exercised in the integration tests.
  normalize <- function(m) {
    vapply(
      m,
      function(s) {
        normalised <- gsub("[^a-z0-9]", "_", tolower(s))
        if (normalised %in% c("importance_sampling", "importancesampling")) {
          normalised <- "imp"
        }
        if (normalised %in% c("bayesian", "mcmc")) {
          normalised <- "bayes"
        }
        match.arg(
          normalised,
          c("foce", "focei", "saem", "gn", "gn_hybrid", "imp", "bayes")
        )
      },
      character(1L),
      USE.NAMES = FALSE
    )
  }
  expect_equal(normalize("bayes"), "bayes")
  expect_equal(normalize("BAYES"), "bayes")
  expect_equal(normalize("bayesian"), "bayes")
  expect_equal(normalize("mcmc"), "bayes")
})
test_that("ferx_fit rejects 'bayes' in a method chain (standalone only)", {
  # Tier 1: chaining bayes would run the MCMC then discard its posterior (the
  # final stage wins), so the wrapper rejects non-standalone bayes R-side.
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "bayes")),
    regexp = "standalone|only method",
    ignore.case = TRUE
  )
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("bayes", "focei")),
    regexp = "standalone|only method",
    ignore.case = TRUE
  )
})
test_that("method = 'bayes' produces a posterior summary on $bayes", {
  # Tier 2/3: a real (short) MCMC run via the engine. Asserts the posterior
  # structure is surfaced and the means settle near the warfarin estimate.
  fit <- warfarin_bayes_fit()
  expect_false(is.null(fit$bayes))
  b <- fit$bayes
  expect_equal(b$n_chains, 2)
  expect_true(is.finite(b$max_rhat))
  expect_gte(length(b$param_names), 4L)
  expect_true(all(b$q025 <= b$median & b$median <= b$q975))
  # Posterior mean of TVCL near the FOCEI / NONMEM-BAYES value (~0.133).
  i_cl <- match("TVCL", b$param_names)
  expect_false(is.na(i_cl))
  expect_gt(b$mean[i_cl], 0.11)
  expect_lt(b$mean[i_cl], 0.16)
  # Bayes reports credible intervals, not a Hessian covariance.
  expect_identical(fit$covariance_status, "not_requested")
})
test_that("ferx_fit enforces `imp` method-chain rules", {
  # Tier 1: surface the engine's `imp`-placement constraints to the R caller
  # *before* the engine round-trip, so the error message references the R
  # argument and avoids spinning up the backend on inputs that will be
  # rejected anyway. Engine-side guards remain as a safety net (and are
  # covered in the ferx-core integration tests).
  ex <- ferx_example("warfarin")
  # `imp` may appear at most once, regardless of mode.
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "imp", "imp")),
    regexp = "at most once",
    ignore.case = TRUE
  )
  # An evaluation-only `imp` (`imp_eval_only = TRUE`, NONMEM EONLY=1) must be the
  # terminal stage; the default estimating `imp` may sit anywhere, so the same
  # chains are accepted without the flag.
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("imp", "focei"),
      settings = list(imp_eval_only = TRUE)),
    regexp = "final stage|must be the final stage",
    ignore.case = TRUE
  )
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "imp", "focei"),
      settings = list(imp_eval_only = TRUE)),
    regexp = "final stage|must be the final stage",
    ignore.case = TRUE
  )
})
test_that("errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, "no_such_data.csv", method = "focei"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})
test_that("fd_hessian_step = 0 is rejected", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, fd_hessian_step = 0),
    regexp = "positive finite",
    ignore.case = TRUE
  )
})
test_that("fd_hessian_step < 0 is rejected", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, fd_hessian_step = -0.01),
    regexp = "positive finite",
    ignore.case = TRUE
  )
})
test_that("fd_hessian_step = Inf is rejected", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, fd_hessian_step = Inf),
    regexp = "positive finite",
    ignore.case = TRUE
  )
})
test_that("fd_hessian_step as character is rejected", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, fd_hessian_step = "0.05"),
    regexp = "positive finite",
    ignore.case = TRUE
  )
})
test_that("fd_hessian_step and settings duplicate does not cause duplicate-key error", {
  ex <- ferx_example("warfarin")
  expect_no_error(
    ferx_fit(ex$model, ex$data,
             fd_hessian_step = 0.05,
             settings = list(fd_hessian_step = 0.1),
             covariance = FALSE, verbose = FALSE)
  )
})
test_that("fd_hessian_step via settings still works when dedicated arg is absent", {
  ex <- ferx_example("warfarin")
  expect_no_error(
    ferx_fit(ex$model, ex$data,
             settings = list(fd_hessian_step = 0.05),
             covariance = FALSE, verbose = FALSE)
  )
})
test_that("outer_ftol / outer_xtol BOBYQA settings are accepted by ferx_fit (ferx-core #469)", {
  # The two derivative-free outer-optimizer stop tolerances flow through the
  # engine's apply_fit_option passthrough (no dedicated arg). This checks the R
  # wrapper forwards them; the numerical effect (TTE frailty omega^2 onto the
  # NONMEM/nlmixr2 minimum) is validated engine-side in tte_convergence_weibull_mixed.
  # Until the ferx-core pin carries #469 these are unknown keys, so skip gracefully;
  # after the pin bump the fit runs end to end.
  ex <- ferx_example("warfarin")
  ran <- tryCatch(
    {
      ferx_fit(ex$model, ex$data,
               settings = list(outer_ftol = 1e-8, outer_xtol = 1e-4),
               covariance = FALSE, verbose = FALSE)
      TRUE
    },
    error = function(e) {
      if (grepl("unknown fit setting", conditionMessage(e), ignore.case = TRUE)) {
        skip("ferx-core pin predates #469 (outer_ftol/outer_xtol) — bump the engine pin")
      }
      stop(e)
    }
  )
  expect_true(ran)
})
test_that("autodiff OFV gradient matches finite-difference gradient within tolerance", {
  skip_if(
    !isTRUE(ferx_rust_autodiff_enabled()),
    "Enzyme autodiff not available — skipping gradient-correctness test"
  )
  skip("Gradient inspection API not yet exposed from ferx-core — implement once available")
})
test_that("$n_threads_used is a positive integer", {
  fit <- warfarin_fit()
  expect_true(is.integer(fit$n_threads_used) || is.numeric(fit$n_threads_used))
  expect_true(fit$n_threads_used >= 1L)
})
test_that("$nlopt_missing_algorithms is a character vector", {
  fit <- warfarin_fit()
  expect_true(is.character(fit$nlopt_missing_algorithms))
})
test_that("$saem_mu_ref_m_step_evals_saved is NULL or numeric for non-SAEM fit", {
  fit <- warfarin_fit()
  val <- fit$saem_mu_ref_m_step_evals_saved
  expect_true(is.null(val) || (is.numeric(val) && length(val) == 1L))
})
test_that("$covariance_n_evals_estimated is NULL or numeric", {
  fit <- warfarin_fit()
  val <- fit$covariance_n_evals_estimated
  expect_true(is.null(val) || (is.numeric(val) && length(val) == 1L))
})
test_that("shrinkage_kappa_by_occ is a data frame with correct structure for IOV fit", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  df  <- fit$shrinkage_kappa_by_occ
  # Non-NULL: IOV model should always produce per-occasion shrinkage
  expect_false(is.null(df), info = "shrinkage_kappa_by_occ should not be NULL for IOV fit")
  expect_s3_class(df, "data.frame")
  # Must contain an 'occ' column
  expect_true("occ" %in% colnames(df), info = "data frame must have occ column")
  # One kappa column per kappa parameter
  kap_cols <- setdiff(colnames(df), "occ")
  expect_equal(length(kap_cols), length(fit$kappa_names))
  if (length(fit$kappa_names) > 0L) {
    expect_equal(sort(kap_cols), sort(fit$kappa_names))
  }
  # Shrinkage values must be finite (can exceed 1 in underdetermined models)
  for (cn in kap_cols) {
    expect_true(all(is.finite(df[[cn]])),
                info = sprintf("shrinkage values for %s must be finite", cn))
  }
})
test_that("shrinkage_kappa_by_occ is NULL for non-IOV fit", {
  fit <- warfarin_fit()
  expect_null(fit$shrinkage_kappa_by_occ)
})
test_that("print.ferx_fit shows per-occasion shrinkage table for IOV fit", {
  skip_on_cran()
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  skip_if(is.null(fit$shrinkage_kappa_by_occ),
          "IOV fit did not produce shrinkage_kappa_by_occ")
  out <- capture.output(suppressWarnings(print(fit)))
  expect_true(any(grepl("Shrinkage by occasion", out, fixed = TRUE)))
  expect_true(any(grepl("Occasion", out, fixed = TRUE)))
})
test_that("SAEM with n_leapfrog > 0 uses HMC proposals on the AD build", {
  skip_if(
    !isTRUE(ferx_rust_autodiff_enabled()),
    "Enzyme autodiff not available — skipping SAEM HMC behavioural test"
  )

  ex  <- ferx_example("warfarin_saem")
  fit <- ferx_fit(
    ex$model, ex$data,
    method   = "saem",
    verbose  = FALSE,
    settings = list(
      n_leapfrog    = 3L,
      n_exploration = 50L,
      n_convergence = 50L,
      omega_burnin  = 20L,
      seed          = 42L
    )
  )

  expect_s3_class(fit, "ferx_fit")
  # The count must be a non-NA integer, and at least one subject must have
  # used an HMC proposal (warfarin_saem is an analytical model, so no ODE
  # fall-back to MH).
  expect_false(is.na(fit$saem_n_subjects_hmc))
  expect_type(fit$saem_n_subjects_hmc, "integer")
  expect_gt(fit$saem_n_subjects_hmc, 0L)
  expect_lte(fit$saem_n_subjects_hmc, fit$n_subjects)
})
test_that("SAEM with n_leapfrog = 0 uses MH (no HMC subjects)", {
  skip_if(
    !isTRUE(ferx_rust_autodiff_enabled()),
    "Enzyme autodiff not available — skipping SAEM HMC behavioural test"
  )

  ex  <- ferx_example("warfarin_saem")
  fit <- ferx_fit(
    ex$model, ex$data,
    method   = "saem",
    verbose  = FALSE,
    settings = list(
      n_leapfrog    = 0L,
      n_exploration = 50L,
      n_convergence = 50L,
      seed          = 42L
    )
  )

  # n_leapfrog = 0 is pure Metropolis-Hastings: no subject uses HMC.
  expect_true(is.na(fit$saem_n_subjects_hmc) || fit$saem_n_subjects_hmc == 0L)
})

# ---- header from test-fit-options-conflicts.R ----
# Tests for the model-file vs. call-time [fit_options] conflict detection
# (#31 / PR #66). Covers the parser ferx:::.ferx_parse_model_fit_options() and the
# warning helper ferx:::.ferx_warn_fit_option_conflicts() in isolation, plus a small
# end-to-end check via ferx_fit() against the bundled warfarin example.

# ---------------------------------------------------------------------------
# .ferx_parse_model_fit_options
# ---------------------------------------------------------------------------

write_fit_options_model <- function(lines) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)",
               "[fit_options]", lines), path)
  path
}








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








test_that("no warning when model_file_opts is empty", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = setNames(character(0), character(0)),
      dedicated_explicit = list(method = "focei"),
      settings_parts     = empty_settings()
    )
  )
})
test_that("no warning when explicit dedicated arg matches model file", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "foce"),
      dedicated_explicit = list(method = "foce"),
      settings_parts     = empty_settings()
    )
  )
})
test_that("warning when explicit dedicated arg disagrees with model file", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "foce"),
      dedicated_explicit = list(method = "focei"),
      settings_parts     = empty_settings()
    ),
    "method = foce.*overrides it with `focei`"
  )
})
test_that("no warning when dedicated arg is NULL (default accepted, not explicit)", {
  # This is the regression test for the spurious-warning bug: a user who does
  # not pass `method` should not be warned that they "overrode" the model file.
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(method = "focei"),
      dedicated_explicit = list(method = NULL),  # user did not pass `method`
      settings_parts     = empty_settings()
    )
  )
})
test_that("comparison is case-insensitive (TRUE vs true does not warn)", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(covariance = "TRUE"),
      dedicated_explicit = list(covariance = "true"),
      settings_parts     = empty_settings()
    )
  )
})
test_that("settings list key conflict warns", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(maxiter = "300"),
      dedicated_explicit = list(),
      settings_parts     = list(keys = "maxiter", values = "500")
    ),
    "maxiter = 300.*overrides it with `500`"
  )
})
test_that("settings list key matching model file does not warn", {
  expect_silent(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(maxiter = "300"),
      dedicated_explicit = list(),
      settings_parts     = list(keys = "maxiter", values = "300")
    )
  )
})
test_that("aliased model-file keys are matched (bloq → bloq_method)", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(bloq = "drop"),
      dedicated_explicit = list(bloq_method = "m3"),
      settings_parts     = empty_settings()
    ),
    "bloq = drop.*overrides it with `m3`"
  )
})
test_that("aliased model-file keys are matched (gradient_method → gradient)", {
  expect_warning(
    ferx:::.ferx_warn_fit_option_conflicts(
      model_file_opts    = c(gradient_method = "ad"),
      dedicated_explicit = list(gradient = "fd"),
      settings_parts     = empty_settings()
    ),
    "gradient_method = ad.*overrides it with `fd`"
  )
})
test_that("ferx_fit() does not warn when defaults are accepted (regression for #66 review)", {
  # The model file specifies method = foce, covariance = true. Calling
  # ferx_fit() without overriding either should NOT warn — match.call() must
  # distinguish "user accepted the default" from "user explicitly passed the
  # default".
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  # Use expect_no_warning() rather than expect_silent() — ferx_fit() may print
  # informational messages (e.g. "Mu-referencing detected for..."), which are
  # unrelated to the conflict-detection logic under test.
  expect_no_warning(
    fit <- suppressMessages(ferx_fit(
      fx$model, fx$data,
      verbose = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ))
  )
  expect_s3_class(fit, "ferx_fit")
  expect_equal(fit$model_file_settings$method, "foce")
})
test_that("mu-reference detection messages ignore non-mu-reference warnings", {
  warnings <- c(
    "mu-ref: ETA_CL, ETA_V",
    "SAEM: individual parameter(s) not mu-referenced: CL. This can strongly affect convergence.",
    "[SAEM] individual parameter(s) not mu-referenced: V",
    "[FOCEI] mu-ref: ETA_KA"
  )

  expect_equal(
    .ferx_mu_ref_detection_names(warnings),
    c("ETA_CL, ETA_V", "ETA_KA")
  )
})
test_that("ferx_fit() warns when an explicit dedicated arg overrides the model file", {
  fx <- make_fast_warfarin()  # model file: method = foce
  on.exit(unlink(fx$dir, recursive = TRUE))

  expect_warning(
    ferx_fit(
      fx$model, fx$data,
      method = "focei",
      verbose = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ),
    "method = foce.*overrides it with `focei`"
  )
})
test_that("ferx_fit() runs the model file's method when none is passed (#558)", {
  # The model file sets `method = foce`. Before #558 the R-side default
  # method = "focei" silently overrode it; now passing no method keeps FOCE.
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressMessages(ferx_fit(
    fx$model, fx$data,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  ))
  expect_equal(fit$method, "FOCE")
})
test_that("ferx_fit() method argument overrides the model file's method (#558)", {
  fx <- make_fast_warfarin()  # model file: method = foce
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressWarnings(suppressMessages(ferx_fit(
    fx$model, fx$data,
    method = "focei",
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  )))
  expect_equal(fit$method, "FOCEI")
})
test_that("ferx_fit() honours the model file's `covariance = false` when none is passed (#558)", {
  # Before extending #558, the R default `covariance = TRUE` ran the covariance
  # step even though the model file disabled it. Passing no `covariance` must now
  # keep the model file's setting.
  fx <- make_fast_warfarin_no_cov()
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressMessages(ferx_fit(
    fx$model, fx$data,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  ))
  expect_equal(fit$covariance_status, "not_requested")
})
test_that("ferx_fit() covariance argument overrides the model file (#558)", {
  fx <- make_fast_warfarin_no_cov()  # model file: covariance = false
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressWarnings(suppressMessages(ferx_fit(
    fx$model, fx$data,
    covariance = TRUE,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  )))
  expect_true(fit$covariance_status %in% c("computed", "failed", "sir_fallback"))
})
test_that("ferx_fit() does not warn when a flag arg is passed explicitly as NULL (#558)", {
  # Regression: naming covariance / sir / mu_referencing / gradient = NULL marks
  # them "explicit" in match.call(), but NULL means "keep the model file value".
  # The empty sentinel must not leak into the conflict comparison and emit a
  # garbled override warning.
  fx <- make_fast_warfarin_no_cov()  # model file: covariance = false, method = foce
  on.exit(unlink(fx$dir, recursive = TRUE))

  expect_no_warning(
    suppressMessages(ferx_fit(
      fx$model, fx$data,
      covariance     = NULL,
      sir            = NULL,
      mu_referencing = NULL,
      gradient       = NULL,
      verbose        = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ))
  )
})
test_that("ferx_fit() records the requested gradient token on a default fit (#558)", {
  # Regression: `gradient` defaults to NULL; assigning that raw NULL to
  # result$gradient dropped the field entirely (and the RUN INFO line). The
  # field must survive as a token.
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  fit <- suppressWarnings(suppressMessages(ferx_fit(
    fx$model, fx$data,
    verbose = FALSE,
    settings = list(max_unconverged_frac = 1.0)
  )))
  expect_false(is.null(fit$gradient))
  expect_true(is.character(fit$gradient) && length(fit$gradient) == 1L)
  expect_true(fit$gradient %in% c("auto", "ad", "fd"))
})
test_that("ferx_fit() errors when sir = TRUE but the model file disables covariance (#558)", {
  # The model file sets covariance = false. Requesting SIR without an explicit
  # covariance arg must error R-side rather than silently skipping SIR.
  fx <- make_fast_warfarin_no_cov()
  on.exit(unlink(fx$dir, recursive = TRUE))

  expect_error(
    ferx_fit(
      fx$model, fx$data,
      sir = TRUE,
      verbose = FALSE,
      settings = list(max_unconverged_frac = 1.0)
    ),
    "requires covariance"
  )
})
test_that("ferx_fit() errors when sir = TRUE and covariance = FALSE explicitly", {
  fx <- make_fast_warfarin()
  on.exit(unlink(fx$dir, recursive = TRUE))

  expect_error(
    ferx_fit(
      fx$model, fx$data,
      sir = TRUE,
      covariance = FALSE,
      verbose = FALSE
    ),
    "requires covariance"
  )
})

# ---- header from test-fit-eigenvalues.R ----
# Tests for the R-layer handling of eigenvalues and condition_number fields
# that are pre-computed by the Rust backend and exposed via the FFI.
#
# The Rust side (ferx-core) is tested for computation correctness
# (fixed-param exclusion, Inf for non-positive eigenvalue) in api.rs
# tests_cov_diagnostics. These tests cover the R-side sentinel-to-NULL
# conversion and the high-condition-number warning injection, exercised
# directly through the production helper .ferx_apply_cov_sentinels().

make_fit <- function(cov_eigenvalues = numeric(0),
                     cov_condition_number = NaN,
                     warnings = character(0)) {
  list(
    cov_eigenvalues      = cov_eigenvalues,
    cov_condition_number = cov_condition_number,
    warnings             = warnings
  )
}

apply_sentinels <- ferx:::.ferx_apply_cov_sentinels








test_that("empty cov_eigenvalues and NaN cov_condition_number become NULL", {
  r <- apply_sentinels(make_fit())
  expect_null(r$eigenvalues)
  expect_null(r$condition_number)
})
test_that("populated eigenvalues and finite condition_number are forwarded as-is", {
  r <- apply_sentinels(make_fit(
    cov_eigenvalues = c(1.5, 0.8, 0.3),
    cov_condition_number = 5.0
  ))
  expect_equal(r$eigenvalues, c(1.5, 0.8, 0.3))
  expect_equal(r$condition_number, 5.0)
})
test_that("Inf condition_number is forwarded as Inf (not NULL)", {
  r <- apply_sentinels(make_fit(
    cov_eigenvalues = c(2.5, -0.5),
    cov_condition_number = Inf
  ))
  expect_true(is.infinite(r$condition_number))
  expect_false(is.null(r$condition_number))
})
test_that("condition_number > 1000 appends a warning", {
  r <- apply_sentinels(make_fit(
    cov_eigenvalues = c(1500.0, 1.5),
    cov_condition_number = 1000.1
  ))
  expect_true(any(grepl("High condition number", r$warnings)))
})
test_that("condition_number <= 1000 does not append a warning", {
  r <- apply_sentinels(make_fit(
    cov_eigenvalues = c(1.2, 1.0),
    cov_condition_number = 1000.0
  ))
  expect_false(any(grepl("High condition number", r$warnings)))
})
test_that("Inf condition_number does not trigger the high-condition-number warning", {
  r <- apply_sentinels(make_fit(
    cov_eigenvalues = c(2.5, -0.5),
    cov_condition_number = Inf
  ))
  expect_false(any(grepl("High condition number", r$warnings)))
})

# ---- header from test-settings-imp.R ----
to_strings <- getFromNamespace(".ferx_settings_to_strings", "ferx")





test_that("importance-sampling settings use imp_* names", {
  parts <- to_strings(list(
    imp_samples = 2000L,
    imp_proposal_df = 4,
    imp_seed = 1L,
    imp_low_ess_threshold = 0.1
  ))

  expect_identical(parts$keys, c(
    "imp_samples",
    "imp_proposal_df",
    "imp_seed",
    "imp_low_ess_threshold"
  ))
  expect_identical(parts$values, c("2000", "4", "1", "0.10000000000000001"))
})
test_that("short-lived is_* settings are not translated by R", {
  parts <- to_strings(list(is_seed = 7L))

  expect_identical(parts$keys, "is_seed")
  expect_identical(parts$values, "7")
})
test_that("duplicated settings are rejected before serialization", {
  expect_error(
    to_strings(list(imp_samples = 1L, imp_samples = 2L)),
    "uniquely-named"
  )
})

# ---- header from test-omega-se.R ----
# Unit tests for internal omega-SE lookup and IMPMAP trace reconstruction.
# Local aliases avoid the ::: operator (undesirable_operator_linter).
omega_se_at <- getFromNamespace(".omega_se_at", "ferx")
reconstruct_impmap_trace <- getFromNamespace(".reconstruct_impmap_trace", "ferx")










test_that(".reconstruct_impmap_trace rebuilds a data.frame column-major", {
  tr <- list(
    flat      = c(1, 2, 3, 4, 5, 6),  # 3 rows x 2 cols, column-major
    n_rows    = 3L,
    n_cols    = 2L,
    col_names = c("CL", "V")
  )
  df <- reconstruct_impmap_trace(tr)
  expect_s3_class(df, "data.frame")
  expect_identical(dim(df), c(3L, 2L))
  expect_identical(colnames(df), c("CL", "V"))
  expect_equal(df$CL, c(1, 2, 3))
  expect_equal(df$V, c(4, 5, 6))
})

# ---- header from test-summary.R ----
# Helpers ---------------------------------------------------------------

make_fake_fit <- function(...) {
  defaults <- list(
    model_name = "test_model",
    data_name = "test_data",
    converged = TRUE,
    method = "focei",
    method_chain = NULL,
    gradient = "auto",
    ofv = -100.0,
    aic = -90.0,
    bic = -80.0,
    n_subjects = 10L,
    n_obs = 100L,
    n_parameters = 5L,
    n_iterations = 50L,
    theta = c(CL = 1.0, V = 10.0),
    se_theta = NULL,
    omega = NULL,
    se_omega = NULL,
    sigma = NULL,
    se_sigma = NULL,
    shrinkage_eta = NULL,
    shrinkage_eps = NULL,
    covariance_status = "not_requested",
    wall_time_secs = 5.0,
    ferx_version = "0.1.0",
    ebe_convergence_warnings = 0L,
    max_unconverged_subjects = NULL,
    total_ebe_fallbacks = 0L,
    call_settings = list(),
    sir_ess = NULL,
    saem_n_subjects_hmc = NULL,
    warnings = character(0)
  )
  structure(modifyList(defaults, list(...)), class = "ferx_fit")
}

# summary.ferx_fit structure -------------------------------------------







# print.ferx_summary output --------------------------------------------

















# gradient_used display ------------------------------------------------





# .ferx_parse_setting_value --------------------------------------------





# call_settings round-trip type fidelity --------------------------------




test_that("summary.ferx_fit returns a ferx_summary with expected fields", {
  fit <- make_fake_fit()
  s <- summary(fit)

  expect_s3_class(s, "ferx_summary")
  expect_named(s, c(
    "model_name", "data_name", "gradient", "gradient_used",
    "method", "method_chain",
    "converged", "ofv", "aic", "bic", "n_subjects", "n_obs",
    "n_parameters", "n_iterations", "theta", "se_theta", "omega",
    "se_omega", "sigma", "se_sigma", "shrinkage_eta", "shrinkage_eps",
    "covariance_status", "eigenvalues", "condition_number",
    "wall_time_secs", "ferx_version",
    "ebe_convergence_warnings", "max_unconverged_subjects",
    "total_ebe_fallbacks", "model_structure", "call_settings",
    "model_file_settings", "sir_ess", "importance_sampling", "bayes",
    "warnings",
    "uses_sde", "dw_statistic", "iwres_lag1_r", "saem_n_subjects_hmc"
  ), ignore.order = TRUE)
})
test_that("summary.ferx_fit passes through call_settings", {
  fit <- make_fake_fit(call_settings = list(optimizer = "slsqp", max_iter = 200L))
  s <- summary(fit)

  expect_equal(s$call_settings$optimizer, "slsqp")
  expect_equal(s$call_settings$max_iter, 200L)
})
test_that("summary.ferx_fit defaults call_settings to empty list when absent", {
  fit <- make_fake_fit()
  fit$call_settings <- NULL  # simulate pre-PR fit object
  s <- summary(fit)

  expect_equal(s$call_settings, list())
})
test_that("summary.ferx_fit passes through warnings", {
  fit <- make_fake_fit(warnings = c("ETA1 may be non-normal", "ETA2 boundary"))
  s <- summary(fit)

  expect_equal(s$warnings, c("ETA1 may be non-normal", "ETA2 boundary"))
})
test_that("summary.ferx_fit passes through sir_ess", {
  fit <- make_fake_fit(sir_ess = 847.2)
  s <- summary(fit)

  expect_equal(s$sir_ess, 847.2)
})
test_that("summary.ferx_fit passes through method_chain", {
  fit <- make_fake_fit(method = "focei", method_chain = c("saem", "focei"))
  s <- summary(fit)

  expect_equal(s$method_chain, c("saem", "focei"))
})
test_that("print.ferx_summary shows method in uppercase for a single method", {
  s <- summary(make_fake_fit(method = "focei", method_chain = NULL))
  out <- capture.output(print(s))

  expect_true(any(grepl("Method:\\s+FOCEI", out)))
})
test_that("print.ferx_summary shows full chain in uppercase when method_chain > 1", {
  s <- summary(make_fake_fit(method = "focei", method_chain = c("saem", "focei")))
  out <- capture.output(print(s))

  expect_true(any(grepl("Method:\\s+SAEM -> FOCEI", out)))
  expect_false(any(grepl("saem|focei", out)))
})
test_that("print.ferx_summary omits Settings block when both settings sources are empty", {
  s <- summary(make_fake_fit(call_settings = list()))
  out <- capture.output(print(s))

  expect_false(any(grepl("^Settings", out)))
})
test_that("print.ferx_summary shows Settings block when call_settings is non-empty", {
  s <- summary(make_fake_fit(call_settings = list(optimizer = "slsqp", max_iter = 200)))
  out <- capture.output(print(s))

  expect_true(any(grepl("^Settings", out)))
  expect_true(any(grepl("optimizer", out)))
  expect_true(any(grepl("slsqp", out)))
  expect_true(any(grepl("max_iter", out)))
})
test_that("print.ferx_summary omits Warnings block when no warnings", {
  s <- summary(make_fake_fit(
    warnings = character(0), ebe_convergence_warnings = 0L, total_ebe_fallbacks = 0L
  ))
  out <- capture.output(print(s))

  expect_false(any(grepl("^Warnings:", out)))
})
test_that("print.ferx_summary shows model warnings", {
  s <- summary(make_fake_fit(warnings = "ETA1 Shapiro-Wilk p=0.03"))
  out <- capture.output(print(s))

  expect_true(any(grepl("Warnings:", out)))
  expect_true(any(grepl("ETA1 Shapiro-Wilk", out)))
})
test_that("print.ferx_summary aggregates EBE convergence warnings", {
  s <- summary(make_fake_fit(
    ebe_convergence_warnings = 3L,
    max_unconverged_subjects = 2L
  ))
  out <- capture.output(print(s))

  expect_true(any(grepl("Warnings:", out)))
  expect_true(any(grepl("3 EBE convergence", out)))
  expect_true(any(grepl("max unconverged subjects: 2", out)))
})
test_that("print.ferx_summary aggregates EBE fallback warnings", {
  s <- summary(make_fake_fit(total_ebe_fallbacks = 5L))
  out <- capture.output(print(s))

  expect_true(any(grepl("5 EBE fallback", out)))
})
test_that("print.ferx_summary omits SIR ESS line when sir_ess is NULL", {
  s <- summary(make_fake_fit(sir_ess = NULL))
  out <- capture.output(print(s))

  expect_false(any(grepl("SIR ESS", out)))
})
test_that("print.ferx_summary shows SIR ESS when sir_ess is present", {
  s <- summary(make_fake_fit(sir_ess = 512.3))
  out <- capture.output(print(s))

  expect_true(any(grepl("SIR ESS.*512\\.3", out)))
})
test_that("print.ferx_summary shows R package version with core version in parens", {
  s <- summary(make_fake_fit(ferx_version = "0.1.0"))
  out <- capture.output(print(s))
  pkg_v <- as.character(utils::packageVersion("ferx"))
  expect_true(any(grepl(
    paste0("ferx v", pkg_v, " \\(core v0\\.1\\.0\\)"), out, fixed = FALSE
  )))
})
test_that("print.ferx_fit shows R package version with core version in parens", {
  om <- matrix(0.04, 1, 1, dimnames = list("ETA_CL", "ETA_CL"))
  fit <- make_fake_fit(
    ferx_version = "0.1.0",
    omega = om,
    sigma = 0.01
  )
  out <- capture.output(suppressWarnings(print(fit)))
  pkg_v <- as.character(utils::packageVersion("ferx"))
  expect_true(any(grepl(
    paste0("ferx v", pkg_v, " \\(core v0\\.1\\.0\\)"), out, fixed = FALSE
  )))
})
test_that("print.ferx_summary wraps output in dashed borders", {
  s <- summary(make_fake_fit())
  out <- capture.output(print(s))

  bars <- grep("^-{10,}$", out)
  expect_gte(length(bars), 2L)
})
test_that("print.ferx_summary shows Structure line when model_structure is non-NULL", {
  s <- summary(make_fake_fit(model_structure = list(
    theta_names = c("TVCL", "TVV"),
    model_type  = "1-cpt oral",
    iiv         = c("ETA_CL", "ETA_V"),
    iov         = character(0),
    residual    = "proportional"
  )))
  out <- capture.output(print(s))

  expect_true(any(grepl("^Structure:", out)))
  expect_true(any(grepl("1-cpt oral", out)))
  expect_true(any(grepl("ETA_CL", out)))
  expect_true(any(grepl("proportional", out)))
})
test_that("print.ferx_summary omits Structure line when model_structure is NULL", {
  s <- summary(make_fake_fit())
  s$model_structure <- NULL
  out <- capture.output(print(s))

  expect_false(any(grepl("^Structure:", out)))
})
test_that(".ferx_short_gradient_label maps engine labels to short tokens", {
  expect_identical(ferx:::.ferx_short_gradient_label("Enzyme AD"), "ad")
  expect_identical(ferx:::.ferx_short_gradient_label("finite differences"), "fd")
  expect_identical(ferx:::.ferx_short_gradient_label("N/A"), "N/A")
  expect_identical(ferx:::.ferx_short_gradient_label(""), NA_character_)
  expect_identical(ferx:::.ferx_short_gradient_label(NULL), NA_character_)
})
test_that("summary.ferx_fit carries gradient_used through", {
  fit <- make_fake_fit(gradient = "auto", gradient_used = "ad")
  s <- summary(fit)
  expect_identical(s$gradient_used, "ad")
})
test_that("print.ferx_summary shows requested -> used arrow when used is known", {
  s <- summary(make_fake_fit(gradient = "auto", gradient_used = "ad"))
  out <- capture.output(print(s))
  expect_true(any(grepl("Gradient:\\s+auto \\(requested\\) -> ad \\(used\\)", out)))
})
test_that("print.ferx_summary falls back to (requested) only when used is NA", {
  s <- summary(make_fake_fit(gradient = "auto", gradient_used = NA_character_))
  out <- capture.output(print(s))
  expect_true(any(grepl("Gradient:\\s+auto \\(requested\\)$", out)))
  expect_false(any(grepl("\\(used\\)", out)))
})
test_that(".ferx_parse_setting_value converts booleans", {
  expect_identical(ferx:::.ferx_parse_setting_value("true"), TRUE)
  expect_identical(ferx:::.ferx_parse_setting_value("false"), FALSE)
})
test_that(".ferx_parse_setting_value converts numerics", {
  expect_equal(ferx:::.ferx_parse_setting_value("200"), 200)
  expect_equal(ferx:::.ferx_parse_setting_value("3.14"), 3.14)
})
test_that(".ferx_parse_setting_value converts null to NA", {
  expect_identical(ferx:::.ferx_parse_setting_value("null"), NA)
})
test_that(".ferx_parse_setting_value leaves strings as character", {
  expect_identical(ferx:::.ferx_parse_setting_value("slsqp"), "slsqp")
})
test_that("call_settings stores typed values after settings are serialised", {
  # Simulate what ferx_fit() does: serialise then parse back.
  parts <- ferx:::.ferx_settings_to_strings(list(
    optimizer = "slsqp",
    max_iter = 200L,
    scale_params = TRUE,
    stagnation_guard = FALSE,
    reconverge_gradient_interval = 5L
  ))
  typed <- setNames(
    lapply(parts$values, ferx:::.ferx_parse_setting_value),
    parts$keys
  )

  expect_identical(typed$optimizer, "slsqp")
  expect_equal(typed$max_iter, 200)
  expect_identical(typed$scale_params, TRUE)
  expect_identical(typed$stagnation_guard, FALSE)
  expect_equal(typed$reconverge_gradient_interval, 5)
})
test_that("reconverge_gradient_interval serialises as a plain integer string", {
  # Integer-valued settings must not pick up scientific notation or decimals
  # when forwarded to the Rust fit-option parser (which expects a usize).
  parts <- ferx:::.ferx_settings_to_strings(list(reconverge_gradient_interval = 10L))
  expect_identical(parts$keys, "reconverge_gradient_interval")
  expect_identical(parts$values, "10")

  # 0 is the documented "off" value and must round-trip, not become "" or NA.
  parts0 <- ferx:::.ferx_settings_to_strings(list(reconverge_gradient_interval = 0L))
  expect_identical(parts0$values, "0")
})

# ---- header from test-bayes.R ----
# Tier-1 tests for the Bayes (method = "bayes") R-side surface: the R-hat
# convergence flag, the print/summary rendering of `$bayes`, and the persist
# wire builder. These exercise the wrapper logic with a mocked posterior so
# they run in milliseconds — the real MCMC engine round-trip is covered by the
# (slow) integration tests in test-fit.R / test-persist.R.

# A minimal but complete `$bayes` posterior summary, shaped exactly as
# `fit_result_to_list()` (Rust) and `ferx_load_fit()` (persist) produce it.
fake_bayes <- function(max_rhat = 1.004) {
  list(
    n_chains          = 2,
    n_warmup          = 300,
    n_draws_per_chain = 600,
    n_divergent       = 0,
    max_rhat          = max_rhat,
    param_names       = c("TVCL", "TVV", "OMEGA_CL", "SIGMA_PROP"),
    mean              = c(0.133, 8.0, 0.09, 0.10),
    sd                = c(0.01, 0.4, 0.02, 0.01),
    q025              = c(0.115, 7.2, 0.05, 0.08),
    median            = c(0.133, 8.0, 0.09, 0.10),
    q975              = c(0.151, 8.8, 0.13, 0.12),
    rhat              = c(1.001, 1.002, 1.004, 1.000),
    ess_bulk          = c(820, 760, 540, 910),
    ess_tail          = c(800, 740, 520, 900),
    mcse              = c(3e-4, 0.01, 8e-4, 3e-4)
  )
}

# .ferx_bayes_rhat_flag ------------------------------------------------------



# print.ferx_fit — BAYES posterior block -------------------------------------




# summary.ferx_fit / print.ferx_summary --------------------------------------



# .fitrx_build_bayes_wire ----------------------------------------------------




test_that(".ferx_bayes_rhat_flag flags only when max R-hat exceeds the threshold", {
  # Threshold is 1.01; below/at it there is no marker.
  expect_identical(ferx:::.ferx_bayes_rhat_flag(1.004), "")
  expect_identical(ferx:::.ferx_bayes_rhat_flag(ferx:::.FERX_BAYES_RHAT_THRESHOLD), "")
  # Above the threshold a "[!]" marker is returned (colour codes may wrap it).
  flag <- ferx:::.ferx_bayes_rhat_flag(1.2)
  expect_true(nzchar(flag))
  expect_true(grepl("\\[!\\]", flag))
})
test_that(".ferx_bayes_rhat_flag returns no marker for non-finite R-hat", {
  # A degenerate sampler can yield NA/Inf R-hat; that must not crash or flag.
  expect_identical(ferx:::.ferx_bayes_rhat_flag(NA_real_), "")
  expect_identical(ferx:::.ferx_bayes_rhat_flag(Inf), "")
})
test_that("print.ferx_fit renders the BAYES posterior summary block", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1), bayes = fake_bayes())
  out <- capture.output(print(fit))

  expect_true(any(grepl("BAYES", out)))
  expect_true(any(grepl("posterior summary", out)))
  # Scalar metadata line.
  meta <- out[grepl("Chains:", out)]
  expect_length(meta, 1L)
  expect_true(grepl("Chains: 2", meta))
  expect_true(grepl("Warmup: 300", meta))
  expect_true(grepl("Draws/chain: 600", meta))
  expect_true(grepl("Divergent: 0", meta))
  expect_true(grepl("Max R-hat: 1.0040", meta))
  # Per-parameter table: header + one row per parameter.
  expect_true(any(grepl("Parameter", out) & grepl("R-hat", out) & grepl("ESS", out)))
  expect_true(any(grepl("TVCL", out)))
  expect_true(any(grepl("SIGMA_PROP", out)))
})
test_that("print.ferx_fit marks an unconverged BAYES fit with [!]", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1), bayes = fake_bayes(max_rhat = 1.25))
  out <- capture.output(print(fit))
  meta <- out[grepl("Max R-hat:", out)]
  expect_length(meta, 1L)
  expect_true(grepl("\\[!\\]", meta))
})
test_that("print.ferx_fit omits the BAYES block for non-Bayes fits", {
  out <- capture.output(print(make_fake_fit(omega = matrix(0.09, 1, 1))))
  expect_false(any(grepl("posterior summary", out)))
})
test_that("summary.ferx_fit carries $bayes and print.ferx_summary shows it", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1), bayes = fake_bayes())
  s <- summary(fit)
  expect_false(is.null(s$bayes))
  expect_identical(s$bayes$param_names, fit$bayes$param_names)

  out <- capture.output(print(s))
  line <- out[grepl("^Bayes:", out)]
  expect_length(line, 1L)
  expect_true(grepl("2 chains", line))
  expect_true(grepl("600 draws/chain", line))
  expect_true(grepl("max R-hat = 1.0040", line))
})
test_that("print.ferx_summary flags an unconverged BAYES fit", {
  s <- summary(make_fake_fit(omega = matrix(0.09, 1, 1), bayes = fake_bayes(max_rhat = 1.3)))
  out <- capture.output(print(s))
  expect_true(any(grepl("^Bayes:", out) & grepl("\\[!\\]", out)))
})

# ---- header from test-sde.R ----
# Tests for SDE (diffusion / EKF) surface on the R side.
# All tests use make_fake_fit() from helper-trace.R — no real fit needed.

# uses_sde flag in print.ferx_fit ----------------------------------------





# uses_sde flag in summary / print.ferx_summary ---------------------------






# DIFF_* theta parameters flow through naturally --------------------------


# ferx_estimates() with DIFF_* theta -------------------------------------


# ferx_save_fit / ferx_load_fit round-trip --------------------------------




test_that("print.ferx_fit shows '+ SDE (EKF)' when uses_sde is TRUE", {
  fit <- make_fake_fit(uses_sde = TRUE, omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  method_line <- out[grepl("Method:", out)]
  expect_length(method_line, 1L)
  expect_match(method_line, "SDE \\(EKF\\)")
})
test_that("print.ferx_fit does not show SDE tag when uses_sde is FALSE", {
  fit <- make_fake_fit(uses_sde = FALSE, omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("SDE", out)))
})
test_that("print.ferx_fit does not show SDE tag when uses_sde is absent", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("SDE", out)))
})
test_that("print.ferx_fit shows SDE tag after method chain when uses_sde TRUE", {
  fit <- make_fake_fit(
    uses_sde     = TRUE,
    method       = "focei",
    method_chain = c("saem", "focei"),
    omega        = matrix(0.09, 1, 1)
  )
  out <- capture.output(print(fit))
  chain_line <- out[grepl("Method:", out)]
  expect_length(chain_line, 1L)
  expect_match(chain_line, "SDE \\(EKF\\)")
})
test_that("summary.ferx_fit carries uses_sde = TRUE into ferx_summary", {
  fit <- make_fake_fit(uses_sde = TRUE)
  s   <- summary(fit)
  expect_true(isTRUE(s$uses_sde))
})
test_that("summary.ferx_fit carries uses_sde = FALSE into ferx_summary", {
  fit <- make_fake_fit(uses_sde = FALSE)
  s   <- summary(fit)
  expect_false(s$uses_sde)
})
test_that("summary.ferx_fit defaults uses_sde to FALSE when field absent", {
  fit <- make_fake_fit()
  s   <- summary(fit)
  expect_false(s$uses_sde)
})
test_that("print.ferx_summary shows '+ SDE (EKF)' when uses_sde TRUE", {
  fit <- make_fake_fit(uses_sde = TRUE)
  out <- capture.output(print(summary(fit)))
  method_line <- out[grepl("Method:", out)]
  expect_length(method_line, 1L)
  expect_match(method_line, "SDE \\(EKF\\)")
})
test_that("print.ferx_summary omits SDE tag when uses_sde FALSE", {
  fit <- make_fake_fit(uses_sde = FALSE)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("SDE", out)))
})
test_that("DIFF_* theta parameters appear in print.ferx_fit theta table", {
  fit <- make_fake_fit(
    uses_sde  = TRUE,
    theta     = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    se_theta  = NULL,
    omega     = matrix(0.09, 1, 1)
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("DIFF_CENTRAL", out)))
})
test_that("ferx_estimates() includes DIFF_* theta rows", {
  fit <- make_fake_fit(
    uses_sde = TRUE,
    theta    = c(TVCL = 5.0, TVV = 50.0, DIFF_CENTRAL = 0.5),
    omega    = matrix(0.09, 1, 1)
  )
  est <- ferx_estimates(fit)
  expect_true(any(est$param == "DIFF_CENTRAL"))
  diff_row <- est[est$param == "DIFF_CENTRAL", ]
  expect_equal(diff_row$estimate, 0.5)
})

# ---- header from test-covtab.R ----
# Covariate table (`fit$covtab`) + types (`fit$covariate_types`), surfaced from
# a `[covariates]` block. See ferx-core #182/#184 for the underlying DSL + table.

# Cached fit on the covariate example so the multiple tests below share one
# estimation (mirrors helper-warfarin-fit.R's memoised pattern).
cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("two_cpt_oral_cov")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method = "focei", verbose = FALSE,
        covariance = FALSE, settings = list(maxiter = 5L)
      )
    }
    fit
  }
})








test_that("covtab is a data frame: ID/TIME/EVID + declared covariates in order, one row per input record", {
  fit <- cov_fit()
  ex  <- ferx_example("two_cpt_oral_cov")

  expect_s3_class(fit$covtab, "data.frame")
  # Exact columns + order: ID, TIME, EVID, then declared covariates in
  # declaration order.
  expect_identical(names(fit$covtab), c("ID", "TIME", "EVID", "WT", "CRCL"))

  # One row per input dataset record (incl. dose / EVID rows), unlike sdtab.
  expect_equal(nrow(fit$covtab), nrow(utils::read.csv(ex$data)))
  expect_gt(nrow(fit$covtab), nrow(fit$sdtab)) # sdtab is observation rows only

  # Column types: ID character, TIME/covariates numeric, EVID integer.
  expect_type(fit$covtab$ID, "character")
  expect_type(fit$covtab$TIME, "double")
  expect_type(fit$covtab$WT, "double")
  expect_type(fit$covtab$EVID, "integer")

  # Dose rows are present (EVID == 1).
  expect_true(any(fit$covtab$EVID == 1L))
})
test_that("covariate_types carries the declared continuous/categorical tag", {
  fit <- cov_fit()
  expect_type(fit$covariate_types, "character")
  expect_identical(
    fit$covariate_types,
    c(WT = "continuous", CRCL = "continuous")
  )
})
test_that("a missing covariate cell surfaces as NA (not NaN) in covtab", {
  ex <- ferx_example("two_cpt_oral_cov")
  # Blank the WT value on the first record, write to a temp dataset, and fit.
  d <- utils::read.csv(ex$data, check.names = FALSE, colClasses = "character")
  d$WT[1] <- ""
  tmp_data <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_data), add = TRUE)
  utils::write.csv(d, tmp_data, row.names = FALSE, quote = FALSE)

  fit <- ferx_fit(
    ex$model, tmp_data,
    method = "focei", verbose = FALSE,
    covariance = FALSE, settings = list(maxiter = 2L)
  )
  expect_true(is.na(fit$covtab$WT[1]))
  expect_false(is.nan(fit$covtab$WT[1])) # specifically NA_real_, not NaN
})
test_that("covtab and covariate_types are NULL when the model declares no [covariates] block", {
  fit <- warfarin_fit() # warfarin example has no [covariates] block
  expect_null(fit$covtab)
  expect_null(fit$covariate_types)
})

# ---- header from test-diagnostics-more.R ----
# Tests for diagnostic functions that operate on a fit's fields and can be
# driven with crafted inputs (no live model fit needed):
#   ferx_estimates(), .ferx_est_row(), ferx_eta_cov(), ferx_cor_matrix(),
#   ferx_get_warnings(), and .ferx_compute_eta_normality().
# make_fake_fit() comes from helper-trace.R.

.est_row       <- getFromNamespace(".ferx_est_row",            "ferx")
.compute_norm  <- getFromNamespace(".ferx_compute_eta_normality", "ferx")

# ---------------------------------------------------------------------------
# ferx_estimates — theta (unnamed), scalar omega, sigma, and IOV kappa rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_eta_cov — message branches plus the correlation path
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# ferx_cor_matrix — non-positive diagonal warning
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# ferx_get_warnings — no-warnings branch and per-severity labels
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# .ferx_compute_eta_normality — NULL / empty / large-N branches
# ---------------------------------------------------------------------------




test_that(".ferx_compute_eta_normality returns NULL for missing or eta-less input", {
  expect_null(.compute_norm(NULL))
  expect_null(.compute_norm("not a df"))
  expect_null(.compute_norm(data.frame(ID = 1:3)))   # only ID, no eta columns
})
test_that(".ferx_compute_eta_normality subsamples and warns above the 5000 cap", {
  set.seed(1)
  big <- data.frame(ID = seq_len(5001L), ETA_CL = rnorm(5001L))
  expect_warning(res <- .compute_norm(big), "exceeds")
  expect_s3_class(res, "data.frame")
  expect_true(is.finite(res$W))
})
test_that(".ferx_compute_eta_normality flags non-normal etas", {
  set.seed(2)
  df <- data.frame(ID = 1:200, ETA_SKEW = rexp(200))  # clearly non-normal
  res <- .compute_norm(df)
  expect_true(res$p_val < 0.05)
  expect_identical(res$flag, "[!]")
})

# ---- header from test-map-estimates.R ----
# Per-subject EBE / individual-parameter outputs.
#
# Pins the contract that ferx_fit() returns three per-subject diagnostic
# tables — `ebe_etas`, `individual_estimates`, and (for IOV models)
# `ebe_kappas` — and that `sdtab` no longer carries ETA columns. Also
# checks the downstream consumers (`eta_normality`, `ferx_eta_cov`) read
# from `ebe_etas`.
#
# Uses the shared `warfarin_fit()` helper (helper-warfarin-fit.R) so the
# real FOCEI fit is reused across test files.







# Warfarin data has no non-NONMEM covariate columns, so a synthetic WT column
# is added to exercise the ferx_eta_cov() happy path. One value per subject
# (constant within subject) so it qualifies as a covariate.
warfarin_dat_with_cov <- function() {
  ex  <- ferx_example("warfarin")
  dat <- read.csv(ex$data)
  ids <- unique(dat$ID)
  set.seed(42L)
  wt_map <- setNames(rnorm(length(ids), mean = 70, sd = 10), ids)
  dat$WT <- wt_map[as.character(dat$ID)]
  dat
}





test_that("fit$ebe_etas has one row per subject with named ETA columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$ebe_etas, "data.frame")
  expect_equal(nrow(fit$ebe_etas), fit$n_subjects)
  expect_true("ID" %in% names(fit$ebe_etas))
  # Warfarin declares ETA_CL, ETA_V, ETA_KA in [parameters].
  expect_setequal(
    setdiff(names(fit$ebe_etas), "ID"),
    c("ETA_CL", "ETA_V", "ETA_KA")
  )
  expect_true(all(vapply(fit$ebe_etas[, -1], is.numeric, logical(1L))))
})
test_that("fit$individual_estimates has one row per subject with named parameter columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$individual_estimates, "data.frame")
  expect_equal(nrow(fit$individual_estimates), fit$n_subjects)
  expect_true("ID" %in% names(fit$individual_estimates))
  # Warfarin declares CL, V, KA in [individual_parameters].
  expect_setequal(
    setdiff(names(fit$individual_estimates), "ID"),
    c("CL", "V", "KA")
  )
  # Values come out of pk_param_fn — must be finite and positive for a
  # log-normal parameterisation.
  for (p in c("CL", "V", "KA")) {
    vals <- fit$individual_estimates[[p]]
    expect_true(all(is.finite(vals)), info = sprintf("%s has non-finite values", p))
    expect_true(all(vals > 0), info = sprintf("%s has non-positive values", p))
  }
})
test_that("ebe_etas and individual_estimates align on ID and row count", {
  fit <- warfarin_fit()
  expect_equal(fit$ebe_etas$ID, fit$individual_estimates$ID)
})
test_that("sdtab no longer contains ETA columns", {
  fit <- warfarin_fit()
  expect_s3_class(fit$sdtab, "data.frame")
  eta_cols <- grep("^ETA", names(fit$sdtab), value = TRUE)
  expect_length(eta_cols, 0L)
})
test_that("eta_normality reads from ebe_etas, one row per ETA", {
  fit <- warfarin_fit()
  expect_s3_class(fit$eta_normality, "data.frame")
  expect_setequal(fit$eta_normality$eta, c("ETA_CL", "ETA_V", "ETA_KA"))
  expect_true(all(c("eta", "W", "p_val", "flag") %in% names(fit$eta_normality)))
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



test_that("print.ferx_fit uses bare eta_names in OMEGA section", {
  fit <- make_fake_fit(
    omega         = matrix(0.09, 1, 1),
    eta_names     = "ETA_CL",
    eta_param_types = "log_normal"
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("ETA_CL", out)))
  expect_false(any(grepl("OMEGA\\(ETA_CL\\)", out)))
})
test_that("print.ferx_fit emits [log scale] tag for log-transformed theta", {
  fit <- make_fake_fit(
    theta            = c(TVCL = log(0.134)),
    se_theta         = c(TVCL = 0.15),
    theta_transforms = "log",
    omega            = matrix(0.07, 1, 1),
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  theta_line <- out[grepl("TVCL", out)]
  expect_true(any(grepl("log scale", theta_line)))
  typical_line <- out[grepl("typical", out)]
  expect_true(length(typical_line) >= 1L)
  expect_true(any(grepl("95%", typical_line)))
})
test_that("print.ferx_fit emits [logit scale] tag for logit-transformed theta", {
  fit <- make_fake_fit(
    theta            = c(THETA_F = log(0.7 / 0.3)),
    se_theta         = c(THETA_F = 0.3),
    theta_transforms = "logit",
    omega            = matrix(0.10, 1, 1),
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  theta_line <- out[grepl("THETA_F", out)]
  expect_true(any(grepl("logit scale", theta_line)))
})
test_that("print.ferx_fit shows [log-normal] label and CV% for log_normal ETA", {
  fit <- make_fake_fit(
    theta           = c(CL = 1),
    omega           = matrix(0.40, 1, 1),
    eta_param_types = "log_normal",
    sigma           = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("log-normal", omega_line)))
  expect_true(any(grepl("CV%", omega_line)))
  expect_true(any(grepl("70\\.1", omega_line)))
})
test_that("print.ferx_fit shows [additive] label and SD for additive ETA", {
  fit <- make_fake_fit(
    theta            = c(TVTLAG = 0.5),
    theta_transforms = "identity",
    omega            = matrix(0.04, 1, 1),
    eta_param_types  = "additive",
    eta_linked_theta = "TVTLAG",
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("additive", omega_line)))
  expect_true(any(grepl("SD = 0\\.2", omega_line)))
})
test_that("print.ferx_fit shows [logit] label and +/-1SD range for logit ETA", {
  logit_theta <- log(0.7 / 0.3)
  fit <- make_fake_fit(
    theta            = c(THETA_F = logit_theta),
    theta_transforms = "logit",
    omega            = matrix(0.10, 1, 1),
    eta_param_types  = "logit",
    eta_linked_theta = "THETA_F",
    sigma            = 0.01
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("logit", omega_line)))
  expect_true(any(grepl("+/-1SD", omega_line, fixed = TRUE)))
})
test_that("print.ferx_fit shows [proportional] label and CV% for proportional sigma", {
  fit <- make_fake_fit(
    theta       = c(CL = 1),
    omega       = matrix(0.07, 1, 1),
    sigma       = 0.01,
    sigma_types = "proportional"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("SIGMA\\(1\\)", out)]
  expect_true(any(grepl("proportional", sigma_line)))
  expect_true(any(grepl("CV%", sigma_line)))
})
test_that("print.ferx_fit shows [additive] label and SD for additive sigma", {
  fit <- make_fake_fit(
    theta       = c(CL = 1),
    omega       = matrix(0.07, 1, 1),
    sigma       = 1.0,
    sigma_types = "additive"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("SIGMA\\(1\\)", out)]
  expect_true(any(grepl("additive", sigma_line)))
  expect_true(any(grepl("SD = 1\\.0", sigma_line)))
})
test_that("print.ferx_fit works without eta_param_types (falls back to log_normal)", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), sigma = 0.01)
  # Remove transform fields to simulate old ferx-core binary
  fit$eta_param_types  <- NULL
  fit$theta_transforms <- NULL
  fit$sigma_types      <- NULL
  expect_no_error(capture.output(print(fit)))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(any(grepl("70\\.1", omega_line)))  # still shows log-normal CV%
})

# ---- header from test-data-selection.R ----
# Tests for [data_selection] / ferx_apply_selection() / ferx_fit(ignore=) support.

# Local alias avoids the ::: operator.
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# .sel_attr()/warfarin_sel_fit() come from helper-selection.R

# ---------------------------------------------------------------------------
# 1. Plain warfarin fit has NULL exclusions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 2. [data_selection] model populates exclusions
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 3. ferx_fit(data, ignore = ...) produces exclusions matching model-file path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 4. ferx_apply_selection() pure-R preview
# ---------------------------------------------------------------------------







# ---------------------------------------------------------------------------
# 5. ferx_data passed directly to ferx_fit()
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6. print.ferx_fit() shows DATA SELECTION block
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 7. ferx_runlog() shows exclusion line
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 8. persist round-trip preserves exclusions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 9. print.ferx_data
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 10. ferx_apply_selection(excluded = TRUE) on a ferx_data
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 11. accept= and ignore_ids= coverage for ferx_apply_selection()
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 12. ferx_apply_selection() data.frame input
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 13. ferx_apply_selection(fit, excluded = TRUE) replays fired conditions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 14. ferx_apply_selection(fit, excluded = TRUE) NULL exclusions path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 15. ferx_fit() type validation stops for ignore/accept/ignore_ids
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# 16. ferx_fit() accept= and ignore_ids= settings_parts code paths
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 17. print.ferx_fit fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 18. ferx_runlog fired_accept and excluded_subject_ids display paths
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 19. .cmp() operator coverage (pure R, via ferx_apply_selection() with inline df)
# ---------------------------------------------------------------------------

# .sel_test_df() comes from helper-selection.R





# ---------------------------------------------------------------------------
# 20. .parse_filter_expr() coverage: character rhs, &&, !found
# ---------------------------------------------------------------------------





test_that("print.ferx_fit shows DATA SELECTION when exclusions present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  out <- capture.output(print(fit))
  expect_true(any(grepl("DATA SELECTION", out, fixed = TRUE)))
  expect_true(any(grepl("Obs excl", out, fixed = TRUE)))
})
test_that("print.ferx_fit does NOT show DATA SELECTION for plain fit", {
  skip_on_cran()
  fit <- warfarin_fit()
  out <- capture.output(print(fit))
  expect_false(any(grepl("DATA SELECTION", out, fixed = TRUE)))
})
test_that("print.ferx_fit shows fired_accept and excluded_subject_ids when present", {
  skip_on_cran()
  fit <- warfarin_sel_fit()
  # Augment the exclusions in-place to exercise the display branches.
  fit$exclusions$fired_accept          <- c("accept: DV > -999")
  fit$exclusions$excluded_subject_ids  <- c("3", "5")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Fired accept", out, fixed = TRUE)))
  expect_true(any(grepl("Subjects excluded", out, fixed = TRUE)))
})

# ---- header from test-derived-output.R ----
# Tests for [derived] and [output] blocks in ferx model files.
#
# These are integration tests that require a real fit; they are skipped on
# CRAN. The `warfarin_derived_fit()` helper in helper-warfarin-fit.R caches
# the result so the fit runs at most once per test session.
#
# Coverage:
#   - [derived] columns (per-row and aggregate) appear in fit$sdtab
#   - [output] columns (individual PK parameters) appear in fit$sdtab
#   - ETA columns are absent from sdtab (regression: ferx-core#188)
#   - fit$covariate_names is populated when the model uses covariates
#   - Derived column values are numerically plausible

# -- [derived] columns present in sdtab ------------------------------------




# -- [output] columns present in sdtab ------------------------------------



# -- ETAs absent from sdtab (regression guard) ----------------------------


# -- ETA values are in ebe_etas, accessible for merging -------------------



# -- Numerical plausibility -----------------------------------------------





test_that("[derived] per-row columns appear in sdtab", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # KE and T_HALF are per-row expressions declared in warfarin_derived.ferx
  expect_true("KE"     %in% names(fit$sdtab),
              label = "KE (per-row [derived]) must be in sdtab")
  expect_true("T_HALF" %in% names(fit$sdtab),
              label = "T_HALF (per-row [derived]) must be in sdtab")
  expect_true(all(is.finite(fit$sdtab$KE)),     label = "KE must be finite")
  expect_true(all(is.finite(fit$sdtab$T_HALF)), label = "T_HALF must be finite")
})
test_that("[derived] aggregate columns appear in sdtab and are subject-constant", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # CMAX and AUC_0_72 are aggregate [derived] columns
  expect_true("CMAX"     %in% names(fit$sdtab))
  expect_true("AUC_0_72" %in% names(fit$sdtab))
  # Aggregates are constant within a subject (repeat the same value)
  for (col in c("CMAX", "AUC_0_72")) {
    by_id <- tapply(fit$sdtab[[col]], fit$sdtab$ID, function(x) length(unique(x)))
    expect_true(all(by_id == 1L),
                label = sprintf("%s must be constant within each subject", col))
  }
})
test_that("[derived] TAFD/TAD auxiliary columns are present and positive", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_true("DAY"      %in% names(fit$sdtab))
  expect_true("TAU_TIME" %in% names(fit$sdtab))
  expect_true(all(fit$sdtab$DAY >= 1L))
})
test_that("[output] individual PK parameters appear in sdtab", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  # warfarin_derived.ferx has [output] CL V KA
  for (col in c("CL", "V", "KA")) {
    expect_true(col %in% names(fit$sdtab),
                label = sprintf("[output] column %s must be in sdtab", col))
    expect_true(all(is.finite(fit$sdtab[[col]])),
                label = sprintf("[output] column %s must be finite", col))
    expect_true(all(fit$sdtab[[col]] > 0),
                label = sprintf("[output] column %s must be positive", col))
  }
})
test_that("[output] PK params are constant within subject (no random noise)", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  for (col in c("CL", "V", "KA")) {
    by_id <- tapply(fit$sdtab[[col]], fit$sdtab$ID, function(x) length(unique(x)))
    expect_true(all(by_id == 1L),
                label = sprintf("%s must be constant within each subject", col))
  }
})
test_that("[derived] fit: sdtab does not contain ETA columns", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  eta_cols <- grep("^ETA_|^ETA[0-9]", names(fit$sdtab), value = TRUE)
  expect_length(eta_cols, 0L)
})
test_that("fit$ebe_etas has ETA columns matching model declarations", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_s3_class(fit$ebe_etas, "data.frame")
  expect_true("ID" %in% names(fit$ebe_etas))
  # warfarin_derived declares ETA_CL, ETA_V, ETA_KA
  for (eta in c("ETA_CL", "ETA_V", "ETA_KA")) {
    expect_true(eta %in% names(fit$ebe_etas),
                label = sprintf("%s must be in ebe_etas", eta))
  }
  expect_equal(nrow(fit$ebe_etas), fit$n_subjects)
})
test_that("ebe_etas can be merged onto sdtab by ID for downstream analysis", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  subj <- fit$sdtab[!duplicated(fit$sdtab$ID), ]
  merged <- merge(subj, fit$ebe_etas, by = "ID")
  expect_equal(nrow(merged), fit$n_subjects)
  expect_true("ETA_CL" %in% names(merged))
})
test_that("T_HALF and KE are numerically consistent: T_HALF == 0.693 / KE", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_equal(fit$sdtab$T_HALF, 0.6931472 / fit$sdtab$KE,
               tolerance = 1e-6,
               label = "T_HALF must equal 0.6931472 / KE within tolerance")
})
test_that("CL / V == KE within tolerance", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_equal(fit$sdtab$CL / fit$sdtab$V, fit$sdtab$KE,
               tolerance = 1e-6,
               label = "KE must equal CL / V")
})
test_that("CMAX is non-negative and >= all IPRED values per subject", {
  skip_on_cran()
  fit <- warfarin_derived_fit()
  expect_true(all(fit$sdtab$CMAX >= 0))
  max_ipred <- tapply(fit$sdtab$IPRED, fit$sdtab$ID, max)
  cmax_val  <- tapply(fit$sdtab$CMAX, fit$sdtab$ID, unique)
  expect_true(all(unlist(cmax_val) >= unlist(max_ipred) - 1e-6),
              label = "CMAX must be >= max(IPRED) for every subject")
})

test_that("ferx_fit inits_from_nca arg validation", {
  expect_error(
    ferx_fit(ferx_example("warfarin")$model,
             ferx_example("warfarin")$data,
             inits_from_nca = NA),
    "TRUE/FALSE"
  )
  expect_error(
    ferx_fit(ferx_example("warfarin")$model,
             ferx_example("warfarin")$data,
             inits_from_nca = "bogus"),
    "should be one of"
  )
})

# ---- header from test-print-fit.R ----
# make_fake_fit() lives in helper-trace.R and is auto-loaded by testthat.

# print.ferx_fit CV% formula --------------------------------------------





# SIGMA display (ferx#59 / ferx-core#57): sigma is on the SD scale, so
# print() must show variance = sigma^2 always and CV% = sigma * 100 only
# for proportional components.





# Block-omega correlations (#60): print uses omega_param_corr when the engine
# provides it (bivariate-lognormal formula), and falls back to the eta-level
# Pearson formula when it is absent.





# NN-aware output (Option E, ferx-core#52): per-weight thetas are hidden from
# the THETA table and a compact NEURAL NETWORKS block summarises them.



# --- init_as_sd annotation tests ---





# STATUS line surfaces failing critical category --------------------------









test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is absent (fallback)", {
  # No eta_param_types field -> defaults to log_normal.
  # omega = 0.40: exact CV% = sqrt(exp(0.40) - 1) * 100 ≈ 70.1 (not 63.2 from sqrt(0.40)*100)
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(length(omega_line) == 1L)
  expect_false(grepl("63\\.2", omega_line))  # old approximate value
  expect_true(grepl("70\\.1", omega_line))   # exact log-normal CV%
})
test_that("print.ferx_fit uses exact log-normal CV% when eta_param_types is 'log_normal'", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "log_normal")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("70\\.1", omega_line))
})
test_that("print.ferx_fit shows SD_logit for logit eta_param_types", {
  fit <- make_fake_fit(omega = matrix(0.40, 1, 1), eta_param_types = "logit")
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("logit", omega_line))
  expect_true(grepl("SD_logit", omega_line))
  # CV% is not meaningful for logit ETAs — must not appear
  expect_false(grepl("CV% =", omega_line))
})
test_that("print.ferx_fit CV% equals zero when omega diagonal is zero", {
  fit <- make_fake_fit(omega = matrix(0, 1, 1))
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("CV% = 0\\.0", omega_line))
})
test_that("print.ferx_fit shows variance + CV% for proportional sigma using its declared name", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = "PROP_ERR",
    sigma_types = "proportional"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("PROP_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.1^2 = 0.01; CV% = 0.1 * 100 = 10.0.
  expect_true(grepl("var = 0\\.010000", sigma_line))
  expect_true(grepl("CV% = 10\\.0", sigma_line))
  # No legacy SIGMA(1) label when a name is supplied.
  expect_false(any(grepl("SIGMA\\(1\\)", out)))
})
test_that("print.ferx_fit shows variance but no CV% for additive sigma", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.5,
    sigma_names = "ADD_ERR",
    sigma_types = "additive"
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("ADD_ERR", out)]
  expect_length(sigma_line, 1L)
  # variance = 0.5^2 = 0.25; no CV% on observation-unit scale.
  expect_true(grepl("var = 0\\.250000", sigma_line))
  expect_false(grepl("CV%", sigma_line))
})
test_that("print.ferx_fit handles combined error: CV% on prop component only", {
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = c(0.2, 0.5),
    sigma_names = c("PROP_ERR", "ADD_ERR"),
    sigma_types = c("proportional", "additive")
  )
  out <- capture.output(print(fit))
  prop_line <- out[grepl("PROP_ERR", out)]
  add_line  <- out[grepl("ADD_ERR",  out)]
  expect_length(prop_line, 1L)
  expect_length(add_line,  1L)
  expect_true(grepl("CV% = 20\\.0", prop_line))
  expect_false(grepl("CV%", add_line))
  expect_true(grepl("var = 0\\.040000", prop_line))
  expect_true(grepl("var = 0\\.250000", add_line))
})
test_that("print.ferx_fit falls back to SIGMA(i) when sigma_names is missing", {
  # sigma_names absent (older Rust binary or unit test that didn't pass them).
  fit <- make_fake_fit(
    omega       = matrix(0.09, 1, 1),
    sigma       = 0.1,
    sigma_names = NULL,
    sigma_types = NULL
  )
  out <- capture.output(print(fit))
  expect_true(any(grepl("SIGMA\\(1\\)", out)))
})
test_that("print.ferx_fit uses omega_param_corr value and 'param corr' label when present", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = pc,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
})
test_that("print.ferx_fit falls back to eta-level corr label when omega_param_corr is NULL", {
  # cov = 0.025, vars = 0.10 -> Pearson corr = 0.025 / 0.10 = 0.25
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega            = om,
    omega_param_corr = NULL,
    eta_param_types  = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("\\(corr = 0\\.2500\\)", corr_line))
  expect_false(grepl("param corr", corr_line))
})
test_that("print.ferx_fit uses omega_iov_param_corr when present for IOV correlations", {
  om     <- matrix(0.10, 1, 1)
  iov    <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  iov_pc <- matrix(c(1.0, 0.5227, 0.5227, 1.0), 2, 2)
  fit <- make_fake_fit(
    omega                 = om,
    omega_iov             = iov,
    omega_iov_param_corr  = iov_pc,
    kappa_names           = c("KAPPA1", "KAPPA2"),
    se_kappa              = NULL,
    shrinkage_kappa       = NULL
  )
  out <- capture.output(print(fit))
  corr_line <- out[grepl("KAPPA2 ~ KAPPA1", out)]
  expect_length(corr_line, 1L)
  expect_true(grepl("param corr = 0\\.5227", corr_line))
})
test_that("print.ferx_fit uses exact log-normal CV% for OMEGA_IOV when kappa_param_types absent", {
  # kappa = 0.20 -> exact CV% = sqrt(exp(0.20) - 1) * 100 ≈ 47.1 (not 44.7 from sqrt(0.20)*100)
  fit <- make_fake_fit(
    omega           = matrix(0.10, 1, 1),
    omega_iov       = matrix(0.20, 1, 1),
    kappa_names     = "KAPPA1",
    se_kappa        = NULL,
    shrinkage_kappa = NULL
  )
  out <- capture.output(print(fit))
  kappa_line <- out[grepl("KAPPA1", out)]
  expect_true(length(kappa_line) >= 1L)
  expect_false(grepl("44\\.7", kappa_line[1]))  # old approximate value
  expect_true(grepl("47\\.1", kappa_line[1]))   # exact log-normal CV%
})
test_that("print.ferx_fit skips NN-weight thetas and emits NEURAL NETWORKS block", {
  # 2 baseline thetas (TVCL, TVV1) + 4 NN weights at offset=2.
  theta <- c(TVCL = 4.0, TVV1 = 40.0, W1 = 0.1, W2 = -0.2, W3 = 0.3, W4 = 0.0)
  nn <- list(list(
    name              = "TYPICAL_PK",
    shape             = c(2L, 2L, 2L),
    hidden_activation = "tanh",
    output_activation = "exp",
    n_weights         = 4L,
    weights_offset    = 2L,
    input_names       = c("WT", "CRCL"),
    output_names      = c("CL", "V1"),
    weights           = c(0.1, -0.2, 0.3, 0.0)
  ))
  fit <- make_fake_fit(
    omega           = matrix(0.10, 1, 1),
    theta           = theta,
    neural_networks = nn
  )
  out <- capture.output(print(fit))

  # Baseline thetas appear in the THETA table.
  expect_true(any(grepl("^TVCL\\s", out)))
  expect_true(any(grepl("^TVV1\\s", out)))
  # NN-weight thetas are skipped from the THETA table (no `W1`/`W2`/... rows).
  theta_block <- out[seq(
    which(grepl("^\\s*THETA\\s*$|--- THETA", out))[1L],
    (which(grepl("NEURAL NETWORKS", out)) - 1L)[1L]
  )]
  expect_false(any(grepl("^W[1-4]\\s", theta_block)))

  # NEURAL NETWORKS block contains the network metadata and weight summary.
  expect_true(any(grepl("NEURAL NETWORKS", out)))
  expect_true(any(grepl("TYPICAL_PK.*shape=\\[2, 2, 2\\].*activation=tanh/exp.*n_weights=4", out)))
  expect_true(any(grepl("inputs:\\s+\\[WT, CRCL\\]", out)))
  expect_true(any(grepl("outputs:\\s+\\[CL, V1\\]", out)))
  expect_true(any(grepl("weights: min -0\\.2000\\s+max 0\\.3000", out)))
})
test_that("print.ferx_fit omits NEURAL NETWORKS block when neural_networks is empty/NULL", {
  fit <- make_fake_fit(omega = matrix(0.10, 1, 1))
  out <- capture.output(print(fit))
  expect_false(any(grepl("NEURAL NETWORKS", out)))
})
test_that("print.ferx_fit annotates omega with [initial specified as SD] when flag is TRUE", {
  fit <- make_fake_fit(
    omega           = matrix(0.09, 1, 1),
    eta_names       = "ETA_CL",
    omega_init_as_sd = TRUE
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("ETA_CL", out)]
  expect_true(length(omega_line) > 0L)
  expect_true(any(grepl("\\[initial specified as SD\\]", omega_line)))
})
test_that("print.ferx_fit does NOT annotate omega when flag is FALSE", {
  fit <- make_fake_fit(
    omega           = matrix(0.09, 1, 1),
    eta_names       = "ETA_CL",
    omega_init_as_sd = FALSE
  )
  out <- capture.output(print(fit))
  expect_false(any(grepl("\\[initial specified as SD\\]", out)))
})
test_that("print.ferx_fit annotates sigma with [initial specified as SD] when flag is TRUE", {
  fit <- make_fake_fit(
    omega            = matrix(0.09, 1, 1),
    sigma            = c(PROP_ERR = 0.01),
    sigma_names      = "PROP_ERR",
    sigma_types      = "proportional",
    sigma_init_as_sd = TRUE
  )
  out <- capture.output(print(fit))
  sigma_line <- out[grepl("PROP_ERR", out)]
  expect_true(length(sigma_line) > 0L)
  expect_true(any(grepl("\\[initial specified as SD\\]", sigma_line)))
})
test_that("print.ferx_fit handles missing omega_init_as_sd gracefully (old fit object)", {
  fit <- make_fake_fit(omega = matrix(0.09, 1, 1), eta_names = "ETA_CL")
  # omega_init_as_sd not set at all — should not error, no annotation
  out <- capture.output(print(fit))
  expect_false(any(grepl("\\[initial specified as SD\\]", out)))
})
test_that("print.ferx_fit appends critical category to STATUS when fit did not converge", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = c("critical", "info"),
      category = c("convergence", "mu_referencing"),
      message  = c("Outer optimization did not converge",
                   "mu-ref: CL, V"),
      source_method = c("", ""),
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_length(status_line, 1L)
  expect_true(grepl("NOT CONVERGED", status_line))
  # Inline category appears between the label and the iteration/time tail.
  expect_true(grepl("\\(convergence\\)", status_line))
})
test_that("print.ferx_fit deduplicates and joins multiple critical categories", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = c("critical", "critical", "critical"),
      category = c("convergence", "covariance_step", "convergence"),
      message  = c("a", "b", "c"),
      source_method = c("", "", ""),
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("\\(convergence, covariance_step\\)", status_line))
})
test_that("print.ferx_fit STATUS has no inline category when there are no critical warnings", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = TRUE,
    warnings_structured = data.frame(
      severity = "info",
      category = "mu_referencing",
      message  = "mu-ref: CL",
      source_method = "",
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("CONVERGED", status_line))
  # No parenthesised category snippet because nothing is critical.
  expect_false(grepl("\\(mu_referencing\\)|\\(convergence\\)", status_line))
})
test_that("print.ferx_fit STATUS has no inline category when non-converged fit has only info warnings", {
  fit <- make_fake_fit(
    omega     = matrix(0.10, 1, 1),
    converged = FALSE,
    warnings_structured = data.frame(
      severity = "info",
      category = "mu_referencing",
      message  = "mu-ref: CL",
      source_method = "",
      stringsAsFactors = FALSE
    )
  )
  out <- capture.output(print(fit))
  status_line <- out[grepl("STATUS", out)]
  expect_true(grepl("NOT CONVERGED", status_line))
  expect_false(grepl("\\(", status_line))
})

# se_omega display (#143): print() shows "SE = ..." for omega diagonal and
# block off-diagonal elements, reading via .omega_se_at(). Block omega carries
# a full lower-triangle se_omega; diagonal omega carries one SE per variance.
test_that("print.ferx_fit shows SE for omega diagonal when se_omega is present", {
  fit <- make_fake_fit(
    omega    = matrix(0.40, 1, 1),
    se_omega = 0.05
  )
  out <- capture.output(print(fit))
  omega_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_length(omega_line, 1L)
  expect_true(grepl("SE = 0\\.050000", omega_line))
})
test_that("print.ferx_fit shows off-diagonal SE for block omega (full lower-triangle se_omega)", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  # Full lower-triangle, column-major: [(1,1), (2,1), (2,2)]
  fit <- make_fake_fit(
    omega    = om,
    se_omega = c(0.011, 0.022, 0.033),
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  # Diagonal SE on the OMEGA(1,1) line.
  diag_line <- out[grepl("OMEGA\\(1,1\\)", out)]
  expect_true(grepl("SE = 0\\.011000", diag_line))
  # Off-diagonal covariance line carries the (2,1) SE.
  cov_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(cov_line, 1L)
  expect_true(grepl("SE = 0\\.022000", cov_line))
})
test_that("print.ferx_fit omits off-diagonal SE when se_omega is diagonal-only", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega    = om,
    se_omega = c(0.011, 0.033),  # length n_eta -> diagonal-only, no cov SE
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  cov_line <- out[grepl("OMEGA\\(2,1\\) : cov =", out)]
  expect_length(cov_line, 1L)
  expect_false(grepl("SE =", cov_line))
})
test_that("print.ferx_fit shows off-diagonal SE with named etas (block omega)", {
  om <- matrix(c(0.10, 0.025, 0.025, 0.10), 2, 2)
  fit <- make_fake_fit(
    omega     = om,
    se_omega  = c(0.011, 0.022, 0.033),
    eta_names = c("ETA_CL", "ETA_V"),
    eta_param_types = c("log_normal", "log_normal")
  )
  out <- capture.output(print(fit))
  # Loop emits the higher-index eta first (i outer, j < i), so the (2,1)
  # pair prints as "ETA_V ~ ETA_CL".
  cov_line <- out[grepl("ETA_V ~ ETA_CL : cov =", out)]
  expect_length(cov_line, 1L)
  expect_true(grepl("SE = 0\\.022000", cov_line))
})

# ---- header from test-sir.R ----
# Tests for ferx_sir() and the path/hash provenance plumbing that backs it.
#
# SIR happy-path tests are gated on `skip_if(is.null(fit$cov_matrix), ...)`
# because in the no-autodiff CI build the warfarin cov step occasionally
# fails to converge with maxiter = 30 (same pattern as test-simulate-uncertainty.R).

sir_cov_skip <- "covariance step did not converge - skipping"









test_that("ferx_fit records model_path, data_path, and hashes on the fit", {
  fit <- warfarin_fit_cov()
  ex <- ferx_example("warfarin")

  expect_true(!is.na(fit$model_path))
  expect_true(!is.na(fit$data_path))
  # Hashes are 64-char lowercase hex (SHA-256).
  expect_match(fit$model_hash, "^[0-9a-f]{64}$")
  expect_match(fit$data_hash, "^[0-9a-f]{64}$")

  # Paths point at the actual on-disk files used for the fit.
  expect_true(file.exists(fit$model_path))
  expect_true(file.exists(fit$data_path))
  expect_identical(normalizePath(fit$model_path),
                   normalizePath(ex$model))
  expect_identical(normalizePath(fit$data_path),
                   normalizePath(ex$data))
})

# ---------------------------------------------------------------------------
# .ferx_parse_model_fit_options
# ---------------------------------------------------------------------------

write_fit_options_model <- function(lines) {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)",
               "[fit_options]", lines), path)
  path
}

test_that("ferx:::.ferx_parse_model_fit_options() returns named char vec of raw values", {
  path <- write_fit_options_model(c("  method = focei", "  maxiter = 300"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "focei", maxiter = "300"))
})
test_that("ferx:::.ferx_parse_model_fit_options() lower-cases keys, preserves value case", {
  path <- write_fit_options_model(c("  Method = FOCE", "  Covariance = TRUE"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(names(opts), c("method", "covariance"))
  expect_equal(unname(opts), c("FOCE", "TRUE"))
})
test_that("ferx:::.ferx_parse_model_fit_options() strips comments via .ferx_extract_blocks()", {
  path <- write_fit_options_model(c(
    "  method = foce  # the FOCE method",
    "  # a full-line comment",
    "  maxiter = 300 // also a comment"
  ))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "foce", maxiter = "300"))
})
test_that("ferx:::.ferx_parse_model_fit_options() preserves embedded `=` in values", {
  path <- write_fit_options_model("  custom = a=b=c")
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(unname(opts["custom"]), "a=b=c")
})
test_that("ferx:::.ferx_parse_model_fit_options() returns empty char vec when block absent", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)"), path)
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_length(opts, 0L)
  expect_type(opts, "character")
})
test_that("ferx:::.ferx_parse_model_fit_options() returns empty char vec when block is empty", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVCL(1.0, 0.001, 100.0)",
               "[fit_options]"), path)
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_length(opts, 0L)
})
test_that("ferx:::.ferx_parse_model_fit_options() ignores malformed lines (no `=`)", {
  path <- write_fit_options_model(c("  method = foce", "  garbage_no_equals"))
  on.exit(unlink(path))

  opts <- ferx:::.ferx_parse_model_fit_options(path)
  expect_equal(opts, c(method = "foce"))
})
