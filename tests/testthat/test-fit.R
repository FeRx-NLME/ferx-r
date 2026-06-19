# Local aliases avoid the ::: operator (undesirable_operator_linter).
ferx_rust_autodiff_enabled <- getFromNamespace("ferx_rust_autodiff_enabled", "ferx")

# Return structure — Tier 1

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

# These tests require the covariance step to have succeeded. With maxiter = 30L
# the outer optimisation may not converge on all machines, so we skip rather
# than fail — a skip here means "covariance step needs more iterations", not
# a bug. A full-convergence run (no maxiter cap) is tested manually / locally.
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

test_that("ferx_fit rejects malformed `imp` method chains", {
  # Tier 1: surface the engine's `imp`-placement constraints to the R caller
  # *before* the engine round-trip, so the error message references the R
  # argument and avoids spinning up the backend on inputs that will be
  # rejected anyway. Engine-side guards remain as a safety net (and are
  # covered in the ferx-core integration tests).
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("imp", "focei")),
    regexp = "final stage|must be the final stage",
    ignore.case = TRUE
  )
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "imp", "focei")),
    regexp = "final stage|must be the final stage",
    ignore.case = TRUE
  )
  expect_error(
    ferx_fit(ex$model, ex$data, method = c("focei", "imp", "imp")),
    regexp = "at most once",
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

# fd_hessian_step argument — Tier 1 (R-side validation, no Rust call)

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

# ferx_check_init — Tier 1

test_that("ferx_check_init returns without error on valid inputs", {
  ex <- ferx_example("warfarin")
  expect_no_error(ferx_check_init(ex$model, ex$data))
})

test_that("ferx_check_init errors on missing model file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_check_init("no_such_model.ferx", ex$data),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

test_that("ferx_check_init errors on missing data file", {
  ex <- ferx_example("warfarin")
  expect_error(
    ferx_check_init(ex$model, "no_such_data.csv"),
    regexp = "not found|No such file|does not exist|is not TRUE",
    ignore.case = TRUE
  )
})

# ferx_cor_matrix — Tier 1, requires covariance = TRUE

test_that("ferx_cor_matrix returns a matrix with same dimnames as $cov_matrix", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  expect_true(is.matrix(cor))
  expect_equal(dimnames(cor), dimnames(fit$cov_matrix))
})

test_that("ferx_cor_matrix diagonal is all 1s", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  expect_true(all(diag(cor) == 1))
})

test_that("ferx_cor_matrix off-diagonal values are in [-1, 1]", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge — skipping")
  cor <- ferx_cor_matrix(fit)
  d <- nrow(cor)
  if (d > 1L) {
    off <- cor[row(cor) != col(cor)]
    expect_true(all(abs(off) <= 1))
  }
})

test_that("ferx_cor_matrix errors gracefully when covariance was FALSE", {
  fit <- warfarin_fit()
  expect_error(ferx_cor_matrix(fit))
})

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

test_that("autodiff OFV gradient matches finite-difference gradient within tolerance", {
  skip_if(
    !isTRUE(ferx_rust_autodiff_enabled()),
    "Enzyme autodiff not available — skipping gradient-correctness test"
  )
  skip("Gradient inspection API not yet exposed from ferx-core — implement once available")
})

# -- Diagnostic fields (Step 3 feature-parity) --------------------------------

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

# -- IOV: shrinkage_kappa_by_occ (Step 1 feature-parity) ---------------------

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

# SAEM HMC proposals — [ENZYME ONLY]
#
# HMC proposals in the SAEM E-step are gated on the Enzyme autodiff build
# (`hmc_step` is `#[cfg(feature = "autodiff")]` in ferx-core). On a stable /
# FERX_NO_AUTODIFF=1 build, `n_leapfrog > 0` is silently ignored and the
# sampler uses Metropolis-Hastings, so `saem_n_subjects_hmc` is NA. The
# guard below makes this test inert on Tier 1 machines; CI with the Enzyme
# toolchain exercises the behavioural assertions.

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
