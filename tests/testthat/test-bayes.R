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

# print.ferx_fit — BAYES posterior block -------------------------------------

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

# summary.ferx_fit / print.ferx_summary --------------------------------------

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

# .fitrx_build_bayes_wire ----------------------------------------------------

test_that(".fitrx_build_bayes_wire returns NULL when there is no posterior", {
  expect_null(ferx:::.fitrx_build_bayes_wire(list(bayes = NULL)))
  expect_null(ferx:::.fitrx_build_bayes_wire(list()))
})

test_that(".fitrx_build_bayes_wire serialises every posterior field", {
  wire <- ferx:::.fitrx_build_bayes_wire(list(bayes = fake_bayes()))
  expect_false(is.null(wire))
  # Scalars survive as plain numerics; param_names become a JSON-friendly list.
  expect_equal(wire$n_chains, 2)
  expect_equal(wire$n_draws_per_chain, 600)
  expect_equal(wire$max_rhat, 1.004)
  expect_equal(unlist(wire$param_names), c("TVCL", "TVV", "OMEGA_CL", "SIGMA_PROP"))
  # Per-parameter vectors round-trip through the unwrap helpers unchanged.
  expect_equal(ferx:::.fitrx_unwrap_opt_num_vec(wire$mean), fake_bayes()$mean)
  expect_equal(ferx:::.fitrx_unwrap_opt_num_vec(wire$ess_bulk), fake_bayes()$ess_bulk)
})
