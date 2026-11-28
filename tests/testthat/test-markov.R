# [markov_model] from R (#788). The package builds ferx-core with its `markov`
# feature, so a continuous-time (CTMM) or discrete-time (DTMM) Markov endpoint
# fits from R as it does from the `ferx` command line. The oracles are NONMEM
# 7.5.1 on the same likelihoods, recorded in ferx-core at the pinned revision
# (tests/markov_nonmem.rs, tests/dtmm_nonmem_anchor.rs), and for the DTMM the
# exact closed-form Laplace objective
# (nonmem_anchor/dtmm_2state_exact_laplace.py). The fixtures under
# fixtures/markov/ are copies of ferx-core's models and data.

markov_fixture <- function(name) test_path("fixtures", "markov", name)

test_that("the shipped build accepts [markov_model] (#788)", {
  expect_true("markov_model" %in% ferx:::ferx_rust_known_blocks())
  out <- capture.output(
    res <- ferx_model_validate(markov_fixture("ctmm_2state.ferx"))
  )
  expect_true(isTRUE(res$ok))
  expect_false("E_BLOCK_FEATURE_DISABLED" %in% res$diagnostics$code)
  expect_false(any(grepl("feature not enabled", out, fixed = TRUE)))
})

test_that("a 2-state CTMM fit matches NONMEM METHOD=0 (#788)", {
  # NONMEM: LQ01 = -0.540113, LQ10 = -1.215660, OFV = 606.98181. The OFV is
  # flat near the optimum, so ferx's thetas agree to ~0.001 (SE ~0.14).
  fit <- ferx_fit(markov_fixture("ctmm_2state.ferx"),
                  markov_fixture("ctmm_2state.csv"), verbose = FALSE)
  expect_s3_class(fit, "ferx_fit")
  expect_true(isTRUE(fit$converged))
  expect_lt(abs(fit$ofv - 606.98181), 0.01)
  expect_named(fit$theta, c("LQ01", "LQ10"))
  expect_lt(max(abs(unname(fit$theta) - c(-0.540113, -1.215660))), 0.01)
  # A fixed-effects model: an empty omega matrix, finite theta SEs.
  expect_true(is.matrix(fit$omega))
  expect_identical(nrow(fit$omega), 0L)
  expect_true(all(is.finite(fit$se_theta)))
})

test_that("a 2-state DTMM at NONMEM's estimates reproduces its OFV (#788)", {
  # maxiter = 0: one evaluation at NONMEM's final estimates. NONMEM's LAPLACE
  # OFV is 1049.4031789; the exact Laplace objective is 1049.4079889, which
  # ferx's sits 2e-4 from (Linux and macOS alike).
  fit <- ferx_fit(markov_fixture("dtmm_2state_at_nonmem.ferx"),
                  markov_fixture("dtmm_2state.csv"), verbose = FALSE)
  expect_lt(abs(fit$ofv - 1049.4079889), 1e-3)
  expect_lt(abs(fit$ofv - 1049.4031789), 0.01)
  expect_named(fit$theta, c("A_AS", "B_DOSE", "A_SA"))
  expect_identical(rownames(fit$omega), c("ETA_AS", "ETA_SA"))

  # Each subject's EBEs, under the declared eta names.
  dat <- read.csv(markov_fixture("dtmm_2state.csv"))
  expect_identical(fit$ebe_etas$ID, as.character(unique(dat$ID)))
  expect_named(fit$ebe_etas, c("ID", "ETA_AS", "ETA_SA"))
  expect_true(all(is.finite(as.matrix(fit$ebe_etas[, -1L]))))

  # The engine builds per-record diagnostics only for the binary and count
  # endpoints (ferx-core categorical::binary_diagnostics); a Markov record has
  # none yet (CTMM rows: ferx-core #820), so sdtab is empty and the
  # per-subject discrete count is 0. When the engine starts emitting them this
  # reddens: check that the rows reach sdtab and survive a .fitrx round trip,
  # then assert them here.
  expect_identical(nrow(fit$sdtab), 0L)
  expect_true(all(fit$subject_summary$n_obs_discrete == 0L))
})

test_that("a Markov fit survives ferx_save_fit()/ferx_load_fit() (#788)", {
  fit <- ferx_fit(markov_fixture("dtmm_2state_at_nonmem.ferx"),
                  markov_fixture("dtmm_2state.csv"), verbose = FALSE)
  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out))
  ferx_save_fit(fit, out)
  back <- ferx_load_fit(out)
  expect_identical(back$ofv, fit$ofv)
  expect_identical(back$theta, fit$theta)
  expect_identical(back$omega, fit$omega)
  expect_identical(back$ebe_etas, fit$ebe_etas)
  expect_identical(back$eta_names, fit$eta_names)
  expect_identical(nrow(back$sdtab), nrow(fit$sdtab))
})
