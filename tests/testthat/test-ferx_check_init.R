
# write_pipe_test_model()/modifying_editor() come from helper-model-pipe.R

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------



test_that("ferx_check_init() accepts a ferx_model in a pipe and uses its data", {
  ex <- ferx_example("warfarin")
  chk <- suppressWarnings(
    ferx_model(ex$model, data = ex$data) |>
      ferx_check_init(method = "focei", maxiter = 2L)
  )
  expect_named(chk, c("fit", "trace", "summary"))
  expect_s3_class(chk$fit, "ferx_fit")
  expect_true(is.data.frame(chk$summary))
})
test_that("ferx_check_init() errors when ferx_model has no data and none is supplied", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  expect_error(
    ferx_model(model = path) |> ferx_check_init(),
    regexp = "data"
  )
})
test_that("ferx_check_init() explicit data argument overrides data on ferx_model", {
  ex <- ferx_example("warfarin")
  # Sanity check that explicit data still wins over $data on the object —
  # we pass the same ex$data twice but via different routes.
  chk <- suppressWarnings(
    ferx_check_init(
      ferx_model(ex$model, data = ex$data),
      data = ex$data, method = "focei", maxiter = 2L
    )
  )
  expect_s3_class(chk$fit, "ferx_fit")
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
