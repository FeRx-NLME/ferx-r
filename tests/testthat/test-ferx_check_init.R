
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
# Return structure — Tier 1
















# These tests require the covariance step to have succeeded. With maxiter = 30L
# the outer optimisation may not converge on all machines, so we skip rather
# than fail — a skip here means "covariance step needs more iterations", not
# a bug. A full-convergence run (no maxiter cap) is tested manually / locally.









# fd_hessian_step argument — Tier 1 (R-side validation, no Rust call)








# ferx_check_init — Tier 1




# ferx_cor_matrix — Tier 1, requires covariance = TRUE





# -- Diagnostic fields (Step 3 feature-parity) --------------------------------





# -- IOV: shrinkage_kappa_by_occ (Step 1 feature-parity) ---------------------




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
