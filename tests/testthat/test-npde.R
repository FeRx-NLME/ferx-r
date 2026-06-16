# Tests for ferx_npde(): post-hoc simulation-based NPDE/NPD from a fit.
#
# The argument-guard tests are pure R and run without the compiled engine: each
# error is raised before the FFI call. The end-to-end tests need the real
# backend and the bundled warfarin example (shared `warfarin_fit()` helper).

# --- argument validation (no engine needed) ----------------------------------

test_that("ferx_npde requires a fit with theta/omega/sigma", {
  expect_error(ferx_npde(list()), "theta, omega, and sigma")
})

# A minimal structurally-valid fit that passes validate_fit_for_params but
# carries no model/data path, so the file guards fire before any FFI call.
.stub_fit <- function() {
  list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2), sigma = 0.05)
}

test_that("ferx_npde rejects non-positive or non-scalar nsim", {
  expect_error(ferx_npde(.stub_fit(), nsim = 0L),  "positive integer")
  expect_error(ferx_npde(.stub_fit(), nsim = -5L), "positive integer")
  expect_error(ferx_npde(.stub_fit(), nsim = c(10L, 20L)), "positive integer")
})

test_that("ferx_npde rejects a non-integer / NA seed", {
  expect_error(ferx_npde(.stub_fit(), seed = c(1L, 2L)), "single integer or NULL")
  expect_error(ferx_npde(.stub_fit(), seed = NA_integer_), "single integer or NULL")
})

test_that("ferx_npde errors when no model/data path is available", {
  expect_error(ferx_npde(.stub_fit()), "model file")
  expect_error(ferx_npde(.stub_fit(), model = "no_such.ferx"), "model file")
})

# --- end-to-end against the real engine --------------------------------------

test_that("ferx_npde returns ID/TIME/NPDE/NPD, one row per observation", {
  ex   <- ferx_example("warfarin")
  npde <- ferx_npde(warfarin_fit(), nsim = 200L, seed = 1L,
                    model = ex$model, data = ex$data)
  expect_s3_class(npde, "data.frame")
  expect_true(all(c("ID", "TIME", "NPDE", "NPD") %in% names(npde)))

  dat   <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(npde), n_obs)
})

test_that("ferx_npde NPD is finite and roughly mean-zero on a sensible fit", {
  ex   <- ferx_example("warfarin")
  npde <- ferx_npde(warfarin_fit(), nsim = 500L, seed = 1L,
                    model = ex$model, data = ex$data)
  expect_true(all(is.finite(npde$NPD)))
  expect_lt(abs(mean(npde$NPD)), 0.3)
})

test_that("ferx_npde is reproducible for a fixed seed", {
  ex <- ferx_example("warfarin")
  a  <- ferx_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  b  <- ferx_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  expect_equal(a$NPDE, b$NPDE)
  expect_equal(a$NPD,  b$NPD)
})

test_that("ferx_npde uses fit$model_path / fit$data_path by default", {
  # The cached fit records the paths it was run on, so no explicit model/data
  # is needed.
  npde <- ferx_npde(warfarin_fit(), nsim = 100L, seed = 1L)
  expect_true(all(c("ID", "TIME", "NPDE", "NPD") %in% names(npde)))
})
