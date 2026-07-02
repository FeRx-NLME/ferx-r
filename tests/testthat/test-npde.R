# Tests for ferx_calc_npde(): post-hoc simulation-based NPDE/NPD from a fit.
#
# The argument-guard tests are pure R and run without the compiled engine: each
# error is raised before the FFI call. The end-to-end tests need the real
# backend and the bundled warfarin example (shared `warfarin_fit()` helper).

# --- argument validation (no engine needed) ----------------------------------

test_that("ferx_calc_npde requires a fit with theta/omega/sigma", {
  expect_error(ferx_calc_npde(list()), "theta, omega, and sigma")
})

# A minimal structurally-valid fit that passes validate_fit_for_params and the
# sdtab gate but carries no model/data path, so the nsim/seed/file guards fire
# before any FFI call.
.stub_fit <- function() {
  list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2), sigma = 0.05,
       sdtab = data.frame(ID = 1, TIME = 1, DV = 5))
}

test_that("ferx_calc_npde rejects non-positive or non-scalar nsim", {
  expect_error(ferx_calc_npde(.stub_fit(), nsim = 0L),  "positive integer")
  expect_error(ferx_calc_npde(.stub_fit(), nsim = -5L), "positive integer")
  expect_error(ferx_calc_npde(.stub_fit(), nsim = c(10L, 20L)), "positive integer")
})

test_that("ferx_calc_npde rejects a non-integer / NA seed", {
  expect_error(ferx_calc_npde(.stub_fit(), seed = c(1L, 2L)), "single integer or NULL")
  expect_error(ferx_calc_npde(.stub_fit(), seed = NA_integer_), "single integer or NULL")
})

test_that("ferx_calc_npde rejects a negative seed (would silently become the default)", {
  expect_error(ferx_calc_npde(.stub_fit(), seed = -5L), "non-negative integer")
})

test_that("ferx_calc_npde errors when no model/data path is available", {
  expect_error(ferx_calc_npde(.stub_fit()), "model file")
  expect_error(ferx_calc_npde(.stub_fit(), model = "no_such.ferx"), "model file")
})

# --- merge helper (pure R) ---------------------------------------------------

.attach <- getFromNamespace(".ferx_attach_npde", "ferx")

test_that(".ferx_attach_npde copies NPDE/NPD positionally and keeps existing cols", {
  sd <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1), DV = c(5, 6, 7))
  np <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1),
                   NPDE = c(0.1, 0.2, 0.3), NPD = c(0.4, 0.5, 0.6))
  out <- .attach(sd, np)
  expect_equal(out$NPDE, c(0.1, 0.2, 0.3))
  expect_equal(out$NPD,  c(0.4, 0.5, 0.6))
  expect_true(all(c("DV", "NPDE", "NPD") %in% names(out)))
})

test_that(".ferx_attach_npde errors (not silently joins) when row counts differ", {
  sd <- data.frame(ID = c(1, 9), TIME = c(1, 9))
  np <- data.frame(ID = 1, TIME = 1, NPDE = 0.1, NPD = 0.2)
  expect_error(.attach(sd, np), "cannot align")
})

test_that(".ferx_attach_npde errors when rows do not line up by ID/TIME", {
  sd <- data.frame(ID = c(2, 1), TIME = c(1, 2))
  np <- data.frame(ID = c(1, 2), TIME = c(2, 1), NPDE = c(0.9, 0.1), NPD = c(0.8, 0.2))
  expect_error(.attach(sd, np), "do not line up")
})

# --- end-to-end against the real engine --------------------------------------

test_that("ferx_calc_npde returns the fit with NPDE/NPD in sdtab, one row per obs", {
  ex  <- ferx_example("warfarin")
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 1L,
                   model = ex$model, data = ex$data)
  expect_s3_class(fit, "ferx_fit")
  expect_true(all(c("NPDE", "NPD") %in% names(fit$sdtab)))

  dat   <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(fit$sdtab), n_obs)
})

test_that("ferx_calc_npde NPD is finite and roughly mean-zero on a sensible fit", {
  ex  <- ferx_example("warfarin")
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 500L, seed = 1L,
                   model = ex$model, data = ex$data)
  expect_true(all(is.finite(fit$sdtab$NPD)))
  expect_lt(abs(mean(fit$sdtab$NPD)), 0.3)
})

test_that("ferx_calc_npde is reproducible for a fixed seed", {
  ex <- ferx_example("warfarin")
  a  <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  b  <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  expect_equal(a$sdtab$NPDE, b$sdtab$NPDE)
  expect_equal(a$sdtab$NPD,  b$sdtab$NPD)
})

test_that("ferx_calc_npde uses fit$model_path / fit$data_path by default", {
  # The cached fit records the paths it was run on, so no explicit model/data
  # is needed.
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 100L, seed = 1L)
  expect_true(all(c("NPDE", "NPD") %in% names(fit$sdtab)))
})

test_that("ferx_calc_npde errors when the fit carries no sdtab", {
  bad <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2),
              sigma = 0.05)  # no sdtab
  expect_error(ferx_calc_npde(bad), "sdtab` is empty")
})

test_that("ferx_calc_npde surfaces a clean error when the engine returns NULL", {
  # The FFI prints its own message and returns NULL on failure; ferx_calc_npde must
  # turn that into a clear R error rather than crash in the alignment step.
  fit <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2),
              sigma = 0.05, sdtab = data.frame(ID = 1, TIME = 1, DV = 5),
              model_path = tempfile(fileext = ".ferx"),
              data_path  = tempfile(fileext = ".csv"))
  file.create(fit$model_path, fit$data_path)
  testthat::local_mocked_bindings(ferx_rust_npde_from_fit = function(...) NULL)
  expect_error(ferx_calc_npde(fit), "engine returned no NPDE table")
})
