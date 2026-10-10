# SIR draws kept with `sir_keep_samples = TRUE` survive a .fitrx round trip
# (#549), and the sir_ci_* intervals keep their row names (#482).
#
# Before #549 `.fitrx_build_sir_wire()` wrote `resamples_packed = NULL`, so a
# reloaded fit could not drive `ferx_simulate_with_uncertainty(method =
# "sir")`, and the loader rebuilt sir_ci_theta / _omega / _sigma without
# dimnames. The oracle for the wire shape is ferx-core's own reader and writer
# (`fitrx_engine_trip()`), not R's loader agreeing with R's writer.

sir_res_cov_skip <- "covariance step did not converge - skipping"

# The SIR settings every test here uses; `sir_seed` fixes the draws.
sir_res_settings <- list(sir_samples = 400L, sir_resamples = 200L,
                         sir_seed = 1L, sir_keep_samples = TRUE)

# Standalone SIR on the cached covariance fit, draws kept.
sir_res_standalone <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      base <- warfarin_fit_cov()
      if (is.null(base$cov_matrix)) return(NULL)
      fit <<- ferx_sir(base, sir_samples = 400L, sir_resamples = 200L,
                       sir_seed = 1L, sir_keep_samples = TRUE,
                       verbose = FALSE)
    }
    fit
  }
})

# In-fit SIR, draws kept: the call validate_fit_for_uncertainty() recommends.
sir_res_infit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                       covariance = TRUE, sir = TRUE,
                       settings = c(list(maxiter = 30L), sir_res_settings))
    }
    fit
  }
})

sir_res_fields <- c("sir_resamples", "sir_resamples_n", "sir_resamples_dim",
                    "sir_ci_theta", "sir_ci_omega", "sir_ci_sigma",
                    "sir_ci_kappa")

sir_res_roundtrip <- function(fit) {
  path <- withr::local_tempfile(fileext = ".fitrx", .local_envir = parent.frame())
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

expect_sir_fields_identical <- function(loaded, fit) {
  for (f in sir_res_fields) {
    expect_identical(loaded[[f]], fit[[f]], info = f)
  }
}

sir_res_simulate <- function(fit) {
  ex <- ferx_example("warfarin")
  ferx_simulate_with_uncertainty(ex$model, ex$data, fit,
                                 n_uncertainty_draws = 5L, n_sim_per_draw = 1L,
                                 method = "sir", seed = 1L)
}

test_that("ferx_sir() draws and intervals survive save -> load (#549, #482)", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  # The fixture really kept draws, and its intervals carry names - otherwise
  # identical() below would pass on two NULLs.
  expect_identical(fit$sir_resamples_n, 200L)
  expect_gt(fit$sir_resamples_dim, 0L)
  expect_identical(rownames(fit$sir_ci_theta), names(fit$theta))
  expect_identical(rownames(fit$sir_ci_omega), fit$eta_names)
  expect_identical(rownames(fit$sir_ci_sigma), fit$sigma_names)

  loaded <- sir_res_roundtrip(fit)
  expect_sir_fields_identical(loaded, fit)
  expect_identical(sir_res_simulate(loaded), sir_res_simulate(fit))
})

test_that("in-fit SIR draws and intervals survive save -> load (#549, #482)", {
  fit <- sir_res_infit()
  # Skip only when the covariance step failed, which SIR needs. A fit that ran
  # but kept no draws is the regression this test is for, so it must fail.
  skip_if(is.null(fit$cov_matrix), sir_res_cov_skip)
  expect_identical(fit$sir_resamples_n, 200L)
  expect_identical(rownames(fit$sir_ci_theta), names(fit$theta))

  loaded <- sir_res_roundtrip(fit)
  expect_sir_fields_identical(loaded, fit)
  expect_identical(sir_res_simulate(loaded), sir_res_simulate(fit))
})

test_that("fit.json carries the draws in ferx-core's resamples_packed shape", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  trip <- fitrx_engine_trip(fit)
  r_wire <- fitrx_fit_json(trip$r_path)$sir$resamples_packed
  core_wire <- fitrx_fit_json(trip$core_path)$sir$resamples_packed

  # One array per resample, each one row of the flat row-major vector.
  expect_length(r_wire, fit$sir_resamples_n)
  expect_identical(lengths(r_wire), rep(fit$sir_resamples_dim, fit$sir_resamples_n))
  m <- matrix(fit$sir_resamples, nrow = fit$sir_resamples_n, byrow = TRUE)
  expect_identical(as.numeric(unlist(r_wire[[2L]])), m[2L, ])
  # ferx-core read the R bundle's draws and wrote back the same values.
  expect_identical(core_wire, r_wire)

  # And R reads the engine-written bundle exactly as it reads its own.
  expect_sir_fields_identical(ferx_load_fit(trip$core_path), fit)
})

test_that("a bundle without resamples_packed still loads, draws NULL", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path, include_data = TRUE)

  # What every R writer before #549 wrote (`null`), and a bundle with no key.
  as_null <- fitrx_edit_bundle(path, function(w) {
    w$sir["resamples_packed"] <- list(NULL)
    w
  })
  no_key <- fitrx_edit_bundle(path, function(w) {
    w$sir$resamples_packed <- NULL
    w
  })
  for (p in c(as_null, no_key)) {
    loaded <- ferx_load_fit(p)
    expect_null(loaded$sir_resamples)
    expect_null(loaded$sir_resamples_n)
    expect_null(loaded$sir_resamples_dim)
    expect_identical(loaded$sir_ci_theta, fit$sir_ci_theta)
    expect_error(sir_res_simulate(loaded), "sir_resamples` is empty")
  }
})

test_that("ferx_save_fit() refuses draws that do not fill n x dim", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  path <- withr::local_tempfile(fileext = ".fitrx")
  short <- fit
  short$sir_resamples <- short$sir_resamples[-1L]
  expect_error(ferx_save_fit(short, path),
               paste0("holds ", length(short$sir_resamples), " values.*edited after the fit"))
  # Draws without their shape are the same edit, not a silent drop (round 1 #3).
  for (field in c("sir_resamples_n", "sir_resamples_dim")) {
    no_shape <- fit
    no_shape[field] <- list(NULL)
    expect_error(ferx_save_fit(no_shape, path), "edited after the fit", info = field)
  }
})

test_that("ferx_save_fit() refuses non-finite draws (round 1 #1)", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  path <- withr::local_tempfile(fileext = ".fitrx")
  for (bad in c(NaN, NA, Inf, -Inf)) {
    f <- fit
    f$sir_resamples[3L] <- bad
    expect_error(ferx_save_fit(f, path), "1 non-finite value", info = format(bad))
  }
})

test_that("ferx_load_fit() refuses resamples_packed rows that are not n x dim numbers", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path, include_data = TRUE)
  corrupt <- function(edit) {
    p <- fitrx_edit_bundle(path, function(w) {
      w$sir$resamples_packed <- edit(w$sir$resamples_packed)
      w
    })
    expect_error(ferx_load_fit(p), "the bundle is corrupt")
  }
  # A shorter row.
  corrupt(function(rp) { rp[[2L]] <- rp[[2L]][-1L]; rp })
  # Every row empty (round 1 #7).
  corrupt(function(rp) lapply(rp, function(r) list()))
  # A `null` cell, as a NaN would be written: same length, one value fewer
  # after unlist() (round 1 #1).
  corrupt(function(rp) { rp[[1L]][3L] <- list(NULL); rp })
  # A text cell.
  corrupt(function(rp) { rp[[1L]][[3L]] <- "x"; rp })
})

test_that("ferx_sir() without sir_keep_samples drops the input fit's draws", {
  fit <- sir_res_standalone()
  skip_if(is.null(fit), sir_res_cov_skip)
  rerun <- ferx_sir(fit, sir_samples = 400L, sir_resamples = 200L,
                    sir_seed = 2L, sir_keep_samples = FALSE, verbose = FALSE)
  expect_false(isTRUE(rerun$sir_settings$keep_samples))
  expect_null(rerun$sir_resamples)
  expect_null(rerun$sir_resamples_n)
  expect_null(rerun$sir_resamples_dim)
})
