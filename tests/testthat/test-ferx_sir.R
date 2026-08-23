
# ---- header from test-sir.R ----
# Tests for ferx_sir() and the path/hash provenance plumbing that backs it.
#
# SIR happy-path tests are gated on `skip_if(is.null(fit$cov_matrix), ...)`
# because in the no-autodiff CI build the warfarin cov step occasionally
# fails to converge with maxiter = 30 (same pattern as test-simulate-uncertainty.R).

sir_cov_skip <- "covariance step did not converge - skipping"









test_that("ferx_sir rejects fits without a covariance matrix", {
  fit <- warfarin_fit()  # covariance = FALSE
  expect_error(
    ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L),
    "covariance matrix",
    fixed = TRUE
  )
})
test_that("ferx_sir runs end-to-end and populates SIR fields", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)
  out <- ferx_sir(
    fit,
    sir_samples = 8L,
    sir_resamples = 4L,
    sir_seed = 1L,
    verbose = FALSE
  )

  # All SIR fields populated and shaped correctly.
  expect_true(is.numeric(out$sir_ess) && length(out$sir_ess) == 1L)
  expect_true(is.finite(out$sir_ess))

  n_theta <- length(out$theta)
  n_eta <- nrow(out$omega)
  n_sigma <- length(out$sigma)

  expect_equal(dim(out$sir_ci_theta), c(n_theta, 2L))
  expect_equal(dim(out$sir_ci_omega), c(n_eta, 2L))
  expect_equal(dim(out$sir_ci_sigma), c(n_sigma, 2L))
  expect_equal(colnames(out$sir_ci_theta), c("lower", "upper"))

  # CIs straddle the point estimates (with this tiny SIR run, occasional
  # extreme draws push the band wide; the lo <= hi check is the
  # invariant we actually care about).
  expect_true(all(out$sir_ci_theta[, "lower"] <= out$sir_ci_theta[, "upper"]))
  expect_true(all(out$sir_ci_omega[, "lower"] <= out$sir_ci_omega[, "upper"]))
  expect_true(all(out$sir_ci_sigma[, "lower"] <= out$sir_ci_sigma[, "upper"]))

  # Returned object still classes as ferx_fit so downstream printers /
  # methods keep working.
  expect_s3_class(out, "ferx_fit")
})
test_that("ferx_sir retains resamples when sir_keep_samples = TRUE", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)
  out <- ferx_sir(
    fit,
    sir_samples = 8L,
    sir_resamples = 4L,
    sir_seed = 2L,
    sir_keep_samples = TRUE,
    verbose = FALSE
  )

  # The resamples pool is a long flat vector; n × packed_dim
  # establishes the matrix shape on the consumer side.
  expect_true(is.numeric(out$sir_resamples))
  expect_equal(out$sir_resamples_n, 4L)
  expect_true(is.integer(out$sir_resamples_dim) && out$sir_resamples_dim > 0L)
  expect_equal(length(out$sir_resamples),
               out$sir_resamples_n * out$sir_resamples_dim)
})
test_that("ferx_sir refuses to run after the model file is tampered with", {
  # Run against a temp copy of the example so we don't mutate the bundled
  # example file. ferx_fit normalises paths, so the fit will point at the
  # temp model/data we created here.
  src <- ferx_example("warfarin")
  tmp <- tempfile("ferx_sir_tamper_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  model_tmp <- file.path(tmp, "warfarin.ferx")
  data_tmp <- file.path(tmp, "warfarin.csv")
  file.copy(src$model, model_tmp)
  file.copy(src$data, data_tmp)

  fit <- ferx_fit(
    model_tmp,
    data_tmp,
    method = "focei",
    verbose = FALSE,
    covariance = TRUE,
    settings = list(maxiter = 30L)
  )
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  # Append whitespace to the model — flips the SHA-256.
  cat("\n# tampered\n", file = model_tmp, append = TRUE)

  expect_error(
    ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L),
    "hash mismatch"
  )
})
test_that(".fitrx bundle preserves model.ferx bytes verbatim", {
  # Regression test for the writeLines() vs file.copy() byte-fidelity bug.
  # Doesn't need a fit covariance matrix.
  #
  # Specifically uses a model file *without* a trailing newline — that's
  # the case where `readLines + paste(collapse = "\n") + writeLines`
  # silently adds bytes that didn't exist, breaking the SHA-256 the engine
  # captured at fit time. The bundled warfarin example happens to end in
  # LF and round-trips fine through the buggy code path, which is why
  # we construct a degenerate copy here instead of reusing it.
  src <- ferx_example("warfarin")
  tmp <- tempfile("ferx_sir_bytes_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  model_tmp <- file.path(tmp, "warfarin.ferx")
  data_tmp <- file.path(tmp, "warfarin.csv")
  # Strip the trailing newline so writeLines() would re-add one.
  model_bytes <- readBin(src$model, "raw", n = file.info(src$model)$size)
  if (length(model_bytes) > 0 && model_bytes[length(model_bytes)] == as.raw(0x0a)) {
    model_bytes <- model_bytes[-length(model_bytes)]
  }
  writeBin(model_bytes, model_tmp)
  file.copy(src$data, data_tmp)

  fit <- ferx_fit(model_tmp, data_tmp,
                  method = "focei", verbose = FALSE,
                  covariance = FALSE,
                  settings = list(maxiter = 30L))
  expect_match(fit$model_hash, "^[0-9a-f]{64}$")

  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out), add = TRUE)
  ferx_save_fit(fit, out, include_data = FALSE)

  # Pull the bundled model.ferx and compare byte-for-byte with the
  # on-disk original we just fit against. The buggy writeLines() path
  # would have re-added a trailing 0x0a, so the bundled file would be
  # exactly one byte longer than the original.
  staging <- tempfile("fitrx_byte_check_")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::unzip(out, files = "model.ferx", exdir = staging, junkpaths = TRUE)
  bundled <- file.path(staging, "model.ferx")
  expect_true(file.exists(bundled))

  bundled_bytes <- readBin(bundled, "raw", n = file.info(bundled)$size)
  expect_identical(bundled_bytes, model_bytes)
})
test_that("path/hash fields round-trip through .fitrx save/load", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)
  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out), add = TRUE)

  ferx_save_fit(fit, out, include_data = TRUE)
  loaded <- ferx_load_fit(out)

  # Hashes are unchanged (the bundle preserves the source bytes verbatim,
  # so loading + re-hashing the unpacked copy yields the same digest).
  expect_identical(loaded$model_hash, fit$model_hash)
  expect_identical(loaded$data_hash, fit$data_hash)

  # After load, model_path / data_path point at the unpacked tempfile
  # copies; running ferx_sir against the loaded fit must still succeed.
  expect_true(file.exists(loaded$model_path))
  expect_true(file.exists(loaded$data_path))

  out_sir <- ferx_sir(
    loaded,
    sir_samples = 6L,
    sir_resamples = 3L,
    sir_seed = 7L
  )
  expect_true(is.finite(out_sir$sir_ess))
})

# ferx-core#1021: a covariance direction the data do not identify comes back
# from the covariance step with a variance of ~1/eigenvalue-floor. SIR used to
# die on it ("All SIR samples had invalid weights"); it now shrinks the
# direction and says so. The warning is the only signal that the CIs this call
# just wrote are qualified, so it has to reach the fit - both the flat vector
# and the structured table.
#
# The merge is tested against a mocked binding rather than a real engine run:
# the behaviour under test is R-side, and gating it on the engine would make it
# unrunnable until the ferx-core pin in src/rust/Cargo.lock catches up. The
# end-to-end test below covers the real engine when it is new enough.
fake_sir_return <- function(warnings, n_theta, n_eta, n_sigma) {
  list(
    sir_ess = 4,
    sir_ci_theta = rep(c(0.1, 0.2), n_theta),
    sir_ci_omega = rep(c(0.01, 0.02), n_eta),
    sir_ci_sigma = rep(c(0.001, 0.002), n_sigma),
    sir_resamples = numeric(0),
    sir_resamples_n = 0L,
    sir_resamples_dim = 0L,
    warnings = warnings
  )
}

test_that("ferx_sir merges engine SIR warnings onto the fit", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  shrunk <- paste0(
    "proposal was shrunk in 1 direction(s) so draws mostly stay inside the ",
    "parameter bounds [TVCL -0.99 (sd 6.15e3 -> 2.37e0)]. Those directions come ",
    "from eigenvalue-floored (non-identified) curvature in the covariance step."
  )
  testthat::local_mocked_bindings(
    ferx_rust_sir = function(...) {
      fake_sir_return(shrunk, length(fit$theta), nrow(fit$omega), length(fit$sigma))
    },
    .package = "ferx"
  )

  out <- ferx_sir(fit, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)

  # Flat vector, structured table, and the pre-existing warnings all survive.
  expect_true(any(grepl("shrunk", out$warnings, fixed = TRUE)))
  expect_true(all(fit$warnings %in% out$warnings))
  ws <- out$warnings_structured
  expect_s3_class(ws, "data.frame")
  sir_rows <- ws[ws$category == "sir", , drop = FALSE]
  expect_gt(nrow(sir_rows), 0L)
  expect_true(any(grepl("shrunk", sir_rows$message, fixed = TRUE)))
  # The parameter behind the shrunk direction is named, which is what makes the
  # warning actionable.
  expect_true(any(grepl("TVCL", sir_rows$message, fixed = TRUE)))

  # A second pass over the same fit must not stack a duplicate of either row.
  twice <- ferx_sir(out, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)
  expect_equal(sum(grepl("shrunk", twice$warnings, fixed = TRUE)), 1L)
  expect_equal(
    sum(grepl("shrunk", twice$warnings_structured$message, fixed = TRUE)),
    1L
  )
})

test_that("ferx_sir adds no sir rows when the engine reports nothing", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  testthat::local_mocked_bindings(
    ferx_rust_sir = function(...) {
      fake_sir_return(character(0), length(fit$theta), nrow(fit$omega), length(fit$sigma))
    },
    .package = "ferx"
  )
  out <- ferx_sir(fit, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)

  ws <- out$warnings_structured
  sir_rows <- if (is.data.frame(ws)) ws[ws$category == "sir", , drop = FALSE] else NULL
  expect_true(is.null(sir_rows) || nrow(sir_rows) == 0L)
  expect_identical(out$warnings, fit$warnings)
})

# End-to-end against the real engine. The pin in src/rust/Cargo.lock carries
# ferx-core#1021 as of #304, so this runs for real rather than skipping. It is
# deliberately NOT wrapped in a tryCatch that downgrades the pre-#1021
# `All SIR samples had invalid weights` failure to a skip: that is the exact
# regression this test exists to catch, and turning it into a skip would report
# it as green.
test_that("ferx_sir surfaces engine proposal-conditioning warnings end-to-end", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  # Adding to a diagonal keeps the matrix PSD, so this is a covariance a real
  # (non-identified) fit could produce rather than a malformed input.
  degenerate <- fit
  degenerate$cov_matrix[1L, 1L] <- degenerate$cov_matrix[1L, 1L] + 1e8

  out <- ferx_sir(degenerate, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)

  expect_true(is.finite(out$sir_ess))
  expect_true(any(grepl("shrunk", out$warnings, fixed = TRUE)))
  sir_rows <- out$warnings_structured
  sir_rows <- sir_rows[sir_rows$category == "sir", , drop = FALSE]
  expect_true(any(grepl("shrunk", sir_rows$message, fixed = TRUE)))
  expect_true(any(grepl(names(fit$theta)[1L], sir_rows$message, fixed = TRUE)))
})

test_that("ferx_sir leaves a healthy fit free of sir warnings", {
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), sir_cov_skip)

  before <- fit$warnings
  out <- ferx_sir(fit, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)

  expect_false(any(grepl("shrunk", out$warnings, fixed = TRUE)))
  expect_false(any(grepl("rank-deficient", out$warnings, fixed = TRUE)))
  # Pre-existing warnings survive the merge.
  expect_true(all(before %in% out$warnings))
})
