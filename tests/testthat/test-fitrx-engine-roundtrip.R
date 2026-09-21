# A `.fitrx` written from R has to load in ferx-core, not only in R (#379).
#
# `ferx_save_fit()` wrote the bundle with `jsonlite::write_json(auto_unbox =
# TRUE)`, which collapses every length-1 vector to a JSON scalar. ferx-core's
# reader declares those fields as sequences, so a one-sigma / one-warning /
# unchained fit - which is most fits - produced a bundle the engine refused
# with `invalid type: string "foce", expected a sequence`. The existing
# `ferx_save_fit()` -> `ferx_load_fit()` tests all passed throughout, because
# R's own reader is happy either way; that is precisely why this survived, so
# the tests below go through the engine instead.

# The bundle under test: one sigma, one (unchained) method, and a covariance
# step, written from a real fit.
fitrx_r_bundle <- function(env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".fitrx", .local_envir = env)
  ferx_save_fit(warfarin_fit_cov(), path)
  path
}

# `//` starts a comment in a `.ferx` file, and `tempdir()` on macOS is
# `/var/folders/.../T//Rtmp...`, so a path pasted into a block is truncated
# before the engine sees it and the error names the truncated path. Collapse
# the run of slashes before writing one into a model file.
fitrx_block_path <- function(path) gsub("/+", "/", normalizePath(path))

# -- The shape itself, no engine required -------------------------------------

test_that("fit.json writes the schema's array fields as JSON arrays", {
  bundle <- fitrx_r_bundle()
  staging <- withr::local_tempdir()
  utils::unzip(bundle, exdir = staging)
  wire <- jsonlite::read_json(file.path(staging, "fit.json"),
                              simplifyVector = FALSE)

  # The field the engine trips on first, and the ones the issue's patch loop
  # reached after it.
  expect_type(wire$method_chain, "list")
  expect_length(wire$method_chain, 1L) # unchained: this is the length-1 case
  expect_type(wire$warnings, "list")
  for (key in c("names", "estimates", "se", "fixed", "types", "init_as_sd")) {
    expect_type(wire$sigma[[key]], "list")
  }
  expect_length(wire$sigma$names, 1L) # the model has exactly one sigma
  for (key in c("names", "estimates", "se", "fixed", "transform")) {
    expect_type(wire$theta[[key]], "list")
  }
  expect_type(wire$omega$matrix$data, "list")
  expect_type(wire$covariance_matrix$data, "list")
  expect_type(wire$cov_eigenvalues, "list")

  # Scalars stay scalars: `auto_unbox` is still on, and a wire scalar written
  # as a one-element array would fail the engine the other way round.
  expect_type(wire$method, "character")
  expect_type(wire$converged, "logical")
  expect_type(wire$ofv, "double")
  expect_type(wire$omega$matrix$rows, "integer")

  manifest <- jsonlite::read_json(file.path(staging, "manifest.json"),
                                  simplifyVector = FALSE)
  expect_type(manifest$entries, "list")
})

test_that("a one-kappa IOV fit writes the nested arrays as arrays", {
  # `shrinkage_kappa_by_occ` is `Vec<Vec<f64>>`: with a single kappa each
  # occasion's row is a length-1 vector, so the *inner* sequence is the one
  # `auto_unbox` collapses. The same shape covers the SIR CI pairs.
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
                  covariance = FALSE, settings = list(maxiter = 5L))
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  staging <- withr::local_tempdir()
  utils::unzip(path, exdir = staging)
  wire <- jsonlite::read_json(file.path(staging, "fit.json"),
                              simplifyVector = FALSE)

  expect_type(wire$iov$kappa_names, "list")
  expect_length(wire$iov$kappa_names, 1L)
  expect_type(wire$iov$kappa_fixed, "list")
  expect_type(wire$iov$shrinkage_kappa, "list")
  expect_type(wire$iov$shrinkage_kappa_by_occ, "list")
  expect_gt(length(wire$iov$shrinkage_kappa_by_occ), 0L)
  for (row in wire$iov$shrinkage_kappa_by_occ) {
    expect_type(row, "list")
    expect_length(row, 1L)
  }
  expect_type(wire$iov$omega_iov$data, "list")

  ferx_bin <- Sys.which("ferx")
  if (nzchar(ferx_bin)) {
    out <- suppressWarnings(system2(ferx_bin, c("summary", path),
                                    stdout = TRUE, stderr = TRUE))
    expect_identical(as.integer(attr(out, "status") %||% 0L), 0L,
                     info = paste(out, collapse = "\n"))
  }
})

test_that("the R -> R round trip still reads the bundle", {
  # The direction that was never broken, kept green next to the one that was.
  bundle <- fitrx_r_bundle()
  fit <- warfarin_fit_cov()
  loaded <- ferx_load_fit(bundle)
  expect_s3_class(loaded, "ferx_fit")
  expect_equal(loaded$theta, fit$theta, tolerance = 1e-10)
  # Unnamed on both sides: `ferx_load_fit()` names the sigma vector from
  # `sigma$names` and a fresh fit does not, which is older drift than this
  # change and not what this test is about.
  expect_equal(unname(loaded$sigma), unname(fit$sigma), tolerance = 1e-10)
  expect_identical(names(loaded$sigma), as.character(fit$sigma_names))
  expect_equal(loaded$ofv, fit$ofv, tolerance = 1e-10)
  expect_identical(as.character(loaded$method_chain %||% loaded$method),
                   as.character(fit$method_chain %||% fit$method))
  expect_identical(loaded$warnings, fit$warnings)
})

# -- Through the engine: [priors] from_fit ------------------------------------
#
# `from_fit` is read by the engine's own `.fitrx` loader while the model file
# is parsed, so this fails on an unreadable bundle without needing the CLI -
# it runs wherever the package's tests run.

test_that("[priors] from_fit reads a bundle written by ferx_save_fit()", {
  ex     <- ferx_example("warfarin")
  bundle <- fitrx_r_bundle()
  model  <- withr::local_tempfile(fileext = ".ferx")
  writeLines(
    c(readLines(ex$model), "",
      "[priors]",
      paste0("  from_fit = ", fitrx_block_path(bundle))),
    model
  )

  res <- ferx_model_validate(model, ex$data)
  expect_true(
    isTRUE(res$ok),
    info = paste(utils::capture.output(print(res$diagnostics)), collapse = "\n")
  )
  # The failure this test exists for arrives as a parse error naming the file.
  msgs <- paste(as.character(res$diagnostics$message %||% character()),
                collapse = " ")
  expect_false(grepl("from_fit", msgs, fixed = TRUE), info = msgs)
  expect_false(grepl("expected a sequence", msgs, fixed = TRUE), info = msgs)

  # ... and the priors are actually applied: a priored fit reports the split
  # objective, which a fit without a `[priors]` block does not.
  fit <- ferx_fit(model, ex$data, method = "foce", verbose = FALSE,
                  covariance = FALSE, settings = list(maxiter = 3L))
  expect_true(is.finite(fit$ofv))
  expect_false(is.null(fit$prior_summary))
  expect_gt(nrow(as.data.frame(fit$prior_summary)), 0L)
})

# -- Through the engine: the `ferx` CLI ---------------------------------------
#
# The reproducer from #379 verbatim. Skipped where no binary is on PATH, which
# is every CI runner - the `from_fit` test above is the one that guards this
# there.

test_that("`ferx summary` loads a bundle written by ferx_save_fit()", {
  ferx_bin <- Sys.which("ferx")
  skip_if(!nzchar(ferx_bin), "no `ferx` binary on PATH")

  bundle <- fitrx_r_bundle()
  # No shell here: `system2()` passes the vector straight to execve, so the
  # path needs no quoting (and quoting it would be part of the file name).
  out <- suppressWarnings(system2(ferx_bin, c("summary", bundle),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status") %||% 0L
  expect_identical(
    as.integer(status), 0L,
    info = paste(out, collapse = "\n")
  )
  expect_false(any(grepl("expected a sequence", out, fixed = TRUE)),
               info = paste(out, collapse = "\n"))
})
