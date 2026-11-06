
# ---- header from test-persist.R ----
# ferx_save_fit() / ferx_load_fit() round-trip and CLI-flag tests.
#
# These exercise the .fitrx bundle format. The cross-tool test (Rust CLI
# writes, R reads) lives in ferx-core; here we only verify the R-side
# writer + reader are inverses, and that the `output` argument of
# ferx_fit() saves a valid bundle.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.fitrx_build_iov_wire   <- getFromNamespace(".fitrx_build_iov_wire",   "ferx")





















test_that("bad path errors cleanly", {
  expect_error(ferx_load_fit("does-not-exist.fitrx"), "does not exist")
  bogus <- tempfile(fileext = ".fitrx")
  writeLines("not a zip", bogus)
  on.exit(unlink(bogus), add = TRUE)
  expect_error(ferx_load_fit(bogus))
})
# ---- header from test-persist-helpers.R ----
# Unit tests for the pure-R serialization helpers in persist.R. These convert
# between the in-memory ferx_fit fields and the cross-language .fitrx wire
# schema, and need no model fit to exercise — so they live in the fast tier.

# ---------------------------------------------------------------------------
# Method name <-> token round-trip
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Covariance status <-> token
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Error model inference from sigma types
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Matrix <-> wire (row-major)
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Confidence-interval <-> wire
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Optional scalar wrappers (NULL/NA/empty collapse to NULL)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Named SE unwrap + subject string IDs
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-subject OFV / N_OBS extraction from sdtab
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# EBEs CSV writer — header-only path when there are no random effects
# ---------------------------------------------------------------------------


test_that(".fitrx_num_or_nan reads a JSON null back as NaN (#508)", {
  expect_identical(ferx:::.fitrx_num_or_nan(NULL), NaN)
  expect_identical(ferx:::.fitrx_num_or_nan(list()), NaN)
  expect_identical(ferx:::.fitrx_num_or_nan(NA_real_), NA_real_)
  expect_identical(ferx:::.fitrx_num_or_nan(12.5), 12.5)
  expect_identical(ferx:::.fitrx_num_or_nan(list(-3)), -3)
})
test_that("a loaded fit with zero-padded IDs runs ferx_covariance like the live one (#468)", {
  fit <- relabelled_fit("warfarin", "padded", function(id) sprintf("%03d", id),
                        covariance = FALSE)
  expect_identical(fit$ebe_etas$ID[1:2], c("001", "002"))
  path <- tempfile(fileext = ".fitrx")
  on.exit(unlink(path), add = TRUE)
  ferx_save_fit(fit, path)
  loaded <- ferx_load_fit(path)
  expect_identical(loaded$ebe_etas$ID, fit$ebe_etas$ID)

  # The bundle carries the packed estimate (ferx-core #1815), so the loaded
  # fit repeats the live step exactly (#511).
  expect_identical(loaded$packed_estimate, fit$packed_estimate)
  live <- ferx_covariance(fit)
  skip_if(is.null(live$cov_matrix), "covariance step did not converge - skipping")
  expect_identical(ferx_covariance(loaded)$se_theta, live$se_theta)
})

# ---- text subject IDs through save / load (#475) ----
# Each key gives subject 1 one awkward label and leaves the others numeric.
# "numeric" relabels nothing, so it pins the unchanged case.
text_id_labels <- list(na = "NA", comma = "Smith, 2019", quote = "q\"x",
                       padded = "007", numeric = "1")

text_id_fit <- function(key) {
  label <- text_id_labels[[key]]
  relabelled_fit("warfarin", paste0("text_", key),
                 function(id) ifelse(id == 1L, label, as.character(id)))
}

test_that("a reloaded fit keeps text IDs and runs ferx_covariance and ferx_sir (#475)", {
  fits <- lapply(names(text_id_labels), text_id_fit)
  names(fits) <- names(text_id_labels)
  skip_if(any(vapply(fits, function(f) is.null(f$cov_matrix), logical(1))),
          "covariance step did not converge - skipping")
  for (key in names(text_id_labels)) {
    fit <- fits[[key]]
    expect_identical(fit$ebe_etas$ID[1], text_id_labels[[key]], info = key)
    path <- withr::local_tempfile(fileext = ".fitrx")
    ferx_save_fit(fit, path)
    loaded <- ferx_load_fit(path)
    expect_identical(loaded$ebe_etas$ID, fit$ebe_etas$ID, info = key)
    expect_identical(loaded$sdtab$ID, fit$sdtab$ID, info = key)

    # The bundle carries the packed estimate (ferx-core #1815), so the loaded
    # fit repeats the live step exactly (#511).
    expect_identical(ferx_covariance(loaded)$se_theta,
                     ferx_covariance(fit)$se_theta, info = key)
    sir <- function(f) ferx_sir(f, sir_samples = 50L, sir_resamples = 20L, sir_seed = 1L)
    expect_identical(sir(loaded)$sir_ci_theta, sir(fit)$sir_ci_theta, info = key)
  }
})

test_that("ferx-core's loader reads an R-written bundle with the same text IDs (#475)", {
  # The oracle is the engine's own load_fit / save_fit, not R reading back
  # what R wrote; its loader also checks predictions.csv against ebes.csv.
  for (key in names(text_id_labels)) {
    fit <- text_id_fit(key)
    trip <- fitrx_engine_trip(fit)
    from_core <- ferx_load_fit(trip$core_path)
    expect_identical(from_core$ebe_etas$ID, fit$ebe_etas$ID, info = key)
    expect_identical(from_core$sdtab$ID, fit$sdtab$ID, info = key)
  }
})

test_that("a bundle without reader settings loads with neither field and still runs (#462)", {
  # As written before #462: a fit that records neither leaves both keys out.
  fit <- rs_legacy(warfarin_fit_cov())
  skip_if(is.null(fit$cov_matrix), "covariance step did not converge - skipping")
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  staging <- withr::local_tempdir()
  utils::unzip(path, exdir = staging)
  wire <- jsonlite::read_json(file.path(staging, "fit.json"), simplifyVector = FALSE)
  expect_false(any(c("reader_settings", "population_fingerprint") %in% names(wire)))

  loaded <- ferx_load_fit(path)
  expect_null(loaded$reader_settings)
  expect_null(loaded$population_fingerprint)
  out <- suppressWarnings(ferx_covariance(loaded))
  expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
})
