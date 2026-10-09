
# ---- header from test-persist.R ----
# ferx_save_fit() / ferx_load_fit() round-trip and CLI-flag tests.
#
# These exercise the .fitrx bundle format. The cross-tool test (Rust CLI
# writes, R reads) lives in ferx-core; here we only verify the R-side
# writer + reader are inverses, and that the `output` argument of
# ferx_fit() saves a valid bundle.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.fitrx_wire_to_fit      <- getFromNamespace(".fitrx_wire_to_fit",      "ferx")
.fitrx_build_iov_wire   <- getFromNamespace(".fitrx_build_iov_wire",   "ferx")





















test_that("bad path errors cleanly", {
  expect_error(ferx_load_fit("does-not-exist.fitrx"), "does not exist")
  bogus <- tempfile(fileext = ".fitrx")
  writeLines("not a zip", bogus)
  on.exit(unlink(bogus), add = TRUE)
  expect_error(ferx_load_fit(bogus))
})
test_that("ferx_load_fit on old .fitrx without init_as_sd produces empty logicals", {
  # Build a fake wire list that omits init_as_sd (pre-PR#57 bundle)
  wire_omega <- list(
    names = list("ETA_CL"),
    matrix = list(list(0.09)),
    fixed = list(FALSE),
    log_transformed = list(TRUE),
    shrinkage = list(0.1),
    se = NULL,
    param_corr = NULL
    # init_as_sd intentionally absent
  )
  wire_sigma <- list(
    names = list("PROP_ERR"),
    estimates = list(0.01),
    fixed = list(FALSE),
    types = list("proportional"),
    se = NULL
    # init_as_sd intentionally absent
  )
  wire <- list(
    method = "focei", method_chain = list("focei"),
    converged = TRUE, ofv = -100, aic = -90, bic = -80,
    n_obs = 10L, n_subjects = 2L, n_parameters = 2L, n_iterations = 5L,
    interaction = TRUE, wall_time_secs = 1.0,
    gradient_method_inner = "", gradient_method_outer = "",
    covariance_status = "not_requested",
    omega = wire_omega, sigma = wire_sigma,
    theta = list(names = list("TVCL"), estimates = list(1.0),
                 fixed = list(FALSE), transform = list("log"), se = NULL)
  )
  result <- .fitrx_wire_to_fit(wire)
  expect_equal(result$omega_init_as_sd, logical(0L))
  expect_equal(result$sigma_init_as_sd, logical(0L))
  expect_equal(result$kappa_init_as_sd, logical(0L))
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


test_that(".fitrx_method_label is the inverse for every known token", {
  # One entry per token ferx-core's `method_to_str()` can write into a .fitrx. The test
  # name has always claimed "every known token" but the list had drifted: imp, impmap and
  # bayes were missing and fell through the `as.character()` passthrough, so a reloaded IMP
  # fit reported its method as lowercase "imp" instead of "IMP". Keep this exhaustive — a
  # missing entry does not error, it silently reports the wrong label.
  expect_identical(ferx:::.fitrx_method_label(NULL), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label(""), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label("foce"), "FOCE")
  expect_identical(ferx:::.fitrx_method_label("focei"), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label("laplace"), "LAPLACE")
  expect_identical(ferx:::.fitrx_method_label("agq"), "AGQ")
  expect_identical(ferx:::.fitrx_method_label("foce_gn"), "FOCE-GN")
  expect_identical(ferx:::.fitrx_method_label("foce_gn_hybrid"), "FOCE-GN-Hybrid")
  expect_identical(ferx:::.fitrx_method_label("saem"), "SAEM")
  expect_identical(ferx:::.fitrx_method_label("imp"), "IMP")
  expect_identical(ferx:::.fitrx_method_label("impmap"), "IMPMAP")
  expect_identical(ferx:::.fitrx_method_label("bayes"), "BAYES")
  expect_identical(ferx:::.fitrx_method_label("other"), "other") # unknown passthrough
})
test_that(".fitrx_covariance_status_label maps tokens back to display form", {
  expect_identical(ferx:::.fitrx_covariance_status_label(NULL), "NotRequested")
  expect_identical(ferx:::.fitrx_covariance_status_label("computed"), "Computed")
  expect_identical(ferx:::.fitrx_covariance_status_label("failed"), "Failed")
  expect_identical(ferx:::.fitrx_covariance_status_label("not_requested"), "NotRequested")
  expect_identical(ferx:::.fitrx_covariance_status_label("sir_fallback"), "SirFallback")
  expect_identical(ferx:::.fitrx_covariance_status_label("other"), "other")
})
test_that(".fitrx_matrix_from_wire rebuilds the matrix and guards bad input", {
  expect_null(ferx:::.fitrx_matrix_from_wire(NULL))
  expect_null(ferx:::.fitrx_matrix_from_wire(list(rows = 0L, cols = 0L, data = numeric())))
  expect_null(ferx:::.fitrx_matrix_from_wire(list(rows = 2L, cols = 3L, data = c(1, 2)))) # length mismatch
})
test_that(".fitrx_unwrap_named_se names values when lengths match", {
  expect_null(ferx:::.fitrx_unwrap_named_se(NULL, NULL))
  expect_identical(
    ferx:::.fitrx_unwrap_named_se(c(1, 2), c("a", "b")),
    stats::setNames(c(1, 2), c("a", "b"))
  )
  # Mismatched name length -> unnamed values
  expect_identical(ferx:::.fitrx_unwrap_named_se(c(1, 2), c("only_one")), c(1, 2))
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

  # A bundle carries no packed estimate (ferx-core#1815), so the live side
  # drops it too (#511).
  fit$packed_estimate <- NULL
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

    # A bundle carries no packed estimate (ferx-core#1815), so the live side
    # drops it too (#511).
    fit$packed_estimate <- NULL
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

test_that("predictions.csv IDs read back as the engine's sdtab number (#475)", {
  # Rust's f64 parse where it succeeds, else the subject's 1-based position.
  expect_identical(
    ferx:::.fitrx_sdtab_id_number(c("NA", "NA", "007", "Smith, 2019", "1.",
                                    "+.5", "1E2", "-inf", "Infinity", "nan",
                                    "0x10", " 7")),
    c(1, 1, 7, 3, 1, 0.5, 100, -Inf, Inf, NaN, 10, 11)
  )
  # A subject with no observations has no rows but keeps its position, so the
  # positions come from every subject's ID (ebes.csv), not the block count.
  expect_identical(
    ferx:::.fitrx_sdtab_id_number(c("B", "B", "C", "7"), c("Z", "B", "C", "7")),
    c(2, 2, 3, 7)
  )
  # A block whose text is not among the subjects: block indices, as without.
  expect_identical(ferx:::.fitrx_sdtab_id_number(c("B", "Q"), c("Z", "B")), c(1, 2))
})

test_that("ferx_load_fit numbers a core-spelled predictions.csv past a zero-obs subject (#475)", {
  # ferx-core's writer puts the text ID in predictions.csv and leaves out a
  # subject with no observations; its sdtab number for `B` is then 2, not 1.
  fake <- structure(
    list(
      sdtab = data.frame(ID = c(2, 2), TIME = c(0.5, 1), DV = c(5, 6),
                         PRED = c(5, 6), IPRED = c(5, 6), CWRES = c(0, 0),
                         IWRES = c(0, 0)),
      ebe_etas = data.frame(ID = c("Z", "B"), ETA_CL = c(0, 0.1),
                            stringsAsFactors = FALSE),
      theta = c(TVCL = 1.0),
      omega = matrix(0.04, 1L, 1L, dimnames = list("ETA_CL", "ETA_CL")),
      eta_names = "ETA_CL",
      sigma = c(prop = 0.05), sigma_names = "prop", sigma_types = "proportional",
      theta_fixed = FALSE, omega_fixed = FALSE, sigma_fixed = FALSE,
      ofv = 0, aic = 2, bic = 4,
      n_obs = 2L, n_subjects = 2L, n_parameters = 1L, n_iterations = 1L,
      method = "FOCEI", method_chain = "FOCEI", converged = TRUE,
      warnings = character(), shrinkage_eta = 0, shrinkage_eps = 0,
      wall_time_secs = 0, model_name = "fake", ferx_version = "0.1.0",
      gradient_method_inner = "Enzyme AD", gradient_method_outer = "N/A",
      covariance_status = "NotRequested", model_source = "model fake\n",
      data_path = NA_character_
    ),
    class = "ferx_fit"
  )
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fake, path)
  staging <- withr::local_tempdir()
  utils::unzip(path, exdir = staging)
  preds <- file.path(staging, "predictions.csv")
  lines <- readLines(preds)
  lines[-1] <- sub("^[^,]*", "B", lines[-1])
  writeLines(lines, preds)
  core_spelled <- withr::local_tempfile(fileext = ".fitrx")
  utils::zip(core_spelled, list.files(staging, full.names = TRUE), flags = "-j -q")
  expect_identical(ferx_load_fit(core_spelled)$sdtab$ID, c(2, 2))
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
