
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
  expect_identical(ferx:::.fitrx_method_label(NULL), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label(""), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label("foce"), "FOCE")
  expect_identical(ferx:::.fitrx_method_label("focei"), "FOCEI")
  expect_identical(ferx:::.fitrx_method_label("foce_gn"), "FOCE-GN")
  expect_identical(ferx:::.fitrx_method_label("foce_gn_hybrid"), "FOCE-GN-Hybrid")
  expect_identical(ferx:::.fitrx_method_label("saem"), "SAEM")
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
