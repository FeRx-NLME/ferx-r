# Unit tests for the pure-R serialization helpers in persist.R. These convert
# between the in-memory ferx_fit fields and the cross-language .fitrx wire
# schema, and need no model fit to exercise — so they live in the fast tier.

# ---------------------------------------------------------------------------
# Method name <-> token round-trip
# ---------------------------------------------------------------------------
test_that(".fitrx_method_to_token normalises labels and defaults to focei", {
  expect_identical(ferx:::.fitrx_method_to_token(NULL), "focei")
  expect_identical(ferx:::.fitrx_method_to_token(""), "focei")
  expect_identical(ferx:::.fitrx_method_to_token("FOCEI"), "focei")   # case-insensitive
  expect_identical(ferx:::.fitrx_method_to_token("FOCE"), "foce")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce_gn"), "foce_gn")
  expect_identical(ferx:::.fitrx_method_to_token("foce-gn-hybrid"), "foce_gn_hybrid")
  expect_identical(ferx:::.fitrx_method_to_token("SAEM"), "saem")
  expect_identical(ferx:::.fitrx_method_to_token("Mystery"), "mystery") # unknown -> tolower
})

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

test_that("method label/token round-trips through the canonical display names", {
  for (lbl in c("FOCE", "FOCEI", "FOCE-GN", "FOCE-GN-Hybrid", "SAEM")) {
    expect_identical(ferx:::.fitrx_method_label(ferx:::.fitrx_method_to_token(lbl)), lbl)
  }
})

# ---------------------------------------------------------------------------
# Covariance status <-> token
# ---------------------------------------------------------------------------
test_that(".fitrx_covariance_status_to_token normalises and defaults", {
  expect_identical(ferx:::.fitrx_covariance_status_to_token(NULL), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token(""), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("Computed"), "computed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("FAILED"), "failed")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("NotRequested"), "not_requested")
  expect_identical(ferx:::.fitrx_covariance_status_to_token("weird"), "weird")
})

test_that(".fitrx_covariance_status_label maps tokens back to display form", {
  expect_identical(ferx:::.fitrx_covariance_status_label(NULL), "NotRequested")
  expect_identical(ferx:::.fitrx_covariance_status_label("computed"), "Computed")
  expect_identical(ferx:::.fitrx_covariance_status_label("failed"), "Failed")
  expect_identical(ferx:::.fitrx_covariance_status_label("not_requested"), "NotRequested")
  expect_identical(ferx:::.fitrx_covariance_status_label("other"), "other")
})

# ---------------------------------------------------------------------------
# Error model inference from sigma types
# ---------------------------------------------------------------------------
test_that(".fitrx_error_model_from_sigma_types classifies sigma vectors", {
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(NULL), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(character(0)), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("proportional", "Proportional")), "proportional")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive")), "additive")
  expect_identical(ferx:::.fitrx_error_model_from_sigma_types(c("additive", "proportional")), "combined")
})

# ---------------------------------------------------------------------------
# Matrix <-> wire (row-major)
# ---------------------------------------------------------------------------
test_that(".fitrx_matrix_to_wire emits row-major data and dims", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE) # [[1,2,3],[4,5,6]]
  w <- ferx:::.fitrx_matrix_to_wire(m)
  expect_identical(w$rows, 2L)
  expect_identical(w$cols, 3L)
  expect_identical(w$data, c(1, 2, 3, 4, 5, 6)) # row-major
})

test_that(".fitrx_matrix_to_wire returns NULL for NULL / non-matrix input", {
  expect_null(ferx:::.fitrx_matrix_to_wire(NULL))
  expect_null(ferx:::.fitrx_matrix_to_wire(c(1, 2, 3)))
})

test_that(".fitrx_matrix_from_wire rebuilds the matrix and guards bad input", {
  expect_null(ferx:::.fitrx_matrix_from_wire(NULL))
  expect_null(ferx:::.fitrx_matrix_from_wire(list(rows = 0L, cols = 0L, data = numeric())))
  expect_null(ferx:::.fitrx_matrix_from_wire(list(rows = 2L, cols = 3L, data = c(1, 2)))) # length mismatch
})

test_that("matrix round-trips through wire form unchanged", {
  m <- matrix(c(0.1, 0.0, 0.0, 0.2), nrow = 2)
  expect_equal(ferx:::.fitrx_matrix_from_wire(ferx:::.fitrx_matrix_to_wire(m)), m)
})

# ---------------------------------------------------------------------------
# Confidence-interval <-> wire
# ---------------------------------------------------------------------------
test_that(".fitrx_ci_to_wire accepts a 2-column matrix", {
  ci <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE) # rows (1,2) (3,4)
  w  <- ferx:::.fitrx_ci_to_wire(ci)
  expect_length(w, 2L)
  expect_identical(w[[1L]], c(1, 2))
  expect_identical(w[[2L]], c(3, 4))
})

test_that(".fitrx_ci_to_wire accepts the legacy flat layout and rejects odd length", {
  expect_equal(ferx:::.fitrx_ci_to_wire(c(1, 2, 3, 4)), list(c(1, 2), c(3, 4)))
  expect_null(ferx:::.fitrx_ci_to_wire(c(1, 2, 3)))   # odd length
  expect_null(ferx:::.fitrx_ci_to_wire(NULL))
})

test_that("CI round-trips matrix -> wire -> matrix", {
  ci <- matrix(c(0.5, 1.5, 2.5, 3.5), ncol = 2, byrow = TRUE)
  back <- ferx:::.fitrx_unwrap_ci(ferx:::.fitrx_ci_to_wire(ci))
  expect_equal(back, ci)
  expect_null(ferx:::.fitrx_unwrap_ci(NULL))
})

# ---------------------------------------------------------------------------
# Optional scalar wrappers (NULL/NA/empty collapse to NULL)
# ---------------------------------------------------------------------------
test_that("optional scalar wrappers collapse missing values to NULL", {
  for (f in list(ferx:::.fitrx_opt_num, ferx:::.fitrx_opt_int, ferx:::.fitrx_opt_chr,
                 ferx:::.fitrx_unwrap_opt_num, ferx:::.fitrx_unwrap_opt_int,
                 ferx:::.fitrx_unwrap_opt_chr)) {
    expect_null(f(NULL))
    expect_null(f(character(0)))
  }
  expect_identical(ferx:::.fitrx_opt_num(2.5), 2.5)
  expect_null(ferx:::.fitrx_opt_num(NA))
  expect_identical(ferx:::.fitrx_opt_int(5L), 5L)
  expect_null(ferx:::.fitrx_opt_int(NA))
  expect_identical(ferx:::.fitrx_opt_chr("x"), "x")
  expect_null(ferx:::.fitrx_opt_chr(""))     # empty string is "missing"
  expect_null(ferx:::.fitrx_opt_chr(NA))
})

test_that("optional vector wrappers keep non-empty vectors and drop empties", {
  expect_identical(ferx:::.fitrx_opt_num_vec(c(1, 2, 3)), c(1, 2, 3))
  expect_null(ferx:::.fitrx_opt_num_vec(numeric(0)))
  expect_null(ferx:::.fitrx_opt_num_vec(NULL))
  expect_identical(ferx:::.fitrx_unwrap_opt_num_vec(list(1, 2)), c(1, 2))
  expect_null(ferx:::.fitrx_unwrap_opt_num_vec(NULL))
})

# ---------------------------------------------------------------------------
# Named SE unwrap + subject string IDs
# ---------------------------------------------------------------------------
test_that(".fitrx_unwrap_named_se names values when lengths match", {
  expect_null(ferx:::.fitrx_unwrap_named_se(NULL, NULL))
  expect_identical(
    ferx:::.fitrx_unwrap_named_se(c(1, 2), c("a", "b")),
    stats::setNames(c(1, 2), c("a", "b"))
  )
  # Mismatched name length -> unnamed values
  expect_identical(ferx:::.fitrx_unwrap_named_se(c(1, 2), c("only_one")), c(1, 2))
})

test_that(".fitrx_subject_string_ids extracts character IDs or NULL", {
  expect_null(ferx:::.fitrx_subject_string_ids(list()))
  expect_null(ferx:::.fitrx_subject_string_ids(list(ebe_etas = data.frame(x = 1))))
  fit <- list(ebe_etas = data.frame(ID = c(1, 2, 10), eta = c(0, 0, 0)))
  expect_identical(ferx:::.fitrx_subject_string_ids(fit), c("1", "2", "10"))
})

# ---------------------------------------------------------------------------
# Per-subject OFV / N_OBS extraction from sdtab
# ---------------------------------------------------------------------------
test_that(".fitrx_per_subject_ofv_nobs pulls one row per subject from sdtab", {
  expect_null(ferx:::.fitrx_per_subject_ofv_nobs(list()))
  expect_null(ferx:::.fitrx_per_subject_ofv_nobs(list(sdtab = data.frame(ID = 1)))) # missing cols
  fit <- list(sdtab = data.frame(
    ID = c("A", "A", "B"),
    EBE_OFV = c(1.5, 1.5, 2.5),
    N_OBS = c(3L, 3L, 4L)
  ))
  res <- ferx:::.fitrx_per_subject_ofv_nobs(fit)
  expect_identical(res$id, c("A", "B"))
  expect_identical(res$ofv, c(1.5, 2.5))
  expect_identical(res$n_obs, c(3L, 4L))
})

# ---------------------------------------------------------------------------
# EBEs CSV writer — header-only path when there are no random effects
# ---------------------------------------------------------------------------
test_that(".fitrx_write_ebes_csv writes a stable header when ebes are absent", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  ferx:::.fitrx_write_ebes_csv(list(ebe_etas = NULL), path)
  expect_identical(readLines(path), "ID,ofv_contribution,n_obs")
})
