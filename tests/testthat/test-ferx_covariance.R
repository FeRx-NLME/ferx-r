# Tests for ferx_covariance() - the standalone covariance step (issue #738),
# the covariance-step analogue of ferx_sir().
#
# Happy-path assertions are gated on `skip_if(is.null(fit$cov_matrix), ...)`
# because in the no-autodiff CI build the warfarin FD covariance step
# occasionally fails to converge at maxiter = 30 (same pattern as
# test-ferx_sir.R).

cov_skip <- "covariance step did not converge - skipping"

test_that(".ferx_fit_interaction follows the last estimating (non-IMP) stage", {
  # fit$interaction is never plumbed to R, so the flag is derived from the
  # method chain. A trailing IMP stage is diagnostic-only and must be skipped,
  # otherwise the covariance step would differentiate the wrong (FOCE) NLL.
  expect_true(ferx:::.ferx_fit_interaction(list(method_chain = "FOCEI")))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = "FOCE")))
  # c("focei", "imp"): interaction TRUE despite the terminal IMP.
  expect_true(ferx:::.ferx_fit_interaction(list(method_chain = c("FOCEI", "IMP"))))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = c("FOCE", "IMP"))))
  # Falls back to fit$method when method_chain is absent.
  expect_true(ferx:::.ferx_fit_interaction(list(method = "FOCEI")))
  # Degenerate / all-IMP chains don't error.
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = character(0))))
  expect_false(ferx:::.ferx_fit_interaction(list(method_chain = "IMP")))
})

test_that("ferx_covariance surfaces cov-step warnings in the structured table", {
  # Covariance-step warnings must reach fit$warnings_structured, not just the
  # flat vector: ferx_get_warnings() and the print tally read the structured
  # table. Force the step to fail (bad FD step) so the engine emits a diagnostic.
  fit <- warfarin_fit()
  out <- ferx_covariance(fit, covariance_method = "r")
  skip_if(is.null(out$cov_matrix), cov_skip)

  # A well-conditioned refit should carry no stale critical condition_number row.
  ws <- out$warnings_structured
  if (is.data.frame(ws) && nrow(ws) > 0L) {
    expect_true(all(c("severity", "category", "message") %in% names(ws)))
  }
  # ferx_get_warnings() reads the structured table; it must run without error
  # and return a data frame on the refreshed fit.
  gw <- ferx_get_warnings(out, as_df = TRUE)
  expect_s3_class(gw, "data.frame")
})

test_that("ferx_covariance validates its inputs", {
  fit <- warfarin_fit()  # covariance = FALSE

  expect_error(ferx_covariance("not a fit"), "ferx_fit object")
  expect_error(
    ferx_covariance(fit, covariance_method = "bhhh"),
    "covariance_method"
  )
  expect_error(
    ferx_covariance(fit, covariance_method = c("r", "s")),
    "single string"
  )
})

test_that("ferx_covariance runs end-to-end and populates covariance fields", {
  fit <- warfarin_fit()  # no inline covariance step
  expect_null(fit$cov_matrix)

  out <- ferx_covariance(fit, verbose = FALSE)
  skip_if(is.null(out$cov_matrix), cov_skip)

  expect_s3_class(out, "ferx_fit")
  expect_equal(out$covariance_status, "computed")

  # Covariance matrix is square, symmetric-ish, and named by parameter.
  d <- nrow(out$cov_matrix)
  expect_equal(ncol(out$cov_matrix), d)
  expect_false(is.null(rownames(out$cov_matrix)))
  expect_identical(rownames(out$cov_matrix), colnames(out$cov_matrix))

  # Standard errors populated and named for theta.
  expect_true(is.numeric(out$se_theta))
  expect_equal(length(out$se_theta), length(out$theta))
  expect_identical(names(out$se_theta), names(out$theta))
  expect_true(all(out$se_theta > 0))

  # Derived correlation matrix refreshed alongside cov_matrix.
  expect_false(is.null(out$cor_matrix))
  expect_equal(dim(out$cor_matrix), c(d, d))
  expect_equal(unname(diag(out$cor_matrix)), rep(1, d), tolerance = 1e-8)

  # Non-covariance fields are untouched.
  expect_identical(out$theta, fit$theta)
  expect_identical(out$omega, fit$omega)
  expect_equal(out$ofv, fit$ofv)
})

test_that("ferx_covariance closely reproduces the inline covariance step", {
  # Re-running the covariance step against a fit that already carries an inline
  # covariance (same point, same starting EBEs) reproduces that fit's own
  # covariance matrix and SEs closely - the numerics route through the same
  # engine covariance step. This is the true apples-to-apples parity check
  # (comparing two independently-converged fits would confound the covariance
  # with tiny differences in the optimum).
  #
  # Parity is close but NOT bit-exact: run_covariance recomputes the inner-loop
  # EBEs (it cold-starts warm_etas = None to mirror the inline step), so the
  # FD-Hessian is evaluated around slightly different etas. ferx-core's own
  # parity test (run_covariance_matches_inline_covariance) documents that
  # "strict sub-1e-6 parity is not achievable through this path" and bounds the
  # gap at abs 1e-4 for a small converged synthetic model. The shared warfarin
  # fit (maxiter = 30) shows a larger but stable gap, so assert an absolute
  # bound that still catches any real regression - a wrong interaction flag or
  # wrong EBEs diverge by orders of magnitude more.
  fit <- warfarin_fit_cov()   # covariance = TRUE
  skip_if(is.null(fit$cov_matrix), cov_skip)

  out <- ferx_covariance(fit)
  skip_if(is.null(out$cov_matrix), cov_skip)

  expect_equal(dim(out$cov_matrix), dim(fit$cov_matrix))
  expect_lt(max(abs(unname(out$cov_matrix) - unname(fit$cov_matrix))), 2e-3)
  # se_theta is the sqrt of the cov diagonal, so the same recompute-path divergence
  # lands slightly larger here. The ferx-core f7d52ec pin (categorical endpoints,
  # #900) grew the stable warfarin gap from <2e-3 to a deterministic ~2.8e-3, so
  # bound se_theta at 5e-3 - still orders of magnitude below any real covariance
  # regression (a wrong interaction flag / wrong EBEs).
  expect_lt(max(abs(unname(out$se_theta) - unname(fit$se_theta))), 5e-3)
})

test_that("ferx_covariance can re-run with a different covariance_method", {
  fit <- warfarin_fit()
  out_r   <- ferx_covariance(fit, covariance_method = "r")
  skip_if(is.null(out_r$cov_matrix), cov_skip)
  out_rsr <- ferx_covariance(fit, covariance_method = "rsr")
  skip_if(is.null(out_rsr$cov_matrix), cov_skip)

  # Both produce a usable matrix; the sandwich generally differs from the
  # plain inverse-Hessian, but at minimum both are the same shape and named.
  expect_equal(dim(out_r$cov_matrix), dim(out_rsr$cov_matrix))
  expect_identical(rownames(out_r$cov_matrix), rownames(out_rsr$cov_matrix))
})

test_that("ferx_covariance refuses to run after the model file is tampered with", {
  src <- ferx_example("warfarin")
  tmp <- tempfile("ferx_cov_tamper_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  model_tmp <- file.path(tmp, "warfarin.ferx")
  data_tmp <- file.path(tmp, "warfarin.csv")
  file.copy(src$model, model_tmp)
  file.copy(src$data, data_tmp)

  fit <- ferx_fit(
    model_tmp, data_tmp,
    method = "focei", verbose = FALSE,
    covariance = FALSE, settings = list(maxiter = 30L)
  )
  expect_match(fit$model_hash, "^[0-9a-f]{64}$")

  # Append whitespace to the model - flips the SHA-256.
  cat("\n# tampered\n", file = model_tmp, append = TRUE)

  expect_error(ferx_covariance(fit), "hash mismatch")
})

test_that("ferx_covariance errors when the fit has no recorded model path", {
  fit <- warfarin_fit()
  fit$model_path <- NULL
  expect_error(ferx_covariance(fit), "no recorded model_path")
})

test_that("ferx_covariance labels a block omega's rows column-major", {
  # The second call site of the cov_matrix labels (#367, #437): `ferx_fit()`
  # and this function name the same matrix, and a fit that skipped the
  # covariance step and picked it up here must arrive at identical dimnames.
  ex  <- ferx_example("warfarin_block_omega")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE,
                                   covariance = FALSE))
  out <- suppressWarnings(ferx_covariance(fit))
  skip_if(is.null(out$cov_matrix), cov_skip)

  d <- diag(out$cov_matrix)
  expect_true(all(
    c("ETA_CL,ETA_CL", "ETA_V,ETA_CL", "ETA_KA,ETA_CL",
      "ETA_V,ETA_V", "ETA_KA,ETA_V", "ETA_KA,ETA_KA") %in% names(d)
  ))
  # The held covariances of the partial block are the zeros; no variance is.
  expect_identical(unname(d[["ETA_KA,ETA_CL"]]), 0)
  expect_identical(unname(d[["ETA_KA,ETA_V"]]), 0)
  expect_gt(d[["ETA_V,ETA_V"]], 0)
  expect_identical(rownames(out$cor_matrix), rownames(out$cov_matrix))
})

# -- IOV models: the kappa segment gets its own labels (#437) ------------------
#
# The engine packs theta / omega / sigma / kappa / [mixture] omega overrides /
# [mixture] sigma overrides / block_sigma rho. R used to label the matrix from
# counts, read every coordinate after theta and before sigma as omega, and so
# came up one name short on any IOV model - and dropped all of them.
#
# The labels are now the glue's, from its own walk of that order, and only
# their total length is checked against the matrix. So the tests below pin the
# order of every segment, one model per neighbouring pair: a block omega ahead
# of a diagonal kappa, a diagonal omega ahead of a block kappa, kappa ahead of
# the mixture overrides, and the mixture overrides ahead of rho. (kappa next to
# rho cannot be built: the engine rejects block_sigma with IOV,
# E_BLOCK_SIGMA_IOV_UNSUPPORTED.)
#
# The two IOV variants are checked at both call sites, ferx_fit()'s covariance
# step and ferx_covariance() on a fit that skipped it, one test per site so a
# site that stops producing a matrix shows as a skip of its own.

expect_cov_labels <- function(x, expected) {
  skip_if(is.null(x$cov_matrix), cov_skip)
  expect_identical(dimnames(x$cov_matrix), list(expected, expected))
  expect_identical(rownames(x$cor_matrix), expected)
  expect_identical(names(x$cov_fixed), expected)
}

iov_variant_fit <- function(edit, covariance) {
  ex  <- ferx_example("warfarin_iov")
  txt <- edit(readLines(ex$model))
  mod <- tempfile(fileext = ".ferx")
  writeLines(txt, mod)
  suppressWarnings(ferx_fit(mod, ex$data, method = "foce", verbose = FALSE,
                            covariance = covariance,
                            settings = list(maxiter = 30L)))
}

block_omega_iov <- function(txt) {
  txt <- sub("^\\s*omega ETA_CL ~ .*$",
             "  block_omega (ETA_CL, ETA_V) = [0.07, 0.01, 0.02]", txt)
  txt <- txt[!grepl("^\\s*omega ETA_V\\s+~", txt)]
  stopifnot(sum(grepl("block_omega", txt)) == 1L,
            !any(grepl("^\\s*omega ETA_V", txt)))
  txt
}

block_kappa_iov <- function(txt) {
  txt <- sub("^\\s*kappa KAPPA_CL ~ .*$",
             "  block_kappa (KAPPA_CL, KAPPA_V) = [0.04, 0.005, 0.02]", txt)
  txt <- sub("TVV  \\* exp\\(ETA_V\\)", "TVV  * exp(ETA_V + KAPPA_V)", txt)
  stopifnot(sum(grepl("block_kappa", txt)) == 1L,
            sum(grepl("exp\\(ETA_V \\+ KAPPA_V\\)", txt)) == 1L)
  txt
}

# Mixed block + diagonal omega packs its full lower triangle (6), then sigma,
# then the one kappa.
block_omega_iov_labels <- c(
  "TVCL", "TVV", "TVKA",
  "ETA_CL,ETA_CL", "ETA_V,ETA_CL", "ETA_KA,ETA_CL",
  "ETA_V,ETA_V", "ETA_KA,ETA_V", "ETA_KA,ETA_KA",
  "PROP_ERR", "KAPPA_CL"
)

block_kappa_iov_labels <- c(
  "TVCL", "TVV", "TVKA", "ETA_CL", "ETA_V", "ETA_KA", "PROP_ERR",
  "KAPPA_CL,KAPPA_CL", "KAPPA_V,KAPPA_CL", "KAPPA_V,KAPPA_V"
)

test_that("ferx_fit(): a block omega with a diagonal kappa labels kappa after sigma", {
  expect_cov_labels(iov_variant_fit(block_omega_iov, covariance = TRUE),
                    block_omega_iov_labels)
})

test_that("ferx_covariance(): a block omega with a diagonal kappa labels kappa after sigma", {
  fit <- iov_variant_fit(block_omega_iov, covariance = FALSE)
  expect_cov_labels(suppressWarnings(ferx_covariance(fit)), block_omega_iov_labels)
})

test_that("ferx_fit(): a diagonal omega with a block kappa labels kappa column-major", {
  expect_cov_labels(iov_variant_fit(block_kappa_iov, covariance = TRUE),
                    block_kappa_iov_labels)
})

# Shared by the label test below and the #473 wrong-size refusal.
block_kappa_cov_free_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) fit <<- iov_variant_fit(block_kappa_iov, covariance = FALSE)
    fit
  }
})

test_that("ferx_covariance(): a diagonal omega with a block kappa labels kappa column-major", {
  fit <- block_kappa_cov_free_fit()
  expect_cov_labels(suppressWarnings(ferx_covariance(fit)), block_kappa_iov_labels)
})

# -- IOV fits: the covariance step runs at the fitted kappa (#473) -------------
#
# The twin of the #465 SIR block in test-ferx_sir.R. The binding used to pass
# the fitted kappa only when `fit$omega_iov` was present and the right length,
# and `None` otherwise, which the engine fills from the model file's *initial*
# kappa: standard errors around the wrong point, with no error. Both bindings
# now build their skeleton through one helper that refuses instead. The fixture
# starts KAPPA_CL ten times above where the data put it, so a fallback to the
# initial value cannot hide.
kappa_far_from_init <- function(txt) {
  out <- sub("kappa KAPPA_CL ~ 0.04", "kappa KAPPA_CL ~ 0.4", txt, fixed = TRUE)
  stopifnot(sum(out != txt) == 1L)
  out
}

warfarin_iov_cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) fit <<- iov_variant_fit(kappa_far_from_init, covariance = TRUE)
    fit
  }
})

test_that("ferx_covariance on an IOV fit reproduces the inline covariance step (#473)", {
  fit <- warfarin_iov_cov_fit()
  skip_if(is.null(fit$cov_matrix), cov_skip)
  expect_lt(fit$omega_iov[1L, 1L], 0.4 / 5)

  out <- suppressWarnings(ferx_covariance(fit))
  skip_if(is.null(out$cov_matrix), cov_skip)
  # The bound is borrowed from the warfarin parity test above, where re-solving
  # the EBEs leaves a small gap; on this fixture the measured gap is 0. Centred
  # on the initial kappa instead it is 0.211 (kappa SE: 0.390).
  expect_lt(max(abs(unname(out$cov_matrix) - unname(fit$cov_matrix))), 2e-3)
  expect_lt(abs(out$se_kappa - fit$se_kappa), 5e-3)
})

test_that("ferx_covariance refuses a kappa fit that has lost its omega_iov (#473)", {
  fit <- warfarin_iov_cov_fit()
  fit$omega_iov <- NULL
  # The prefix is the shared builder's only per-binding input: anchoring it
  # catches the two call sites swapping labels.
  expect_error(ferx_covariance(fit), "^ferx_covariance: .*carries no omega_iov",
               info = "ferx_covariance side of the shared skeleton (#473)")
})

test_that("ferx_covariance refuses a kappa matrix of the wrong size (#473)", {
  fit <- block_kappa_cov_free_fit()
  expect_identical(dim(fit$omega_iov), c(2L, 2L))
  fit$omega_iov <- fit$omega_iov[1L, 1L, drop = FALSE]
  expect_error(ferx_covariance(fit),
               "omega_iov dim 1 does not match model (2 expected)", fixed = TRUE,
               info = "ferx_covariance side of the shared skeleton (#473)")
})

# A two-class [mixture] model with an omega(2) and a sigma(2) override, plus
# either a kappa (on the warfarin_iov data) or a block_sigma correlation (on
# the warfarin data). Checked through ferx_fit() only: ferx_covariance() does
# not compute a matrix for a [mixture] fit today (its Hessian comes out flat in
# the class-2 theta and the mixing probability), so there is nothing to label.
mixture_fit <- function(kappa = FALSE, rho = FALSE) {
  cl_eta <- if (kappa) "ETA_CL + KAPPA_CL" else "ETA_CL"
  mod <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL1(0.1, 0.001, 10.0)",
    "  theta TVCL2(0.2, 0.001, 10.0)",
    "  theta TVV(8.0, 0.1, 500.0)",
    "  theta TVKA(1.0, 0.01, 50.0)",
    "  theta P1(0.5, 0.01, 0.99)",
    "  omega ETA_CL ~ 0.07",
    "  omega ETA_V ~ 0.02",
    if (kappa) "  kappa KAPPA_CL ~ 0.04",
    if (rho) "  block_sigma (PROP_ERR, ADD_ERR) = [0.01, 0.005, 0.10]"
    else "  sigma PROP_ERR ~ 0.01",
    "[mixture]",
    "  nsub = 2",
    "  p(1) = P1",
    "  omega(2) ETA_CL ~ 0.10",
    "  sigma(2) PROP_ERR ~ 0.02",
    "[individual_parameters]",
    sprintf("  CL = if (MIXNUM == 1) TVCL1 * exp(%s) else TVCL2 * exp(%s)",
            cl_eta, cl_eta),
    "  V  = TVV * exp(ETA_V)",
    "  KA = TVKA",
    "[structural_model]",
    "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
    "[error_model]",
    if (rho) "  DV ~ combined(PROP_ERR, ADD_ERR)" else "  DV ~ proportional(PROP_ERR)",
    if (kappa) c("[fit_options]", "  iov_column = OCC")
  ), mod)
  data <- ferx_example(if (kappa) "warfarin_iov" else "warfarin")$data
  suppressWarnings(ferx_fit(mod, data, method = "foce", verbose = FALSE,
                            covariance = TRUE, settings = list(maxiter = 30L)))
}

test_that("kappa packs ahead of the [mixture] overrides, which carry their class", {
  expect_cov_labels(mixture_fit(kappa = TRUE), c(
    "TVCL1", "TVCL2", "TVV", "TVKA", "P1", "ETA_CL", "ETA_V", "PROP_ERR",
    "KAPPA_CL", "ETA_CL (class 2)", "PROP_ERR (class 2)"
  ))
})

test_that("the [mixture] overrides pack ahead of the block_sigma correlation", {
  expect_cov_labels(mixture_fit(rho = TRUE), c(
    "TVCL1", "TVCL2", "TVV", "TVKA", "P1", "ETA_CL", "ETA_V",
    "PROP_ERR", "ADD_ERR", "ETA_CL (class 2)", "PROP_ERR (class 2)",
    "ADD_ERR ~ PROP_ERR"
  ))
})

# ---- subject IDs on the standalone skeleton (#468) ----
# The skeleton FitResult carries the fit's own subject IDs, which the engine
# checks by position against the population it re-reads from the data. Before
# #468 it invented `1..n`, so every other ID set was refused.

test_that("ferx_covariance gives the same answer on any subject labels (#468)", {
  base <- relabelled_fit("warfarin", "identity", identity)
  gappy <- relabelled_fit("warfarin", "gappy", relabel_gappy)
  expect_identical(gappy$ebe_etas$ID, as.character(relabel_gappy(1:10)))

  out_base <- ferx_covariance(base)
  skip_if(is.null(out_base$cov_matrix), cov_skip)
  out_gappy <- ferx_covariance(gappy)
  expect_identical(out_gappy$se_theta, out_base$se_theta)
  expect_identical(out_gappy$cov_matrix, out_base$cov_matrix)
})

test_that("a fit without random effects takes its IDs from individual_estimates (#468)", {
  base <- relabelled_fit("one_cpt_iv_pooled", "identity", identity)
  shifted <- relabelled_fit("one_cpt_iv_pooled", "plus100", function(id) id + 100L)
  # n_eta = 0: no EBE table, so the IDs come from the per-subject table.
  expect_null(shifted$ebe_etas)
  expect_identical(shifted$individual_estimates$ID[1:2], c("101", "102"))

  out_base <- ferx_covariance(base)
  skip_if(is.null(out_base$cov_matrix), cov_skip)
  expect_identical(ferx_covariance(shifted)$se_theta, out_base$se_theta)

  skip_if(is.null(base$cov_matrix) || is.null(shifted$cov_matrix), cov_skip)
  sir <- function(f) ferx_sir(f, sir_samples = 50L, sir_resamples = 20L, sir_seed = 1L)
  sir_base <- sir(base)
  sir_shifted <- sir(shifted)
  expect_identical(sir_shifted$sir_ess, sir_base$sir_ess)
  expect_identical(sir_shifted$sir_ci_theta, sir_base$sir_ci_theta)
})

test_that("the IDs reach the engine's check in fit order, not sorted or re-read (#468)", {
  gappy <- relabelled_fit("warfarin", "gappy", relabel_gappy)
  rotated <- gappy
  ids <- gappy$ebe_etas$ID
  rotated$ebe_etas$ID <- c(ids[-1], ids[1])
  expect_error(
    ferx_covariance(rotated),
    "subject 1 of the population is `110`, but the fit's is `113`",
    fixed = TRUE
  )
})

test_that("ferx_covariance hands the binding the fit's IDs verbatim (#468)", {
  skip_if_not_installed("mockery")
  gappy <- relabelled_fit("warfarin", "gappy", relabel_gappy)
  cap <- capture_binding_args()
  mockery::stub(ferx_covariance, "ferx_rust_covariance", cap$fake)
  expect_error(ferx_covariance(gappy), "captured")
  args <- cap$seen()
  expect_identical(args$subject_ids, gappy$ebe_etas$ID)
  expect_null(args$n_subjects)
})

test_that(".ferx_fit_subject_ids reads ebe_etas, then individual_estimates (#468)", {
  subject_ids <- getFromNamespace(".ferx_fit_subject_ids", "ferx")
  both <- list(
    ebe_etas = data.frame(ID = c("PT1", "PT2"), ETA_CL = c(0, 0)),
    individual_estimates = data.frame(ID = c("X", "Y"))
  )
  expect_identical(subject_ids(both, "ferx_sir"), c("PT1", "PT2"))
  pooled <- list(ebe_etas = NULL, individual_estimates = data.frame(ID = c("7", "9")))
  expect_identical(subject_ids(pooled, "ferx_sir"), c("7", "9"))
  # A numeric ID column (hand-built list) still crosses as text.
  numeric_ids <- list(individual_estimates = data.frame(ID = c(7, 9)))
  expect_identical(subject_ids(numeric_ids, "ferx_sir"), c("7", "9"))
})

test_that(".ferx_fit_subject_ids refuses a fit with no IDs, naming both fields (#468)", {
  subject_ids <- getFromNamespace(".ferx_fit_subject_ids", "ferx")
  none <- list(ebe_etas = NULL, individual_estimates = NULL, n_subjects = 10L)
  err <- tryCatch(subject_ids(none, "ferx_covariance"), error = conditionMessage)
  expect_match(err, "^ferx_covariance: the fit carries no subject IDs")
  expect_match(err, "fit$ebe_etas$ID", fixed = TRUE)
  expect_match(err, "fit$individual_estimates$ID", fixed = TRUE)
  expect_match(err, "cannot be matched to the data", fixed = TRUE)
  expect_match(err, "Re-fit via ferx_fit(model, data).", fixed = TRUE)
  # A recorded subject count is not an ID set: it must not be used, or named.
  expect_no_match(err, "n_subjects", fixed = TRUE)
  expect_no_match(err, "SIR", fixed = TRUE)
  expect_no_match(err, "data file", fixed = TRUE)
  sir_err <- tryCatch(subject_ids(none, "ferx_sir"), error = conditionMessage)
  expect_match(sir_err, "^ferx_sir: ")
  expect_no_match(sir_err, "covariance", fixed = TRUE)
})

test_that(".ferx_fit_subject_ids refuses NA IDs, counting them (#468)", {
  subject_ids <- getFromNamespace(".ferx_fit_subject_ids", "ferx")
  with_na <- list(ebe_etas = data.frame(ID = c("1", NA, NA), ETA_CL = 0))
  err <- tryCatch(subject_ids(with_na, "ferx_covariance"), error = conditionMessage)
  expect_match(err, "^ferx_covariance: 2 of the 3 subject IDs in fit\\$ebe_etas\\$ID are NA\\.")
  expect_match(err, "matched to the data subject by subject", fixed = TRUE)
  expect_match(err, "every one must be present", fixed = TRUE)
  expect_no_match(err, "covariance = TRUE", fixed = TRUE)
})

test_that(".ferx_fit_subject_ids never borrows IDs for EBE rows without an ID column (#468)", {
  subject_ids <- getFromNamespace(".ferx_fit_subject_ids", "ferx")
  # The warm-start is built from these rows; another table's IDs carry no
  # guarantee of the same order.
  no_id <- list(ebe_etas = data.frame(ETA_CL = c(0.1, 0.2)),
                individual_estimates = data.frame(ID = c("B", "A")))
  err <- tryCatch(subject_ids(no_id, "ferx_sir"), error = conditionMessage)
  expect_match(err, "^ferx_sir: fit\\$ebe_etas has no ID column")
  expect_match(err, "cannot be matched to the data's subjects", fixed = TRUE)
})

# --- #462: the fit's own data selection -----------------------------------
#
# A record-only selection (`ferx_fit(ignore =, accept =, ignore_ids =)`) used to
# be lost on the way to the standalone step, which re-read the data with the
# model file's selection and scored other rows than the fit. See
# helper-reader-settings.R for the oracle and the margins.

test_that("ferx_covariance on an `ignore =` fit reproduces the in-fit step (#462)", {
  fit <- rs_fit(ignore = rs_ignore)
  skip_if(is.null(fit$cov_matrix), cov_skip)
  expect_type(fit$reader_settings, "character")
  expect_type(fit$population_fingerprint, "character")
  out <- ferx_covariance(fit)
  expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
})

test_that("ferx_covariance on an `accept =` fit reproduces the in-fit step (#462)", {
  fit <- rs_fit(accept = "TIME < 120")
  # The settings the fit reads with are the ones it records, so a selection
  # lost from them is lost from the fit too: assert the fit applied it.
  expect_identical(fit$exclusions$fired_accept, "accept: TIME < 120")
  skip_if(is.null(fit$cov_matrix), cov_skip)
  out <- ferx_covariance(fit)
  expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
})

test_that("ferx_covariance on an `ignore_ids =` fit runs and reproduces the in-fit step (#462)", {
  # Refused before #462: "the population has 10 subjects but the fit has 7".
  fit <- rs_fit(ignore_ids = 1:3)
  expect_identical(nrow(fit$ebe_etas), 7L)
  skip_if(is.null(fit$cov_matrix), cov_skip)
  out <- ferx_covariance(fit)
  expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
})

test_that("the engine checks the re-read population against the fit's fingerprint (#462)", {
  # The only test that sees the fingerprint go missing: every number above
  # still matches with the settings alone. Editing the recorded selection makes
  # the re-read keep a record the fit dropped.
  fit <- rs_fit(ignore = rs_ignore)
  edited <- sub("DV < 1.0", "DV < 0.9", fit$reader_settings, fixed = TRUE)
  expect_false(identical(edited, fit$reader_settings))
  fit$reader_settings <- edited
  expect_error(ferx_covariance(fit), "is not the one the fit was given",
               fixed = TRUE)
})

test_that("a malformed fit$reader_settings is refused by name, not a panic (#462)", {
  fit <- rs_fit(ignore = rs_ignore)
  fit$reader_settings <- "{not json"
  expect_error(ferx_covariance(fit), "ferx_covariance: `fit$reader_settings`",
               fixed = TRUE)
  fit <- rs_fit(ignore = rs_ignore)
  fit$population_fingerprint <- "{not json"
  expect_error(ferx_covariance(fit), "ferx_covariance: `fit$population_fingerprint`",
               fixed = TRUE)
})

test_that("a fit records a `settings =` iov_column in its reader settings (#462)", {
  ex <- ferx_example("warfarin_iov")
  data <- withr::local_tempfile(fileext = ".csv")
  rows <- utils::read.csv(ex$data)
  rows$VISIT <- rows$OCC
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  fit <- suppressWarnings(ferx_fit(ex$model, data, method = "foce", verbose = FALSE,
                                   covariance = FALSE,
                                   settings = list(maxiter = 2L, iov_column = "VISIT")))
  rs <- jsonlite::fromJSON(fit$reader_settings)
  expect_identical(rs$iov_column, "VISIT")
})

# -- The legacy guard: a fit that records no selection ----------------------

test_that("a legacy fit with a record-only `ignore =` is refused, naming the clause (#462)", {
  fit <- rs_legacy(rs_fit(ignore = rs_ignore))
  e <- tryCatch(ferx_covariance(fit), error = function(e) e)
  expect_s3_class(e, "error")
  msg <- conditionMessage(e)
  expect_match(msg, "`ignore: EVID == 0 && DV < 1.0`", fixed = TRUE)
  expect_match(msg, "predates the record of the data selection", fixed = TRUE)
  expect_match(msg, "ferx_fit(..., covariance = TRUE)", fixed = TRUE)
  expect_no_match(msg, "edited", fixed = TRUE)
  expect_no_match(msg, "reads the file differently", fixed = TRUE)
})

test_that("a legacy fit with `ignore_ids =` is refused, naming the subjects (#462)", {
  fit <- rs_legacy(rs_fit(ignore_ids = 1:3))
  e <- tryCatch(ferx_covariance(fit), error = function(e) e)
  msg <- conditionMessage(e)
  expect_match(msg, "subject(s) 1, 2, 3 (from `ignore_ids` / `ignore_subjects`)",
               fixed = TRUE)
  expect_no_match(msg, "`ignore_subjects: 1`", fixed = TRUE)
})

test_that("a legacy fit whose selection the model file states still runs (#462)", {
  fit <- rs_legacy(warfarin_sel_fit())
  expect_gt(length(fit$exclusions$fired_ignore), 0L)
  out <- ferx_covariance(fit)
  expect_true(is.numeric(out$se_theta))
})

test_that("a legacy fit whose record-only clause fired nothing still runs (#462)", {
  fit <- rs_legacy(rs_fit(ignore = "DV < -1"))
  expect_length(fit$exclusions$fired_ignore, 0L)
  out <- ferx_covariance(fit)
  expect_true(is.numeric(out$se_theta))
})

test_that("a legacy fit with no exclusion record runs as before (#462)", {
  fit <- rs_legacy(warfarin_fit_cov())
  expect_null(fit$exclusions)
  skip_if(is.null(fit$cov_matrix), cov_skip)
  out <- ferx_covariance(fit)
  # M1: no selection, unchanged.
  expect_lt(rs_rel(out$se_theta, fit$se_theta), 1e-8)
})

test_that("a legacy fit on an edited model file hears the hash mismatch, not the guard (#462)", {
  fit <- rs_legacy(rs_fit(ignore = rs_ignore))
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(c(readLines(fit$model_path), "# edited"), model)
  fit$model_path <- model
  e <- tryCatch(ferx_covariance(fit), error = function(e) e)
  msg <- conditionMessage(e)
  expect_match(msg, "hash mismatch", fixed = TRUE)
  expect_no_match(msg, "predates the record", fixed = TRUE)
})

test_that("the guard labels the model file's clauses through the engine's parse (#462)", {
  # A quoted spelling in the file reads as the same clause the reader logs;
  # gluing "ignore: " onto the raw text would keep the quotes and refuse.
  model <- withr::local_tempfile(fileext = ".ferx")
  ex <- ferx_example("warfarin_data_selection")
  lines <- readLines(ex$model)
  at <- grep("ignore = DV < 1.0", lines, fixed = TRUE)
  expect_length(at, 1L)
  lines[at] <- '  ignore = "DV < 1.0"'
  writeLines(lines, model)
  unstated <- ferx:::ferx_rust_unstated_selection(
    "probe", model, "", c("ignore: DV < 1.0", "ignore: TIME > 100", "ignore_subjects: 4"),
    "accept: TIME < 120", ""
  )
  # Only the clauses the file does not state, the file's own dropped.
  expect_identical(unstated,
                   c("ignore: TIME > 100", "ignore_subjects: 4", "accept: TIME < 120"))
})

test_that("the guard compares a call-time iov_column with the file's, case-insensitively (#526)", {
  ex <- ferx_example("warfarin_iov")
  unstated <- function(col) {
    ferx:::ferx_rust_unstated_selection("probe", ex$model, "", character(), character(), col)
  }
  # The file says `iov_column = OCC`.
  expect_identical(unstated("VISIT"), "iov_column: VISIT")
  expect_identical(unstated("occ"), character())
  expect_identical(unstated(""), character())
})

test_that("a legacy fit read with a `settings =` iov_column is refused, naming it (#526)", {
  fit <- rs_legacy(rs_iov_fit())
  expect_identical(fit$call_settings$iov_column, "VISIT")
  e <- tryCatch(suppressWarnings(ferx_covariance(fit)), error = function(e) e)
  expect_s3_class(e, "error")
  msg <- conditionMessage(e)
  expect_match(msg, "ferx_covariance: this fit predates", fixed = TRUE)
  expect_match(msg, '`settings = list(iov_column = "VISIT")`', fixed = TRUE)
  expect_match(msg, "ferx_fit(..., covariance = TRUE)", fixed = TRUE)
  # The same fit with its settings recorded runs.
  out <- suppressWarnings(ferx_covariance(rs_iov_fit()))
  expect_true(is.numeric(out$se_theta))
})
