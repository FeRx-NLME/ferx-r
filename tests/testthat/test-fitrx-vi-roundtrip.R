# A VI fit through ferx_save_fit() / ferx_load_fit() (#508). VI reports a
# bound, not an objective, so under the default `vi_final_ofv = none` its
# OFV, AIC and BIC are NaN by design; the bundle has to bring them back as
# NaN, the method as the engine's "VI" label, and the VI result itself
# (`fit$vi`, #507) in ferx-core's own `vi` slot (ferx-core #2058).

vi_fit <- function(final_ofv = NULL, example = "warfarin") {
  ex <- ferx_example(example)
  src <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", src)
  # Keep the example's own options (an IOV model's `iov_column`), not its
  # method or covariance step.
  opts <- src[-seq_len(at)]
  opts <- opts[!grepl("^\\s*(method|covariance)\\s*=", opts)]
  src <- c(
    src[seq_len(at - 1L)],
    "[fit_options]",
    opts,
    "  method     = vi",
    "  vi_iters   = 30",
    "  vi_seed    = 42",
    if (!is.null(final_ofv)) paste0("  vi_final_ofv = ", final_ofv),
    "  covariance = false"
  )
  model <- tempfile(fileext = ".ferx")
  writeLines(src, model)
  suppressWarnings(ferx_fit(model, ex$data, verbose = FALSE))
}

test_that("a VI fit carries its VI result as fit$vi (#507)", {
  fit <- vi_fit()
  v <- fit$vi
  expect_type(v, "list")
  for (field in c("neg_two_elbo", "data_term", "kl_term",
                  "elbo_tightness_ratio")) {
    expect_true(is.double(v[[field]]) && length(v[[field]]) == 1L,
                info = field)
  }
  expect_true(is.finite(v$neg_two_elbo))
  expect_identical(v$n_iterations, as.numeric(length(v$elbo_trace)))
  expect_type(v$converged, "logical")
  expect_true(is.character(v$family) && nzchar(v$family))
  expect_true(is.character(v$kl) && nzchar(v$kl))
  # Unregularized: the engine stores the objective trace as empty, meaning
  # "the bound itself".
  expect_identical(v$objective_trace, v$elbo_trace)
  expect_identical(v$superseded_by, NA_character_)

  # Posterior means keyed by subject, named as the output convention asks.
  expect_s3_class(v$eta_means, "data.frame")
  expect_identical(names(v$eta_means), c("ID", fit$eta_names))
  expect_identical(v$eta_means$ID, fit$subject_ids)
  expect_identical(names(v$eta_covs), fit$subject_ids)
  cov1 <- v$eta_covs[[1L]]
  expect_identical(dimnames(cov1), list(fit$eta_names, fit$eta_names))
  expect_true(all(diag(cov1) > 0))
  # No IOV: one empty occasion block per subject.
  expect_length(v$kappa_means, length(fit$subject_ids))
  expect_true(all(vapply(v$kappa_means, nrow, 1L) == 0L))
})

test_that("a VI fit round-trips fit$vi, method \"VI\" and its NaN OFV (#508)", {
  fit <- vi_fit()
  # The premise: a live VI fit is labelled by the engine and has no OFV.
  expect_identical(fit$method, "VI")
  expect_true(is.nan(fit$ofv))
  expect_false(is.null(fit$vi))

  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)

  # NaN is written as JSON null (JSON has no NaN), not as a string.
  wire <- fitrx_fit_json(path)
  expect_identical(wire$method, "vi")
  expect_null(wire$ofv)
  # The VI result is in the engine's slot and layout, keyed by subject ID.
  expect_false(is.null(wire$vi))
  expect_null(wire$r_extras$vi)
  expect_identical(unlist(wire$vi$eta_names), fit$eta_names)
  expect_identical(vapply(wire$vi$eta_posterior, `[[`, "", "id"),
                   fit$subject_ids)
  expect_null(wire$vi$objective_trace) # omitted when it is the bound

  back <- ferx_load_fit(path)
  expect_identical(back$method, "VI")
  expect_identical(back$method_chain, fit$method_chain)
  for (field in c("ofv", "ofv_data", "aic", "bic")) {
    expect_length(back[[field]], 1L)
    expect_identical(is.nan(back[[field]]), is.nan(fit[[field]]), info = field)
  }
  expect_true(is.nan(back$ofv))
  expect_identical(back$vi, fit$vi)

  # What a VI fit carries besides: the estimates, and the posterior means as
  # the EBEs.
  expect_equal(back$theta, fit$theta, tolerance = 1e-12)
  expect_equal(back$omega, fit$omega, tolerance = 1e-12)
  # (A reloaded sigma carries names a live one does not; not VI-specific.)
  expect_equal(unname(back$sigma), unname(fit$sigma), tolerance = 1e-12)
  eta_cols <- setdiff(names(fit$ebe_etas), "ID")
  expect_equal(back$ebe_etas[eta_cols], fit$ebe_etas[eta_cols], tolerance = 1e-12)

  # And it survives a second trip unchanged.
  path2 <- tempfile(fileext = ".fitrx")
  ferx_save_fit(back, path2)
  again <- ferx_load_fit(path2)
  expect_identical(again$method, "VI")
  expect_true(is.nan(again$ofv))
  expect_identical(again$vi, fit$vi)
})

test_that("ferx-core loads an R VI bundle and R loads core's (#508)", {
  fit <- vi_fit()
  trip <- fitrx_engine_trip(fit)
  # The engine's writer puts back the `vi` section it read, byte for byte as
  # data: R writes the engine's layout.
  r_vi <- fitrx_fit_json(trip$r_path)$vi
  core_vi <- fitrx_fit_json(trip$core_path)$vi
  expect_false(is.null(core_vi))
  expect_identical(core_vi, r_vi)
  # And R reads the engine's bundle back to the same fit$vi.
  expect_identical(ferx_load_fit(trip$core_path)$vi, fit$vi)
})

test_that("a distinct objective trace and a superseding stage round-trip (#508)", {
  # A live warfarin fit takes the sentinel branches only (objective trace ==
  # bound, nothing superseded); set both the other way, as a regularized
  # `[covariate_nn]` fit and a `methods = [vi, focei]` chain would, so a
  # writer that always omitted either field fails here.
  fit <- vi_fit()
  fit$vi$objective_trace <- fit$vi$elbo_trace + 1.5
  fit$vi$superseded_by <- "FOCEI"

  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  wire <- fitrx_fit_json(path)$vi
  expect_equal(unlist(wire$objective_trace), fit$vi$objective_trace,
               tolerance = 0)
  expect_identical(wire$superseded_by, "FOCEI")
  expect_identical(ferx_load_fit(path)$vi, fit$vi)

  # The engine reads both and writes them back.
  trip <- fitrx_engine_trip(fit)
  expect_identical(fitrx_fit_json(trip$core_path)$vi,
                   fitrx_fit_json(trip$r_path)$vi)
  expect_identical(ferx_load_fit(trip$core_path)$vi, fit$vi)
})

test_that("an IOV VI fit round-trips its per-occasion kappa means (#508)", {
  fit <- vi_fit(example = "warfarin_iov")
  km <- fit$vi$kappa_means
  expect_identical(names(km), fit$subject_ids)
  expect_true(all(vapply(km, nrow, 1L) > 0L))
  expect_identical(colnames(km[[1L]]), fit$kappa_names)

  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  expect_identical(ferx_load_fit(path)$vi, fit$vi)
  trip <- fitrx_engine_trip(fit)
  expect_identical(ferx_load_fit(trip$core_path)$vi, fit$vi)
})

test_that("a VI posterior filed under the wrong subject is refused (#508)", {
  fit <- vi_fit()
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  bad <- fitrx_edit_bundle(path, function(w) {
    w$vi$eta_posterior[[2L]]$id <- "not-a-subject"
    w
  })
  expect_error(ferx_load_fit(bad), "vi.eta_posterior\\[1\\] is subject")
})

test_that("a non-VI fit writes no vi section and reloads without one (#508)", {
  fit <- suppressWarnings(warfarin_fit())
  expect_null(fit$vi)
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  expect_false("vi" %in% names(fitrx_fit_json(path)))
  expect_null(ferx_load_fit(path)$vi)
})

test_that("a VI fit with a Laplace OFV keeps the finite OFV (#508)", {
  fit <- vi_fit("laplace")
  expect_identical(fit$method, "VI")
  expect_true(is.finite(fit$ofv))
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  back <- ferx_load_fit(path)
  expect_identical(back$method, "VI")
  expect_identical(back$ofv, fit$ofv)
  expect_identical(back$aic, fit$aic)
  expect_identical(back$vi, fit$vi)
})
