# Tests for parameter priors (ferx-core #254) as they surface on a ferx_fit,
# and for the two standalone paths that rebuild a skeleton FitResult from the
# fields R holds - ferx_sir() and ferx_covariance() (ferx-r #366).
#
# The engine applies a prior declared inline in the `.ferx` file; nothing in the
# R API switches it on. So these tests write a priored copy of the bundled
# warfarin model to a temp dir and fit that.

prior_cov_skip <- "covariance step did not converge - skipping"

# A copy of the bundled warfarin model with a prior on TVCL, plus the matching
# data path. Written once per test file run.
priored_warfarin_paths <- local({
  paths <- NULL
  function() {
    if (is.null(paths)) {
      ex <- ferx_example("warfarin")
      txt <- readLines(ex$model, warn = FALSE)
      i <- grep("theta TVCL(0.134, 0.001, 10.0)", txt, fixed = TRUE)
      stopifnot(length(i) == 1L)
      txt[i] <- "  theta TVCL(0.134, 0.001, 10.0) prior(0.15, rse = 10%)"
      model <- file.path(tempdir(), "warfarin_prior.ferx")
      writeLines(txt, model)
      paths <<- list(model = model, data = ex$data)
    }
    paths
  }
})

# Cached FOCEI fit of that priored model, with the covariance step on (both the
# SIR and the covariance test need the matrix).
priored_warfarin_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      p <- priored_warfarin_paths()
      fit <<- ferx_fit(p$model, p$data,
                       method = "focei", verbose = FALSE,
                       covariance = TRUE, settings = list(maxiter = 30L))
    }
    fit
  }
})

test_that("an unpriored fit reports the whole objective as the data half", {
  fit <- warfarin_fit()
  expect_equal(fit$ofv_data, fit$ofv)
  expect_equal(fit$ofv_prior, 0)
  expect_null(fit$prior_summary)
  # AIC is computed from the data half, which here is the whole objective.
  expect_equal(fit$aic, fit$ofv_data + 2 * fit$n_parameters)
})

test_that("a priored fit splits its objective and reports each priored parameter", {
  fit <- priored_warfarin_fit()

  # fit$ofv is the penalized total; the two halves add up to it.
  expect_true(fit$ofv_prior > 0)
  expect_equal(fit$ofv, fit$ofv_data + fit$ofv_prior)
  # AIC / BIC stay on the data half - a penalized objective is not a log
  # likelihood, so ranking a priored fit against an unpriored one needs ofv_data.
  expect_equal(fit$aic, fit$ofv_data + 2 * fit$n_parameters)

  ps <- fit$prior_summary
  expect_s3_class(ps, "data.frame")
  expect_equal(nrow(ps), 1L)
  expect_equal(
    names(ps),
    c("name", "prior_value", "estimate", "shift_in_prior_sds", "penalty",
      "family", "prior_lower_95", "prior_upper_95")
  )
  expect_equal(ps$name, "TVCL")
  expect_equal(ps$prior_value, 0.15)
  # A positive declared lower bound packs the coordinate on the log scale, so
  # the realised prior family is lognormal.
  expect_equal(ps$family, "lognormal")
  expect_equal(ps$estimate, unname(fit$theta[["TVCL"]]))
  # The per-parameter penalties sum to ofv_prior, and each is the square of the
  # reported shift in prior SDs.
  expect_equal(sum(ps$penalty), fit$ofv_prior)
  expect_equal(ps$penalty, ps$shift_in_prior_sds^2)
  expect_true(ps$prior_lower_95 < ps$prior_value)
  expect_true(ps$prior_upper_95 > ps$prior_value)
})

test_that("print.ferx_fit annotates the OFV line only for a priored fit", {
  expect_output(print(priored_warfarin_fit()), "prior penalty")
  out <- utils::capture.output(print(warfarin_fit()))
  expect_false(any(grepl("prior penalty", out, fixed = TRUE)))
})

test_that("the prior split survives a .fitrx round-trip", {
  fit <- priored_warfarin_fit()
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  ferx_save_fit(fit, f)
  back <- ferx_load_fit(f)

  expect_equal(back$ofv, fit$ofv)
  expect_equal(back$ofv_data, fit$ofv_data)
  expect_equal(back$ofv_prior, fit$ofv_prior)
  expect_equal(back$prior_summary, fit$prior_summary)
})

test_that("an unpriored bundle round-trips to the unpriored split", {
  fit <- warfarin_fit()
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  ferx_save_fit(fit, f)
  back <- ferx_load_fit(f)

  # The three fields are omitted from the wire for an unpriored fit (matching
  # ferx-core's own writer); the loader reconstructs the unpriored split.
  expect_equal(back$ofv_data, fit$ofv)
  expect_equal(back$ofv_prior, 0)
  expect_null(back$prior_summary)
})

test_that(".ferx_ofv_prior reads 0 off a fit that predates the field", {
  # An older .fitrx bundle carries no ofv_prior at all; 0 is the truth there,
  # because no prior could have been applied when it was written.
  expect_equal(ferx:::.ferx_ofv_prior(list(ofv = 10)), 0)
  expect_equal(ferx:::.ferx_ofv_prior(list(ofv_prior = NA_real_)), 0)
  expect_equal(ferx:::.ferx_ofv_prior(list(ofv_prior = 2.5)), 2.5)
})

test_that("ferx_sir takes the data half of a priored objective as its reference", {
  fit <- priored_warfarin_fit()
  skip_if(is.null(fit$cov_matrix), prior_cov_skip)

  # The engine's SIR reference objective is `ofv - ofv_prior` and it adds the
  # prior penalty back itself. So handing it the penalized total beside the
  # matching prior half must be identical to handing it the data half beside a
  # zero prior half - the two spellings of the same reference. If ferx_sir()
  # stopped passing ofv_prior (or the engine started reading `ofv` directly),
  # the first call would run against a doubly-penalized reference and these
  # would diverge.
  as_data_half <- fit
  as_data_half$ofv <- fit$ofv_data
  as_data_half$ofv_prior <- 0

  a <- ferx_sir(fit, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)
  b <- ferx_sir(as_data_half, sir_samples = 8L, sir_resamples = 4L, sir_seed = 1L)

  expect_equal(a$sir_ess, b$sir_ess)
  expect_equal(a$sir_ci_theta, b$sir_ci_theta)
  expect_equal(a$sir_ci_omega, b$sir_ci_omega)
  expect_equal(a$sir_ci_sigma, b$sir_ci_sigma)
})

test_that("ferx_covariance re-derives the prior curvature from the model file", {
  fit <- priored_warfarin_fit()
  skip_if(is.null(fit$cov_matrix), prior_cov_skip)

  # The standalone covariance step is handed a skeleton FitResult that carries
  # no prior information beyond the model path, and re-derives the penalty's
  # curvature from the model's own `prior(...)` declaration. So its standard
  # errors must match the ones the inline covariance step produced. Not
  # bit-exact: the standalone path cold-starts the EBEs (see PR #253).
  out <- ferx_covariance(fit)
  skip_if(is.null(out$cov_matrix), prior_cov_skip)
  expect_equal(out$se_theta, fit$se_theta, tolerance = 1e-3)

  # And the prior really does enter that curvature: a prior at rse = 10% is
  # tighter than the data's own information on TVCL, so the priored SE is the
  # smaller of the two. Without this the check above would pass with both
  # paths ignoring the prior.
  plain <- warfarin_fit_cov()
  skip_if(is.null(plain$cov_matrix), prior_cov_skip)
  expect_true(out$se_theta[["TVCL"]] < plain$se_theta[["TVCL"]])
})
