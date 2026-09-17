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

# Capture the arguments a binding is called with, then abort. The SIR and
# covariance steps are expensive and their output cannot show what reference
# objective they were handed, so the forwarding is asserted on the call itself.
capture_binding_args <- function() {
  seen <- NULL
  list(
    fake = function(...) {
      seen <<- list(...)
      stop("captured", call. = FALSE)
    },
    seen = function() seen
  )
}

test_that("ferx_sir forwards the fit's own prior half to the binding", {
  skip_if_not_installed("mockery")
  fit <- priored_warfarin_fit()

  # This is the assertion that actually pins ferx-r #366. Nothing in SIR's
  # *output* can: the engine takes its reference objective as
  # `ofv - ofv_prior` and adds the prior penalty back itself, so passing a
  # penalized `ofv` beside a zero prior half shifts `ofv_hat` by a constant -
  # and that constant cancels in the normalized importance weights, leaving the
  # weights, the ESS and every interval bit-identical. The double count is only
  # visible at the call boundary, so that is where it is checked.
  cap <- capture_binding_args()
  mockery::stub(ferx_sir, "ferx_rust_sir", cap$fake)
  expect_error(ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L), "captured")

  args <- cap$seen()
  expect_true(args$ofv_prior > 0)
  expect_equal(args$ofv_prior, fit$ofv_prior)
  # `ofv` still crosses as the penalized total - the engine subtracts the half
  # itself. Handing it the data half here instead would double-subtract.
  expect_equal(args$ofv, fit$ofv)
})

test_that("ferx_covariance forwards the fit's own prior half to the binding", {
  skip_if_not_installed("mockery")
  fit <- priored_warfarin_fit()

  cap <- capture_binding_args()
  mockery::stub(ferx_covariance, "ferx_rust_covariance", cap$fake)
  expect_error(ferx_covariance(fit), "captured")

  args <- cap$seen()
  expect_equal(args$ofv_prior, fit$ofv_prior)
  expect_equal(args$ofv, fit$ofv)
})

test_that("an unpriored fit forwards a zero prior half", {
  skip_if_not_installed("mockery")
  fit <- warfarin_fit_cov()
  skip_if(is.null(fit$cov_matrix), prior_cov_skip)

  cap <- capture_binding_args()
  mockery::stub(ferx_sir, "ferx_rust_sir", cap$fake)
  expect_error(ferx_sir(fit, sir_samples = 4L, sir_resamples = 2L), "captured")
  expect_equal(cap$seen()$ofv_prior, 0)
})

test_that("SIR output is unchanged by which spelling of the reference it gets", {
  fit <- priored_warfarin_fit()
  skip_if(is.null(fit$cov_matrix), prior_cov_skip)

  # This pins the *cancellation*, not the forwarding (the test above does that).
  # `(penalized ofv, matching prior half)` and `(data half, zero prior half)`
  # name the same reference objective, so SIR must not be able to tell them
  # apart. It is worth keeping because the cancellation is what makes the
  # pre-#366 double count harmless: if a future engine used `ofv_hat` as
  # anything but a normalized difference, this would start failing and say so.
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

test_that("a pre-#366 bundle of a priored fit recovers its split from the AIC", {
  fit <- priored_warfarin_fit()
  f <- tempfile(fileext = ".fitrx")
  on.exit(unlink(f), add = TRUE)
  ferx_save_fit(fit, f)

  # Reproduce the object shape ferx_save_fit() wrote before this PR: a priored
  # fit whose bundle carries the penalized `ofv` and no split at all. Reading
  # that as `ofv_data = ofv` would relabel the penalized objective as the
  # likelihood and put ferx_sir() back on the #366 double count.
  staging <- file.path(tempdir(), "legacy_fitrx")
  unlink(staging, recursive = TRUE)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::unzip(f, exdir = staging, junkpaths = TRUE)
  fj <- file.path(staging, "fit.json")
  w <- jsonlite::fromJSON(fj, simplifyVector = FALSE)
  w$ofv_data <- NULL
  w$ofv_prior <- NULL
  w$prior_summary <- NULL
  jsonlite::write_json(w, fj, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null", na = "null")
  legacy <- tempfile(fileext = ".fitrx")
  on.exit(unlink(legacy), add = TRUE)
  owd <- setwd(staging)
  utils::zip(legacy, list.files(staging), flags = "-q")
  setwd(owd)

  expect_warning(back <- ferx_load_fit(legacy), "recovered from the stored AIC")

  expect_equal(back$ofv, fit$ofv)
  expect_equal(back$ofv_data, fit$ofv_data)
  expect_equal(back$ofv_prior, fit$ofv_prior)
  # The per-parameter report was never written, so it cannot come back.
  expect_null(back$prior_summary)
  # The invariant the stored AIC was computed under holds again.
  expect_equal(back$aic, back$ofv_data + 2 * back$n_parameters)
  # The OFV line still reports the split, without a parameter count.
  expect_output(print(back), "prior: data OFV")

  # And re-saving that loaded fit keeps the split, rather than dropping it for
  # want of a prior_summary.
  again <- tempfile(fileext = ".fitrx")
  on.exit(unlink(again), add = TRUE)
  ferx_save_fit(back, again)
  # Loading that is not a guess, so it does not warn a second time.
  expect_no_warning(back2 <- ferx_load_fit(again))
  expect_equal(back2$ofv_prior, fit$ofv_prior)
  expect_equal(back2$ofv_data, fit$ofv_data)
})

test_that(".fitrx_recover_ofv_split leaves an unpriored or unusable wire alone", {
  rec <- ferx:::.fitrx_recover_ofv_split
  # Unpriored: ofv == aic - 2k exactly, so the prior half is reported as an
  # exact 0 rather than a rounding artefact.
  none <- function(ofv) list(ofv_data = ofv, ofv_prior = 0, recovered = FALSE)
  expect_equal(rec(list(ofv = -300, aic = -290, n_parameters = 5)), none(-300))
  # Priored: the AIC identity recovers the half.
  expect_equal(rec(list(ofv = -279.1978, aic = -266.1263, n_parameters = 7)),
               list(ofv_data = -280.1263, ofv_prior = 0.9285, recovered = TRUE))
  # A negative recovered half means the identity does not hold for this file;
  # a penalty is a sum of squares, so fall back rather than report nonsense.
  expect_equal(rec(list(ofv = -300, aic = -280, n_parameters = 5)), none(-300))
  # Missing or non-finite ingredients fall back too.
  expect_equal(rec(list(ofv = -300)), none(-300))
  expect_equal(rec(list(ofv = -300, aic = NA_real_, n_parameters = 5)), none(-300))
  expect_equal(rec(list(ofv = NA_real_, aic = -290, n_parameters = 5)),
               none(NA_real_))
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
