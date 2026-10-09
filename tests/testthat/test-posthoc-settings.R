# ferx_covariance(fit) and ferx_sir(fit) repeat the in-fit steps (#511, #472).
#
# A fit records how it was scored: `fit$scoring_settings` (the inner-loop and
# ODE settings of the stage that produced the estimates, ferx-core #1805),
# `fit$sir_settings` (what its SIR ran under, ferx-core #1758) and, in memory,
# `fit$packed_estimate` (the optimizer's exact packed vector). The standalone
# steps hand all three back, so with no arguments they are identical() to the
# in-fit steps, and an explicit argument overrides by editing the record.
#
# Every row asserts its premise first - dropping the field it is about moves the
# result on that fixture - so a setting that does not move its fixture cannot
# pass by accident (#511's `inner_tol = 1e-5` is the engine default and
# `mu_referencing = FALSE` is inert on warfarin; both rows moved, plan section 0b).
#
# Engine: the R glue is built FERX_NO_AUTODIFF=1 in the measurements quoted
# here; ferx-core's inner gradients on these fixtures are analytic either way.

ps_scratch <- function(ext, env = parent.frame()) {
  withr::local_tempfile(fileext = ext, .local_envir = env)
}

# SIR fits draw 300, keep 100, seed 1 - the in-fit run each S4/S5 row repeats.
ps_sir_settings <- list(sir_samples = 300L, sir_resamples = 100L, sir_seed = 1L)

# Fits are shared between tests, so each is made once per file run.
ps_cache <- new.env(parent = emptyenv())
ps_fit <- function(key, example, settings = list(), sir = FALSE, ...) {
  if (is.null(ps_cache[[key]])) {
    ex <- if (is.list(example)) example else ferx_example(example)
    if (sir) settings <- c(settings, ps_sir_settings)
    ps_cache[[key]] <- suppressWarnings(ferx_fit(
      ex$model, ex$data, covariance = TRUE, verbose = FALSE, sir = sir,
      settings = if (length(settings)) settings else NULL, ...
    ))
  }
  fit <- ps_cache[[key]]
  skip_if(is.null(fit$cov_matrix), paste(key, "covariance step did not converge"))
  fit
}

ps_round_trip <- function(fit, env = parent.frame()) {
  path <- ps_scratch(".fitrx", env)
  unlink(path)
  ferx_save_fit(fit, path)
  ferx_load_fit(path)
}

ps_sir_fields <- function(fit) {
  fit[c("sir_ess", "sir_ci_theta", "sir_ci_omega", "sir_ci_sigma", "sir_ci_kappa")]
}

ps_scoring_keys <- c(
  "inner_maxiter", "inner_tol", "inner_restarts", "mu_referencing", "n_agq",
  "inner_optimizer", "ebe_warm_start", "ode_reltol", "ode_abstol",
  "ode_max_steps", "ode_method", "ode_stiff_abort_after", "ode_auto_switch"
)

test_that("S11: a fit shows its scoring record, and its SIR record only after SIR", {
  fit <- ps_fit("w_default", "warfarin")
  expect_identical(names(fit$scoring_settings), ps_scoring_keys)
  expect_identical(fit$scoring_settings$inner_optimizer, "auto")
  expect_true(is.numeric(fit$packed_estimate) && length(fit$packed_estimate) > 0L)
  expect_null(fit$sir_settings)

  sir <- ps_fit("w_sir_df3", "warfarin", list(sir_df = 3), sir = TRUE)
  expect_identical(
    names(sir$sir_settings),
    c("samples", "resamples", "seed", "df", "scale", "keep_samples", "scoring")
  )
  expect_identical(names(sir$sir_settings$scoring), ps_scoring_keys)
  expect_identical(sir$sir_settings[c("samples", "resamples", "seed", "df")],
                   list(samples = 300L, resamples = 100L, seed = 1, df = 3))
})

test_that("S1: on a default fit ferx_covariance(fit) is identical() to the in-fit step", {
  # The original 1e-7 case (#511): warfarin_iov SE(TVKA) 0.76168056 vs
  # 0.76168066 was the packed estimate, not the settings (plan section 0a).
  fit <- ps_fit("iov_default", "warfarin_iov")
  no_packed <- fit
  no_packed$packed_estimate <- NULL
  expect_false(identical(ferx_covariance(no_packed)$cov_matrix, fit$cov_matrix))

  cov <- ferx_covariance(fit)
  expect_identical(cov$cov_matrix, fit$cov_matrix)
  expect_identical(cov$se_theta, fit$se_theta)
})

# Fits made without SIR, so `run_covariance`'s fallback to the SIR record's
# scoring half cannot stand in for the record under test.
ps_cov_rows <- list(
  list(key = "w_lbfgs", example = "warfarin", setting = "inner_optimizer",
       settings = list(inner_optimizer = "lbfgs"), value = "lbfgs"),
  list(key = "w_tol3", example = "warfarin", setting = "inner_tol",
       settings = list(inner_tol = 1e-3), value = 1e-3),
  list(key = "iov_mu", example = "warfarin_iov", setting = "mu_referencing",
       args = list(mu_referencing = FALSE), value = FALSE),
  list(key = "ode_rtol3", example = "warfarin_ode", setting = "ode_reltol",
       settings = list(ode_reltol = 1e-3), value = 1e-3)
)

ps_cov_fit <- function(row) {
  do.call(ps_fit, c(list(row$key, row$example, row$settings %||% list()), row$args))
}

test_that("S2: ferx_covariance(fit) follows the fit's scoring record", {
  for (row in ps_cov_rows) {
    fit <- ps_cov_fit(row)
    # The record holds the setting as the stage ran it.
    expect_identical(fit$scoring_settings[[row$setting]], row$value, info = row$key)
    # Premise: on this fixture the setting moves the covariance.
    no_record <- fit
    no_record$scoring_settings <- NULL
    expect_false(identical(ferx_covariance(no_record)$cov_matrix, fit$cov_matrix),
                 info = row$key)
    cov <- ferx_covariance(fit)
    expect_identical(cov$cov_matrix, fit$cov_matrix, info = row$key)
    expect_identical(cov$se_theta, fit$se_theta, info = row$key)
  }
})

# In-fit SIR at 300 / 100 / seed 1, each under one non-default setting.
ps_sir_rows <- list(
  list(key = "w_sir_df3", example = "warfarin", settings = list(sir_df = 3)),
  list(key = "w_sir_maxit5", example = "warfarin", settings = list(inner_maxiter = 5L)),
  list(key = "iov_sir_natural", example = "warfarin_iov",
       settings = list(sir_scale = "natural")),
  list(key = "ode_sir_rtol3", example = "warfarin_ode",
       settings = list(ode_reltol = 1e-3))
)

ps_sir_fit <- function(row) {
  fit <- ps_fit(row$key, row$example, row$settings, sir = TRUE)
  skip_if(is.null(fit$sir_ess), paste(row$key, "in-fit SIR did not run"))
  fit
}

test_that("S4: ferx_sir(fit) with no arguments repeats the in-fit SIR", {
  for (row in ps_sir_rows) {
    fit <- ps_sir_fit(row)
    # Premise: the same draws without the record score differently.
    no_record <- fit
    no_record$sir_settings <- NULL
    expect_false(identical(
      ps_sir_fields(ferx_sir(no_record, 300L, 100L, sir_seed = 1L)),
      ps_sir_fields(fit)
    ), info = row$key)
    expect_identical(ps_sir_fields(ferx_sir(fit)), ps_sir_fields(fit), info = row$key)
  }
})

test_that("S6: an explicit argument overrides the fit's record, even at its default value", {
  # The engine resolves a record by value, so an explicit default passed as an
  # option would lose to it (plan section 0c); the argument edits the record.
  fit <- ps_cov_fit(ps_cov_rows[[3L]])  # warfarin_iov, recorded mu_referencing = FALSE
  no_record <- fit
  no_record$scoring_settings <- NULL
  over <- ferx_covariance(fit, mu_referencing = TRUE)$cov_matrix
  expect_identical(over, ferx_covariance(no_record)$cov_matrix)
  expect_false(identical(over, fit$cov_matrix))
  expect_identical(ferx_covariance(fit, mu_referencing = FALSE)$cov_matrix,
                   fit$cov_matrix)

  nat <- ps_sir_fit(ps_sir_rows[[3L]])  # warfarin_iov, recorded sir_scale = "natural"
  packed <- ferx_sir(nat, sir_scale = "packed")
  ref <- nat
  ref$sir_settings <- NULL
  expect_identical(ps_sir_fields(packed),
                   ps_sir_fields(ferx_sir(ref, 300L, 100L, sir_seed = 1L)))
  expect_false(identical(ps_sir_fields(packed), ps_sir_fields(nat)))
  # The run records what it did; the settings not passed stay the fit's.
  expect_identical(packed$sir_settings$scale, "packed")
  expect_identical(packed$sir_settings[c("samples", "resamples", "seed")],
                   nat$sir_settings[c("samples", "resamples", "seed")])
})

test_that("S7: a fit with no record is scored as it always was, and says nothing", {
  bare <- function(fit) {
    fit$scoring_settings <- NULL
    fit$sir_settings <- NULL
    fit$packed_estimate <- NULL
    fit
  }
  f0 <- bare(ps_cov_fit(ps_cov_rows[[3L]]))
  expect_no_warning(cov <- ferx_covariance(f0))
  # An explicit argument on a record-less fit edits the engine's default record:
  # at the default value it is the same run.
  expect_identical(cov$cov_matrix,
                   ferx_covariance(f0, mu_referencing = TRUE)$cov_matrix)

  s0 <- bare(ps_sir_fit(ps_sir_rows[[1L]]))
  expect_no_warning(sir <- ferx_sir(s0))
  expect_identical(sir$sir_settings[c("samples", "resamples", "seed", "df", "scale")],
                   list(samples = 1000L, resamples = 250L, seed = 12345, df = 5,
                        scale = "packed"))
  expect_identical(ps_sir_fields(sir),
                   ps_sir_fields(ferx_sir(s0, 1000L, 250L, sir_scale = "packed")))
})

test_that("S8: a fit whose estimates were edited does not keep its packed estimate", {
  fit <- ps_fit("w_default", "warfarin")
  edited <- fit
  edited$theta[1L] <- edited$theta[1L] * 1.01
  edited_np <- edited
  edited_np$packed_estimate <- NULL
  # With the packed vector the step would differentiate at the unedited point,
  # which is the in-fit covariance.
  cov <- ferx_covariance(edited)$cov_matrix
  expect_identical(cov, ferx_covariance(edited_np)$cov_matrix)
  expect_false(identical(cov, fit$cov_matrix))
})

test_that("S9: a [mixture] per-class override fit is scored in memory, refused after a reload", {
  ex <- list(model = test_path("fixtures", "mixture_override.ferx"),
             data = test_path("fixtures", "mixture_iv.csv"))
  fit <- ps_fit("mix_override", ex, method = "focei")
  refusal <- "does not store their fitted values"
  no_packed <- fit
  no_packed$packed_estimate <- NULL
  expect_error(ferx_covariance(no_packed), refusal, fixed = TRUE)

  cov <- ferx_covariance(fit)
  expect_identical(cov$cov_matrix, fit$cov_matrix)
  expect_identical(cov$se_theta, fit$se_theta)

  # .fitrx has no key for the packed estimate (ferx-core#1815).
  loaded <- ps_round_trip(fit)
  expect_null(loaded$packed_estimate)
  expect_error(ferx_covariance(loaded), refusal, fixed = TRUE)
})
