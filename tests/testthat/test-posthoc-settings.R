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
