# Each standard analytical PK model ships with an ODE-form sibling that is its
# hand transcription (amount-based states, concentration read out via
# `[scaling] obs_scale = V`/`V1`). The two must produce identical population
# predictions. We predict both forms on the *analytical* model's dataset so the
# dosing/observation grid is identical, and compare PRED.
#
# We ALSO compare the fit objective (OFV). PRED agrees at eta=0 to RK45
# tolerance, but the FOCE OFV exercises the ODE solve at the conditional modes
# and amplifies solver error, so the ODE examples set a tight `ode_reltol` in
# their `[fit_options]` (see inst/examples/models/*_ode.ferx). The OFV check is
# a *tolerance band* (per-pair `ofv_band` below), not exact equality: on a
# finite-difference (no-autodiff) build the inner EBE optimisation leaves a
# residual that is not removed by tightening the solver (the ODE form can even
# converge to a slightly better OFV). Exact OFV equality needs exact gradients
# (an Enzyme/autodiff build); revisit the bands once that lands or the optimiser
# improves. The band still guards against gross regressions -- e.g. the ODE
# tolerance plumbing breaking, which would push two_cpt_oral_cov / three_cpt_iv
# back to 10-15 OFV units.
#
# The exhaustive PRED equivalence check across every dosing mode (bolus,
# infusion, multi-dose, steady state, lag, bioavailability) lives in ferx-core's
# tests/analytical_ode_equivalence.rs; this guards the *shipped* example files.

# analytical name -> ODE sibling, fit method, and OFV tolerance band.
ode_pairs <- list(
  list(an = "one_cpt_iv",       ode = "one_cpt_iv_ode",       method = "foce",  ofv_band = 0.5),
  list(an = "warfarin",         ode = "warfarin_ode",         method = "foce",  ofv_band = 0.5),
  list(an = "two_cpt_iv",       ode = "two_cpt_iv_ode",       method = "foce",  ofv_band = 0.5),
  list(an = "two_cpt_oral_cov", ode = "two_cpt_oral_cov_ode", method = "focei", ofv_band = 0.5),
  list(an = "three_cpt_iv",     ode = "three_cpt_iv_ode",     method = "foce",  ofv_band = 1.0),
  list(an = "three_cpt_oral",   ode = "three_cpt_oral_ode",   method = "foce",  ofv_band = 8.0),
  list(an = "bioavailability",  ode = "bioavailability_ode",  method = "focei", ofv_band = 0.5)
)

for (pair in ode_pairs) {
  analytical <- pair$an
  ode        <- pair$ode

  test_that(sprintf("ODE form '%s' matches analytical '%s' (PRED)", ode, analytical), {
    an_ex  <- ferx_example(analytical)
    ode_ex <- ferx_example(ode)

    # Use the analytical model's dataset for both so the comparison is
    # apples-to-apples (same doses, same observation times).
    data <- an_ex$data

    an_pred  <- ferx_predict(an_ex$model,  data)
    ode_pred <- ferx_predict(ode_ex$model, data)

    expect_equal(nrow(ode_pred), nrow(an_pred))
    expect_equal(ode_pred$TIME, an_pred$TIME)
    expect_true(all(is.finite(ode_pred$PRED)))

    # RK45 reproduces the analytical closed form to solver tolerance.
    expect_equal(ode_pred$PRED, an_pred$PRED, tolerance = 1e-3)
  })
}

# OFV equivalence (tolerance band). Evaluated at the shared initial parameters
# (settings = list(maxiter = 0)) so this is a deterministic, fit-free check of
# the likelihood surface; the ODE side uses the tight `ode_reltol` baked into
# its [fit_options]. Slower than the PRED check (the tight ODE solve runs the
# inner EBE optimisation), so it is gated behind a normal run.
ofv_at_init <- function(model_path, data_path, method) {
  # covariance = FALSE and maxiter = 0 intentionally override the model file's
  # [fit_options]; suppress the resulting (expected) override warnings.
  fit <- suppressWarnings(ferx_fit(
    model_path, data_path,
    method = method, covariance = FALSE, verbose = FALSE,
    settings = list(maxiter = 0)
  ))
  fit$ofv
}

for (pair in ode_pairs) {
  analytical <- pair$an
  ode        <- pair$ode
  method     <- pair$method
  band       <- pair$ofv_band

  test_that(sprintf("ODE form '%s' OFV is within %.2g of analytical '%s'",
                    ode, band, analytical), {
    skip_on_cran()
    an_ex  <- ferx_example(analytical)
    ode_ex <- ferx_example(ode)
    data   <- an_ex$data

    an_ofv  <- ofv_at_init(an_ex$model,  data, method)
    ode_ofv <- ofv_at_init(ode_ex$model, data, method)

    expect_true(is.finite(an_ofv) && is.finite(ode_ofv))
    expect_lt(abs(ode_ofv - an_ofv), band)
  })
}
