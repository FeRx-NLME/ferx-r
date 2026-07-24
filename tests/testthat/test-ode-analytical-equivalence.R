# Each standard analytical PK model ships with an ODE-form sibling that is its
# hand transcription (amount-based states, concentration read out via
# `[scaling] obs_scale = V`/`V1`). The two must produce identical population
# predictions. We predict both forms on the *analytical* model's dataset so the
# dosing/observation grid is identical, and compare PRED.
#
# We ALSO compare the fit objective (OFV). PRED agrees at eta=0 to RK45
# tolerance, but the FOCE OFV exercises the ODE solve at the per-subject
# conditional modes, so the ODE examples set a tight `ode_reltol` in their
# `[fit_options]` (see inst/examples/models/*_ode.ferx). Both forms compute the
# same objective with *analytic* sensitivities -- closed-form for the analytical
# PK path, an RK45-integrated sensitivity walk for the ODE path (there is no
# finite-difference gradient here; the Enzyme/autodiff path was retired in
# ferx-core #381) -- so once the inner EBE solve lands on the same mode the two
# OFVs agree to solver round-off (~1e-6 here). The check is a per-pair `ofv_band`
# rather than exact equality only to tolerate that RK45 round-off; each band is
# set just above the observed residual (set FERX_OFV_MEASURE=1 to print it). It
# still guards against gross regressions -- e.g. the ODE tolerance plumbing
# breaking, which would push these pairs back to 10-15 OFV units.
#
# The exhaustive PRED equivalence check across every dosing mode (bolus,
# infusion, multi-dose, steady state, lag, bioavailability) lives in ferx-core's
# tests/analytical_ode_equivalence.rs; this guards the *shipped* example files.

# analytical name -> ODE sibling, fit method, and OFV tolerance band.
ode_pairs <- list(
  list(an = "one_cpt_iv",       ode = "one_cpt_iv_ode",       method = "foce",  ofv_band = 0.5),
  list(an = "warfarin",         ode = "warfarin_ode",         method = "foce",  ofv_band = 0.5),
  list(an = "two_cpt_iv",       ode = "two_cpt_iv_ode",       method = "foce",  ofv_band = 0.5),
  # two_cpt_oral_cov once showed a ~0.91 (Linux-release CI) init-point OFV gap:
  # the analytical and ODE inner-EBE solves could land on different modes of a
  # multimodal per-subject objective (the same failure as three_cpt_iv below).
  # ferx-core #378 makes the exact (analytical) path keep the better inner-EBE
  # estimate, matching the ODE path, so the gap collapses to solver round-off
  # (|diff| ~0 measured). Tightened 1.0 -> 0.5 to match the other pairs.
  list(an = "two_cpt_oral_cov", ode = "two_cpt_oral_cov_ode", method = "focei", ofv_band = 0.5),
  # three_cpt_iv's FOCE *marginal* OFV (Sheiner-Beal linearised marginal built on
  # the per-subject inner EBE eta-hat solve -- NOT a Laplace log|H| term) once
  # diverged from its analytical sibling by ~18 OFV units (Linux-release CI; ~0.13
  # macOS-dev) -- platform-sensitive because it hinged on which mode each form's
  # inner solve reached. Root cause (ferx-core #378): on a subject whose
  # individual objective is multimodal (id 14, 6 IIV etas / 10 obs) the closed-form
  # inner BFGS found the better mode but stalled just above the tightened inner_tol
  # (#330); the exact-path recovery then discarded that near-global partial for a
  # worse cold Nelder-Mead basin, while the ODE form kept the good mode. #378 keeps
  # the lower-objective of {BFGS partial, NM} on the non-FREM exact path (matching
  # the ODE path's #555 guard), so both forms condition to the same mode and the
  # marginals agree to round-off (|diff| ~1e-6 measured). Back to a tight guard
  # alongside the other IV pairs.
  list(an = "three_cpt_iv",     ode = "three_cpt_iv_ode",     method = "foce",  ofv_band = 0.5),
  # three_cpt_oral is the highest-dimensional pair (depot + central + 2
  # peripherals, plus KA). Its wide band was attributed to "FD inner-EBE
  # residual", but both forms use analytic sensitivities; the real cause was the
  # same multimodal inner-EBE mode split fixed in ferx-core #378. With that fix
  # the residual collapses to round-off (|diff| ~0 measured), so this tightens
  # 8.0 -> 0.5 like the rest.
  list(an = "three_cpt_oral",   ode = "three_cpt_oral_ode",   method = "foce",  ofv_band = 0.5),
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
# (settings = list(maxiter = 0)): no outer optimisation runs, so it is a
# deterministic check of the likelihood surface at theta0 -- though the FOCE OFV
# still solves the inner EBE modes per subject. That inner solve (at the tight
# `ode_reltol` baked into the ODE side's [fit_options]) is why it is slower than
# the PRED check and why the comparison is a band rather than exact equality.
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

    # Measurement scaffold: set FERX_OFV_MEASURE=1 to print the actual residual
    # per pair without failing. Use this after a ferx-core pin bump (e.g. the
    # inner-EBE accuracy work in ferx-core #337/#289/#354) to see how far each
    # band can be tightened, then set `ofv_band` just above the observed values.
    # Lines are tagged "[OFV-MEASURE]" so they are easy to grep out of test logs.
    if (nzchar(Sys.getenv("FERX_OFV_MEASURE"))) {
      message(sprintf(
        "[OFV-MEASURE] %-18s an=%.6f ode=%.6f |diff|=%.6f band=%.2g headroom=%.6f",
        ode, an_ofv, ode_ofv, abs(ode_ofv - an_ofv), band,
        band - abs(ode_ofv - an_ofv)
      ))
    }

    expect_true(is.finite(an_ofv) && is.finite(ode_ofv))
    expect_lt(abs(ode_ofv - an_ofv), band)
  })
}

# `ode_template NAME(...)` generates the standard disposition ODE from the named
# model (ferx-core #363). The generated form must predict identically to BOTH
# the analytical `pk two_cpt_oral` and the hand-written ODE sibling
# `two_cpt_oral_cov_ode` -- it desugars to exactly the latter. This guards the
# shipped `*_ode_template` example; the exhaustive per-dosing-mode check lives
# in ferx-core's tests/ode_template_equivalence.rs.
test_that("ode_template form matches analytical and hand-ODE (PRED)", {
  an_ex   <- ferx_example("two_cpt_oral_cov")
  ode_ex  <- ferx_example("two_cpt_oral_cov_ode")
  tmpl_ex <- ferx_example("two_cpt_oral_cov_ode_template")
  data    <- an_ex$data

  an_pred   <- ferx_predict(an_ex$model,  data)
  ode_pred  <- ferx_predict(ode_ex$model, data)
  tmpl_pred <- ferx_predict(tmpl_ex$model, data)

  # `ode_template` is ferx-core parser syntax (FeRx-NLME/ferx-core#363). On an
  # engine that predates it the model fails to parse ("No PK model found in
  # [structural_model] block"); ferx_predict() does NOT error in that case -- it
  # prints the parse error and returns a 0-row result -- so detect the empty
  # result and skip. The test activates automatically once src/rust/Cargo.lock
  # is bumped to a ferx-core that supports ode_template (after #363 merges).
  skip_if(is.null(tmpl_pred) || nrow(tmpl_pred) == 0L,
          "ode_template not supported by pinned ferx-core (model did not parse)")

  expect_equal(nrow(tmpl_pred), nrow(an_pred))
  expect_true(all(is.finite(tmpl_pred$PRED)))
  expect_equal(tmpl_pred$PRED, an_pred$PRED,  tolerance = 1e-3)
  expect_equal(tmpl_pred$PRED, ode_pred$PRED, tolerance = 1e-3)
})
