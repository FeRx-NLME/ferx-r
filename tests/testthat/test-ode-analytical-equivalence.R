# Each standard analytical PK model ships with an ODE-form sibling that is its
# hand transcription (amount-based states, concentration read out via
# `[scaling] obs_scale = V`/`V1`). The two must produce identical population
# predictions. We predict both forms on the *analytical* model's dataset so the
# dosing/observation grid is identical, and compare PRED. The ODE path uses the
# RK45 solver (reltol 1e-4), so agreement is asserted to solver tolerance.
#
# The exhaustive equivalence check across every dosing mode (bolus, infusion,
# multi-dose, steady state, lag, bioavailability) lives in ferx-core's
# tests/analytical_ode_equivalence.rs; this guards the *shipped* example files.

# analytical example name -> its ODE-form sibling
ode_pairs <- list(
  c("one_cpt_iv",       "one_cpt_iv_ode"),
  c("warfarin",         "warfarin_ode"),
  c("two_cpt_iv",       "two_cpt_iv_ode"),
  c("two_cpt_oral_cov", "two_cpt_oral_cov_ode"),
  c("three_cpt_iv",     "three_cpt_iv_ode"),
  c("three_cpt_oral",   "three_cpt_oral_ode"),
  c("bioavailability",  "bioavailability_ode")
)

for (pair in ode_pairs) {
  analytical <- pair[[1]]
  ode        <- pair[[2]]

  test_that(sprintf("ODE form '%s' matches analytical '%s'", ode, analytical), {
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
