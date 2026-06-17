library(ferx)

# Two-compartment oral PK with covariates, written with `ode_template`:
# ferx generates the standard disposition ODE (depot/central/periph, the
# micro-constant RHS, and obs_scale = V1) from the named model, so you do not
# hand-write the [odes] block. It takes the same parameters as the analytical
# `pk two_cpt_oral(...)`, including `ka`.
ex <- ferx_example("two_cpt_oral_cov_ode_template")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# The generated disposition is identical to the analytical model and to the
# hand-written ODE sibling -- same predictions.
an  <- ferx_predict(ferx_example("two_cpt_oral_cov")$model, ex$data)
gen <- ferx_predict(ex$model, ex$data)
stopifnot(isTRUE(all.equal(gen$PRED, an$PRED, tolerance = 1e-3)))

# To add absorption, re-declare the depot in [odes] (overrides only that
# compartment); see ferx-core examples/transit_ode_template.ferx for the
# `transit(...)` override.
ferx_model_inspect(fit)
