library(ferx)

# One-compartment model with Weibull absorption, via the built-in
# weibull(td, beta) input rate fed straight into the central compartment. Unlike a
# transit chain, weibull() models the entire absorption delay (no first-order ka);
# the dose feeds the Weibull density over time (not as a bolus). The shape `beta`
# selects the profile: > 1 a delayed peak, = 1 first-order (ka = 1/Td), < 1 fast
# early uptake.
#
# Paired with the ferx-core NONMEM absorption anchor data (igd_oral.csv). That data
# was simulated from a transit model, so the Weibull fit is mildly mis-specified --
# this is a syntax / workflow demo, not a parameter-recovery demo. See ferx-core
# nonmem_anchor/weibull_absorption.ctl for the NONMEM cross-check.
ex <- ferx_example("weibull_absorption")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the weibull() input rate.
ferx_model_inspect(fit)
