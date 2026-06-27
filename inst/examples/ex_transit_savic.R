library(ferx)

# One-compartment oral ODE model with Savic transit-compartment absorption via
# the built-in `transit(n, mtt)` input rate. Unlike `ex_transit_2cpt.R` -- which
# hand-codes three explicit transit ODE states with a fixed integer count -- here
# the number of transit compartments N is a single continuous, estimable
# parameter (KTR = (N+1)/MTT). The oral dose feeds the transit() appearance rate
# over time and is not also added as a bolus.
#
# Paired with the transit-simulated transit_oral.csv anchor data, so this is a
# correct-specification parameter-recovery demo (true TVCL=5, TVV=50, TVKA=1.0,
# TVMTT=1.0, TVN=3).
ex <- ferx_example("transit_savic")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the model structure ferx parsed, including the transit() input rate.
ferx_model_inspect(fit)
