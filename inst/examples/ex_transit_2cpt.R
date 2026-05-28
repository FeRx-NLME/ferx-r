library(ferx)

# Two-compartment ODE model with 3-transit-compartment absorption,
# allometric scaling on CL/V1/V2 (reference WT = 70 kg), and a lag time on
# the first transit compartment. KTR = (n+1) / MTT = 4 / MTT (n = 3 transits).
ex <- ferx_example("transit_2cpt")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the model structure ferx parsed from the ODE block.
ferx_model_inspect(fit)
