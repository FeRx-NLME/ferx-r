## Oral bioavailability with logit-normal F (ODE implementation)
#
# ODE equivalent of ex10_bioavailability.R — same model, same data, same true
# parameters. Useful for verifying that the ODE and analytical paths agree.
#
# The depot compartment carries oral dose in mg (amount); the central
# compartment is in concentration (mg/L). Bioavailability F enters the
# depot→central transfer: d/dt(central) = F * KA * depot / V - CL/V * central.

library(ferx)

ex <- ferx_example("bioavailability_ode")

fit <- ferx_fit(
  model      = ex$model,
  data       = ex$data,
  method     = "focei",
  covariance = TRUE
)

print(fit)
