## Oral bioavailability with logit-normal F (analytical one-compartment)
#
# Bioavailability F is constrained to (0,1) via logit transform.
# THETA_F is specified directly on the (0,1) scale (e.g. 0.70 = 70%);
# ferx applies logit(THETA_F) internally so the optimiser works on an
# unbounded scale while THETA_F stays interpretable.
#
# True simulation parameters:
#   TVCL=5, TVV=50, TVKA=1.5, F=0.70, omega_CL=0.09, omega_F=0.10,
#   sigma_prop=0.15

library(ferx)

ex <- ferx_example("bioavailability")

fit <- ferx_fit(
  model      = ex$model,
  data       = ex$data,
  method     = "focei",
  covariance = TRUE
)

print(fit)
