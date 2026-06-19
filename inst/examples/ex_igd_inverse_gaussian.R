library(ferx)

# One-compartment model with Freijer & Post inverse-Gaussian absorption, via the
# built-in igd(mat, cv2) input rate fed straight into the central compartment.
# Unlike a transit chain, igd() models the entire absorption delay (no first-order
# ka); the dose feeds the IG density over time (not as a bolus).
#
# Paired with the ferx-core NONMEM igd anchor data (igd_oral.csv). That data was
# simulated from a transit model, so the IG fit is mildly mis-specified -- this is
# a syntax / workflow demo, not a parameter-recovery demo. See ferx-core
# nonmem_anchor/freijer_ig.ctl for the NONMEM cross-check.
ex <- ferx_example("igd_inverse_gaussian")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the igd() input rate.
ferx_model_inspect(fit)
