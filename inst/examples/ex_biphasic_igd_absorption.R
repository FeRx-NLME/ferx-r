library(ferx)

# Biphasic inverse-Gaussian (Freijer & Post) absorption, via the shared
# input-rate pathway-fraction multiplier (ferx-core #388): two igd() pathways --
# a fast (short MAT) and a slow (long MAT) one -- feed the central compartment,
# split by a declared dose fraction:
#
#   d/dt(central) = FR1*igd(mat=MAT1, cv2=CV2_1) + FR2*igd(mat=MAT2, cv2=CV2_2) - CL/V*central
#
# The multiplier must be a single declared individual parameter (not an expression
# like (1-FR)), so the complement is declared explicitly (FR2 = 1 - FR1); the two
# pathways together deliver exactly the dose.
#
# Paired with the ferx-core NONMEM absorption anchor data (igd_oral.csv). That data
# was simulated from a transit model, so the biphasic fit is mildly mis-specified --
# this is a syntax / workflow demo, not a parameter-recovery demo. The matched,
# recovering cross-check is ferx-core nonmem_anchor/freijer_biphasic_ig.ctl (ferx
# FOCEI objective vs NONMEM #OBJV to ~1e-5).
ex <- ferx_example("biphasic_igd_absorption")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the two fraction-weighted igd()
# input-rate terms on the central compartment.
ferx_model_inspect(fit)
