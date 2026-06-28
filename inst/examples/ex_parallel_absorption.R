library(ferx)

# Parallel (dual first-order) absorption via first_order() composition
# (ferx-core #505): two first_order() pathways -- a fast (large KA1) and a slow
# (small KA2) one -- feed the central compartment, split by a declared dose
# fraction:
#
#   d/dt(central) = FR1*first_order(ka=KA1) + FR2*first_order(ka=KA2) - CL/V*central
#
# The multiplier must be a single declared individual parameter (not an expression
# like (1-FR)), so the complement is declared explicitly (FR2 = 1 - FR1); the two
# pathways together deliver exactly the dose, and the model keeps exact analytic
# FOCEI gradients (both pathways are smooth exp-only densities).
#
# Paired with the ferx-core NONMEM absorption anchor data (igd_oral.csv). That data
# was simulated from a transit model, so the parallel fit is mildly mis-specified --
# a syntax / workflow demo, not a parameter-recovery demo. The matched, recovering
# cross-check is ferx-core nonmem_anchor/parallel_first_order.ctl (ferx FOCEI
# objective vs NONMEM #OBJV to ~1e-5).
ex <- ferx_example("parallel_absorption")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the two fraction-weighted
# first_order() input-rate terms on the central compartment.
ferx_model_inspect(fit)
