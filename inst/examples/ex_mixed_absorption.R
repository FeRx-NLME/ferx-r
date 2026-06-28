library(ferx)

# Mixed (zero-order + first-order) absorption via zero_order()/first_order()
# composition (ferx-core #505): a zero-order (constant-rate) input over a modeled
# duration DUR plus a first-order input at rate KA feed the central compartment,
# split by a declared dose fraction:
#
#   d/dt(central) = FZO1*first_order(ka=KA) + FZO*zero_order(dur=DUR) - CL/V*central
#
# The multiplier must be a single declared individual parameter (not an expression
# like (1-FZO)), so the first-order fraction is declared explicitly (FZO1 = 1 - FZO);
# the two pathways together deliver exactly the dose. The zero-order duration/fraction
# are differentiated by finite differences (the moving-boundary cutoff), the
# first-order pathway analytically.
#
# Paired with the ferx-core NONMEM absorption anchor data (igd_oral.csv). That data
# was simulated from a transit model, so the mixed fit is mildly mis-specified -- a
# syntax / workflow demo, not a parameter-recovery demo. The matched, recovering
# cross-check is ferx-core nonmem_anchor/mixed_zero_first.ctl (ferx FOCEI objective
# vs NONMEM #OBJV to ~1e-4).
ex <- ferx_example("mixed_absorption")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the fraction-weighted zero_order()
# and first_order() input-rate terms on the central compartment.
ferx_model_inspect(fit)
