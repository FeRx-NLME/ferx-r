library(ferx)

# Per-route absorption lag: immediate-release + delayed-release via the optional
# per-route lag= argument (ferx-core #856). Two first_order() pathways feed the
# central compartment, split by a declared dose fraction, where the delayed-release
# pathway carries its OWN onset lag:
#
#   d/dt(central) = FR1*first_order(ka=KA1)                # IR: absorbs at once
#                 + FR2*first_order(ka=KA2, lag=LAG2)      # DR: onset delayed by LAG2
#                 - CL/V*central
#
# lag= is a universal optional argument on every input-rate function; it delays that
# one pathway on top of any compartment lagtime (a bare lagtime would shift BOTH
# routes together), so parallel / mixed routes can switch on at different times --
# the classic modified-release picture. The multiplier (FR*) and lag= must each be a
# single declared individual parameter, so the complement is declared explicitly
# (FR2 = 1 - FR1). A model carrying a per-route lag is fit over finite differences
# (predictions are exact); the analytic per-route onset saltation is a ferx-core
# follow-up (#859).
#
# Paired with the ferx-core NONMEM absorption anchor data (igd_oral.csv). That data
# was simulated from a transit model, so this fit is mildly mis-specified -- a syntax
# / workflow demo, not a parameter-recovery demo. The matched, recovering cross-check
# is ferx-core nonmem_anchor/per_route_lag.ctl (ferx FOCEI objective vs NONMEM
# #OBJV = -882.357 to ~1e-6).
ex <- ferx_example("per_route_lag_absorption")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the ODE structure ferx parsed, including the two fraction-weighted
# first_order() input-rate terms on central and the per-route lag on the second.
ferx_model_inspect(fit)
