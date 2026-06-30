library(ferx)

# One-compartment model with Savic transit-compartment absorption as a fast
# ANALYTIC closed form -- `pk one_cpt_transit(cl, v, n, mtt)` (ferx-core #611).
# The dose's absorption time is Gamma(N+1, KTR), KTR = (N+1)/MTT, fed straight
# into central; the whole model is closed form (no ODE solve), with exact Dual2
# FOCE/FOCEI sensitivities and a continuous, estimable N. With N = 0 it reduces
# to first-order (Bateman) oral absorption.
#
# Contrast `ex_transit_savic.R`, which models the same Savic transit chain on the
# numerical ODE path AND adds an explicit first-order `ka` step after the chain.
# Here there is no `ka`: the Gamma absorption time IS the entire delay.
#
# Paired with the ferx-core transit anchor data (transit_oral.csv). Because that
# data was simulated from a transit + first-order-`ka` model (see
# ex_transit_savic.R), this analytic fit is mildly mis-specified -- a syntax /
# workflow demo, not a parameter-recovery demo: CL and V recover, while the Gamma
# absorption shape (MTT, N) absorbs the omitted `ka` step.
ex <- ferx_example("one_cpt_transit")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the model structure ferx parsed (the analytic one_cpt_transit model).
ferx_model_inspect(fit)
