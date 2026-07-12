library(ferx)

# Two-compartment oral model with analytic Savic transit-compartment absorption,
# via the closed-form `pk two_cpt_transit(cl, v1, q, v2, n, mtt)` (ferx-core #634).
# The same Gamma(N+1, KTR) transit closed form as one_cpt_transit, superposed
# bi-exponentially onto a two-compartment disposition, with exact tolerance-free
# FOCE/FOCEI gradients and a continuous, estimable N (no ODE solve). The analytic
# counterpart to the ODE `transit_2cpt` example.
ex <- ferx_example("two_cpt_transit")

fit <- ferx_fit(ex$model, ex$data)

# The data is simulated from this same model, so the estimates recover the six
# simulation truths (CL 5, V1 50, Q 10, V2 100, MTT 1, N 3). As with two_cpt_ig,
# the 2-cpt analytic fit may print STATUS: NOT CONVERGED while sitting at the
# optimum -- a known artifact of the default outer optimizer (ferx-core #751),
# not a transit-specific issue.
print(fit)

# Inspect the model structure ferx parsed (the analytic two_cpt_transit model).
ferx_model_inspect(fit)
