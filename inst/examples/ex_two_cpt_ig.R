library(ferx)

# Two-compartment oral model with analytic inverse-Gaussian (Freijer & Post)
# absorption, via the closed-form `pk two_cpt_ig(...)` (ferx-core #790). The same
# exponential-tilting IG closed form as one_cpt_ig, convolved onto a
# two-compartment disposition, with exact tolerance-free FOCE/FOCEI gradients.
ex <- ferx_example("two_cpt_ig")

fit <- ferx_fit(ex$model, ex$data)

# Note: this fit prints STATUS: NOT CONVERGED, yet the estimates sit at the
# optimum -- they recover the six simulation truths (CL 5, V1 40, Q 10, V2 60,
# MAT 2, CV2 0.3) with tight standard errors. The flag is a known artifact of the
# default outer optimizer (ferx-core #751), not an IG-specific issue.
print(fit)

# Inspect the model structure ferx parsed.
ferx_model_inspect(fit)
