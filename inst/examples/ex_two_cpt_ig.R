library(ferx)

# Two-compartment oral model with analytic inverse-Gaussian (Freijer & Post)
# absorption, via the closed-form `pk two_cpt_ig(...)` (ferx-core #790). The same
# exponential-tilting IG closed form as one_cpt_ig, convolved onto a
# two-compartment disposition, with exact tolerance-free FOCE/FOCEI gradients.
ex <- ferx_example("two_cpt_ig")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the model structure ferx parsed.
ferx_model_inspect(fit)
