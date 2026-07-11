library(ferx)

# One-compartment oral model with analytic inverse-Gaussian (Freijer & Post)
# absorption, via the closed-form `pk one_cpt_ig(...)` (ferx-core #790). It gives
# exact FOCE/FOCEI gradients that do not depend on any ODE-solver tolerance, and
# the same uniform pk-line interface as the analytic transit models. Reaches the
# same absorption density as igd_inverse_gaussian, which uses the ODE igd().
ex <- ferx_example("one_cpt_ig")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Inspect the model structure ferx parsed.
ferx_model_inspect(fit)
