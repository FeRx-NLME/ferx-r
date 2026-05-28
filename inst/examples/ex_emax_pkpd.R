library(ferx)

# Simultaneous PK/PD: oral 1-cpt PK + biophase Emax PD, fit jointly with a
# SEPARATE residual error model per endpoint via the per-CMT [error_model]
# block (CMT=2 proportional, CMT=3 additive). Multi-endpoint fitting is
# ODE-only and requires gradient = fd (set in the .ferx file).
ex <- ferx_example("emax_pkpd")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Per-endpoint residual labels carry through to fit$sigma_types and the
# residual summary reported by ferx_model_inspect().
fit$sigma
fit$sigma_types
ferx_model_inspect(fit)$residual  # "per-CMT (CMT2=proportional, CMT3=additive)"
