library(ferx)

# Demonstrates [derived] for computing PD quantities post-fit.
# PD parameters (MEC, EMAX, EC50) are fixed at assumed values -- they
# do not enter the PK likelihood and cannot be estimated from PK data.
# The model is fit as a standard 1-cpt oral PK model; [derived] then
# computes EFF (Emax response), time-above-MEC (TAM_TAU, TAM_D1), and
# a TAD-gated Cmax (CMAX_INTERVAL) for each subject.

ex <- ferx_example("warfarin_derived_pkpd")
ferx_model_show(ex$model)

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# PK-only parameters are estimated; PD thetas are fixed (SE = 0):
fit$estimates

# PD derived quantities are in sdtab:
names(fit$sdtab)

# Per-subject Emax response and time-above-MEC:
unique(fit$sdtab[, c("ID", "CMAX_INTERVAL", "TAM_TAU", "TAM_D1")])

# EFF varies by time point (depends on IPRED at that time):
fit$sdtab[fit$sdtab$ID == 1, c("TIME", "IPRED", "EFF")]
