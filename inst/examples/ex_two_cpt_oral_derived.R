library(ferx)

# Two-compartment oral PK model with [derived] columns.
# Computes micro-rate constants (K10, K12, K21), subject-level CMAX
# and TMAX, a week-1 AUC (AUC_0_168), and a periodic 24-h AUC (AUC_TAU).
# The data alias is two_cpt_oral_cov (WT and CRCL covariates), so the
# dataset has covariate columns that this model does not use.

ex <- ferx_example("two_cpt_oral_derived")
ferx_model_show(ex$model)

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Derived micro-rate constants and PK summary stats in sdtab:
names(fit$sdtab)

# K10, K12, K21 are per-row (same value per subject; depend on EBE CL/V):
head(unique(fit$sdtab[, c("ID", "K10", "K12", "K21")]))

# CMAX, TMAX are aggregate (one value per subject):
head(unique(fit$sdtab[, c("ID", "CMAX", "TMAX", "AUC_0_168")]))

# AUC_TAU is periodic: one value per 24-h window
fit$sdtab[fit$sdtab$ID == 1, c("TIME", "TAFD", "AUC_TAU")]
