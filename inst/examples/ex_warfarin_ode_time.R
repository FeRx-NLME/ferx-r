library(ferx)

# Demonstrates TIME, T, TAFD, and TAD available inside [odes] RHS
# expressions. The model uses a standard 1-cpt oral warfarin ODE but
# adds two accumulator compartments:
#   AUC_D1  -- accumulates central/V only in the first 24 h (TIME < 24)
#   TAM_INT -- accumulates time when central/V > 0.5 AND TAD < 24 h
#
# [derived] then computes KE, T_HALF, CMAX, and AUC_72 using the
# standard post-fit expressions (not the ODE accumulators).
# [scaling] obs_scale = V maps the central amount to a concentration.

ex <- ferx_example("warfarin_ode_time")
ferx_model_show(ex$model)

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Derived columns in sdtab:
names(fit$sdtab)

# KE, T_HALF vary per subject (proportional to EBE CL/V):
unique(fit$sdtab[, c("ID", "KE", "T_HALF")])

# CMAX is the subject-level peak IPRED:
unique(fit$sdtab[, c("ID", "CMAX")])

# AUC_72 from [derived] integral -- uses IPRED at obs times with step=0.2h:
unique(fit$sdtab[, c("ID", "AUC_72")])
