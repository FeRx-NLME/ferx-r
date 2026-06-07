library(ferx)

# Demonstrates ADDL-column multiple dosing.
# A single dose row (ADDL=6, II=24) expands to 7 daily doses at read time;
# TAD resets to 0 at each expanded dose automatically.
# The [derived] block shows trough accumulation (CMIN_TAU) and periodic AUC.

ex <- ferx_example("warfarin_addl")
ferx_model_show(ex$model)

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# --- Derived columns from sdtab -----------------------------------------
# Per-subject scalar derived quantities (same value broadcast to every row):
head(fit$sdtab[, c("ID", "TIME", "IPRED", "KE", "T_HALF", "CMAX", "CMIN_TAU", "AUC_TAU")])

# CMIN_TAU uses TAD < 1e-10 to locate the post-dose trough observation
# (TAD resets to 0 at each ADDL-expanded dose).  Compare troughs across
# subjects to see inter-individual variability in accumulation:
unique(fit$sdtab[, c("ID", "CMIN_TAU", "CMAX", "AUC_TAU")])

# DAY column: floor(TAFD / 24) + 1, so each observation carries its dosing day
table(fit$sdtab$DAY)
