library(ferx)

# Demonstrates the [derived] and [output] blocks.
# [derived] adds computed columns to sdtab after the fit.
# [output] echoes individual PK parameters (CL, V, KA) that are not
# in the mandatory sdtab minimum but are useful for post-hoc analysis.

ex <- ferx_example("warfarin_derived")
ferx_model_show(ex$model)

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# All derived columns appear in sdtab:
names(fit$sdtab)
# Includes: KE, T_HALF, DAY, TAU_TIME, CMAX, TMAX, CTROUGH, CMAX_D1,
#           CMAX_D14, AUC_0_72, AUC_TAU, AUC_DV_72, CL, V, KA

# Per-row derived columns vary across time points (same value per subject):
head(fit$sdtab[, c("ID", "TIME", "KE", "T_HALF")], 12)

# Aggregate columns are constant per subject (same value on every row):
unique(fit$sdtab[, c("ID", "CMAX", "TMAX", "AUC_0_72")])

# AUC_TAU is a periodic (windowed) integral: one value per 24-h interval.
# It changes across intervals (different accumulated AUC per window).
fit$sdtab[fit$sdtab$ID == 1, c("TIME", "AUC_TAU")]

# CTROUGH uses TAD < 1e-10 to safely select the "trough" (dose time).
# Floating-point residuals from modular arithmetic prevent exact TAD == 0.
unique(fit$sdtab[, c("ID", "CTROUGH")])
