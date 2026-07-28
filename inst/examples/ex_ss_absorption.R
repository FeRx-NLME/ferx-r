library(ferx)

# Steady-state dosing (SS=1) into a built-in absorption compartment
# (ferx-core #719, gap 1). A one-compartment oral model whose absorption is the
# built-in first_order(ka) input-rate kernel; the dose is a NONMEM steady-state
# record (SS=1, II=8), so the engine equilibrates the periodic dosing through
# the absorption kernel instead of needing an explicit run-in. See
# ?ferx_example and the model header for details.
ex   <- ferx_example("ss_absorption")
data <- read.csv(ex$data)

# Population prediction at the (NONMEM-matched) typical values. The dataset's DV
# column carries the NONMEM 7.6.0 population prediction (ADVAN2 TRANS2), so
# ferx_predict() should reproduce it closely -- a cross-engine check that the
# steady-state absorption equilibration is correct.
pred <- ferx_predict(ex$model, ex$data)

obs <- data[data$EVID == 0, c("TIME", "DV")]
obs$DV <- as.numeric(obs$DV)   # DV is character: dose rows use "." as a placeholder
cmp <- merge(obs, pred[, c("TIME", "PRED")], by = "TIME")
cmp$rel_err <- abs(cmp$PRED - cmp$DV) / cmp$DV
print(cmp)
cat(sprintf("max relative error vs NONMEM PRED: %.2e\n", max(cmp$rel_err)))

# A fit works the same way -- raise omega ETA_CL and add subjects to estimate
# from steady-state data:
#   fit <- ferx_fit(ex$model, ex$data)
