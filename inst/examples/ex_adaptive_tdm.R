library(ferx)

# Adaptive (feedback) dosing -- a vancomycin-style TDM trough titration (epic
# #391). The dosing regimen is NOT in the data: the model file's
# [adaptive_dosing] block titrates the next dose every 12 h to drive the measured
# (assay-noised) pre-dose trough into a 10-20 mg/L window. The base subjects are
# dose-free; the controller supplies every dose at run time. See
# ?ferx_simulate_adaptive.
ex <- ferx_example("adaptive_tdm")

res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 20L, seed = 1L)

# Four tables: concentration trajectories, the realized dose ledger, the
# per-decision log (including holds), and the per-subject outcome metrics.
str(res, max.level = 1)

# For the first subject/replicate: the trough the controller observed at each
# decision (SIGNAL) and the rule it fired, plus the realized doses.
head(res$decisions[res$decisions$SIM == 1 & res$decisions$ID == "1", ])
head(res$doses[res$doses$SIM == 1 & res$doses$ID == "1", ])

# Median observed trough by decision across replicates -- it should climb from
# sub-therapeutic into the 10-20 mg/L target band as the dose is titrated up.
agg <- aggregate(SIGNAL ~ DECISION, data = res$decisions, FUN = median)
print(agg)

# Per-subject outcome metrics (one row per subject x replicate): cumulative
# dose, realized dose-change counts, and -- because this model declares
# target_window = [10, 20] -- the fraction of monitored troughs that landed in
# the 10-20 mg/L band (PCT_TIME_IN_WINDOW).
head(res$metrics[, c("ID", "SIM", "CUM_DOSE", "N_INCREASES", "N_DECREASES",
                     "PCT_TIME_IN_WINDOW")])

# Mean attainment across all simulated subjects/replicates.
cat(sprintf("Mean time in 10-20 mg/L window: %.1f%%\n",
            100 * mean(res$metrics$PCT_TIME_IN_WINDOW, na.rm = TRUE)))
