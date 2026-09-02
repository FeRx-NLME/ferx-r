library(ferx)

# Adaptive (feedback) dosing with a PRE-SCHEDULED LOADING DOSE (ferx-core #702 /
# #929). Unlike the dose-free adaptive_tdm example, here the subject data carries
# a 1500 mg vancomycin loading dose (an EVID = 1 row at t = 0). The model's
# [adaptive_dosing] block does not replace it -- it AUGMENTS it, titrating the
# daily q24h maintenance dose on the measured trough to hold it in a 10-15 mg/L
# window. The day-1 decision therefore sees the decayed loading dose (~5.6 mg/L),
# not 0. See ?ferx_simulate_adaptive.
ex <- ferx_example("adaptive_vanco_loading")

res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 10L, seed = 1L)

# Four tables: concentration trajectories, the realized dose ledger, the
# per-decision log, and the per-subject outcome metrics.
str(res, max.level = 1)

# The trajectory for the first subject/replicate: concentration starts high from
# the loading dose and is held in-window as the maintenance dose is titrated.
head(res$trajectories[res$trajectories$SIM == 1, c("TIME", "IPRED")], 8)

# The controller's decisions: the trough it observed each day and the rule fired.
res$decisions[res$decisions$SIM == 1,
              c("DECISION", "TIME", "SIGNAL", "OUTCOME")]

# IMPORTANT ledger caveat: the dose ledger holds CONTROLLER-issued doses only.
# The 1500 mg loading dose is a base (pre-scheduled) dose, so it is absent from
# res$doses and from metrics$CUM_DOSE. To get the full realized regimen, add the
# base dose(s) from the input data back in.
loading_dose <- 1500
cum_controller <- res$metrics$CUM_DOSE[res$metrics$SIM == 1]
cat(sprintf(
  "Controller-only cumulative dose: %.0f mg\nFull realized (incl. loading): %.0f mg\n",
  cum_controller, cum_controller + loading_dose
))

# The loading dose does not appear in the ledger:
cat(sprintf("Loading dose (1500 mg) present in $doses? %s\n",
            any(res$doses$AMT == loading_dose)))

# Attainment DOES account for the base regimen (target_window = [10, 15]):
cat(sprintf("Mean time in 10-15 mg/L window: %.1f%%\n",
            100 * mean(res$metrics$PCT_TIME_IN_WINDOW, na.rm = TRUE)))
