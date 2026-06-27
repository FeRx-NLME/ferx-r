library(ferx)

# Standalone time-to-event (TTE) model with a Gompertz hazard via the
# [event_model] block, in the compact TTE-only form. The Gompertz hazard
#   h(t) = alpha * exp(gamma * t)
# has an exponentially increasing baseline. Between-subject variability is placed
# on the growth rate gamma (not the intercept alpha) because the two are collinear
# in the population likelihood; BSV on gamma gives clean recovery.
#
# DV codes the censoring status (0 = right-censored, 1 = exact event), with
# administrative censoring at t = 80 h.
#
# Reference (simulated): TVALPHA = 0.002 /h, TVGAMMA = 0.05 /h, omega(ETA_GAMMA) = 0.04.
ex <- ferx_example("tte_gompertz")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Model-predicted survival S(t); `times` spans the 80 h observation window.
surv <- ferx_predict_survival(ex$model, ex$data, times = seq(0, 80, by = 5), fit = fit)
head(surv)
