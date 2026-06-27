library(ferx)

# Standalone time-to-event (TTE) model with an exponential hazard, registered via
# the `[event_model]` block. Observations on CMT 2 are routed to a parametric
# survival likelihood instead of the Gaussian residual model. Event times follow
# Exp(LAMBDA) with LAMBDA = TVLAMBDA * exp(ETA_LAMBDA); DV codes the censoring
# status (0 = right-censored, 1 = exact event), with administrative censoring at
# t = 30 h.
#
# Reference (simulated): TVLAMBDA = 0.05 /h, omega(ETA_LAMBDA) = 0.09.
ex <- ferx_example("tte_exponential")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Model-predicted survival function S(t) = exp(-H(t)); for competing-risks models
# this also returns the cause-specific cumulative incidence (cif, survival_all).
# `times` spans the 30 h observation window (administrative censoring horizon).
surv <- ferx_predict_survival(ex$model, ex$data, times = seq(0, 30, by = 2), fit = fit)
head(surv)
