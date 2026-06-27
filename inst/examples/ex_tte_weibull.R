library(ferx)

# Standalone time-to-event (TTE) model with a Weibull hazard, written in the
# compact TTE-only form: when [event_model] is the sole endpoint, the
# [structural_model], [error_model], and [individual_parameters] blocks are all
# omitted. (This compact syntax requires a recent ferx-r; see ex_tte_exponential.R
# for the placeholder-PK-block form that older builds need.)
#
# Event times on CMT 2 follow Weibull(scale, shape) with
#   scale = TVSCALE * exp(ETA_SCALE)   (characteristic time; larger = later events)
#   shape = TVSHAPE                     (> 1 = increasing hazard over time)
# DV codes the censoring status (0 = right-censored, 1 = exact event), with
# administrative censoring at t = 60 h.
#
# Reference (simulated): TVSCALE = 20 h, TVSHAPE = 1.5, omega(ETA_SCALE) = 0.04.
ex <- ferx_example("tte_weibull")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Model-predicted survival function S(t) = exp(-H(t)).
# `times` spans the 60 h observation window (censoring horizon).
surv <- ferx_predict_survival(ex$model, ex$data, times = seq(0, 60, by = 5), fit = fit)
head(surv)
