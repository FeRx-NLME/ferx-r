library(ferx)

# Competing-risks time-to-event (TTE): each subject is at risk of two distinct
# event types, modelled as separate cause-specific exponential hazards on their
# own compartments via two named [event_model] blocks:
#   CMT 2 -- cause A, rate LAMBDA_A = TVLAMBDA_A * exp(ETA_F)
#   CMT 3 -- cause B, rate LAMBDA_B = TVLAMBDA_B * exp(ETA_F)
# A shared frailty ETA_F links the two hazards. The observed event is whichever
# cause occurs first; the subject is right-censored for the other cause at that
# same time. Data: one row per cause CMT per subject (DV=1 observed, DV=0 censored).
#
# Reference (simulated): TVLAMBDA_A = 0.10 /h, TVLAMBDA_B = 0.06 /h,
# omega(ETA_F) = 0.25, administrative censoring at t = 14 h. The cause-specific
# rates recover tightly; the shared frailty variance is weakly identified under
# FOCEI-Laplace for TTE (prefer SAEM/IMP for the variance component).
ex <- ferx_example("tte_competing_risks")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Per-cause cumulative incidence (cif) alongside all-cause survival (survival_all);
# by construction sum_k cif_k(t) + survival_all(t) = 1.
surv <- ferx_predict_survival(ex$model, ex$data, times = seq(0, 14, by = 1), fit = fit)
head(surv)
