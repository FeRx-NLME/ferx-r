library(ferx)

# Joint PK-TTE: a time-to-event endpoint whose hazard depends on the drug
# concentration, estimated jointly with the PK. The hazard expression in the
# [event_model] block,
#   hazard = H0 * exp(BETA * (central / V))
# references the ODE PK state; ferx appends it to the ODE system as a
# cumulative-hazard state (d/dt(__chz) = hazard) and fits PK + TTE together by
# FOCEI, sharing the CL random effect (ferx-core #564).
#
# Layout: oral dose on CMT 1, PK observed on CMT 2 (central), the event on CMT 3
# (DV = 1 exact event, DV = 0 right-censored). Reference (simulated) values:
# TVCL = 1, TVV = 10, TVKA = 1, with a drug-driven hazard H0 * exp(BETA * Cc).
# H0 and BETA trade off along a flat ridge (a single dose gives a narrow exposure
# range), so the exposure-hazard *curve* is better identified than the pair; the
# PK block recovers cleanly. (The bundled dataset is a small syntax demo.)
ex <- ferx_example("pktte_joint")

fit <- ferx_fit(ex$model, ex$data, method = "focei")
print(fit)

# Survival predictions read S(t), H(t), h(t) off the integrated cumulative-hazard
# compartment for the event endpoint (CMT 3), on a user-supplied time grid.
surv <- ferx_predict_survival(ex$model, ex$data, times = seq(0, 24, by = 2), fit = fit)
head(surv)

# Simulate from the fit. A drug-driven hazard needs a finite `horizon`: each
# subject's event time is sampled off the integrated cumulative hazard, and a
# subject with no event by `horizon` is right-censored there. On the event rows
# (CMT 3) TIME is the sampled event or censoring time and OBSERVED is 1 for an
# event, 0 for censored.
sim <- ferx_simulate(ex$model, ex$data, n_sim = 20L, seed = 1L, fit = fit,
                     horizon = 24)
events <- sim[!is.na(sim$OBSERVED), ]
head(events)
mean(events$OBSERVED)  # share of simulated subjects (all 20 replicates) with an event by t = 24
