library(ferx)

# Analytic Savic transit-compartment absorption WITH inter-occasion variability
# (IOV) -- `pk one_cpt_transit(cl, v, n, mtt)` + a `kappa` random effect on CL
# (ferx-core #719, #1560).
#
# The analytic absorption closed forms (one_cpt_transit / two_cpt_transit /
# one_cpt_ig / two_cpt_ig) used to REJECT IOV: a per-dose closed-form
# superposition cannot carry a dose's drug across an occasion boundary where the
# disposition changes. ferx-core #719 rerouted a subject with IOV to the model's
# exact `transit()` ODE twin (d/dt(central) = transit(n, mtt) - (CL/V) * central).
# Since ferx-core #1560 the subject stays analytic: an exact event-driven walk
# carries the delivered amount across each record where CL changes, with
# closed-form predictions and FOCE/FOCEI gradients and no ODE solve. The twin is
# left as the fallback for the flip-flop regime, SS and infusion doses, which
# this dataset does not have. So the model file needs no ODE rewrite: it is the
# same `pk one_cpt_transit(...)` line plus `kappa KAPPA_CL` and
# `iov_column = OCC`.
#
# IOV is diagonal (Option A): each occasion (OCC column) draws its own KAPPA_CL,
# so CL_ik = TVCL * exp(KAPPA_CL_ik). The bundled data is an 8-subject subset of
# the ferx-core transit + IOV NONMEM anchor -- three doses per subject 48 h apart
# with near-complete washout -- simulated from this model, so the fit recovers the
# data-generating parameters (true TVCL=9, TVV=30, TVMTT=1, TVN=3, omega_V=0.09,
# kappa_CL=0.04, sigma=0.1). The fit, covariance step included, takes well
# under a second.
ex <- ferx_example("one_cpt_transit_iov")

fit <- ferx_fit(ex$model, ex$data)

# The printout shows the between-subject omega (ETA_V) and the between-occasion
# IOV kappa (KAPPA_CL) on their own lines (bare declared names).
print(fit)

# Point estimates, including the recovered IOV variance for KAPPA_CL.
fit$estimates

# Inspect the parsed model.
ferx_model_inspect(fit)
