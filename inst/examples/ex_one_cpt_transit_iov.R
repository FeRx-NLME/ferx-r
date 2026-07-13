library(ferx)

# Analytic Savic transit-compartment absorption WITH inter-occasion variability
# (IOV) -- `pk one_cpt_transit(cl, v, n, mtt)` + a `kappa` random effect on CL
# (ferx-core #719).
#
# The analytic absorption closed forms (one_cpt_transit / two_cpt_transit /
# one_cpt_ig / two_cpt_ig) used to REJECT IOV: a per-dose closed-form
# superposition cannot carry a dose's drug across an occasion boundary where the
# disposition changes. ferx now handles it transparently -- a subject with IOV is
# rerouted, per subject, to the model's exact `transit()` ODE twin
# (d/dt(central) = transit(n, mtt) - (CL/V) * central), which integrates the
# cross-occasion carryover exactly. So the model file needs no ODE rewrite: it is
# the same `pk one_cpt_transit(...)` line plus `kappa KAPPA_CL` and
# `iov_column = OCC`.
#
# IOV is diagonal (Option A): each occasion (OCC column) draws its own KAPPA_CL,
# so CL_ik = TVCL * exp(KAPPA_CL_ik). The bundled data is an 8-subject subset of
# the ferx-core transit + IOV NONMEM anchor -- three doses per subject 48 h apart
# with near-complete washout -- simulated from this model, so the fit recovers the
# data-generating parameters (true TVCL=9, TVV=30, TVMTT=1, TVN=3, omega_V=0.09,
# kappa_CL=0.04, sigma=0.1). The subset is kept small because each IOV subject
# integrates an ODE twin, so a larger dataset is markedly slower (ferx-core #812).
ex <- ferx_example("one_cpt_transit_iov")

fit <- ferx_fit(ex$model, ex$data)

# The printout shows the between-subject omega (ETA_V) and the between-occasion
# IOV kappa (KAPPA_CL) on their own lines (bare declared names).
print(fit)

# Point estimates, including the recovered IOV variance for KAPPA_CL.
fit$estimates

# Inspect the parsed model -- note ferx rerouting the analytic one_cpt_transit to
# its transit() ODE twin because the model carries IOV.
ferx_model_inspect(fit)
