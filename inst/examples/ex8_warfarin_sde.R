## Warfarin — one-compartment ODE model with SDE process noise (EKF)
##
## Demonstrates the [diffusion] block: DIFF_CENTRAL is a within-subject
## system-noise variance fitted alongside the structural PK parameters.
## The Extended Kalman Filter (EKF) integrates the covariance ODE alongside
## the state ODE and inflates the observation variance at each time point by
## the predicted P[central], giving a better-calibrated likelihood when
## IWRES shows autocorrelation.
##
## Performance note: SDE models use finite differences (FD) for gradients.
## This makes each function evaluation significantly slower than a comparable
## analytical PK model.
##
## Key output differences vs. the standard warfarin fit:
##   fit$uses_sde               # TRUE
##   fit$theta["DIFF_CENTRAL"]  # fitted diffusion variance (variance units)

library(ferx)

ex <- ferx_example("warfarin_sde")

## Fit — gradient_method = fd is set in the model file (required for SDE)
fit <- ferx_fit(ex$model, ex$data)

## uses_sde flag and DIFF_CENTRAL theta
print(fit$uses_sde)
print(fit$theta["DIFF_CENTRAL"])

## Full parameter table (DIFF_CENTRAL appears alongside structural thetas)
print(ferx_estimates(fit))

## Standard diagnostics — IWRES and CWRES remain available
print(head(fit$sdtab[, c("ID", "TIME", "DV", "IPRED", "CWRES", "IWRES")]))

## .fitrx round-trip — uses_sde and DIFF_CENTRAL survive serialisation
tmp <- tempfile(fileext = ".fitrx")
ferx_save_fit(fit, tmp)
fit2 <- ferx_load_fit(tmp)
print(fit2$uses_sde)
print(fit2$theta["DIFF_CENTRAL"])
