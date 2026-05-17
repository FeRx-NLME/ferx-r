## Warfarin — one-compartment ODE model with SDE process noise (EKF)
##
## Demonstrates the [diffusion] block: DIFF_CENTRAL is a within-subject
## system-noise variance fitted alongside the structural PK parameters.
## The Extended Kalman Filter (EKF) integrates the covariance ODE alongside
## the state ODE and inflates the observation variance at each time point by
## the predicted P[central], giving a better-calibrated likelihood when
## IWRES shows autocorrelation.
##
## Key output differences vs. the standard warfarin fit:
##   fit$uses_sde          # TRUE
##   fit$theta["DIFF_CENTRAL"]  # fitted diffusion variance (variance units)

library(ferx)

ex <- ferx_example("warfarin_sde")

## Fit — autodiff is automatically forced to finite differences for SDE models
fit <- ferx_fit(ex$model, ex$data)

## uses_sde flag and DIFF_CENTRAL theta
fit$uses_sde
fit$theta["DIFF_CENTRAL"]

## Full parameter table (DIFF_CENTRAL appears alongside structural thetas)
ferx_estimates(fit)

## Standard diagnostics — IWRES and CWRES remain available
head(fit$sdtab[, c("ID", "TIME", "DV", "IPRED", "CWRES", "IWRES")])
