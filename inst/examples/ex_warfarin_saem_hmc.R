library(ferx)

# HMC proposals for the SAEM E-step.
# Requires ferx built with Enzyme autodiff (AD build).
# When FERX_NO_AUTODIFF=1 is set, n_leapfrog is silently ignored and
# standard Metropolis-Hastings proposals are used instead.

ex <- ferx_example("warfarin_saem")

fit <- ferx_fit(
  ex$model,
  ex$data,
  method   = "saem",
  settings = list(n_leapfrog = 3L)
)
print(fit)

# Number of subjects whose E-step used HMC proposals.
# NULL when MH was used (no AD build or n_leapfrog = 0).
fit$saem_n_subjects_hmc
