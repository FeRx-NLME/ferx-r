library(ferx)

# HMC proposals for the SAEM E-step.

ex <- ferx_example("warfarin_saem")

fit <- ferx_fit(
  ex$model,
  ex$data,
  method   = "saem",
  settings = list(
    n_leapfrog     = 3L,
    omega_burnin   = 30L,  # hold Omega fixed for the first 30 iterations
    adapt_interval = 10L   # adapt MH/HMC proposal scale every 10 iterations
  )
)
print(fit)

# Number of subjects whose E-step used HMC proposals.
# NULL when MH was used (for example, n_leapfrog = 0).
fit$saem_n_subjects_hmc
