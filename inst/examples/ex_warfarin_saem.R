library(ferx)

# Basic SAEM estimation. SAEM is the stochastic-approximation EM alternative
# to FOCE/FOCEI - more robust on models with many random effects or
# correlated omega, at the cost of stochastic convergence (set `seed` for
# reproducibility).
#
# See ex_warfarin_saem_hmc.R for HMC proposals (n_leapfrog) and the
# omega_burnin / adapt_interval SAEM tuning knobs.
ex <- ferx_example("warfarin_saem")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# SAEM-specific diagnostic: how many subjects used HMC proposals.
# NULL for plain MH (for example, n_leapfrog = 0).
fit$saem_n_subjects_hmc
