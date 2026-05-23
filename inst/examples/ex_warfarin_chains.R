library(ferx)

ex <- ferx_example("warfarin")

# Method chain: SAEM warm-up followed by FOCEI refinement.
# Each stage is seeded with the previous stage's converged parameters.
fit_chain <- ferx_fit(ex$model, ex$data, method = c("saem", "focei"))
print(fit_chain)
fit_chain$method_chain  # c("saem", "focei")

# Terminal IMP stage: run importance sampling after FOCEI to get a more
# accurate marginal log-likelihood. "imp" must be the last entry.
fit_imp <- ferx_fit(ex$model, ex$data, method = c("focei", "imp"))
print(fit_imp)
fit_imp$importance_sampling  # list: minus2_log_likelihood, mc_standard_error, ess_min, ...
