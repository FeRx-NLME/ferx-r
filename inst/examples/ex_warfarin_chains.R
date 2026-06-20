library(ferx)

ex <- ferx_example("warfarin")

# Method chain: SAEM warm-up followed by FOCEI refinement.
# Each stage is seeded with the previous stage's converged parameters.
fit_chain <- ferx_fit(ex$model, ex$data, method = c("saem", "focei"))
print(fit_chain)
fit_chain$method_chain  # c("saem", "focei")

# Terminal IMP stage: run importance sampling after FOCEI to get a more
# accurate marginal log-likelihood. "imp" must be the last entry. Tuned
# via settings: imp_samples (per-subject draws), imp_proposal_df (t-proposal
# degrees of freedom; smaller => heavier tails), imp_seed, and
# imp_low_ess_threshold (per-subject ESS fraction below which the subject
# is flagged in low_ess_subject_*).
fit_imp <- ferx_fit(
  ex$model,
  ex$data,
  method   = c("focei", "imp"),
  settings = list(
    imp_samples           = 2000L,
    imp_proposal_df       = 4,
    imp_seed              = 1L,
    imp_low_ess_threshold = 0.10   # fraction of imp_samples; flag subjects below
  )
)
print(fit_imp)

# IS marginal log-likelihood + Monte-Carlo SE; comparable across models in
# the same nesting family. mc_standard_error scales as 1 / sqrt(imp_samples).
fit_imp$importance_sampling$minus2_log_likelihood
fit_imp$importance_sampling$mc_standard_error

# Per-subject ESS diagnostics (fraction of imp_samples). ess_min flags the
# worst-conditioned subject; low_ess_subject_ids lists every subject whose
# ESS fraction fell below imp_low_ess_threshold.
fit_imp$importance_sampling$ess_min
fit_imp$importance_sampling$ess_median
fit_imp$importance_sampling$low_ess_subject_ids
fit_imp$importance_sampling$low_ess_subject_frac

# For IOV models, kappa_treatment reports whether kappa was marginalised
# ("marginalized"), fixed at its mode ("fixed_at_mode"), or absent
# ("not_applicable"). Warfarin has no IOV here.
fit_imp$importance_sampling$kappa_treatment
