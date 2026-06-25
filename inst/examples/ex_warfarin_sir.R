library(ferx)

# Sampling Importance Resampling (SIR) for non-asymptotic parameter
# uncertainty. SIR draws candidate parameter vectors from a multivariate
# normal proposal centred on the MLE with the covariance step's variance,
# then resamples them weighted by the actual likelihood. The resulting
# 95% CIs are robust to non-Gaussian likelihood shapes around the MLE -
# more reliable than the standard +/- 1.96 * SE intervals when the surface
# is skewed.
#
# Requires covariance = TRUE; SIR is enabled via sir = TRUE.
#
# We use the warfarin_iov example here because its FOCE covariance step
# converges cleanly on this small finite-difference-gradient fit. The same
# pattern works on any model where the covariance step succeeds.
ex <- ferx_example("warfarin_iov")

fit <- ferx_fit(
  ex$model,
  ex$data,
  sir      = TRUE,
  settings = list(
    sir_samples   = 2000L,   # candidate draws from the asymptotic proposal
    sir_resamples = 500L,    # resamples weighted by the true likelihood
    sir_seed      = 1L
  )
)
print(fit)

# SIR effective sample size. Values close to sir_resamples indicate a good
# proposal/likelihood match. Low ESS means the asymptotic covariance is a
# poor approximation of the true uncertainty - inspect the per-parameter
# CIs vs the +/- 1.96 * SE intervals to see which parameters are skewed.
fit$sir_ess

# 95% CIs from SIR. Compare against the asymptotic SE * 1.96 intervals.
fit$sir_ci_theta
fit$sir_ci_omega
fit$sir_ci_sigma
