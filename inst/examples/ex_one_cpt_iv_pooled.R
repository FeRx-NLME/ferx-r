library(ferx)

# Fixed-effects-only (naive-pooled) fit -- no random effects at all
# (ferx-core #989).
#
# The model is one_cpt_iv.ferx with both `omega` declarations and both
# `exp(ETA_*)` terms removed. With n_eta = 0 there is no inner empirical-Bayes
# problem and no `log|Omega|` term, so FOCE/FOCEI collapse to plain maximum
# likelihood: every subject shares one CL and one V, and `sigma` alone carries
# the spread. Useful for naive-pooled analyses, single-subject fits, and reduced
# models where a variance component has collapsed to the boundary and should be
# removed rather than fixed near zero.
#
# An [error_model] and its `sigma` are still required; only Omega is optional.
ex <- ferx_example("one_cpt_iv_pooled")

fit <- ferx_fit(ex$model, ex$data, method = "focei")
print(fit)

# Omega is a 0x0 matrix and there are no ETA columns anywhere. Nothing
# downstream invents a random effect.
stopifnot(
  nrow(fit$omega) == 0L,
  ncol(fit$omega) == 0L,
  is.null(fit$ebe_etas),
  length(fit$shrinkage_eta) == 0L
)

# With no random effects the individual prediction IS the population
# prediction, so the two residual flavours coincide exactly.
sd <- fit$sdtab
stopifnot(
  isTRUE(all.equal(sd$PRED,  sd$IPRED, tolerance = 0)),
  isTRUE(all.equal(sd$CWRES, sd$IWRES, tolerance = 0))
)

# Anchored against a NONMEM 7.6.0 run of the same model written as
# `$OMEGA 0 FIX` on a degenerate ETA (METHOD=0 SIGDIGITS=4). Note that simply
# *omitting* `$OMEGA` is not equivalent -- NM-TRAN then infers a single-subject
# analysis and drops the ID grouping.
cat(sprintf(
  "OFV  %.4f   (NONMEM -269.6370)\nTVCL %.4f   (NONMEM   4.84070)\nTVV  %.4f   (NONMEM  52.8324)\n",
  fit$ofv, fit$theta[["TVCL"]], fit$theta[["TVV"]]
))

# Standard errors need `covariance = "rsr"` to compare against NONMEM, whose
# $COVARIANCE default is the RSR sandwich. ferx defaults to "r" (the inverse
# Hessian). On a naive-pooled model the difference is not cosmetic: the model
# ignores within-subject correlation by construction, so the sandwich runs about
# twice as wide, and comparing the two defaults looks like a factor-of-two bug.
#
# This call is also the regression for ferx-r #290 -- ferx_covariance() used to
# reject any fit with no ETA columns outright.
fit_rsr <- ferx_covariance(fit, covariance_method = "rsr")
cat(sprintf(
  "SE TVCL %.6f (NONMEM 0.166298)   SE TVV %.6f (NONMEM 1.76021)\n",
  fit_rsr$se_theta[["TVCL"]], fit_rsr$se_theta[["TVV"]]
))
