library(ferx)

ex <- ferx_example("warfarin_iov_saem")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

fit$omega_iov           # IOV omega matrix (KAPPA_CL variance)
fit$ebe_kappas          # per-subject, per-occasion kappa EBEs
fit$saem_n_subjects_hmc # NULL (MH proposals; set n_leapfrog > 0 for HMC)
