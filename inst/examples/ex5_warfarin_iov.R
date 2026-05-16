library(ferx)
ex <- ferx_example("warfarin_iov")

fit <- ferx_fit(ex$model, ex$data)
fit

# IOV-specific output
fit$omega_iov       # IOV omega matrix (KAPPA variance)
fit$ebe_kappas      # per-subject, per-occasion kappa EBEs
