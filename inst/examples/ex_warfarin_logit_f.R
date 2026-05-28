library(ferx)

# Logit-normal bioavailability:
#   F_i = inv_logit(logit(THETA_F) + ETA_F_i)
# THETA_F is specified directly on the (0, 1) scale in the [parameters] block;
# ferx maps it to the logit scale internally so the optimiser works on an
# unbounded parameter. omega ETA_F is BSV variance on the logit scale.
ex <- ferx_example("warfarin_logit_f")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# THETA_F estimate stays on the natural (0, 1) scale; ETA_F variance is on
# the logit scale and does not have a CV% interpretation.
fit$theta["THETA_F"]
fit$ebe_etas
