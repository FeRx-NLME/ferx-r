library(ferx)

# Additive ETA on lag time: TLAG = TVTLAG + ETA_TLAG (normal, not log-normal).
# Use when the parameter has a natural additive variability structure rather
# than a multiplicative one. The omega entry for ETA_TLAG is the variance in
# the parameter's own units (here, h^2).
ex <- ferx_example("warfarin_additive_eta")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# ETA_TLAG variance is additive (no log-normal CV interpretation); check the
# per-subject EBEs to see the additive shifts around TVTLAG.
fit$ebe_etas
