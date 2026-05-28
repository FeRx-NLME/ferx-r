library(ferx)

# [scaling] block (Form A): obs_scale converts model predictions before the
# residual is computed.
#
# The dataset records AMT in micrograms (100 mg = 100000 ug) while DV is in
# mg/L. obs_scale = 1000 divides every prediction by 1000 so residuals - and
# the sigma estimate - are on the DV (mg/L) scale.
#
# PK parameter values match the standard warfarin example; the fit should
# converge to the same estimates as ferx_example("warfarin").
ex <- ferx_example("warfarin_scaled")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Compare against the unscaled warfarin reference - theta should agree to
# within numerical noise.
ref <- ferx_fit(ferx_example("warfarin")$model, ferx_example("warfarin")$data)
rbind(scaled = fit$theta, reference = ref$theta)
