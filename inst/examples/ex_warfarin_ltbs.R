library(ferx)

# Log-transform-both-sides (LTBS) error model.
#
# `log(DV) ~ additive(SIGMA_LOG)` tells the engine to log-transform DV and
# predictions, then fit additive (normal) error on the log scale.
# Under LTBS the sigma estimate is the residual SD on the log scale and
# approximates the CV on the natural scale.
ex <- ferx_example("warfarin_ltbs")

fit <- ferx_fit(ex$model, ex$data, method = "focei")
print(fit)

# Inspect the parsed model: residual reported as "additive (log-transformed)".
ferx_model_inspect(fit)

# Sigma is reported as `SIGMA_LOG`. The `.fitrx` and `print()` surfaces label
# the residual as "additive (log-transformed)".
fit$sigma
fit$sigma_types
