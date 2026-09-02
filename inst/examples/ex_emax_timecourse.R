library(ferx)

# Compartment-free Emax response-time model -- the NONMEM `$PRED` equivalent.
# There are no doses, compartments, or ODEs: the assignments in
# [structural_model] calculate the prediction directly for each observation.
ex <- ferx_example("emax_timecourse")

fit <- ferx_fit(ex$model, ex$data)
print(fit)

# The synthetic dataset was generated with TVE0 = 10, TVEMAX = 6, TVET50 = 2,
# omega(E0) = 0.09, and additive sigma = 0.5. The fit is also anchored against
# an equivalent NONMEM `$PRED` model in ferx-core.
fit$theta
fit$omega
fit$sigma

# Population predictions require no dosing columns; TIME and the model
# parameters are all that the structural equation needs.
pred <- ferx_predict(ex$model, ex$data, fit = fit)
head(pred)
