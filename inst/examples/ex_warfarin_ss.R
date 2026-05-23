library(ferx)

ex <- ferx_example("warfarin_ss")

# The dose row in warfarin_ss.csv has SS=1 and II=24.
# The engine computes steady-state initial conditions automatically.
fit <- ferx_fit(ex$model, ex$data, method = "foce")
print(fit)

# Estimates should be close to the standard warfarin example (same true PK),
# though not identical because the dataset is a different simulated sample.
fit$theta
