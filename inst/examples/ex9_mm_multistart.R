## Multi-start Michaelis-Menten elimination
#
# The Vmax/Km pair is weakly identifiable from sparse PK data: the OFV
# surface contains a narrow ridge and a single-start run often settles on a
# local minimum. Starting values are intentionally set 3x above the true
# values (TVVMAX=12 vs true≈3.5, TVKM=20 vs true≈5.5).
#
# n_starts = 8 runs 8 parallel fits from perturbed starting points and
# returns the converged result with the lowest OFV. On an 8-core machine
# wall-clock time is approximately the same as a single run.

library(ferx)

ex <- ferx_example("mm_multistart")

# Single-start — may land on a local minimum from the inflated starting values
fit_single <- ferx_fit(
  model = ex$model,
  data  = ex$data,
  method = "focei"
)

# Multi-start — 8 parallel runs with wider perturbation for a ridge surface
fit_multi <- ferx_fit(
  model  = ex$model,
  data   = ex$data,
  method = "focei",
  settings = list(
    n_starts         = 8L,
    start_sigma      = 0.5,
    multi_start_seed = 42L
  )
)

cat("Single-start OFV:", fit_single$ofv, "\n")
cat("Multi-start OFV: ", fit_multi$ofv,  "\n")

print(fit_multi)
