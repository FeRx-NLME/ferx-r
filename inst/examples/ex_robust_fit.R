library(ferx)

# Robustness knobs for awkward datasets and rough starting values.
ex <- ferx_example("warfarin")

# --- NCA-derived starting values --------------------------------------------
# inits_from_nca overrides the model file's initial theta with values
# estimated from the data via noncompartmental analysis. Most useful when
# the model file's defaults are far from the truth (e.g. multi-start with
# trust_region, or method = "gn").
#
# Strategies:
#   "nca"        : pure NCA per-subject summaries
#   "nca_sweep"  : NCA + a coarse parameter sweep (default for TRUE)
#   "nca_ebe"    : NCA + a brief EBE pass to refine
fit_nca <- ferx_fit(ex$model, ex$data, inits_from_nca = "nca_sweep")
print(fit_nca)

# Inspect the NCA-derived starting values without fitting:
inits <- ferx_inits_from_nca(ex$model, ex$data, method = "nca_sweep")
inits

# --- Tolerant convergence on awkward subjects -------------------------------
# By default ferx_fit() errors if any subject's inner EBE loop fails to
# converge. For large datasets with a handful of sparse / outlier subjects,
# `max_unconverged_frac` flags the fit converged when at most that fraction
# of (non-sparse) subjects fails; `min_obs_for_convergence_check` excludes
# subjects with fewer than N observations from the count entirely.
fit_tolerant <- ferx_fit(
  ex$model, ex$data,
  settings = list(
    max_unconverged_frac          = 0.10,   # allow up to 10% inner-loop failures
    min_obs_for_convergence_check = 2L      # ignore subjects with < 2 obs
  )
)
print(fit_tolerant)
