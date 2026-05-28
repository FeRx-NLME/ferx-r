library(ferx)

# Cookbook of outer-optimizer choices and tuning knobs forwarded via the
# `settings = list(...)` argument to ferx_fit(). All snippets fit the same
# warfarin example so estimates are comparable across optimizers.
#
# Pick an optimizer based on the model surface:
#   - bobyqa (default): derivative-free, robust, slowest convergence
#   - slsqp / lbfgs / nlopt_lbfgs / mma: gradient-based NLopt optimizers
#   - bfgs: hand-rolled quasi-Newton
#   - trust_region: Steihaug CG inside a trust region; good for ill-conditioned
ex <- ferx_example("warfarin")

# --- Outer optimizer selection ----------------------------------------------
fit_bobyqa <- ferx_fit(ex$model, ex$data, settings = list(optimizer = "bobyqa"))
fit_slsqp  <- ferx_fit(ex$model, ex$data, settings = list(optimizer = "slsqp"))
fit_lbfgs  <- ferx_fit(ex$model, ex$data, settings = list(optimizer = "lbfgs"))
fit_mma    <- ferx_fit(ex$model, ex$data, settings = list(optimizer = "mma"))
fit_bfgs   <- ferx_fit(ex$model, ex$data, settings = list(optimizer = "bfgs"))
fit_tr     <- ferx_fit(
  ex$model, ex$data,
  settings = list(optimizer = "trust_region", steihaug_max_iters = 100L)
)

# Compare OFV across optimizers (all should land on the same minimum up to
# tolerance; differences indicate convergence issues for that optimizer).
sapply(
  list(bobyqa = fit_bobyqa, slsqp = fit_slsqp, lbfgs = fit_lbfgs,
       mma = fit_mma, bfgs = fit_bfgs, trust_region = fit_tr),
  function(f) f$ofv
)

# --- Global search before local refinement ----------------------------------
# Useful when initial values are far from the optimum; runs NLopt's global
# optimizer (DIRECT) for `global_maxeval` evaluations, then refines locally.
fit_global <- ferx_fit(
  ex$model, ex$data,
  settings = list(global_search = TRUE, global_maxeval = 1000L)
)

# --- Iteration caps and inner-loop tolerance --------------------------------
# Tighten the per-subject EBE loop for difficult problems, or loosen it for
# faster fits during exploratory model building.
fit_tight <- ferx_fit(
  ex$model, ex$data,
  settings = list(maxiter = 500L, inner_maxiter = 100L, inner_tol = 1e-6)
)

# --- Stagnation guard (NLopt optimizers only) -------------------------------
# Default TRUE: short-circuit NLopt outer loop once OFV improvement falls
# below 1e-3 over a rolling window. Set FALSE to let NLopt run to xtol/ftol
# or maxiter. Useful when debugging slow-but-real OFV improvements.
fit_no_guard <- ferx_fit(
  ex$model, ex$data,
  settings = list(optimizer = "slsqp", stagnation_guard = FALSE)
)

# --- Reconverged-gradient schedule (FOCE/FOCEI only, non-IOV) ---------------
# 0 (default): fixed-EBE gradient (cheap). 1: reconverge every gradient eval
# (~5-6x cost, recovers the full surface). N: every N-th. Use when slsqp
# stalls above the bobyqa optimum on ill-conditioned non-IOV fits.
fit_reconv <- ferx_fit(
  ex$model, ex$data,
  settings = list(optimizer = "slsqp", reconverge_gradient_interval = 1L)
)

# --- Gauss-Newton damping ---------------------------------------------------
# gn_lambda is the Levenberg-Marquardt damping factor for method = "gn".
fit_gn <- ferx_fit(
  ex$model, ex$data,
  method   = "gn",
  settings = list(gn_lambda = 1e-3)
)

# --- Multi-start (covered in detail by ex9_mm_multistart.R) -----------------
# Runs n_starts independent fits from perturbed initial values in parallel
# (rayon) and returns the converged result with the lowest OFV.
fit_multi <- ferx_fit(
  ex$model, ex$data,
  settings = list(n_starts = 4L, start_sigma = 0.3, multi_start_seed = 123L)
)
fit_multi$ofv
