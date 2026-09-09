library(ferx)

# Allometric body-size scaling, in its two uses: as a model transform (the
# scaled model, nothing fitted) and as a tool (both fits, side by side).

ex <- ferx_example("two_cpt_oral_base")

# --- 1. The transform ------------------------------------------------------
# `(WT/70)^0.75` on every clearance the template line binds, `(WT/70)^1.0` on
# every volume, written as `[covariate_model]` relations. Nothing is fitted, so
# this is a step you can put in the middle of a hand-built workflow.
scaled <- ferx_allometry(ex$model, ex$data, fit = FALSE)
scaled
scaled$scalings
cat(scaled$model)

# The product is a model file: fit it, or edit it further.
fit_scaled <- ferx_fit(scaled$model_path, ex$data)
fit_scaled

# --- 2. The tool -----------------------------------------------------------
# The base and the scaled model fitted side by side, which is what makes the
# scaling's cost visible: the dOFV, and each fit's convergence and strictness
# verdict.
res <- ferx_allometry(ex$model, ex$data,
                      covariate = "WT", reference = 70,
                      retries = 2,
                      directory = file.path(tempdir(), "allometry-run"))
res
res$comparison
res$dofv

# --- 3. Estimated exponents ------------------------------------------------
# The same relations, but with the exponents estimated from the convention's
# values instead of fixed at them.
est <- ferx_allometry(ex$model, ex$data, estimate = TRUE, lower = 0, upper = 2,
                      retries = 2)
est$scalings
est$comparison
