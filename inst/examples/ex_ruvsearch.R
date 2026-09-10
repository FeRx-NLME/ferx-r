library(ferx)

# Residual-error model search (Pharmpy `ruvsearch`) on the analytic Savic
# transit example: is a plain proportional error enough, or do the residuals
# ask for a power, a combined, or a time-varying magnitude?
#
# This one fits models, so it takes minutes rather than seconds.

ex <- ferx_example("one_cpt_transit")

# --- 1. The inline form -----------------------------------------------------
# There is no search space to state: the candidates are the four residual-error
# families, and `skip` is what narrows them. The run knobs mirror
# ferx_bootstrap()'s.
run_dir <- file.path(tempdir(), "ruvsearch-inline")
res <- ferx_ruvsearch(
  model     = ex$model,
  data      = ex$data,
  p_value   = 0.05,
  max_iter  = 2,
  retries   = 1,
  directory = run_dir,
  progress  = TRUE
)
res

# The step table is the engine's own: every candidate of every iteration, with
# the likelihood-ratio test, the termination status and the strictness verdict
# on the same row. A form the gate excluded is a row saying why, not an
# absence.
res$steps[, c("iteration", "feature", "family", "ofv", "dofv", "p_value",
              "significant", "selected", "converged", "passed", "failures")]

# Which form won, and at what p-value - a table, not prose.
res$steps[!is.na(res$steps$selected) & res$steps$selected,
          c("iteration", "feature", "family", "dofv", "p_value")]

# summary() adds the forms that were not selected, each with its reason.
summary(res)

# The winner is a fitted model, not a promise of one.
res$fit
res$final_model_path
ferx_model_get_section(res$final_model_path, "error_model")

# Every fitted model's text, named by candidate id - so a form the search
# rejected can still be read, or refitted, without re-running anything.
names(res$model_text)
cat(res$model_text[["power-1"]])

# --- 2. Pharmpy's default level ---------------------------------------------
# transit_oral.csv is simulated with a plain proportional error, so at the
# default p = 0.001 nothing is accepted and the search returns the model it
# started from. "The error model was already right" is a result, and the table
# still carries every p-value that says so.
strict <- ferx_ruvsearch(
  model     = ex$model,
  data      = ex$data,
  max_iter  = 1,
  retries   = 0,
  directory = file.path(tempdir(), "ruvsearch-strict"),
  progress  = FALSE
)
strict$final_features
strict$steps[, c("feature", "dofv", "p_value", "significant")]

# --- 3. The reproducible form -----------------------------------------------
# The same search, stated in the bundled file. Both forms build the same
# configuration and hand it to the same loader.
res2 <- ferx_ruvsearch(config = ex$search,
                       directory = file.path(tempdir(), "ruvsearch-file"))
res2$final_features

# --- 4. Reading a run back --------------------------------------------------
# The table the run wrote, read back with the engine's own column list - the
# same call works on a run produced by `ferx ruvsearch` on the command line.
tab <- ferx_search_results(run_dir, type = "steps")
attr(tab, "tool")
tab[, c("iteration", "candidate", "feature", "ofv", "p_value", "selected")]

# --- 5. Narrowing the families ----------------------------------------------
# `skip` leaves a family out, so the iteration below tests power and combined
# only. A family can also drop out on its own: on a `foce` model IIV_on_RUV is
# not tested at all - it needs eta-epsilon interaction - and the run says so in
# `$notes` rather than reporting a candidate that was never fitted. This model
# is `focei`, so nothing is skipped that way here.
narrow <- ferx_ruvsearch(
  model     = ex$model,
  data      = ex$data,
  skip      = c("time_varying", "IIV_on_RUV"),
  p_value   = 0.05,
  max_iter  = 1,
  retries   = 0,
  directory = file.path(tempdir(), "ruvsearch-narrow"),
  progress  = FALSE
)
narrow$options$skip
narrow$steps$feature
narrow$notes
