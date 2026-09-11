library(ferx)

# Inter-occasion variability search (Pharmpy `iovsearch`) on the warfarin IOV
# example: which parameters earn a kappa, and do the etas beside them still
# earn theirs?
#
# This one fits models, so it takes minutes rather than seconds.

ex <- ferx_example("warfarin_iov")

# The base model reads its occasions from `iov_column = OCC`. That is what
# makes the search possible: without it the dataset's occasions are never read,
# and the tool says so rather than searching nothing.
ferx_model_get_section(ex$model, "fit_options")

# --- 1. The inline form -----------------------------------------------------
# There is no space to state unless you want one: the candidates default to
# every parameter carrying a free eta, which is Pharmpy's own default.
run_dir <- file.path(tempdir(), "iovsearch-inline")
res <- ferx_iovsearch(
  model     = ex$model,
  data      = ex$data,
  retries   = 0,
  directory = run_dir,
  progress  = TRUE
)
res

# --- 2. Two steps, not one table --------------------------------------------
# The full-IOV model - a kappa on every candidate parameter - is fitted first.
# Step 1 removes the kappas that do not earn their parameters; step 2, taking
# step 1's winner as its parent, removes the etas of the parameters that kept
# one.
res$models[, c("id", "parent", "step", "structure", "criterion", "rank",
               "converged", "passed")]
res$steps

# --- 3. The labels ----------------------------------------------------------
# The engine's `description` is written in parameter names; `structure` beside
# it is the same structure in the model's own declared random-effect names,
# which is the convention every ferx output surface follows.
res$models[, c("id", "description", "structure", "eta_labels", "kappa_labels")]

res$final_etas
res$final_kappas

# --- 4. The winning model ---------------------------------------------------
res$fit
res$final_model_path
ferx_model_get_section(res$final_model_path, "parameters")
ferx_model_get_section(res$final_model_path, "individual_parameters")

# summary() adds the structures that were not selected, each with its reason.
summary(res)

# --- 5. How the kappas are declared -----------------------------------------
# `same-as-iiv` - the default - blocks the kappas as their etas are blocked.
# `disjoint` gives each its own line, `joint` one block over all of them, and
# `explicit` takes the blocks from `groups`.
disjoint <- ferx_iovsearch(
  model        = ex$model,
  data         = ex$data,
  distribution = "disjoint",
  retries      = 0,
  directory    = file.path(tempdir(), "iovsearch-disjoint"),
  progress     = FALSE
)
disjoint$distribution
disjoint$final_kappas

# --- 6. The reproducible form -----------------------------------------------
# The same search, stated in the bundled file. Both forms build the same
# configuration and hand it to the same loader.
res2 <- ferx_iovsearch(config = ex$search,
                       directory = file.path(tempdir(), "iovsearch-file"))
res2$final_structure

# --- 7. Reading a run back --------------------------------------------------
# The table the run wrote, read back with the engine's own column list. Three
# tools write a `models.csv`; the header says which one did.
tab <- ferx_search_results(run_dir, type = "models")
attr(tab, "tool")
tab[, c("id", "step", "description", "etas", "kappas", "criterion", "rank")]

head(res$candidates)
