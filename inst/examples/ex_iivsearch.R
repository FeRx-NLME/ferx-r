library(ferx)

# Variability-structure search (Pharmpy `iivsearch`) on the correlated-random-
# effects warfarin example: which of CL, V and KA need an eta at all, and which
# of the survivors are worth correlating?
#
# This one fits models, so it takes minutes rather than seconds.

ex <- ferx_example("warfarin_block_omega")

# --- 1. The inline form -----------------------------------------------------
# The space is MFL, quoted verbatim. `IIV?` makes each PK parameter's eta
# exploratory; `COVARIANCE?(IIV, *)` offers a block over whichever of them
# survive. The run knobs mirror ferx_bootstrap()'s.
run_dir <- file.path(tempdir(), "iivsearch-inline")
res <- ferx_iivsearch(
  model        = ex$model,
  data         = ex$data,
  search_space = "IIV?(@PK, exp); COVARIANCE?(IIV, *)",
  retries      = 0,
  directory    = run_dir,
  progress     = TRUE
)
res

# --- 2. Two stages, not one table -------------------------------------------
# An iivsearch is two searches in sequence: how many etas, then which of them
# are blocked. Every model row carries the stage it was fitted in, so the two
# decisions can be read apart.
res$models[, c("id", "step", "step_kind", "structure", "criterion", "rank",
               "starts", "converged", "passed")]

# `$steps` is the engine's own per-stage ranking: the parent each candidate was
# compared with, and where it placed. `d_criterion` here is positive-is-better
# (the improvement over the parent); on `$models` it is the signed difference,
# so negative is better there.
res$steps

# The block stage starts from the eta stage's winner, not from the input.
res$steps[res$steps$kind == "block_structure", ]

# --- 3. The labels ----------------------------------------------------------
# The engine's `description` is written in parameter names - Pharmpy's
# spelling, and what models.csv carries. `structure` beside it is the same
# structure in the model's own declared eta names, which is the convention
# every ferx output surface follows.
res$models[, c("id", "description", "structure", "eta_labels", "block_labels")]

# The winner, said out loud: which parameters kept an eta, and which ended up
# correlated.
res$final_etas
res$final_blocks

# --- 4. The winning model ---------------------------------------------------
# The winner is a fitted model, not a promise of one, and its omega carries the
# declared eta names as dimnames.
res$fit
res$final_model_path
ferx_model_get_section(res$final_model_path, "parameters")
res$fit$omega

# Every fitted model's text, named by model id - so a structure the search
# rejected can still be read, or refitted, without re-running anything.
names(res$model_text)

# summary() adds the structures that were not selected, each with its reason.
summary(res)

# --- 5. One stage only ------------------------------------------------------
# `correlation_algorithm = "skip"` decides the number of etas and stops;
# `algorithm = "skip"` does the opposite, taking the input's etas as given and
# deciding only which of them are blocked.
blocks_only <- ferx_iivsearch(
  model                 = ex$model,
  data                  = ex$data,
  search_space          = "IIV(@PK, exp); COVARIANCE?(IIV, *)",
  algorithm             = "skip",
  correlation_algorithm = "top_down_exhaustive",
  retries               = 0,
  directory             = file.path(tempdir(), "iivsearch-blocks-only"),
  progress              = FALSE
)
unique(blocks_only$models$step_kind)
blocks_only$final_structure

# --- 6. The reproducible form -----------------------------------------------
# The same search, stated in the bundled file. Both forms build the same
# configuration and hand it to the same loader.
res2 <- ferx_iivsearch(config = ex$search,
                       directory = file.path(tempdir(), "iivsearch-file"))
res2$final_structure

# --- 7. Reading a run back --------------------------------------------------
# The table the run wrote, read back with the engine's own column list - the
# same call works on a run produced by `ferx iivsearch` on the command line.
tab <- ferx_search_results(run_dir, type = "models")
attr(tab, "tool")
tab[, c("id", "step", "description", "etas", "blocks", "criterion", "rank")]

# The runner's candidate table is a different table of the same run, stacked
# across the directories the run wrote.
head(res$candidates)
