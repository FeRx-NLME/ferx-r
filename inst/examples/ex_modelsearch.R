library(ferx)

# Structural PK model search (Pharmpy `modelsearch`) on the warfarin example:
# is one compartment enough, and is there an absorption delay?
#
# This one fits models, so it takes minutes rather than seconds.

ex <- ferx_example("warfarin")

# --- 1. Look before you search ---------------------------------------------
# What the space expands to, and whether the engine can build all of it. Both
# are answered without fitting anything.
space <- ferx_search_space("ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
                           model = ex$model, data = ex$data)
print(space)
ferx_search_coverage(space)

# --- 2. The inline form ----------------------------------------------------
# No file to author: `search_space` is MFL text, quoted verbatim, and the run
# knobs mirror ferx_bootstrap()'s.
run_dir <- file.path(tempdir(), "modelsearch-inline")
res <- ferx_modelsearch(
  model        = ex$model,
  data         = ex$data,
  search_space = "ABSORPTION(FO); PERIPHERALS(0..1); LAGTIME([OFF,ON])",
  algorithm    = "reduced_stepwise",
  rank         = "bic",
  retries      = 2,
  directory    = run_dir,
  progress     = TRUE
)
res

# The model table is the engine's own: every candidate, with its termination
# status and strictness verdict beside the criterion. A model the gate excluded
# is a row saying why, not an absence.
res$models[, c("id", "parent", "layer", "structure", "ofv", "criterion",
               "d_criterion", "rank", "converged", "passed", "failures")]

# The winner is a fitted model, not a promise of one.
res$fit
res$final_model_path

# Every candidate's model text, named by id - so the model ranked second can be
# read, or refitted, without re-running the search.
runner_up <- res$models$id[!is.na(res$models$rank) & res$models$rank == 2]
if (length(runner_up) == 1) cat(res$model_text[[runner_up]])

# --- 3. The reproducible form ----------------------------------------------
# The same search, stated in the bundled file. Both forms build the same
# configuration and hand it to the same loader.
res2 <- ferx_modelsearch(config = ex$search,
                         directory = file.path(tempdir(), "modelsearch-file"))
summary(res2)

# --- 4. Reading a run back -------------------------------------------------
# The table the run wrote, read back with the engine's own column list - the
# same call works on a run produced by `ferx modelsearch` on the command line.
tab <- ferx_search_results(run_dir, type = "models")
tab[order(tab$criterion), c("id", "absorption", "peripherals", "lagtime",
                            "criterion", "rank")]

# --- 5. Resume -------------------------------------------------------------
# A run with a `directory` journals every candidate, so an interrupted search
# picks up where it stopped instead of refitting what it already knows.
res3 <- ferx_modelsearch(config = ex$search,
                         directory = file.path(tempdir(), "modelsearch-file"),
                         resume = TRUE)
res3$models$id
