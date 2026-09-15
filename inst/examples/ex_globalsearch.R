library(ferx)

# Global model search: structure and covariates decided together, on one grid.
#
# The starting model is one-compartment, covariate-free, on a dataset that
# carries WT and CRCL and was simulated from two compartments - so every axis
# of the grid is a real decision.
#
# This one fits models, so it takes minutes rather than seconds.
#
# Not to be confused with `ferx_fit(settings = list(global_search = TRUE))`,
# which is a global *optimizer* phase inside the estimation of one model. This
# is a global search over *models*.

ex <- ferx_example("two_cpt_oral_global")

# --- 1. Look before you search ---------------------------------------------
# What the space expands to, and whether the engine can build all of it. Both
# are answered without fitting anything.
space <- ferx_search_space(
  "PERIPHERALS(0..1); LAGTIME([OFF,ON])
   COVARIATE?(CL, WT, pow); COVARIATE?(CL, CRCL, pow)",
  model = ex$model, data = ex$data
)
print(space)
ferx_search_coverage(space)

# --- 2. The reproducible form ----------------------------------------------
# The bundled file enumerates the grid exhaustively - sixteen points, which is
# small enough to fit every one of them and so is the reference the genetic
# algorithm below can be checked against.
run_dir <- file.path(tempdir(), "globalsearch-exhaustive")
res <- ferx_globalsearch(config = ex$search, directory = run_dir, progress = TRUE)
res

# The model table is the engine's own: one row per grid point evaluated, with
# the criterion, the charges the search added to it, and the fitness they sum
# to. A genome the strictness gate refused is a row carrying its verdict and
# its gate charge - never an absence.
res$models[, c("id", "genome", "ofv", "criterion", "charge_non_influential",
               "charge_gate", "charge_crash", "fitness", "rank", "passed")]

# Which axes there were, and what each could take.
res$axes
res$space_size

# Grid points cost fewer fits than they look: two genomes that render to the
# same model are fitted once, and the second row says which row it duplicates.
c(evaluated = nrow(res$models), fitted = res$n_fitted)
res$models[!is.na(res$models$duplicate_of),
           c("id", "genome", "duplicate_of", "non_influential")]

# The winner is a fitted model, not a promise of one.
res$fit
res$final_model_path

# Every evaluated model's text, named by id - so the model the table ranked
# second can be read, or refitted, without re-running the search.
runner_up <- res$models$id[!is.na(res$models$rank) & res$models$rank == 2]
if (length(runner_up) == 1) cat(res$model_text[[runner_up]])

# --- 3. The genetic algorithm ----------------------------------------------
# The same grid, searched rather than enumerated. On sixteen points the GA
# should find the same winner while fitting fewer models; on a grid of
# thousands, finding it at all is the point.
ga_dir <- file.path(tempdir(), "globalsearch-ga")
res_ga <- ferx_globalsearch(
  model        = ex$model,
  data         = ex$data,
  search_space = c("PERIPHERALS(0..1); LAGTIME([OFF,ON])",
                   "COVARIATE?(CL, WT, pow); COVARIATE?(CL, CRCL, pow)"),
  algorithm    = "ga",
  ga           = list(population_size = 8, generations = 3, seed = 20250914),
  directory    = ga_dir,
  progress     = TRUE
)
res_ga$generations
c(exhaustive = res$final_fitness, ga = res_ga$final_fitness)

# --- 4. What the fitness charged -------------------------------------------
# `[rank] type` defaults to pyDarwin's penalized fitness here, where every
# other search tool defaults to a BIC. The effective schedule is on the
# object, so the fitness column can be read against the charges that built it.
res$criterion
res$penalties

# The three search-level charges apply under *any* criterion - rank on the
# mixed BIC and a gated candidate still pays the gate charge.
res_bic <- ferx_globalsearch(
  model        = ex$model,
  data         = ex$data,
  search_space = "PERIPHERALS(0..1); COVARIATE?(CL, WT, pow)",
  algorithm    = "exhaustive",
  rank         = "bic",
  penalties    = list(gate = 250),
  directory    = file.path(tempdir(), "globalsearch-bic")
)
res_bic$criterion
res_bic$models[, c("id", "genome", "criterion", "charge_gate", "fitness")]

# --- 5. Reading a run back -------------------------------------------------
# The table the run wrote, read back with the engine's own column list - the
# same call works on a run produced by `ferx globalsearch` on the command
# line.
tab <- ferx_search_results(run_dir, type = "models")
tab[order(tab$fitness), c("id", "genome", "criterion", "fitness", "rank")]

# --- 6. Resume -------------------------------------------------------------
# A run with a `directory` journals every candidate, so an interrupted search
# picks up where it stopped instead of refitting what it already knows.
res2 <- ferx_globalsearch(config = ex$search, directory = run_dir, resume = TRUE)
sum(res2$models$reused)
