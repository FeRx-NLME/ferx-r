library(ferx)

# Stepwise covariate modelling (PsN `scm`, Pharmpy `covsearch`) on the
# two-compartment oral example, run both ways: inline for a one-off, and from
# the bundled `.ferxsearch` file for a run someone else has to reproduce.
#
# This one fits models, so it takes minutes rather than seconds.

ex <- ferx_example("two_cpt_oral_base")

# --- 1. Look before you search ---------------------------------------------
# What the space expands to on *this* model, and whether the engine can express
# all of it. Both are answered without fitting anything.
space <- ferx_search_space("COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])",
                           model = ex$model, data = ex$data)
print(space)
ferx_search_coverage(space)

# --- 2. The inline form ----------------------------------------------------
# No file to author: `search_space` is MFL text, quoted verbatim, and the run
# knobs mirror ferx_bootstrap()'s.
run_dir <- file.path(tempdir(), "covsearch-inline")
res <- ferx_covsearch(
  model        = ex$model,
  data         = ex$data,
  search_space = "COVARIATE?(CL, WT, [pow, lin]); COVARIATE?(V1, WT, pow)",
  algorithm    = "scm-forward-then-backward",
  p_forward    = 0.01,
  p_backward   = 0.001,
  retries      = 2,
  directory    = run_dir,
  progress     = TRUE
)
res

# The step table is the engine's own: every candidate of every step, with its
# termination status and strictness verdict beside the dOFV. A candidate the
# gate excluded is a row saying why, not an absence.
res$steps[, c("step", "phase", "parameter", "covariate", "form",
              "dofv", "p_value", "converged", "passed", "failures")]

# The candidates that were actually selected, and the relations they left.
res$steps[res$steps$selected, ]
res$included

# The winner is a fitted model, not a promise of one.
res$fit
res$final_model_path

# Everything the runner journalled, including duplicates and reused fits.
head(res$candidates)

# --- 3. The reproducible form ----------------------------------------------
# The same search, stated in the bundled file. Both forms build the same
# configuration and hand it to the same loader.
res2 <- ferx_covsearch(config = ex$search,
                       directory = file.path(tempdir(), "covsearch-file"))
summary(res2)

# --- 4. Resume ------------------------------------------------------------
# A run with a `directory` journals every candidate, so an interrupted search
# picks up where it stopped instead of refitting what it already knows.
res3 <- ferx_covsearch(config = ex$search,
                       directory = file.path(tempdir(), "covsearch-file"),
                       resume = TRUE)
res3$steps$candidate
