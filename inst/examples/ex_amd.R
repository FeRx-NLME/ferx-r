library(ferx)

# Automatic model development (Pharmpy `amd`): the whole pipeline over one
# model and one search space - structural, variability, residual error,
# inter-occasion variability, allometry, covariates - each step starting from
# the model the previous one selected.
#
# This one fits a lot of models, so it takes minutes rather than seconds.

ex <- ferx_example("warfarin_amd")

# --- 1. What would run, before paying for it --------------------------------
# The plan is the engine's own: which step runs where, and why one will not.
# Nothing is fitted, so this is the cheap way to check a space before starting.
ferx_amd_plan(config = ex$search)

# The occasion, allometry and covariate steps are skipped here, each for a
# reason the plan spells out: warfarin.csv has no occasion column and no
# covariates, so there is nothing for them to decide.
plan <- ferx_amd_plan(config = ex$search)
plan[!is.na(plan$skipped), c("step", "skipped")]

# --- 2. The reproducible form -----------------------------------------------
# The .ferxsearch file is the artifact: the space, the criterion, the
# strictness gate, the per-tool sections and the [amd] options in one file.
run_dir <- file.path(tempdir(), "amd-run")
res <- ferx_amd(config = ex$search, directory = run_dir, progress = TRUE)
res

# --- 3. The step table ------------------------------------------------------
# One row per planned step, skipped ones included. `status` says whether it
# ran, `reason` says why not, and `converged` / `passed` are the termination
# status and the strictness verdict of the model the step selected - a winner
# that stalled is visible as such rather than hidden behind its rank.
res$steps[, c("index", "step", "status", "criterion", "value_before",
              "value_after", "d_value", "converged", "passed", "selected")]

# Why the steps that did not run did not run.
res$steps[res$steps$status != "ran", c("step", "status", "reason")]

# --- 4. Every candidate of every step ---------------------------------------
# The table the pipeline is audited on. `tool` is the step's tool, plus
# `start` for the pipeline's own first fit and `retries` for a perturbed-
# restart pass over a model a step had already selected.
head(res$candidates[, c("step", "tool", "id", "description", "criterion",
                        "value", "d_ofv", "rank", "passed", "selected")], 20)
table(res$candidates$tool)

# `criterion` names what `value` is on, which is not always the step's own: a
# ruvsearch pre-screen candidate is fitted to the parent's CWRES, so its number
# is on that scale and must not be compared with a data OFV.
unique(res$candidates[, c("step", "tool", "criterion")])

# Candidates the strictness gate excluded, with the gate's own words.
res$candidates[!res$candidates$passed, c("step", "tool", "id", "failures")]

# --- 5. The per-step view ---------------------------------------------------
# Each step that ran, with its directory and its own slice of the candidate
# table. The tool's fuller record (models.csv, every candidate's model text)
# stays in that directory.
names(res$tools)
res$tools[[1]]$selected
head(res$tools[[1]]$candidates[, c("id", "description", "value", "rank",
                                   "selected")])
list.files(run_dir)

# --- 6. The model the pipeline ended on -------------------------------------
# A fitted model, seeded with its own estimates - it reads and refits at the
# fit beside it.
res$fit
res$final_model_path
ferx_model_get_section(res$final_model_path, "parameters")
res$final_ofv
res$d_ofv           # final - start; negative is an improvement

# The engine's own report, as `ferx amd` prints it: the pipeline, every step's
# candidates, and the final estimates with their standard errors.
cat(res$summary_text)

# summary() is the same material as an R object: the step table, then each
# step's candidates, then everything the gate excluded.
summary(res)

# --- 7. The inline form, and a shorter strategy ------------------------------
# `SIR` runs structural, IIV and residual only. The inline form takes the space
# as MFL text; only the file can carry a per-tool section.
short <- ferx_amd(
  model        = ex$model,
  data         = ex$data,
  search_space = "ABSORPTION(FO); PERIPHERALS(0..1); IIV?(@PK, exp)",
  strategy     = "SIR",
  retries_on   = "final",
  retries      = 0,
  directory    = file.path(tempdir(), "amd-sir"),
  progress     = FALSE
)
short$steps[, c("index", "step", "status", "d_value", "passed")]

# The step order is the strategy's, asserted on the run rather than assumed.
short$steps$step

# --- 8. Reading a run back --------------------------------------------------
# The two tables are written beside the run, so a pipeline driven from the
# command line (`ferx amd config.ferxsearch`) reads back the same way.
read.csv(res$steps_csv)[, c("step", "status", "criterion", "value_after",
                            "selected")]
nrow(read.csv(res$candidates_csv))
