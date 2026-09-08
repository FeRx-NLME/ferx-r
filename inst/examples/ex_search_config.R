library(ferx)

# Demonstrates the model-space search surface: validating a `.ferxsearch`
# configuration, expanding an MFL search space against a base model, and
# reading a coverage table. None of this runs a fit -- the point is to find out
# what the search would explore, and what it cannot express, *before*
# committing to a run.

ex <- ferx_example("two_cpt_oral_cov")

# --- 1. Load and validate the bundled configuration ------------------------
# The engine's loader is the validation: an unknown section, an unparseable
# space, a coverage gap or an unimplemented rank type is an error here.
cfg <- ferx_search_config(ex$search)
print(cfg)

cfg$base                  # resolved relative to the config file's directory
cfg$mfl                   # the [space] source, verbatim
cfg$space                 # one row per feature
cfg$rank$type             # "bic"
cfg$strictness            # the effective gate (file keys over engine defaults)
cfg$strictness_set        # which of them the file stated
cfg$tools                 # tool sections kept for their tool: "covsearch"

# --- 2. Expand the space against the base model ----------------------------
# `@IIV` and `@CONTINUOUS` mean nothing until they are resolved against a
# model and its dataset. This is what the search would actually explore.
sp <- ferx_search_space(cfg)
print(sp)

# The ground covariate effects: parameter x covariate x effect form.
effects <- attr(sp, "covariate_effects")
head(effects)
nrow(effects)             # 5 parameters x 2 covariates x 2 forms = 20

# Parsing without a model leaves the symbols as written.
ferx_search_space("COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])")

# --- 3. Coverage: what the engine can express ------------------------------
# As a table rather than an aborted run, so an unsupported feature is a row.
ferx_search_coverage(sp)

# Michaelis-Menten elimination has no candidate the search can build, and the
# reason says why. Putting this in a `.ferxsearch` is an error at load.
ferx_search_coverage("ELIMINATION([FO, MM])")

try(ferx_search_config(local({
  bad <- tempfile(fileext = ".ferxsearch")
  writeLines(c(
    'base = "model.ferx"',
    "[space]",
    'mfl = "ELIMINATION([FO, MM])"'
  ), bad)
  bad
})))

# --- 4. Read a finished run's candidate table ------------------------------
# `ferx_search_results()` types the engine's `candidates.csv`: logical
# converged / passed / reused, numeric criterion / ofv / seconds, and NA for
# the empty cells. A cancelled run's `candidates.partial.csv` reads the same
# way. (The runner itself arrives with covsearch; there is no run directory to
# read yet.)
#
#   res <- ferx_search_results("search-run-1")
#   res[res$passed, c("id", "features", "criterion")]
#   res[!res$passed, c("id", "failures")]
