# FeRx-NLME/ferx-r#389's measurement: does a refused call still leak the parsed
# model and the population read from the dataset?
#
# `Rf_error()` longjmps out of Rust without unwinding. While the glue raised
# from inside each entry point, everything that body was holding was skipped by
# every destructor - about 10 KB per refusal on the bundled examples. #388's
# restructure moved the raise into `entry()`, which runs after the body's frame
# has returned, so the refused loop should now grow like the accepted one.
#
# What this script does NOT measure: the *arguments*. extendr protects every
# argument SEXP in a frame outside the body (ferx-r#394). The loops below pass
# a model path, a data path, `n_sim` and a seed - all bytes - so their numbers
# say "the body's locals are gone", not "nothing leaks". The argument half is
# tests/testthat/test-glue-refusal-args.R, which measures Vcells and is exact
# enough to run in CI.
#
# Not a testthat test on purpose: 3000 calls take ~13 minutes, and RSS is noisy
# enough that a threshold on it would be a coin flip in CI (the control loop
# below fell 7 MB over the same 2500 calls on the machine this was written on).
# Run it by hand when the glue's raise path changes.
#
# Usage, from the repo root, once per loop and once per build:
#   Rscript tools/measure-refusal-rss.R refused [n]
#   Rscript tools/measure-refusal-rss.R control [n]
#
# `R_LIBS` picks the build under test; each loop needs its own fresh process.

args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args) >= 1) args[1] else "refused"
n <- if (length(args) >= 2) as.integer(args[2]) else 3000L
stopifnot(which %in% c("refused", "control"))

suppressPackageStartupMessages(library(ferx))
ex <- ferx_example("pktte_joint")

# gc() before every reading, RSS from `ps`, as the issue measured it.
rss_mb <- function() {
  gc()
  as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
}

# The two loops parse the same model and read the same dataset. Only one of
# them then refuses: a joint PK-TTE model cannot be simulated without a finite
# administrative horizon, and the refusal arrives with `parsed` and the
# population live.
one <- if (which == "refused") {
  function() {
    tryCatch(ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L),
             error = function(e) NULL)
  }
} else {
  function() invisible(ferx_predict(ex$model, ex$data))
}

for (i in seq_len(n)) {
  suppressWarnings(one())
  if (i %% 500L == 0L) cat(sprintf("%s i=%d rss=%.1f\n", which, i, rss_mb()))
}
