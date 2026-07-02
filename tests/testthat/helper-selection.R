# Shared fixtures for the data-selection tests (ferx_apply_selection() and
# its call sites in ferx_fit()/ferx_runlog()/ferx_save_fit()), split out of
# the old test-data-selection.R when R/selection.R was split one-function-
# per-file. testthat sources helper-*.R once before any test-*.R file, so a
# single copy (and a single cached fit) here is shared by all of them.

# ferx_data objects are data.frame subclasses; attributes are NOT accessible
# via `$` (which looks at columns). Use attr() for metadata.
.sel_attr <- function(x, name) attr(x, name, exact = TRUE)

# Cached fit for warfarin_data_selection (FD gradient, no Enzyme required).
warfarin_sel_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin_data_selection")
      fit <<- ferx_fit(ex$model, ex$data, verbose = FALSE,
                       settings = list(maxiter = 30L))
    }
    fit
  }
})

# Tiny inline dataset: 4 obs records, predictable values.
.sel_test_df <- function() {
  data.frame(
    ID   = c(1L, 1L, 2L, 2L),
    TIME = c(0, 1, 0, 1),
    DV   = c(0.5, 2.0, 3.0, 5.5),
    EVID = c(0L, 0L, 0L, 0L),
    AMT  = c(0, 0, 0, 0),
    MDV  = c(0L, 0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
}
