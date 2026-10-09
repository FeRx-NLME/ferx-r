# Cached fits of a bundled example with its ID column relabelled (#468).
#
# The same data under different subject labels must give the same answer, so a
# fit on relabelled IDs and a fit on the bundled ones are twins: `ferx_sir()` /
# `ferx_covariance()` run standalone on each must agree exactly. Both arms are
# standalone runs on one platform, so the comparison is `identical()`; nothing
# here compares a standalone run against the in-fit one (which differs at
# ~3e-12 on Linux, #467).
#
# `relabel` maps the bundled ID column to the new labels; `key` names the cache
# entry. The rewritten CSV lives for the whole test run, since the fit records
# its path and hash and the standalone steps re-read it. The cache is keyed on
# `key`, not on `relabel`, so each key's relabelling is defined once, below.
relabelled_fit <- local({
  cache <- list()
  function(example, key, relabel, covariance = TRUE, maxiter = 30L) {
    id <- paste(example, key, covariance, sep = "/")
    if (is.null(cache[[id]])) {
      ex <- ferx_example(example)
      d <- utils::read.csv(ex$data, stringsAsFactors = FALSE)
      d$ID <- relabel(d$ID)
      path <- tempfile(paste0("relabel_", key, "_"), fileext = ".csv")
      # Text IDs are quoted (`"` doubled), so `Smith, 2019` stays one cell.
      utils::write.csv(d, path, row.names = FALSE,
                       quote = if (is.character(d$ID)) match("ID", names(d)) else FALSE)
      cache[[id]] <<- ferx_fit(
        ex$model, path,
        method = "focei", verbose = FALSE,
        covariance = covariance, settings = list(maxiter = maxiter)
      )
    }
    cache[[id]]
  }
})

# The "gappy" key's labels: non-1-based, with gaps.
relabel_gappy <- function(id) 3L * id + 107L

# Capture the arguments a binding is called with, then abort. The SIR and
# covariance steps are expensive and their output cannot show everything they
# were handed, so the forwarding is asserted on the call itself.
capture_binding_args <- function() {
  seen <- NULL
  list(
    fake = function(...) {
      seen <<- list(...)
      stop("captured", call. = FALSE)
    },
    seen = function() seen
  )
}
