# Package-internal state: holds the trace path of the most recent ferx_fit()
# call (when optimizer_trace = TRUE), so ferx_trace() can be called with no
# argument to inspect the last run.
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

# Coerce a fit's recorded file hash into the string the Rust bindings expect,
# where "" means "skip the integrity check". Three states are handled:
#   1. Non-empty hex string -> forwarded as-is; Rust enforces equality.
#   2. NULL (older binary that never populated the field) -> "".
#   3. NA_character_ (sha256_file failed at fit time; api.rs `.ok()` turned the
#      Err into None, stored R-side as NA) -> "". This coercion is essential:
#      `%||%` only catches NULL, and NA at the FFI boundary stringifies to "NA",
#      which compares unequal to every real digest and would trigger a spurious
#      "hash mismatch" error. Shared by ferx_sir() and ferx_covariance().
.ferx_hash_arg <- function(x) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x[[1L]]))) {
    ""
  } else {
    as.character(x)
  }
}

.ferx_state <- new.env(parent = emptyenv())
.ferx_state$last_trace_path  <- NULL
.ferx_state$last_trace_time  <- NULL
.ferx_state$last_trace_model <- NULL

# The worker-pool width a standalone ferx_covariance() / ferx_sir() step
# runs on (#527), as the `threads` the glue takes: 0 is the engine default,
# n > 0 pins. NULL takes the width the fit ran on, `fit$n_threads_used`, so
# a fit pinned to 2 threads is not refined on every core; a fit that records
# no width (NULL / NA, e.g. a hand-built list) gets the engine default.
.ferx_step_threads <- function(threads, fit) {
  if (is.null(threads)) {
    used <- fit$n_threads_used
    if (is.numeric(used) && length(used) == 1L && is.finite(used) && used >= 1) {
      return(as.integer(used))
    }
    return(0L)
  }
  # The range before the whole-number test, as in ferx_sir() (#536).
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
        threads < 0 || threads > .Machine$integer.max ||
        threads != floor(threads)) {
    stop("`threads` must be NULL or a non-negative integer scalar")
  }
  as.integer(threads)
}
