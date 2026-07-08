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

# TRUE when a fit's estimation used FOCEI (interaction), FALSE for FOCE. The
# engine never plumbs `fit$interaction` to R (fit_result_to_list() omits it, so
# `isTRUE(fit$interaction)` is always FALSE), so we derive it from the method
# chain. Mirrors the engine's rule: the interaction flag follows the last
# *estimating* stage, and a trailing IMP (importance-sampling) stage is
# diagnostic-only and skipped (ferx-core: chain.iter().rev().find(|m| **m != Imp)).
# So c("focei", "imp") is interaction = TRUE.
.ferx_fit_interaction <- function(fit) {
  mchain <- toupper(as.character(fit$method_chain %||% fit$method))
  if (length(mchain) == 0L) return(FALSE)
  estimating <- mchain[mchain != "IMP"]
  last_method <- if (length(estimating) > 0L) {
    estimating[length(estimating)]
  } else {
    mchain[length(mchain)]
  }
  identical(last_method, "FOCEI")
}

.ferx_state <- new.env(parent = emptyenv())
.ferx_state$last_trace_path  <- NULL
.ferx_state$last_trace_time  <- NULL
.ferx_state$last_trace_model <- NULL
