# Package-internal state: holds the trace path of the most recent ferx_fit()
# call (when optimizer_trace = TRUE), so ferx_trace() can be called with no
# argument to inspect the last run.
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

.ferx_state <- new.env(parent = emptyenv())
.ferx_state$last_trace_path  <- NULL
.ferx_state$last_trace_time  <- NULL
.ferx_state$last_trace_model <- NULL
