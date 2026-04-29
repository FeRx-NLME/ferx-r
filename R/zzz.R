# Package-internal state: holds the trace path of the most recent ferx_fit()
# call (when optimizer_trace = TRUE), so ferx_trace() can be called with no
# argument to inspect the last run.
.ferx_state <- new.env(parent = emptyenv())
.ferx_state$last_trace_path  <- NULL
.ferx_state$last_trace_time  <- NULL
.ferx_state$last_trace_model <- NULL

.onAttach <- function(libname, pkgname) {
  enabled <- tryCatch(ferx_rust_autodiff_enabled(), error = function(e) NA)
  if (isFALSE(enabled)) {
    packageStartupMessage(
      "ferx: built WITHOUT autodiff. Gradient-based fits will be unavailable ",
      "or fall back to finite differences. Rebuild without FERX_NO_AUTODIFF=1 ",
      "(requires the Enzyme Rust toolchain) to enable autodiff."
    )
  }
}
