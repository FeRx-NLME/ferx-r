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
      "ferx: built WITHOUT autodiff (FERX_NO_AUTODIFF=1).\n",
      "All estimation methods (FOCE/FOCEI/SAEM/IMP) work with finite-difference gradients.\n",
      "Limitations vs. the full Enzyme build:\n",
      "  - gradient = ad in model files falls back to fd (slower, may be less accurate)\n",
      "  - SAEM HMC proposals (n_leapfrog > 0) are unavailable; MH is used instead\n",
      "  - saem_n_subjects_hmc is always NA\n",
      "To enable autodiff: reinstall without FERX_NO_AUTODIFF=1 using the Enzyme toolchain.\n",
      "See: https://ferx-nlme.github.io/ferx-core/installation.html"
    )
  }
}
