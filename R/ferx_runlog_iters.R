#' Print the full optimizer iteration history for a ferx fit
#'
#' Displays (or returns) the complete per-iteration table from the optimizer
#' trace CSV.  Unlike the "Iteration history" section in
#' \code{\link{ferx_runlog}}, this function never truncates: every iteration
#' row is shown without omission.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} with
#'   \code{optimizer_trace = TRUE}, or a path to a trace CSV file.
#' @param verbose Logical.  When \code{TRUE} (the default) the table is
#'   printed.  When \code{FALSE} it is returned invisibly as a character string.
#'
#' @return A character string (invisibly when \code{verbose = TRUE}).
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "focei",
#'                 covariance = FALSE, optimizer_trace = TRUE)
#' ferx_runlog_iters(fit)
#' }
#'
#' @seealso \code{\link{ferx_runlog}}, \code{\link{ferx_trace}},
#'   \code{\link{plot.ferx_fit}}
#' @family diagnostics
#' @export
ferx_runlog_iters <- function(fit, verbose = TRUE) {
  tr <- if (is.character(fit)) {
    ferx_trace(fit)
  } else if (inherits(fit, "ferx_fit")) {
    if (is.null(fit$trace_path) || is.na(fit$trace_path) ||
        !nzchar(fit$trace_path) || !file.exists(fit$trace_path)) {
      stop("No trace file found. Re-run with optimizer_trace = TRUE.")
    }
    ferx_trace(fit)
  } else {
    stop("`fit` must be a ferx_fit object or a path to a trace CSV")
  }

  hdr  <- c(strrep("=", 60), "ITERATION HISTORY", strrep("-", 60))
  body <- .runlog_iter_table(tr, truncate = FALSE)
  foot <- c(sprintf("  Total: %d iteration(s)", nrow(tr)), strrep("=", 60))
  out  <- paste(c(hdr, body, foot), collapse = "\n")
  if (isTRUE(verbose)) cat(out, "\n")
  invisible(out)
}
