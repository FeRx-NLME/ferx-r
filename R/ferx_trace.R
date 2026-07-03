#' Read the optimizer trace from a ferx fit
#'
#' Reads the per-iteration trace CSV written when \code{optimizer_trace = TRUE}
#' was passed to \code{\link{ferx_fit}}. Returns a tidy data frame with one row
#' per optimizer iteration (or per OFV evaluation for NLopt-based methods).
#'
#' Called with no argument, returns the trace from the most recent
#' \code{ferx_fit()} call in the current R session and prints a one-line
#' message indicating when that fit was started. Named \code{ferx_trace()}
#' rather than \code{trace()} to avoid masking \code{base::trace()}.
#'
#' @param fit Optional. A \code{ferx_fit} object returned by
#'   \code{\link{ferx_fit}}, or a character string giving the path to a trace
#'   CSV file written by ferx. If omitted, uses the trace from the last
#'   \code{ferx_fit()} call in this session.
#'
#' @return A data frame with columns:
#'   \item{iter}{Iteration / evaluation index (integer)}
#'   \item{method}{Estimation method: "foce", "focei", "gn", "gn_hybrid", "saem"}
#'   \item{phase}{Sub-phase label (empty string for single-phase methods)}
#'   \item{ofv}{Objective function value at this iteration}
#'   \item{wall_ms}{Wall-clock milliseconds elapsed since the start of the fit}
#'   \item{grad_norm}{L2 norm of the gradient (\code{NA} for gradient-free optimizers)}
#'   \item{step_norm}{L2 norm of the parameter step (\code{NA} when unavailable)}
#'   \item{inner_iter_count}{\code{NA} (reserved for future inner-loop counts)}
#'   \item{optimizer}{Optimizer name string, e.g. "slsqp", "bobyqa", "bfgs"}
#'   \item{lm_lambda}{Levenberg-Marquardt damping factor (Gauss-Newton only)}
#'   \item{ofv_delta}{OFV change from the previous accepted step (Gauss-Newton only)}
#'   \item{step_accepted}{1 if the GN step was accepted, 0 if rejected (\code{NA} otherwise)}
#'   \item{cond_nll}{Conditional negative log-likelihood (SAEM only)}
#'   \item{gamma}{SAEM step-size schedule gamma}
#'   \item{mh_accept_rate}{Metropolis-Hastings acceptance rate (SAEM only)}
#'
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE,
#'                 optimizer_trace = TRUE)
#' tr  <- ferx_trace(fit)
#' head(tr)
#'
#' @family diagnostics
#' @export
ferx_trace <- function(fit) {
  if (missing(fit) || is.null(fit)) {
    path <- .ferx_state$last_trace_path
    if (is.null(path)) {
      stop("No trace from a previous ferx_fit() call found in this session. ",
           "Run ferx_fit(..., optimizer_trace = TRUE) first, or pass a fit ",
           "object / trace path explicitly.")
    }
    started   <- .ferx_state$last_trace_time
    model_lbl <- .ferx_state$last_trace_model
    msg <- sprintf("Last ferx_fit() started at %s",
                   format(started, "%Y-%m-%d %H:%M:%S"))
    if (!is.null(model_lbl)) {
      msg <- paste0(msg, " (model: ", basename(model_lbl), ")")
    }
    message(msg)
  } else if (is.character(fit)) {
    path <- fit
  } else if (inherits(fit, "ferx_fit")) {
    # Prefer the trace already stored on the fit object - it survives a
    # save/load round-trip and doesn't depend on a temp file still existing.
    # Exact-match access (`[[`, not `$`) - a fake/manual fit object with
    # only `trace_path` set would otherwise silently partial-match `$trace`
    # to `trace_path` (a string, not a data frame).
    if (!is.null(fit[["trace"]])) {
      return(fit[["trace"]])
    }
    path <- fit$trace_path
    if (is.null(path)) {
      stop("No trace path found in fit object. ",
           "Re-run with optimizer_trace = TRUE.")
    }
  } else {
    stop("`fit` must be a ferx_fit object or a path to a trace CSV")
  }

  if (!file.exists(path)) {
    stop("Trace file not found: ", path)
  }

  .ferx_read_trace_csv(path)
}

# Internal: read + normalize a trace CSV. Shared between ferx_trace() (reading
# from a live trace_path) and ferx_fit() (populating fit$trace right after
# the fit completes, while the temp file still exists).
.ferx_read_trace_csv <- function(path) {
  tr <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  na_cols <- c("grad_norm", "step_norm", "inner_iter_count",
               "lm_lambda", "ofv_delta", "step_accepted",
               "cond_nll", "gamma", "mh_accept_rate")
  for (col in na_cols) {
    if (col %in% names(tr)) {
      v <- tr[[col]]
      if (is.character(v)) {
        tr[[col]] <- suppressWarnings(as.numeric(v))
      }
    }
  }

  tr
}
