#' Plot optimizer trace from a ferx fit
#'
#' Produces a two-panel diagnostic plot from the per-iteration trace written
#' when \code{optimizer_trace = TRUE} was passed to \code{\link{ferx_fit}}.
#' The top panel shows OFV over iterations; the bottom panel shows a
#' method-specific convergence metric (gradient norm for gradient-based
#' methods, MH accept rate for SAEM, LM lambda for Gauss-Newton). Phase
#' boundaries are drawn as vertical dashed lines.
#'
#' @param fit A \code{ferx_fit} object or path to a trace CSV file (see
#'   \code{\link{ferx_trace}}).
#' @param log_ofv Logical; plot OFV on a log scale relative to the final value
#'   (\eqn{OFV - OFV_{final}}). Default \code{FALSE}.
#'
#' @return Invisibly returns the trace data frame (from
#'   \code{\link{ferx_trace}}).
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
ferx_plot_trace <- function(fit, log_ofv = FALSE) {
  tr <- ferx_trace(fit)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2L, 1L), mar = c(3, 4, 2, 1))

  methods_present <- unique(tr$method)
  is_saem   <- any(methods_present == "saem")
  is_gn     <- any(grepl("^gn", methods_present))

  # --- panel 1: OFV ---
  ofv_vals <- tr$ofv
  y_label  <- "OFV"
  if (isTRUE(log_ofv)) {
    ofv_vals <- ofv_vals - min(ofv_vals, na.rm = TRUE)
    y_label  <- "OFV - min(OFV)"
  }
  plot(tr$iter, ofv_vals,
       type = "l", col = "steelblue", lwd = 1.5,
       xlab = "", ylab = y_label,
       main = "Optimizer trace - OFV")

  # phase boundaries (vertical lines)
  .add_phase_lines(tr)

  # --- panel 2: method-specific metric ---
  if (is_saem) {
    metric   <- tr$mh_accept_rate
    m_label  <- "MH accept rate"
    m_color  <- "darkorange"
  } else if (is_gn) {
    metric   <- tr$lm_lambda
    m_label  <- "LM lambda"
    m_color  <- "darkgreen"
  } else {
    metric   <- tr$grad_norm
    m_label  <- "Gradient norm"
    m_color  <- "firebrick"
  }

  has_metric <- !all(is.na(metric))
  if (has_metric) {
    plot(tr$iter, metric,
         type = "l", col = m_color, lwd = 1.5,
         xlab = "Iteration", ylab = m_label,
         main = paste("Optimizer trace -", m_label))
    .add_phase_lines(tr)
  } else {
    plot.new()
    mtext("(metric not available for this method)", side = 3)
  }

  invisible(tr)
}

# Internal: add vertical dashed lines at phase-transition rows
.add_phase_lines <- function(tr) {
  if (!"phase" %in% names(tr)) return(invisible())
  phase_vals <- tr$phase
  if (all(phase_vals == "" | is.na(phase_vals))) return(invisible())
  # Find rows where the phase label changes
  rle_phases <- rle(phase_vals)
  boundaries <- cumsum(rle_phases$lengths)
  boundaries <- boundaries[-length(boundaries)]  # not the final boundary
  if (length(boundaries) > 0L) {
    abline(v = tr$iter[boundaries], lty = 2L, col = "grey50")
  }
  invisible()
}
