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
#' \dontrun{
#' fit <- ferx_fit("warfarin.ferx", "warfarin.csv", optimizer_trace = TRUE)
#' head(ferx_trace(fit))
#' head(ferx_trace())              # same trace, recalled from last run
#' }
#'
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
#' \dontrun{
#' fit <- ferx_fit("warfarin.ferx", "warfarin.csv", optimizer_trace = TRUE)
#' ferx_plot_trace(fit)
#' }
#'
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

#' Quick convergence check with an optimizer trace
#'
#' Runs a short pilot fit (5 or 20 outer iterations depending on the method)
#' with \code{optimizer_trace = TRUE} so you can inspect how the optimizer
#' moves from the initial parameter values before committing to a full fit.
#' Useful for diagnosing poor starting values, ill-scaled parameters, or
#' structural model issues.
#'
#' @param model Path to a \code{.ferx} model file, or a \code{ferx_model}
#'   object created by \code{\link{ferx_model}}. When a \code{ferx_model} is
#'   passed and \code{data} is not supplied, the data path on the object is
#'   used.
#' @param data  Path to a NONMEM-format CSV file. Optional when \code{model}
#'   is a \code{ferx_model} that already carries a data path.
#' @param method Estimation method string, passed to \code{\link{ferx_fit}}.
#'   Default \code{"focei"}.
#' @param maxiter Maximum iterations for the pilot fit. Default: 5 for
#'   gradient-based methods, 20 for SAEM (to get a meaningful accept-rate
#'   trajectory).
#' @param ... Additional arguments forwarded verbatim to
#'   \code{\link{ferx_fit}} (e.g. \code{threads}, \code{settings}).
#'
#' @return A named list with:
#'   \item{fit}{The \code{ferx_fit} object from the pilot run}
#'   \item{trace}{The trace data frame (from \code{\link{ferx_trace}})}
#'   \item{summary}{A one-row data frame with \code{n_iter}, \code{ofv_start},
#'     \code{ofv_end}, \code{ofv_drop}, and \code{converged}}
#'
#' @examples
#' \dontrun{
#' # Plain paths
#' chk <- ferx_check_init("warfarin.ferx", "warfarin.csv")
#' ferx_plot_trace(chk$fit)
#' chk$summary
#'
#' # Pipe from a ferx_model — data is picked up from the object
#' ex <- ferx_example("warfarin")
#' ferx_model(ex$model, data = ex$data) |>
#'   ferx_check_init(method = "focei")
#' }
#'
#' @export
ferx_check_init <- function(model, data = NULL, method = "focei", maxiter = NULL, ...) {
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model <- model$model
  }
  if (is.null(data)) {
    stop(
      "No data path supplied. Either pass `data` to ferx_check_init() ",
      "or include it in ferx_model()."
    )
  }
  if (is.null(maxiter)) {
    maxiter <- if (tolower(method) == "saem") 20L else 5L
  }
  dots <- list(...)
  # maxiter flows through settings in ferx_fit (not a dedicated arg); merge
  # with any user-supplied settings, letting maxiter win.
  user_settings <- dots$settings
  dots$settings <- c(list(maxiter = as.integer(maxiter)), user_settings)
  dots$settings[["maxiter"]] <- as.integer(maxiter)
  dots$settings <- dots$settings[!duplicated(names(dots$settings))]

  fit <- do.call(ferx_fit, c(
    list(model           = model,
         data            = data,
         method          = method,
         covariance      = FALSE,
         verbose         = FALSE,
         optimizer_trace = TRUE),
    dots
  ))
  tr <- ferx_trace(fit)
  ofv_start <- tr$ofv[1L]
  ofv_end   <- tr$ofv[nrow(tr)]
  summary_df <- data.frame(
    n_iter     = nrow(tr),
    ofv_start  = ofv_start,
    ofv_end    = ofv_end,
    ofv_drop   = ofv_start - ofv_end,
    converged  = isTRUE(fit$converged),
    stringsAsFactors = FALSE
  )
  list(fit = fit, trace = tr, summary = summary_df)
}
