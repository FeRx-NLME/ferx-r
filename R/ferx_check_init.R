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
#'   is a \code{ferx_model} that already carries a data path, or when the model
#'   file declares a \code{[data]} block (\code{path = ...}).
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
#' ex  <- ferx_example("warfarin")
#'
#' # Check initial values before committing to a full FOCEI fit
#' chk <- ferx_check_init(ex$model, ex$data)
#'
#' # OFV drop confirms the gradient is pointing in the right direction
#' chk$summary   # ofv_start, ofv_end, ofv_drop
#'
#' # Inspect per-iteration trace to spot early divergence or slow descent
#' plot(chk$fit)
#'
#' # Pipe style: data flows in from the ferx_model object
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_check_init(method = "focei")
#' }
#'
#' @family fitting
#' @export
ferx_check_init <- function(model, data = NULL, method = "focei", maxiter = NULL, ...) {
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model <- model$model
  }
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
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
