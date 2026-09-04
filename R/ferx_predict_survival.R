#' Survival-function predictions for time-to-event endpoints
#'
#' For a model with one or more \code{[event_model]} (TTE) endpoints, evaluates
#' the survival function \eqn{S(t)}, cumulative hazard \eqn{H(t)} and hazard
#' \eqn{h(t)} on a user-supplied time grid, for every subject and every TTE CMT.
#' Population-level (\eqn{\eta = 0}) unless a \code{fit} is supplied, in which
#' case its fitted \code{theta} is used. Each subject also carries its median and
#' mean survival time (analytic where available).
#'
#' For a model with multiple TTE CMTs (competing risks), each row also carries the
#' cause-specific cumulative incidence \code{cif} (\eqn{F_k(t)}) and the all-cause
#' survival \code{survival_all} (\eqn{S_\mathrm{all}(t)}), with
#' \eqn{\sum_k F_k(t) + S_\mathrm{all}(t) = 1}. For a single endpoint
#' \code{cif} \eqn{= 1 -} \code{survival} and \code{survival_all} \eqn{=}
#' \code{survival}.
#'
#' @param model Path to a .ferx model file containing at least one
#'   \code{[event_model]} block.
#' @param data Path to a NONMEM-format CSV. When omitted, the model file's
#'   \code{[data]} block (\code{path = ...}) is used.
#' @param times Numeric vector of times at which to evaluate \eqn{S(t)},
#'   \eqn{H(t)}, \eqn{h(t)}.
#' @param fit Optional \code{ferx_fit} result. When provided, predictions use
#'   \code{fit$theta} instead of the model file's initial theta.
#'
#' @return A data.frame with columns: ID, CMT, TIME, survival, cum_hazard,
#'   hazard, cif, survival_all, median_survival, mean_survival (one row per
#'   subject \eqn{\times} TTE CMT \eqn{\times} time).
#'
#' @examples
#' \dontrun{
#' # `model` must contain at least one [event_model] (TTE) block.
#' preds <- ferx_predict_survival(
#'   "tte_model.ferx", "tte_data.csv",
#'   times = seq(0, 24, by = 2)
#' )
#' head(preds)
#' }
#'
#' @family simulation
#' @export
ferx_predict_survival <- function(model, data = NULL, times, fit = NULL) {
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  stopifnot(
    file.exists(model),
    file.exists(data),
    is.numeric(times),
    length(times) > 0,
    all(is.finite(times))
  )
  times <- as.numeric(times)

  if (is.null(fit)) {
    return(ferx_rust_predict_survival(
      model_path = normalizePath(model),
      data_path = normalizePath(data),
      times = times
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_predict_survival_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    times = times,
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma,
    omega_iov_flat = fit_pieces$omega_iov_flat,
    omega_iov_dim = fit_pieces$omega_iov_dim,
    residual_rho = fit_pieces$residual_rho
  )
}
