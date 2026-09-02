#' Population predictions from a NLME model
#'
#' Computes population-level predictions (eta = 0) for all subjects.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV. When omitted, the model file's
#'   \code{[data]} block (\code{path = ...}) is used.
#'   The \code{DV} column may be left empty (\code{.} / \code{NA}) on the
#'   sampling rows - \code{PRED} is what this function produces and the DV is
#'   never read, so an empty cell means "predict here". Rows marked
#'   \code{MDV = 1} are excluded, as always.
#' @param fit Optional \code{ferx_fit} result. When provided, predictions use
#'   \code{fit$theta} instead of the model file's initial estimate for theta.
#'
#' @return A data.frame with columns: ID, TIME, PRED
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' preds <- ferx_predict(ex$model, ex$data, fit = fit)
#' head(preds)
#'
#' @family simulation
#' @export
ferx_predict <- function(model, data = NULL, fit = NULL) {
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  stopifnot(file.exists(model), file.exists(data))

  if (is.null(fit)) {
    return(ferx_rust_predict(
      model_path = normalizePath(model),
      data_path = normalizePath(data)
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_predict_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma,
    omega_iov_flat = fit_pieces$omega_iov_flat,
    omega_iov_dim = fit_pieces$omega_iov_dim
  )
}
