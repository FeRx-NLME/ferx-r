#' Structured diagnostic flags from a fit result
#'
#' Returns a named list of diagnostic data frames covering IWRES autocorrelation
#' and parameter shrinkage. Designed to be inspected programmatically or used
#' as the basis for custom plots.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return A named list with elements:
#'   \describe{
#'     \item{autocorrelation}{1-row data frame: \code{dw_statistic},
#'       \code{lag1_r}, \code{flag} (character interpretation of the DW value).
#'       \code{NULL} when no IWRES autocorrelation data is available.}
#'     \item{shrinkage}{Data frame with columns \code{param}, \code{type}
#'       (\code{"eta"} or \code{"eps"}), \code{shrinkage} (proportion 0-1),
#'       \code{shrinkage_pct} (percentage). \code{NULL} when no shrinkage data
#'       is available.}
#'   }
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' diag <- check_diagnostics(fit)
#' diag$autocorrelation
#' diag$shrinkage
#' }
#' @family diagnostics
#' @export
check_diagnostics <- function(fit) {
  # Autocorrelation
  dw  <- fit$dw_statistic
  r   <- fit$iwres_lag1_r
  autocorr <- if (!is.null(dw) && !is.na(dw)) {
    data.frame(dw_statistic = dw, lag1_r = r %||% NA_real_,
               flag = .dw_label(dw),
               stringsAsFactors = FALSE)
  } else NULL

  # Shrinkage
  rows <- list()
  if (!is.null(fit$shrinkage_eta)) {
    eta_lbls <- if (!is.null(fit$eta_names) && length(fit$eta_names) == length(fit$shrinkage_eta))
      fit$eta_names else sprintf("ETA%d", seq_along(fit$shrinkage_eta))
    for (i in seq_along(fit$shrinkage_eta)) {
      sh <- fit$shrinkage_eta[i]
      if (!is.na(sh)) {
        rows[[length(rows) + 1L]] <- data.frame(
          param = eta_lbls[i], type = "eta",
          shrinkage = sh, shrinkage_pct = sh * 100,
          stringsAsFactors = FALSE)
      }
    }
  }
  eps_vec <- fit$shrinkage_eps
  if (!is.null(eps_vec)) {
    eps_vec <- eps_vec[!is.na(eps_vec)]
    if (length(eps_vec) > 0L) {
      eps_lbls <- if (!is.null(fit$sigma_names) && length(fit$sigma_names) == length(fit$shrinkage_eps))
        fit$sigma_names[!is.na(fit$shrinkage_eps)]
      else if (length(eps_vec) == 1L)
        "EPS"
      else
        sprintf("EPS%d", seq_along(eps_vec))
      for (i in seq_along(eps_vec)) {
        rows[[length(rows) + 1L]] <- data.frame(
          param = eps_lbls[i], type = "eps",
          shrinkage = eps_vec[i], shrinkage_pct = eps_vec[i] * 100,
          stringsAsFactors = FALSE)
      }
    }
  }
  shrinkage_df <- if (length(rows) > 0L) do.call(rbind, rows) else NULL

  list(autocorrelation = autocorr, shrinkage = shrinkage_df)
}
