#' Correlation matrix of estimated parameters
#'
#' Converts the parameter covariance matrix (already stored on the fit object)
#' into a correlation matrix, making off-diagonal structure immediately visible.
#' A correlation close to \eqn{\pm 1} between two parameters flags a structural
#' identifiability problem in the model.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @return Named correlation matrix (invisibly). Printed to the console.
#' @seealso \code{\link{ferx_estimates}} for SEs and \%RSE.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = TRUE)
#' if (!is.null(fit$cov_matrix)) ferx_cor_matrix(fit)
#' @family diagnostics
#' @export
ferx_cor_matrix <- function(fit) {
  if (is.null(fit$cov_matrix)) {
    stop(
      "Covariance matrix not available. ",
      "Run ferx_fit() with covariance = TRUE and check fit$covariance_status."
    )
  }
  se <- sqrt(diag(fit$cov_matrix))
  if (any(se <= 0, na.rm = TRUE)) {
    warning("One or more diagonal elements are non-positive; ",
            "correlation matrix may not be meaningful.")
    se[se <= 0] <- NA_real_
  }
  cor_mat <- fit$cov_matrix / outer(se, se)
  # Clip to [-1, 1] for numerical noise on the diagonal
  cor_mat[cor_mat >  1] <-  1
  cor_mat[cor_mat < -1] <- -1
  # A correlation matrix has an exactly-unit diagonal by definition. Force it:
  # cov[i,i] / sqrt(cov[i,i])^2 can drift to 0.9999999999999999, which the `> 1`
  # clip above does not catch. Leave NA where a non-positive variance was
  # flagged above (those rows/cols are already NA and warned).
  ok <- !is.na(diag(cor_mat))
  diag(cor_mat)[ok] <- 1
  print(round(cor_mat, 3))
  invisible(cor_mat)
}
