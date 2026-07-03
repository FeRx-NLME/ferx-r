# Correlation matrix of estimated parameters, derived from `cov_matrix`.
# Stored on the fit object as `fit$cor_matrix` (NULL when the covariance step
# was not run or failed). A correlation close to +/-1 between two parameters
# flags a structural identifiability problem in the model. (Formerly the
# exported ferx_cor_matrix(fit); see issue #226.)
.ferx_compute_cor_matrix <- function(cov_matrix) {
  if (is.null(cov_matrix)) return(NULL)
  se <- sqrt(diag(cov_matrix))
  if (any(se <= 0, na.rm = TRUE)) {
    warning("One or more diagonal elements are non-positive; ",
            "correlation matrix may not be meaningful.")
    se[se <= 0] <- NA_real_
  }
  cor_mat <- cov_matrix / outer(se, se)
  # Clip to [-1, 1] for numerical noise on the diagonal
  cor_mat[cor_mat >  1] <-  1
  cor_mat[cor_mat < -1] <- -1
  # A correlation matrix has an exactly-unit diagonal by definition. Force it:
  # cov[i,i] / sqrt(cov[i,i])^2 can drift to 0.9999999999999999, which the `> 1`
  # clip above does not catch. Leave NA where a non-positive variance was
  # flagged above (those rows/cols are already NA and warned).
  ok <- !is.na(diag(cor_mat))
  diag(cor_mat)[ok] <- 1
  cor_mat
}
