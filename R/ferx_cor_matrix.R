# Correlation matrix of estimated parameters, derived from `cov_matrix`.
# Stored on the fit object as `fit$cor_matrix` (NULL when the covariance step
# was not run or failed). A correlation close to +/-1 between two parameters
# flags a structural identifiability problem in the model. (Formerly the
# exported ferx_cor_matrix(fit); see issue #226.)
#
# `fixed` is `fit$cov_fixed`: the engine's mask of the coordinates the
# optimizer held (FIX, and the structural zeros of a mixed block Omega). The
# covariance step gives a held coordinate an all-zero row and column, so it has
# no correlation to report: its row and column are NA, and that is not a
# warning - every MBMA fit holds `sigma ... FIX` (#424). The warning is kept
# for what it was meant for, a non-positive variance on a parameter that WAS
# estimated.
.ferx_compute_cor_matrix <- function(cov_matrix, fixed = NULL) {
  if (is.null(cov_matrix)) return(NULL)
  v <- diag(cov_matrix)
  held <- .ferx_cov_held(cov_matrix, fixed)
  nms <- rownames(cov_matrix)
  label <- function(i) {
    paste(if (is.null(nms)) paste0("#", which(i)) else nms[i], collapse = ", ")
  }
  # The covariance step writes a held coordinate as an all-zero row. A held
  # row that is not zero means the mask and the matrix disagree; it stays NA
  # (the mask is the engine's word on what was estimated), but not silently.
  odd <- held & !.ferx_cov_zero_rows(cov_matrix)
  if (any(odd)) {
    warning("cov_matrix has non-zero entries for held parameter(s) ", label(odd),
            "; their correlations are NA.", call. = FALSE)
  }
  bad <- !held & !is.na(v) & v <= 0
  if (any(bad)) {
    warning("Non-positive variance for estimated parameter(s) ", label(bad),
            " in cov_matrix; their correlations are NA.", call. = FALSE)
  }
  # Test before sqrt(): a negative variance would otherwise come back NaN with
  # R's own "NaNs produced" warning and slip past the check above.
  se <- rep(NA_real_, length(v))
  ok <- !held & !bad & !is.na(v)
  se[ok] <- sqrt(v[ok])
  cor_mat <- cov_matrix / outer(se, se)
  # Clip to [-1, 1] for numerical noise on the diagonal
  cor_mat[!is.na(cor_mat) & cor_mat >  1] <-  1
  cor_mat[!is.na(cor_mat) & cor_mat < -1] <- -1
  # A correlation matrix has an exactly-unit diagonal by definition. Force it:
  # cov[i,i] / sqrt(cov[i,i])^2 can drift to 0.9999999999999999, which the `> 1`
  # clip above does not catch. Held and non-positive rows stay NA.
  diag(cor_mat)[ok] <- 1
  cor_mat
}

# Which rows of `cov_matrix` the engine held. The engine's own mask when the
# fit carries one of the right length. Without one - a bundle written before
# #424, or by a writer that does not record it - a row that is exactly zero
# is read as held, since that is how the covariance step writes a
# coordinate it did not estimate; an estimated parameter does not come out
# with an exactly-zero row.
.ferx_cov_held <- function(cov_matrix, fixed = NULL) {
  if (length(fixed) == nrow(cov_matrix)) return(as.logical(unname(fixed)) %in% TRUE)
  .ferx_cov_zero_rows(cov_matrix)
}

# Rows of `cov_matrix` that are exactly zero. Rows only: a covariance matrix
# is symmetric, so the matching column is zero too.
.ferx_cov_zero_rows <- function(cov_matrix) {
  unname(rowSums(!is.na(cov_matrix) & cov_matrix == 0) == ncol(cov_matrix))
}

# `fit$cov_fixed` from the engine's mask: a logical vector named like
# the rows of `cov_matrix`, or NULL when there is no matrix or the mask does
# not cover it.
.ferx_cov_fixed_named <- function(mask, cov_matrix) {
  if (is.null(cov_matrix)) return(NULL)
  mask <- as.logical(unlist(mask, use.names = FALSE))
  if (length(mask) != nrow(cov_matrix)) return(NULL)
  names(mask) <- rownames(cov_matrix)
  mask
}
