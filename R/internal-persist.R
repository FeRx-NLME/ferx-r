# Write-path (ferx_save_fit) and read-path (ferx_load_fit) scalar/vector
# NULL-or-value coercion helpers. Kept side by side (rather than split across
# ferx_save_fit.R / ferx_load_fit.R) since they are mirror images of each
# other and a change to one direction's coercion rule usually needs checking
# against the other. Not merged into single functions: the two directions
# have different semantics (write coerces R values for JSON serialisation,
# read coerces JSON-decoded values back, including its own NA/length quirks),
# so the small differences between e.g. .fitrx_opt_num and
# .fitrx_unwrap_opt_num are intentional, not accidental drift.

.fitrx_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.numeric(x)
  if (is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
  v
}

.fitrx_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
  v
}

.fitrx_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.character(x)
  if (is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.character(x)
  if (length(v) == 0L || is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(x)
  if (length(v) == 0L) return(NULL)
  v
}

.fitrx_unwrap_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(unlist(x, use.names = FALSE))
  if (length(v) == 0L) return(NULL)
  v
}

# Fields derived from other fields already present on `result` (correlation
# matrix, tidy estimates table, eta-covariate correlations). Shared between
# ferx_fit() and ferx_load_fit() - the two fit-construction entry points -
# so they can't drift on how a derived field is computed. Requires
# `result$data_path` to already be set (normalised, in ferx_fit()'s case).
.ferx_populate_derived_fields <- function(result) {
  result$cor_matrix <- .ferx_compute_cor_matrix(result$cov_matrix)
  result$estimates  <- .ferx_compute_estimates(result)
  result$eta_cov    <- .ferx_compute_eta_cov(result$ebe_etas, result$data_path)
  result
}
