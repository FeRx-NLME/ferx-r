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

# Read back a JSON array that may contain nulls - a Rust `Vec<Option<T>>` such
# as the per-kappa `weight =` source text and its typical value. `unlist()`
# would silently *drop* the nulls and shift every later element onto the wrong
# kappa, so map them to NA element-wise instead. NULL when absent or empty.
.fitrx_unwrap_nullable_vec <- function(x, na) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  if (!is.list(x)) return(c(x, na)[seq_along(x)])
  vapply(x, function(el) {
    if (is.null(el) || length(el) == 0L) na else c(el, na)[[1L]]
  }, na)
}

.fitrx_unwrap_opt_chr_vec <- function(x) .fitrx_unwrap_nullable_vec(x, NA_character_)

.fitrx_unwrap_nullable_num_vec <- function(x) .fitrx_unwrap_nullable_vec(x, NA_real_)

.fitrx_opt_lgl <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.logical(x))
  if (length(v) != 1L || is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_lgl <- function(x) {
  if (is.null(x)) return(NA)
  v <- suppressWarnings(as.logical(x))
  if (length(v) != 1L) return(NA)
  v
}

# BicInputs <-> wire object (ferx-core #1177). The engine's struct is all
# integers plus one flag, and it reads the object back with `serde(default)`,
# so an absent field on either side degrades to zero rather than erroring.
.FITRX_BIC_INPUT_COUNTS <- c(
  "n_obs", "theta_random", "theta_fixed", "omega", "kappa", "sigma"
)

.fitrx_bic_inputs_to_wire <- function(x) {
  if (is.null(x)) return(NULL)
  out <- lapply(
    stats::setNames(nm = .FITRX_BIC_INPUT_COUNTS),
    function(k) as.integer(x[[k]] %||% 0L)
  )
  out$sigma_random <- isTRUE(x$sigma_random)
  out
}

.fitrx_bic_inputs_from_wire <- function(w) {
  if (is.null(w)) return(NULL)
  out <- lapply(
    stats::setNames(nm = .FITRX_BIC_INPUT_COUNTS),
    function(k) as.integer(w[[k]] %||% 0L)
  )
  out$sigma_random <- isTRUE(w$sigma_random)
  out
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
  # `ferx_fit()` gets the BIC variants straight from the engine; a fit rebuilt
  # from a `.fitrx` bundle carries only the tally they are computed from.
  result <- .ferx_populate_bic_variants(result)
  result
}
