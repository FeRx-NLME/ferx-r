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

# A Rust `Option<bool>`. Both directions treat NA as absent: a tri-state flag
# reaches R as NA when the engine recorded no verdict, and writing NA back
# would claim a verdict of "false" to any reader that coerces it.
.fitrx_opt_lgl <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.logical(x))
  if (length(v) != 1L || is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_lgl <- function(x) {
  if (is.null(x)) return(NA)
  v <- suppressWarnings(as.logical(unlist(x, use.names = FALSE)))
  if (length(v) != 1L || is.na(v)) return(NA)
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
  result$eta_cov    <- .ferx_compute_eta_cov(result$ebe_etas, result$data_path,
                                           result$model_path)
  result
}

# -- Parameter priors (ferx-core #254) ---------------------------------------
#
# `PriorSummary` on the wire is a plain array of objects with eight required
# fields (no serde defaults), so the write path emits every column and the read
# path tolerates a missing one only by filling NA. The pair is written *only*
# for a priored fit, matching ferx-core's own writer: absent means "this bundle
# predates priors", and the loader then reconstructs `ofv_data = ofv`, which is
# right for every such file because no prior could have been applied.

.fitrx_prior_summary_to_wire <- function(ps) {
  if (!is.data.frame(ps) || nrow(ps) == 0L) return(NULL)
  lapply(seq_len(nrow(ps)), function(i) {
    list(
      name               = as.character(ps$name[[i]]),
      prior_value        = as.numeric(ps$prior_value[[i]]),
      estimate           = as.numeric(ps$estimate[[i]]),
      shift_in_prior_sds = as.numeric(ps$shift_in_prior_sds[[i]]),
      penalty            = as.numeric(ps$penalty[[i]]),
      family             = as.character(ps$family[[i]]),
      prior_lower_95     = as.numeric(ps$prior_lower_95[[i]]),
      prior_upper_95     = as.numeric(ps$prior_upper_95[[i]])
    )
  })
}

.fitrx_prior_summary_from_wire <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  num <- function(el, k) {
    v <- suppressWarnings(as.numeric(el[[k]] %||% NA_real_))
    if (length(v) != 1L) NA_real_ else v
  }
  chr <- function(el, k) {
    v <- as.character(el[[k]] %||% NA_character_)
    if (length(v) != 1L) NA_character_ else v
  }
  data.frame(
    name               = vapply(x, chr, character(1L), "name"),
    prior_value        = vapply(x, num, numeric(1L), "prior_value"),
    estimate           = vapply(x, num, numeric(1L), "estimate"),
    shift_in_prior_sds = vapply(x, num, numeric(1L), "shift_in_prior_sds"),
    penalty            = vapply(x, num, numeric(1L), "penalty"),
    family             = vapply(x, chr, character(1L), "family"),
    prior_lower_95     = vapply(x, num, numeric(1L), "prior_lower_95"),
    prior_upper_95     = vapply(x, num, numeric(1L), "prior_upper_95"),
    stringsAsFactors   = FALSE
  )
}

# Recover the objective's data / prior split from a bundle that does not carry
# it, returning list(ofv_data =, ofv_prior =).
#
# ferx-core's own loader reads a missing split as "this file predates priors",
# and for a bundle *it* wrote that is sound: once the feature existed, its
# writer always emitted the split for a priored fit. It is NOT sound for a
# bundle this package wrote. `ferx_save_fit()` shipped before ferx-r #366 and
# emitted neither half, while R could already fit a priored model - so an
# R-written bundle of a priored fit carries a penalized `ofv` and nothing to
# say so. Reading that as the data half relabels the penalized objective as the
# likelihood, breaks the `aic == ofv_data + 2k` invariant the stored AIC was
# computed under, and hands `ferx_sir()` back the #366 double count.
#
# The split is recoverable from what such a bundle does carry: the engine
# computes `aic = ofv_data + 2 * n_parameters` for every fit, priored or not
# (ferx-core api/fit.rs), and both `aic` and `n_parameters` are on the wire. So
# `ofv_data = aic - 2 * n_parameters` and `ofv_prior = ofv - ofv_data`.
#
# Guarded, because the identity is only as good as the two fields it reads. A
# penalty is a sum of squares, so a recovered prior half below zero means the
# identity does not hold for this file (a hand-edited bundle, a writer that
# computed AIC differently) and the unpriored reading is restored. A recovered
# half within rounding distance of zero is taken as exactly zero, which is what
# the overwhelmingly common unpriored bundle should report.
.fitrx_recover_ofv_split <- function(w) {
  ofv <- suppressWarnings(as.numeric(w$ofv %||% NA_real_))
  unpriored <- list(ofv_data = ofv, ofv_prior = 0, recovered = FALSE)
  if (length(ofv) != 1L || !is.finite(ofv)) return(unpriored)

  aic <- suppressWarnings(as.numeric(w$aic %||% NA_real_))
  k <- suppressWarnings(as.numeric(w$n_parameters %||% NA_real_))
  if (length(aic) != 1L || !is.finite(aic) ||
        length(k) != 1L || !is.finite(k) || k < 0) {
    return(unpriored)
  }

  ofv_data <- aic - 2 * k
  ofv_prior <- ofv - ofv_data
  if (!is.finite(ofv_prior)) return(unpriored)
  # The JSON carries full f64 precision, so the identity round-trips to within
  # rounding; anything below this is noise, not a prior.
  tol <- 1e-9 * max(1, abs(ofv))
  if (ofv_prior < -tol) return(unpriored)
  if (ofv_prior <= tol) return(unpriored)

  list(ofv_data = ofv_data, ofv_prior = ofv_prior, recovered = TRUE)
}
