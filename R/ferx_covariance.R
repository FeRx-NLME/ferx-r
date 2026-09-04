#' Run the covariance step against an existing fit
#'
#' Run the finite-difference-Hessian covariance step against a fit that was
#' produced earlier - the covariance-step analogue of [ferx_sir()]. Useful
#' when the original fit was expensive and you want to add standard errors
#' without re-estimating, when you want to re-run the step with a different
#' `covariance_method` (e.g. the `rsr` sandwich), or when working with a fit
#' loaded from a `.fitrx` bundle.
#'
#' `ferx_covariance()` reconstructs the fitted parameters from the fit,
#' re-runs the inner loop (seeded from the per-subject empirical Bayes ETAs) to
#' rebuild the covariance-step inputs, and calls the same covariance
#' step [ferx_fit()] runs inline. The result closely matches fitting with
#' `covariance = TRUE` (the same engine step; agreement is close but not
#' bit-exact, since the standalone re-reads the data and cold-starts the inner
#' EBE loop). The returned fit is the input with
#' `cov_matrix`, `cor_matrix`, `se_theta` / `se_omega` / `se_sigma` /
#' `se_kappa`, `covariance_status`, `eigenvalues`, and `condition_number`
#' refreshed.
#'
#' A covariance step that runs but cannot produce a usable matrix (a
#' non-positive-definite or structurally-unusable FD Hessian) is **not** an
#' error: the returned fit reports `covariance_status = "failed"` with a
#' diagnostic appended to `warnings`, mirroring [ferx_fit()]. An error is
#' reserved for input problems (missing / hash-mismatched model or data).
#'
#' ## Integrity check
#'
#' Like [ferx_sir()], `ferx_covariance()` re-reads the model and data files
#' recorded on the fit (`fit$model_path`, `fit$data_path`) and verifies their
#' SHA-256 hashes (`fit$model_hash`, `fit$data_hash`); if either file changed
#' since the fit, the call is a **hard error**. Running the covariance step
#' against a modified model or dataset would be meaningless.
#'
#' Edge case: if hashing failed at fit time the corresponding `fit$*_hash` is
#' `NA` and the integrity check silently passes on that side. Check
#' `!is.na(fit$model_hash) && !is.na(fit$data_hash)` before relying on the
#' protection.
#'
#' @param fit A `ferx_fit` object produced by [ferx_fit()] or
#'   [ferx_load_fit()].
#' @param covariance_method Covariance estimator: `"r"` (inverse Hessian, the
#'   default, NONMEM `MATRIX=R`), `"s"` (inverse cross-product,
#'   `MATRIX=S`), or `"rsr"` (the Huber-White sandwich, `MATRIX=RSR`).
#'   `"hessian"` / `"cross_product"` / `"sandwich"` are accepted as aliases.
#' @param mu_referencing Use mu-referencing for the inner-loop warm restart.
#'   Default `TRUE` (the engine default). Set to match the setting used for
#'   the original fit.
#' @param verbose When `TRUE`, the engine prints progress to stderr.
#'   Default `FALSE`.
#'
#' @return The input `fit`, with `cov_matrix`, `cor_matrix`, `se_theta`,
#'   `se_omega`, `se_sigma`, `se_kappa`, `covariance_status`, `eigenvalues`,
#'   `condition_number`, `max_abs_correlation` (the largest absolute
#'   off-diagonal parameter correlation, `NA` when the step produced no
#'   matrix - the value [check_strictness()] gates on), and the `se` column of
#'   `residual_correlations` refreshed.
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#'
#' # Fit without the covariance step, then add it post-hoc.
#' fit <- ferx_fit(ex$model, ex$data, covariance = FALSE)
#' fit <- ferx_covariance(fit)
#' fit$se_theta       # standard errors now populated
#' fit$cov_matrix     # parameter covariance matrix
#'
#' # Re-run with the robust sandwich estimator instead.
#' fit_rsr <- ferx_covariance(fit, covariance_method = "rsr")
#' fit_rsr$se_theta
#' }
#'
#' @seealso [ferx_fit()] for the inline covariance step (`covariance = TRUE`),
#'   [ferx_sir()] for the SIR uncertainty step.
#' @family fitting
#' @export
ferx_covariance <- function(fit,
                            covariance_method = "r",
                            mu_referencing = TRUE,
                            verbose = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }

  # Normalise + validate covariance_method up-front so a typo surfaces as a
  # clean R condition rather than the engine's "unknown value" further in.
  if (!is.character(covariance_method) || length(covariance_method) != 1L ||
        is.na(covariance_method)) {
    stop("`covariance_method` must be a single string (\"r\", \"s\", or \"rsr\").")
  }
  cov_method <- tolower(covariance_method)
  valid_methods <- c("r", "hessian", "s", "cross_product", "rsr", "sandwich")
  if (!cov_method %in% valid_methods) {
    stop("`covariance_method` must be one of \"r\"/\"hessian\", ",
         "\"s\"/\"cross_product\", or \"rsr\"/\"sandwich\" (got: ",
         covariance_method, ").")
  }

  model_path <- fit$model_path
  if (is.null(model_path) || is.na(model_path) || !nzchar(model_path)) {
    stop(
      "ferx_covariance: fit has no recorded model_path; cannot locate the ",
      "model file. Re-fit via ferx_fit(model, data) so the path is recorded."
    )
  }
  if (!file.exists(model_path)) {
    stop("ferx_covariance: model file not found at ", model_path)
  }

  data_path <- fit$data_path
  if (is.null(data_path) || is.na(data_path) || !nzchar(data_path)) {
    stop(
      "ferx_covariance: fit has no recorded data_path; cannot locate the ",
      "data file. Re-fit via ferx_fit(model, data) so the path is recorded."
    )
  }
  if (!file.exists(data_path)) {
    stop("ferx_covariance: data file not found at ", data_path)
  }

  # Build the flat eta_hats matrix (warm-start for the inner loop). `ebe_etas`
  # is a data frame: first column ID, subsequent columns one per ETA. Strip ID
  # and flatten row-major so each subject contributes `n_eta` contiguous
  # values. Identical handling to ferx_sir().
  ebes <- fit$ebe_etas
  n_eta <- nrow(fit$omega)
  if (is.null(n_eta)) n_eta <- 0L

  if (n_eta == 0L) {
    # Fixed-effects-only (naive-pooled) fit - ferx-core #989. There is no inner
    # empirical-Bayes problem, so there are no EBEs to warm-start from and
    # `fit$ebe_etas` is NULL rather than an empty data frame. The covariance
    # step still applies: it is the theta/sigma block that carries the standard
    # errors, and that block is exactly what a naive-pooled user needs, since
    # NONMEM's `$COVARIANCE` default (`rsr`) is the estimator that accounts for
    # the within-subject correlation such a model deliberately ignores.
    eta_hats_flat <- numeric(0)
    # `nrow(eta_mat)` is the usual subject count, but there is no eta matrix
    # here; the engine still needs the count to size its warm-start scaffold.
    n_subj <- as.integer(fit$n_subjects %||% 0L)
    if (n_subj <= 0L) {
      stop(
        "ferx_covariance: fit has no random effects and no recorded ",
        "n_subjects, so the subject count cannot be determined. Re-fit via ",
        "ferx_fit(model, data)."
      )
    }
  } else {
    if (is.null(ebes) || nrow(ebes) == 0L) {
      stop(
        "ferx_covariance: fit$ebe_etas is empty, but fit$omega is ", n_eta, "x",
        n_eta, " so this fit should carry per-subject EBEs. Cannot warm-start ",
        "the inner loop."
      )
    }
    eta_cols <- setdiff(names(ebes), c("ID", "ofv_contribution", "n_obs"))
    if (length(eta_cols) != n_eta) {
      stop(
        "ferx_covariance: fit$ebe_etas has ", length(eta_cols), " ETA columns (",
        paste(eta_cols, collapse = ", "),
        ") but fit$omega is ", n_eta, "x", n_eta, ". ",
        "The EBE table and the omega matrix must agree on n_eta - was ",
        "this fit object hand-edited or assembled from incompatible parts?"
      )
    }
    eta_mat <- as.matrix(ebes[, eta_cols, drop = FALSE])
    storage.mode(eta_mat) <- "double"
    eta_hats_flat <- as.numeric(t(eta_mat))  # row-major
    n_subj <- nrow(eta_mat)
  }

  # Hash plumbing, identical to ferx_sir(): non-empty hex string forwards and
  # Rust enforces equality; NULL (older binary) or NA (hashing failed at fit
  # time) both coerce to "" so Rust skips the check on that side. See
  # `.ferx_hash_arg()` in zzz.R.
  model_hash_arg <- .ferx_hash_arg(fit$model_hash)
  data_hash_arg <- .ferx_hash_arg(fit$data_hash)
  if (!nzchar(model_hash_arg) || !nzchar(data_hash_arg)) {
    warning(
      "ferx_covariance: one or more file hashes are missing on the fit; ",
      "the integrity check will be skipped for the affected file(s). ",
      "Re-fit (or load a newer .fitrx) to enable hash verification."
    )
  }

  # FOCEI (interaction) vs FOCE controls the inner-loop NLL the covariance step
  # differentiates, so it must match the original fit. Derived from the method
  # chain because `fit$interaction` is not plumbed to R - see .ferx_fit_interaction().
  interaction <- .ferx_fit_interaction(fit)

  omega_flat <- as.numeric(t(fit$omega))

  # IOV omega (empty when no IOV). Forwarding the fitted matrix keeps IOV kappa
  # standard errors on the estimated scale; when absent the engine falls back
  # to the model-file init.
  if (!is.null(fit$omega_iov) && is.matrix(fit$omega_iov) && nrow(fit$omega_iov) > 0L) {
    omega_iov_flat <- as.numeric(t(fit$omega_iov))
    omega_iov_dim <- nrow(fit$omega_iov)
  } else {
    omega_iov_flat <- numeric(0)
    omega_iov_dim <- 0L
  }

  raw <- ferx_rust_covariance(
    model_path = model_path,
    data_path = data_path,
    model_hash = model_hash_arg,
    data_hash = data_hash_arg,
    ofv = as.numeric(fit$ofv),
    interaction = interaction,
    theta = as.numeric(fit$theta),
    omega_flat = omega_flat,
    omega_dim = nrow(fit$omega),
    sigma = as.numeric(fit$sigma),
    omega_iov_flat = omega_iov_flat,
    omega_iov_dim = as.integer(omega_iov_dim),
    residual_rho = .ferx_residual_rho_vec(fit),
    eta_hats_flat = eta_hats_flat,
    n_subjects = n_subj,
    covariance_method = cov_method,
    mu_referencing = isTRUE(mu_referencing),
    verbose = isTRUE(verbose)
  )

  # All error paths inside `ferx_rust_covariance` throw an R condition (via
  # `throw_r_error`) and never return `NULL`, so we don't test for that here.

  # --- Merge the refreshed covariance fields onto the fit --------------------
  # The reshaping mirrors the covariance block in ferx_fit()'s post-processing,
  # reading parameter names off the (already-finalised) fit object so a fit and
  # a run_covariance-refreshed fit carry identical shapes.

  fit$covariance_status <- raw$covariance_status %||% "not_requested"
  # The correlation gate check_strictness() applies is read off the matrix this
  # call just produced, so it is refreshed with it (ferx-core #1177); NaN when
  # the step still produced no matrix.
  fit$max_abs_correlation <- if (!is.null(raw$max_abs_correlation) &&
                                   is.finite(raw$max_abs_correlation)) {
    as.numeric(raw$max_abs_correlation)
  } else {
    NA_real_
  }

  theta_names <- names(fit$theta)
  n_theta <- length(fit$theta)
  n_sigma <- length(fit$sigma)
  sig_nms <- if (!is.null(fit$sigma_names) && length(fit$sigma_names) == n_sigma) {
    fit$sigma_names
  } else {
    names(fit$sigma)
  }
  eta_nms <- if (!is.null(fit$eta_names) && length(fit$eta_names) == n_eta) {
    fit$eta_names
  } else {
    NULL
  }

  # The engine packs the `block_sigma` correlations last, after sigma, so their
  # coordinates come out of the count before the rest can be read as omega -
  # otherwise the sigma rows are labelled at the wrong offset.
  rho_nms <- .ferx_residual_corr_labels(fit)
  n_rho <- length(rho_nms)

  d <- raw$cov_matrix_dim %||% 0L
  if (!is.null(raw$cov_matrix) && length(raw$cov_matrix) > 0L && d > 0L) {
    m <- matrix(raw$cov_matrix, nrow = d, ncol = d, byrow = TRUE)
    n_omega_packed <- d - n_theta - n_sigma - n_rho
    omega_names <- if (n_omega_packed == n_eta) {
      if (!is.null(eta_nms)) eta_nms else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
    } else {
      # Block lower-triangle: L(i,j) for i >= j, column-major.
      nm <- character(max(n_omega_packed, 0L))
      k <- 0L
      for (i in seq_len(n_eta)) {
        for (j in seq_len(i)) {
          k <- k + 1L
          nm[k] <- if (!is.null(eta_nms)) sprintf("%s,%s", eta_nms[i], eta_nms[j]) else sprintf("OMEGA(%d,%d)", i, j)
        }
      }
      nm
    }
    pnames <- c(
      theta_names,
      if (n_omega_packed > 0L) omega_names else character(0L),
      if (n_sigma > 0L) (if (!is.null(sig_nms)) sig_nms else paste0("SIGMA(", seq_len(n_sigma), ")")) else character(0L),
      rho_nms
    )
    if (length(pnames) == d) rownames(m) <- colnames(m) <- pnames
    fit$cov_matrix <- m
  } else {
    fit$cov_matrix <- NULL
  }

  # Standard errors: name theta SEs, drop empties to NULL (consistent with
  # ferx_fit()). se_omega / se_sigma stay as bare vectors - the printers index
  # them positionally via helpers.
  se_theta <- raw$se_theta
  if (length(se_theta) > 0L) {
    if (length(se_theta) == length(theta_names)) names(se_theta) <- theta_names
    fit$se_theta <- se_theta
  } else {
    fit$se_theta <- NULL
  }
  fit$se_omega <- if (length(raw$se_omega) > 0L) raw$se_omega else NULL
  fit$se_sigma <- if (length(raw$se_sigma) > 0L) raw$se_sigma else NULL
  # `block_sigma` correlation SEs come out of the same step; write them back
  # onto the frame so a refreshed fit carries uncertainty for every estimated
  # parameter, not just theta/omega/sigma.
  if (n_rho > 0L && is.data.frame(fit$residual_correlations)) {
    se_rho <- as.numeric(raw$se_residual_correlation %||% numeric())
    fit$residual_correlations$se <-
      if (length(se_rho) == n_rho) se_rho else rep(NA_real_, n_rho)
  }

  # se_kappa: name per the IOV convention (diagonal -> kappa_names, block ->
  # lower-triangle labels). NULL when no IOV or the step produced none.
  d_iov <- if (!is.null(fit$omega_iov) && is.matrix(fit$omega_iov)) nrow(fit$omega_iov) else 0L
  if (length(raw$se_kappa) == 0L || d_iov == 0L) {
    fit$se_kappa <- NULL
  } else {
    se_kappa <- raw$se_kappa
    n_tri <- d_iov * (d_iov + 1L) / 2L
    if (length(se_kappa) == d_iov && length(fit$kappa_names) == d_iov) {
      names(se_kappa) <- fit$kappa_names
    } else if (length(se_kappa) == n_tri && length(fit$kappa_names) == d_iov) {
      tri_names <- character(n_tri)
      idx <- 1L
      for (j in seq_len(d_iov)) {
        for (i in j:d_iov) {
          tri_names[idx] <- if (i == j) {
            fit$kappa_names[i]
          } else {
            paste0("COV_", fit$kappa_names[j], "_", fit$kappa_names[i])
          }
          idx <- idx + 1L
        }
      }
      names(se_kappa) <- tri_names
    }
    fit$se_kappa <- se_kappa
  }

  # Eigenvalues / condition number. Stage the raw FFI values on the wire-named
  # fields and let the shared .ferx_apply_cov_sentinels() do the conversion:
  # it sets fit$eigenvalues / fit$condition_number, clears the (now stale) wire
  # fields, and appends the flat high-condition-number warning - exactly what
  # ferx_fit() does, so the two paths can't drift and a loaded fit's old
  # cov_eigenvalues / cov_condition_number don't survive alongside a new matrix.
  fit$cov_eigenvalues <- if (length(raw$cov_eigenvalues) == 0L) NULL else raw$cov_eigenvalues
  fit$cov_condition_number <- raw$cov_condition_number
  fit <- .ferx_apply_cov_sentinels(fit)

  # Covariance-step warnings must reach the STRUCTURED table, not just the flat
  # vector: ferx_get_warnings() and the print/summary tally read
  # fit$warnings_structured and only fall back to fit$warnings when it is
  # absent (never, on a real fit). Mirror ferx_sir()'s handling.
  if (length(raw$warnings) > 0L) {
    fit$warnings <- unique(c(fit$warnings, raw$warnings))
  }
  ws <- fit$warnings_structured
  if (!is.data.frame(ws)) ws <- NULL
  # Drop any stale condition_number row so a re-run that is now well-conditioned
  # (or has a different condition number) doesn't leave a contradictory flag;
  # .ferx_assemble_structured_warnings() re-adds the current one below.
  if (!is.null(ws) && nrow(ws) > 0L) {
    ws <- ws[ws$category != "condition_number", , drop = FALSE]
  }
  # Assembled rows: the current condition_number flag + ETA-normality (both
  # derived from the fit). The engine's covariance-step warnings arrive as flat
  # strings (the binding returns no severity/category), so fold them in as
  # covariance-category rows.
  assembled <- .ferx_assemble_structured_warnings(raw, fit)
  cov_rows <- if (length(raw$warnings) > 0L) {
    data.frame(
      severity      = "warning",
      category      = "covariance",
      message       = as.character(raw$warnings),
      source_method = "",
      stringsAsFactors = FALSE
    )
  } else {
    assembled[0, , drop = FALSE]
  }
  merged <- rbind(ws %||% assembled[0, , drop = FALSE], assembled, cov_rows)
  fit$warnings_structured <- if (nrow(merged) > 0L) unique(merged) else NULL

  # Refresh derived fields (cor_matrix, estimates table) from the new cov_matrix
  # so they can't drift from the fitting path.
  fit <- .ferx_populate_derived_fields(fit)

  fit
}
