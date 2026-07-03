# ----- load -----------------------------------------------------------------

#' Load a fit from a \code{.fitrx} bundle
#'
#' Reads a portable \code{.fitrx} archive written by \code{\link{ferx_save_fit}}
#' (or by the ferx-core Rust CLI's \code{--output} flag) and reconstructs a
#' \code{ferx_fit}-classed list. The returned object is structurally
#' compatible with one returned by \code{\link{ferx_fit}} for the fields
#' covered by the cross-language schema, plus any R-specific fields the
#' writer preserved under \code{r_extras}.
#'
#' To call \code{\link{ferx_predict}} on a loaded fit, the embedded model
#' source can be re-parsed through the existing pipeline.
#'
#' @param path Path to a \code{.fitrx} file.
#' @return A \code{ferx_fit}-classed list.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' tmp <- tempfile(fileext = ".rds")
#' ferx_save_fit(fit, tmp)
#' fit2 <- ferx_load_fit(tmp)
#' identical(fit$theta, fit2$theta)
#' @family utilities
#' @export
#' @seealso \code{\link{ferx_save_fit}}.
ferx_load_fit <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("`path` must be a single character path to a .fitrx file.")
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }

  staging <- tempfile("fitrx_load_")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  # `junkpaths = TRUE` flattens every entry into `staging` using only its
  # basename. The .fitrx layout is intentionally flat - no subdirectories -
  # so this is a no-op for well-formed archives and a defence against Zip
  # Slip (entries with `../` or absolute paths) for malformed ones.
  files <- tryCatch(
    utils::unzip(path, exdir = staging, junkpaths = TRUE, overwrite = TRUE),
    error = function(e) stop("Failed to read .fitrx archive '", path, "': ", conditionMessage(e))
  )
  if (length(files) == 0L) {
    stop("'", path, "' is not a valid .fitrx archive (no entries).")
  }

  manifest_path <- file.path(staging, "manifest.json")
  if (!file.exists(manifest_path)) {
    stop("'", path, "' is missing manifest.json")
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  if (!identical(as.character(manifest$format_version), FITRX_FORMAT_VERSION)) {
    stop(sprintf(
      "Unsupported .fitrx format_version %s; this reader expects %s",
      as.character(manifest$format_version),
      FITRX_FORMAT_VERSION
    ))
  }

  wire <- jsonlite::read_json(file.path(staging, "fit.json"), simplifyVector = FALSE)
  result <- .fitrx_wire_to_fit(wire)

  # Subject-level tables
  ebes_path <- file.path(staging, "ebes.csv")
  if (file.exists(ebes_path)) {
    result$ebe_etas <- utils::read.csv(ebes_path, stringsAsFactors = FALSE, check.names = FALSE)
    # Preserve the live-fit type of the ID column (character); read.csv would
    # otherwise infer integer for numeric-looking IDs, diverging from
    # cond_dist$data$ID (and covtab$ID) below.
    result$ebe_etas$ID <- as.character(result$ebe_etas$ID)
  }
  preds_path <- file.path(staging, "predictions.csv")
  if (file.exists(preds_path)) {
    result$sdtab <- utils::read.csv(preds_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  covtab_path <- file.path(staging, "covtab.csv")
  if (file.exists(covtab_path)) {
    result$covtab <- utils::read.csv(covtab_path, stringsAsFactors = FALSE, check.names = FALSE)
    # Preserve the live-fit type of the ID column (character); read.csv would
    # otherwise infer integer for numeric-looking IDs.
    result$covtab$ID <- as.character(result$covtab$ID)
  }
  kappa_path <- file.path(staging, "ebes_kappa.csv")
  if (file.exists(kappa_path)) {
    result$ebe_kappas <- utils::read.csv(kappa_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    result$ebe_kappas <- NULL
  }
  conddist_path <- file.path(staging, "conddist.csv")
  if (file.exists(conddist_path)) {
    cd_df <- utils::read.csv(conddist_path, stringsAsFactors = FALSE, check.names = FALSE)
    cd_df$ID <- as.character(cd_df$ID)
    result$cond_dist <- list(
      data      = cd_df,
      # `conddist.csv` carries no distribution-based shrinkage / provenance
      # fields, so recompute shrinkage from the loaded mean and omega (same
      # formula as ferx-core's own conddist.csv reader, ferx-core #675).
      # nsamp/burnin aren't recoverable from the CSV alone; NA (not 0L) so
      # "unknown after reload" isn't confused with "zero draws retained".
      shrinkage = .fitrx_cond_dist_shrinkage(cd_df, result$omega, result$eta_names),
      nsamp     = NA_integer_,
      burnin    = NA_integer_
    )
  } else {
    result$cond_dist <- NULL
  }
  trace_csv_path <- file.path(staging, "trace.csv")
  result$trace <- if (file.exists(trace_csv_path)) {
    .ferx_read_trace_csv(trace_csv_path)
  } else {
    NULL
  }

  # Model source + optional data CSV.
  #
  # Both files are bundled into the .fitrx archive. Persist each to a
  # tempfile (outside the staging dir which gets cleaned on exit), and
  # *override* the wire-provided `model_path` / `data_path` with the
  # local copy. That way `ferx_sir(fit)` and `ferx_predict(fit)` work on
  # the loading machine even if the original-machine paths don't exist
  # - and the stored `model_hash` / `data_hash` still match, because
  # the bundle is a verbatim copy of the source bytes.
  model_staging <- file.path(staging, "model.ferx")
  if (file.exists(model_staging)) {
    result$model_source <- paste(readLines(model_staging, warn = FALSE), collapse = "\n")
    result$model_text   <- result$model_source  # ferx_runlog() reads model_text
    persisted_model <- tempfile("fitrx_model_", fileext = ".ferx")
    file.copy(model_staging, persisted_model, overwrite = TRUE)
    result$model_path <- persisted_model
  } else if (is.null(result$model_path) || !nzchar(result$model_path) ||
             !file.exists(result$model_path)) {
    # `model.ferx` should always be present in a well-formed bundle, but
    # be defensive: if it's missing AND the wire path doesn't resolve
    # locally, clear `model_path` to NA so `ferx_sir()`'s "no recorded
    # model path" branch fires with a clear message instead of the
    # generic "file not found".
    result$model_path <- NA_character_
  }
  data_staging <- file.path(staging, "data.csv")
  if (file.exists(data_staging)) {
    # Copy out of the staging dir so it survives the on.exit() cleanup.
    persisted_data <- tempfile("fitrx_data_", fileext = ".csv")
    file.copy(data_staging, persisted_data, overwrite = TRUE)
    result$data_path <- persisted_data
  } else if (is.null(result$data_path) || !nzchar(result$data_path)) {
    # No bundled data and the wire didn't carry a path either.
    result$data_path <- NA_character_
  }
  # If the wire carried model_path / data_path and neither was bundled
  # but the path resolves locally, leave it as-is - `ferx_sir()` will
  # use it directly (and the hash check still applies).

  # Warnings - fit.json is the source of truth, warnings.txt is a mirror.
  if (is.null(result$warnings)) result$warnings <- character()

  # Derived fields computed at fit-time by ferx_fit() are not part of the
  # cross-language .fitrx schema; recompute them here (via the same helper
  # ferx_fit() uses) so a loaded fit is indistinguishable from a
  # freshly-fitted one. eta_cov additionally requires the dataset to still
  # be readable from data_path (bundled via `include_data`, or resolving
  # locally).
  result <- .ferx_populate_derived_fields(result)

  class(result) <- "ferx_fit"
  result
}

.fitrx_wire_to_fit <- function(w) {
  out <- list(
    method = .fitrx_method_label(w$method),
    method_chain = vapply(as.character(w$method_chain), .fitrx_method_label, character(1L), USE.NAMES = FALSE),
    converged = isTRUE(w$converged),
    ofv = as.numeric(w$ofv),
    aic = as.numeric(w$aic),
    bic = as.numeric(w$bic),
    n_obs = as.integer(w$n_obs),
    n_subjects = as.integer(w$n_subjects),
    n_parameters = as.integer(w$n_parameters),
    n_iterations = as.integer(w$n_iterations),
    interaction = isTRUE(w$interaction),
    wall_time_secs = as.numeric(w$wall_time_secs),
    n_threads_used = as.integer(w$n_threads_used %||% NA_integer_),
    uses_ode_solver = isTRUE(w$uses_ode_solver),
    uses_sde = isTRUE(w$uses_sde),
    gradient_method_inner = as.character(w$gradient_method_inner %||% ""),
    gradient_method_outer = as.character(w$gradient_method_outer %||% ""),
    nlopt_missing_algorithms = unlist(w$nlopt_missing_algorithms %||% list(), use.names = FALSE),
    covariance_status = .fitrx_covariance_status_label(w$covariance_status),
    covariance_n_evals_estimated = .fitrx_unwrap_opt_num(w$covariance_n_evals_estimated),
    trace_path = .fitrx_unwrap_opt_chr(w$trace_path),
    ebe_convergence_warnings = as.integer(w$ebe_convergence_warnings %||% 0L),
    max_unconverged_subjects = as.integer(w$max_unconverged_subjects %||% 0L),
    total_ebe_fallbacks = as.integer(w$total_ebe_fallbacks %||% 0L),
    warnings = unlist(w$warnings %||% list(), use.names = FALSE),
    saem_mu_ref_m_step_evals_saved = .fitrx_unwrap_opt_num(w$saem_mu_ref_m_step_evals_saved),

    theta = stats::setNames(
      as.numeric(unlist(w$theta$estimates, use.names = FALSE)),
      unlist(w$theta$names, use.names = FALSE)
    ),
    theta_fixed = as.logical(unlist(w$theta$fixed %||% list(), use.names = FALSE)),
    theta_transforms = stats::setNames(
      as.character(unlist(w$theta$transform, use.names = FALSE)),
      unlist(w$theta$names, use.names = FALSE)
    ),
    se_theta = .fitrx_unwrap_named_se(w$theta$se, w$theta$names),

    eta_names = unlist(w$omega$names, use.names = FALSE),
    omega = .fitrx_named_omega(w$omega),
    omega_fixed = as.logical(unlist(w$omega$fixed %||% list(), use.names = FALSE)),
    eta_log_transformed = as.logical(unlist(w$omega$log_transformed %||% list(), use.names = FALSE)),
    omega_param_corr = .fitrx_matrix_from_wire(w$omega$param_corr),
    shrinkage_eta = as.numeric(unlist(w$omega$shrinkage %||% list(), use.names = FALSE)),
    se_omega = .fitrx_unwrap_opt_num_vec(w$omega$se),
    omega_init_as_sd = as.logical(unlist(w$omega$init_as_sd %||% list(), use.names = FALSE)),

    sigma_names = unlist(w$sigma$names, use.names = FALSE),
    sigma = stats::setNames(
      as.numeric(unlist(w$sigma$estimates, use.names = FALSE)),
      unlist(w$sigma$names, use.names = FALSE)
    ),
    sigma_fixed = as.logical(unlist(w$sigma$fixed %||% list(), use.names = FALSE)),
    sigma_types = stats::setNames(
      as.character(unlist(w$sigma$types, use.names = FALSE)),
      unlist(w$sigma$names, use.names = FALSE)
    ),
    se_sigma = .fitrx_unwrap_opt_num_vec(w$sigma$se),
    sigma_init_as_sd = as.logical(unlist(w$sigma$init_as_sd %||% list(), use.names = FALSE)),

    shrinkage_eps = as.numeric(w$shrinkage_eps),
    dw_statistic = .fitrx_unwrap_opt_num(w$dw_statistic),
    iwres_lag1_r = .fitrx_unwrap_opt_num(w$iwres_lag1_r),
    cov_matrix = .fitrx_matrix_from_wire(w$covariance_matrix),
    cov_eigenvalues = .fitrx_unwrap_opt_num_vec(w$cov_eigenvalues),
    cov_condition_number = .fitrx_unwrap_opt_num(w$cov_condition_number),

    model_name = as.character(w$model_name %||% ""),
    ferx_version = as.character(w$ferx_version %||% ""),

    model_path = .fitrx_unwrap_opt_chr(w$model_path),
    data_path = .fitrx_unwrap_opt_chr(w$data_path),
    model_hash = .fitrx_unwrap_opt_chr(w$model_hash),
    data_hash = .fitrx_unwrap_opt_chr(w$data_hash),

    input_columns   = as.character(unlist(w$input_columns   %||% list(), use.names = FALSE)),
    covariate_names = as.character(unlist(w$covariate_names %||% list(), use.names = FALSE))
  )

  # Data-selection exclusion summary (NULL for older .fitrx files without the field).
  if (!is.null(w$exclusions)) {
    we <- w$exclusions
    out$exclusions <- list(
      n_records_total      = as.integer(we$n_records_total  %||% 0L),
      n_obs_excluded       = as.integer(we$n_obs_excluded   %||% 0L),
      n_dose_excluded      = as.integer(we$n_dose_excluded  %||% 0L),
      n_other_excluded     = as.integer(we$n_other_excluded %||% 0L),
      excluded_subject_ids = as.character(unlist(we$excluded_subject_ids %||% list(), use.names = FALSE)),
      fired_ignore         = as.character(unlist(we$fired_ignore          %||% list(), use.names = FALSE)),
      fired_accept         = as.character(unlist(we$fired_accept          %||% list(), use.names = FALSE))
    )
  } else {
    out$exclusions <- NULL
  }

  # covariate_types: rebuild the named character vector (NULL when absent).
  ct <- w$covariate_types
  if (!is.null(ct) && length(ct) > 0L) {
    out$covariate_types <- unlist(ct, use.names = TRUE)
  }

  # eta_param_info ? parallel R vectors
  epi <- w$eta_param_info
  if (!is.null(epi) && length(epi) > 0L) {
    out$eta_param_types <- vapply(epi, function(x) as.character(x$param_type %||% ""), character(1L))
    out$eta_linked_theta <- vapply(epi, function(x) as.character(x$linked_theta %||% ""), character(1L))
  } else {
    out$eta_param_types <- character()
    out$eta_linked_theta <- character()
  }

  # SIR
  if (!is.null(w$sir)) {
    out$sir_ess <- .fitrx_unwrap_opt_num(w$sir$ess)
    out$sir_ci_theta <- .fitrx_unwrap_ci(w$sir$ci_theta)
    out$sir_ci_omega <- .fitrx_unwrap_ci(w$sir$ci_omega)
    out$sir_ci_sigma <- .fitrx_unwrap_ci(w$sir$ci_sigma)
  }

  # Bayes posterior summary
  if (!is.null(w$bayes)) {
    wb <- w$bayes
    out$bayes <- list(
      n_chains = .fitrx_unwrap_opt_num(wb$n_chains),
      n_warmup = .fitrx_unwrap_opt_num(wb$n_warmup),
      n_draws_per_chain = .fitrx_unwrap_opt_num(wb$n_draws_per_chain),
      n_divergent = .fitrx_unwrap_opt_num(wb$n_divergent),
      max_rhat = .fitrx_unwrap_opt_num(wb$max_rhat),
      param_names = as.character(unlist(wb$param_names %||% list(), use.names = FALSE)),
      mean = .fitrx_unwrap_opt_num_vec(wb$mean),
      sd = .fitrx_unwrap_opt_num_vec(wb$sd),
      q025 = .fitrx_unwrap_opt_num_vec(wb$q025),
      median = .fitrx_unwrap_opt_num_vec(wb$median),
      q975 = .fitrx_unwrap_opt_num_vec(wb$q975),
      rhat = .fitrx_unwrap_opt_num_vec(wb$rhat),
      ess_bulk = .fitrx_unwrap_opt_num_vec(wb$ess_bulk),
      ess_tail = .fitrx_unwrap_opt_num_vec(wb$ess_tail),
      mcse = .fitrx_unwrap_opt_num_vec(wb$mcse)
    )
  }

  # IOV
  if (!is.null(w$iov)) {
    out$omega_iov <- .fitrx_matrix_from_wire(w$iov$omega_iov)
    out$kappa_names <- unlist(w$iov$kappa_names %||% list(), use.names = FALSE)
    out$kappa_fixed <- as.logical(unlist(w$iov$kappa_fixed %||% list(), use.names = FALSE))
    out$se_kappa <- .fitrx_unwrap_opt_num_vec(w$iov$se_kappa)
    out$shrinkage_kappa <- as.numeric(unlist(w$iov$shrinkage_kappa %||% list(), use.names = FALSE))
    raw_by_occ <- w$iov$shrinkage_kappa_by_occ
    if (!is.null(raw_by_occ) && length(raw_by_occ) > 0L) {
      kn <- out$kappa_names
      occ_col <- seq_along(raw_by_occ)
      kappa_mat <- do.call(rbind, lapply(raw_by_occ, as.numeric))
      df_by_occ <- as.data.frame(cbind(occ = occ_col, kappa_mat))
      colnames(df_by_occ) <- c("occ", kn)
      out$shrinkage_kappa_by_occ <- df_by_occ
    } else {
      out$shrinkage_kappa_by_occ <- NULL
    }
    out$omega_iov_param_corr <- .fitrx_matrix_from_wire(w$iov$omega_iov_param_corr)
    out$kappa_init_as_sd <- as.logical(unlist(w$iov$kappa_init_as_sd %||% list(), use.names = FALSE))
  } else {
    out$kappa_names <- character()
    out$kappa_fixed <- logical()
    out$se_kappa <- NULL
    out$shrinkage_kappa <- numeric()
    out$shrinkage_kappa_by_occ <- NULL
    out$omega_iov <- NULL
    out$omega_iov_param_corr <- NULL
    out$kappa_init_as_sd <- logical()
  }

  # R extras
  extras <- w$r_extras
  if (!is.null(extras)) {
    for (key in names(extras)) {
      out[[key]] <- extras[[key]]
    }
  }
  out
}

.fitrx_method_label <- function(token) {
  if (is.null(token) || !nzchar(token)) return("FOCEI")
  switch(token,
    "foce" = "FOCE",
    "focei" = "FOCEI",
    "foce_gn" = "FOCE-GN",
    "foce_gn_hybrid" = "FOCE-GN-Hybrid",
    "saem" = "SAEM",
    as.character(token)
  )
}

.fitrx_covariance_status_label <- function(token) {
  if (is.null(token)) return("NotRequested")
  switch(as.character(token),
    "computed" = "Computed",
    "failed" = "Failed",
    "not_requested" = "NotRequested",
    "sir_fallback" = "SirFallback",
    as.character(token)
  )
}

.fitrx_matrix_from_wire <- function(w) {
  if (is.null(w)) return(NULL)
  rows <- as.integer(w$rows %||% 0L)
  cols <- as.integer(w$cols %||% 0L)
  data <- as.numeric(unlist(w$data %||% list(), use.names = FALSE))
  if (rows == 0L || cols == 0L || length(data) != rows * cols) return(NULL)
  matrix(data, nrow = rows, ncol = cols, byrow = TRUE)
}

.fitrx_unwrap_ci <- function(wire_ci) {
  if (is.null(wire_ci) || length(wire_ci) == 0L) return(NULL)
  mat <- do.call(rbind, lapply(wire_ci, function(p) as.numeric(unlist(p, use.names = FALSE))))
  if (is.null(mat) || ncol(mat) != 2L) return(NULL)
  mat
}

.fitrx_named_omega <- function(omega_wire) {
  om <- .fitrx_matrix_from_wire(omega_wire$matrix)
  if (is.null(om)) return(om)
  en <- unlist(omega_wire$names, use.names = FALSE)
  d  <- nrow(om)
  nms <- if (!is.null(en) && length(en) == d) en else paste0("OMEGA(", seq_len(d), ",", seq_len(d), ")")
  rownames(om) <- colnames(om) <- nms
  om
}

# Distribution-based eta-shrinkage recomputed from a loaded conddist.csv:
# `1 - SD_over_subjects(cond_mean[., j]) / sqrt(Omega_jj)`, parallel to
# `eta_names` (same formula as ferx-core's own conddist.csv reader,
# ferx-core #675/src/io/fitrx.rs). NA when fewer than two subjects or Omega_jj
# is at/below ferx-core's own near-singular floor (`SAEM_OMEGA_DIAG_FLOOR`,
# ferx-core/src/estimation/saem.rs) -- matched here so a near-singular omega
# doesn't silently produce a finite value on this path when ferx-core's own
# reader would report NaN.
.FITRX_SAEM_OMEGA_DIAG_FLOOR <- 1e-6
.fitrx_cond_dist_shrinkage <- function(cd_df, omega, eta_names) {
  if (is.null(eta_names) || length(eta_names) == 0L) return(numeric(0))
  vapply(eta_names, function(nm) {
    vals <- cd_df$COND_MEAN[cd_df$ETA == nm]
    omega_jj <- if (!is.null(omega) && nm %in% rownames(omega)) omega[nm, nm] else NA_real_
    if (length(vals) < 2L || is.na(omega_jj) ||
        omega_jj < .FITRX_SAEM_OMEGA_DIAG_FLOOR) {
      return(NA_real_)
    }
    1 - stats::sd(vals) / sqrt(omega_jj)
  }, numeric(1), USE.NAMES = FALSE)
}

.fitrx_unwrap_named_se <- function(se_wire, names_wire) {
  if (is.null(se_wire)) return(NULL)
  v <- as.numeric(unlist(se_wire, use.names = FALSE))
  if (length(v) == 0L) return(NULL)
  nm <- unlist(names_wire, use.names = FALSE)
  if (length(nm) == length(v)) stats::setNames(v, nm) else v
}
