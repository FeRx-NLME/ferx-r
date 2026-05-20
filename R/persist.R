# Portable .fitrx fit-bundle save/load.
#
# A .fitrx file is a zip archive containing JSON for scalars/vectors/matrices
# and CSV for per-subject/per-observation tables, plus the verbatim .ferx
# model source. The format is shared with ferx-core (Rust) and documented in
# `docs/src/file-formats/fitrx.md` of that crate. R-only fields that have no
# place in the cross-language schema (call settings, R-side derived labels,
# etc.) live under a nested `r_extras` key in `fit.json` and are silently
# ignored by non-R readers.

FITRX_FORMAT_VERSION <- "1"

# ----- save -----------------------------------------------------------------

#' Save a ferx fit to a portable \code{.fitrx} bundle
#'
#' Writes a single, self-describing zip archive that any language with stdlib
#' JSON + CSV + zip readers can open (R, Python, Julia, Rust). The archive
#' includes parameter estimates, per-subject EBEs, per-observation
#' predictions, and the verbatim \code{.ferx} model source.
#'
#' The schema is shared with the ferx-core Rust crate; see its
#' \code{docs/src/file-formats/fitrx.md} for the full field reference.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @param output Path to write. By convention the file extension is
#'   \code{.fitrx}; the function does not enforce this.
#' @param include_data Logical. When \code{TRUE}, the input NONMEM CSV used
#'   to fit the model is embedded verbatim inside the archive (as
#'   \code{data.csv}). Requires that the file is still accessible at the
#'   path captured at fit time (\code{fit$data_path}). Default \code{FALSE}.
#' @return Invisibly returns \code{fit}, so the function can be used inside a
#'   pipe without breaking the chain:
#'   \preformatted{
#'   fit |> ferx_save_fit("run1.fitrx") |> ferx_estimates()
#'   }
#' @export
#' @seealso \code{\link{ferx_load_fit}}, \code{\link{ferx_fit}} (and its
#'   \code{output} argument for save-during-fit).
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' tmp <- tempfile(fileext = ".fitrx")
#' ferx_save_fit(fit, tmp)
#' fit2 <- ferx_load_fit(tmp)
#' identical(fit$theta, fit2$theta)
ferx_save_fit <- function(fit, output, include_data = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (returned by ferx_fit()).")
  }
  if (!is.character(output) || length(output) != 1L || is.na(output)) {
    stop("`output` must be a single character path.")
  }
  if (!is.logical(include_data) || length(include_data) != 1L || is.na(include_data)) {
    stop("`include_data` must be TRUE or FALSE.")
  }

  staging <- tempfile("fitrx_")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  entries <- character()

  # fit.json
  fit_json_path <- file.path(staging, "fit.json")
  .fitrx_write_fit_json(fit, fit_json_path)
  entries <- c(entries, "fit.json")

  # ebes.csv
  .fitrx_write_ebes_csv(fit, file.path(staging, "ebes.csv"))
  entries <- c(entries, "ebes.csv")

  # ebes_kappa.csv — present only when the model uses IOV
  if (.fitrx_has_kappa(fit)) {
    .fitrx_write_ebes_kappa_csv(fit, file.path(staging, "ebes_kappa.csv"))
    entries <- c(entries, "ebes_kappa.csv")
  }

  # predictions.csv
  .fitrx_write_predictions_csv(fit, file.path(staging, "predictions.csv"))
  entries <- c(entries, "predictions.csv")

  # model.ferx
  #
  # Prefer `file.copy(fit$model_path, ...)` so the bundled bytes match
  # the original file exactly — including line endings, BOM, and the
  # presence/absence of a trailing newline. `fit$model_source` was
  # produced by `readLines + paste(collapse = "\n")` and writing it back
  # via `writeLines` loses that fidelity, which causes
  # `sha256(staging/model.ferx)` to drift from the stored
  # `fit$model_hash`. That drift would then trip the hash-mismatch
  # check inside `ferx_sir(loaded_fit)` after a save/load round-trip on
  # any model whose bytes don't round-trip through readLines/writeLines
  # (CRLF endings, missing trailing newline, etc.).
  #
  # Fall back to `writeLines(fit$model_source)` only when no `model_path`
  # is available (e.g. a hand-constructed fit) — in that case the bundle
  # carries the best representation we have.
  model_src_path <- fit$model_path
  if (!is.null(model_src_path) && !is.na(model_src_path) &&
      nzchar(model_src_path) && file.exists(model_src_path)) {
    file.copy(model_src_path, file.path(staging, "model.ferx"), overwrite = TRUE)
  } else {
    writeLines(fit$model_source %||% "", file.path(staging, "model.ferx"))
  }
  entries <- c(entries, "model.ferx")

  # warnings.txt
  warns <- fit$warnings
  if (is.null(warns)) warns <- character()
  writeLines(warns, file.path(staging, "warnings.txt"))
  entries <- c(entries, "warnings.txt")

  # data.csv (optional)
  if (isTRUE(include_data)) {
    data_path <- fit$data_path
    if (is.null(data_path) || is.na(data_path) || !file.exists(data_path)) {
      stop("include_data = TRUE but fit$data_path is missing or does not exist.")
    }
    file.copy(data_path, file.path(staging, "data.csv"), overwrite = TRUE)
    entries <- c(entries, "data.csv")
  }

  # manifest.json — last so it can list every entry
  manifest <- list(
    format_version = FITRX_FORMAT_VERSION,
    ferx_version = fit$ferx_version %||% "",
    model_name = fit$model_name %||% "",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    entries = c("manifest.json", entries)
  )
  jsonlite::write_json(
    manifest,
    file.path(staging, "manifest.json"),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  # Zip everything into `output`. utils::zip needs relative paths from cwd.
  out_abs <- normalizePath(dirname(output), mustWork = FALSE)
  if (!dir.exists(out_abs)) dir.create(out_abs, recursive = TRUE)
  out_abs <- file.path(out_abs, basename(output))

  files_to_zip <- c("manifest.json", entries)
  old_wd <- getwd()
  setwd(staging)
  zip_status <- tryCatch(
    utils::zip(out_abs, files_to_zip, flags = "-q9X"),
    finally = setwd(old_wd)
  )
  if (!is.null(zip_status) && is.numeric(zip_status) && zip_status != 0L) {
    stop(sprintf("zip command failed with status %d while writing %s", zip_status, output))
  }

  invisible(fit)
}

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
  # basename. The .fitrx layout is intentionally flat — no subdirectories —
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
  }
  preds_path <- file.path(staging, "predictions.csv")
  if (file.exists(preds_path)) {
    result$sdtab <- utils::read.csv(preds_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  kappa_path <- file.path(staging, "ebes_kappa.csv")
  if (file.exists(kappa_path)) {
    result$ebe_kappas <- utils::read.csv(kappa_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    result$ebe_kappas <- NULL
  }

  # Model source + optional data CSV.
  #
  # Both files are bundled into the .fitrx archive. Persist each to a
  # tempfile (outside the staging dir which gets cleaned on exit), and
  # *override* the wire-provided `model_path` / `data_path` with the
  # local copy. That way `ferx_sir(fit)` and `ferx_predict(fit)` work on
  # the loading machine even if the original-machine paths don't exist
  # — and the stored `model_hash` / `data_hash` still match, because
  # the bundle is a verbatim copy of the source bytes.
  model_staging <- file.path(staging, "model.ferx")
  if (file.exists(model_staging)) {
    result$model_source <- paste(readLines(model_staging, warn = FALSE), collapse = "\n")
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
  # but the path resolves locally, leave it as-is — `ferx_sir()` will
  # use it directly (and the hash check still applies).

  # Warnings — fit.json is the source of truth, warnings.txt is a mirror.
  if (is.null(result$warnings)) result$warnings <- character()

  class(result) <- "ferx_fit"
  result
}

# ----- internals ------------------------------------------------------------

# Build the on-disk fit.json from the in-memory ferx_fit list. The schema
# matches the Rust FitWire layout (see ferx-core/src/io/fitrx.rs); R-only
# fields land under `r_extras` and are ignored by non-R readers.
.fitrx_write_fit_json <- function(fit, path) {
  wire <- list(
    method = .fitrx_method_to_token(fit$method),
    method_chain = vapply(
      as.character(fit$method_chain %||% fit$method),
      .fitrx_method_to_token, character(1L), USE.NAMES = FALSE
    ),
    converged = isTRUE(fit$converged),
    ofv = as.numeric(fit$ofv),
    aic = as.numeric(fit$aic),
    bic = as.numeric(fit$bic),
    n_obs = as.integer(fit$n_obs),
    n_subjects = as.integer(fit$n_subjects),
    n_parameters = as.integer(fit$n_parameters),
    n_iterations = as.integer(fit$n_iterations),
    interaction = isTRUE(fit$interaction),
    wall_time_secs = as.numeric(fit$wall_time_secs %||% NA_real_),
    n_threads_used = as.integer(fit$n_threads_used %||% NA_integer_),
    uses_ode_solver = isTRUE(fit$uses_ode_solver),
    uses_sde = isTRUE(fit$uses_sde),
    gradient_method_inner = as.character(fit$gradient_method_inner %||% ""),
    gradient_method_outer = as.character(fit$gradient_method_outer %||% ""),
    nlopt_missing_algorithms = as.character(fit$nlopt_missing_algorithms %||% character()),
    covariance_status = .fitrx_covariance_status_to_token(fit$covariance_status),
    covariance_n_evals_estimated = .fitrx_opt_int(fit$covariance_n_evals_estimated),
    trace_path = .fitrx_opt_chr(fit$trace_path),
    ebe_convergence_warnings = as.integer(fit$ebe_convergence_warnings %||% 0L),
    max_unconverged_subjects = as.integer(fit$max_unconverged_subjects %||% 0L),
    total_ebe_fallbacks = as.integer(fit$total_ebe_fallbacks %||% 0L),
    warnings = as.character(fit$warnings %||% character()),
    saem_mu_ref_m_step_evals_saved = .fitrx_opt_int(fit$saem_mu_ref_m_step_evals_saved),

    theta = list(
      names = as.character(names(fit$theta) %||% character()),
      estimates = as.numeric(fit$theta),
      se = .fitrx_opt_num_vec(fit$se_theta),
      fixed = as.logical(fit$theta_fixed %||% rep(FALSE, length(fit$theta))),
      transform = as.character(fit$theta_transforms %||% rep("identity", length(fit$theta)))
    ),
    omega = list(
      names = as.character(fit$eta_names %||% character()),
      matrix = .fitrx_matrix_to_wire(fit$omega),
      se = .fitrx_opt_num_vec(fit$se_omega),
      fixed = as.logical(fit$omega_fixed %||% rep(FALSE, NROW(fit$omega))),
      log_transformed = as.logical(
        fit$eta_log_transformed %||% rep(TRUE, NROW(fit$omega))
      ),
      param_corr = .fitrx_matrix_to_wire(fit$omega_param_corr),
      shrinkage = as.numeric(fit$shrinkage_eta %||% numeric())
    ),
    sigma = list(
      names = as.character(fit$sigma_names %||% character()),
      estimates = as.numeric(fit$sigma %||% numeric()),
      se = .fitrx_opt_num_vec(fit$se_sigma),
      fixed = as.logical(fit$sigma_fixed %||% rep(FALSE, length(fit$sigma %||% numeric()))),
      types = as.character(fit$sigma_types %||% rep("proportional", length(fit$sigma %||% numeric())))
    ),
    error_model = .fitrx_error_model_from_sigma_types(fit$sigma_types),
    shrinkage_eps = as.numeric(fit$shrinkage_eps %||% NA_real_),
    dw_statistic = .fitrx_opt_num(fit$dw_statistic),
    iwres_lag1_r = .fitrx_opt_num(fit$iwres_lag1_r),
    covariance_matrix = .fitrx_matrix_to_wire(fit$cov_matrix),
    cov_eigenvalues = .fitrx_opt_num_vec(fit$cov_eigenvalues),
    cov_condition_number = .fitrx_opt_num(fit$cov_condition_number),

    sir = .fitrx_build_sir_wire(fit),
    iov = .fitrx_build_iov_wire(fit),

    eta_param_info = .fitrx_build_eta_param_info(fit),
    model_name = as.character(fit$model_name %||% ""),
    ferx_version = as.character(fit$ferx_version %||% ""),

    # Source-file provenance (matches the ferx-core FitWire layout). Only
    # written when present so older readers don't trip on unknown nulls.
    model_path = .fitrx_opt_chr(fit$model_path),
    data_path = .fitrx_opt_chr(fit$data_path),
    model_hash = .fitrx_opt_chr(fit$model_hash),
    data_hash = .fitrx_opt_chr(fit$data_hash),

    r_extras = .fitrx_collect_r_extras(fit)
  )

  jsonlite::write_json(
    wire,
    path,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA, # full f64 precision
    null = "null",
    na = "null"
  )
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
    covariance_n_evals_estimated = .fitrx_unwrap_opt_int(w$covariance_n_evals_estimated),
    trace_path = .fitrx_unwrap_opt_chr(w$trace_path),
    ebe_convergence_warnings = as.integer(w$ebe_convergence_warnings %||% 0L),
    max_unconverged_subjects = as.integer(w$max_unconverged_subjects %||% 0L),
    total_ebe_fallbacks = as.integer(w$total_ebe_fallbacks %||% 0L),
    warnings = unlist(w$warnings %||% list(), use.names = FALSE),
    saem_mu_ref_m_step_evals_saved = .fitrx_unwrap_opt_int(w$saem_mu_ref_m_step_evals_saved),

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
    data_hash = .fitrx_unwrap_opt_chr(w$data_hash)
  )

  # eta_param_info → parallel R vectors
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

  # IOV
  if (!is.null(w$iov)) {
    out$omega_iov <- .fitrx_matrix_from_wire(w$iov$omega_iov)
    out$kappa_names <- unlist(w$iov$kappa_names %||% list(), use.names = FALSE)
    out$kappa_fixed <- as.logical(unlist(w$iov$kappa_fixed %||% list(), use.names = FALSE))
    out$se_kappa <- .fitrx_unwrap_opt_num_vec(w$iov$se_kappa)
    out$shrinkage_kappa <- as.numeric(unlist(w$iov$shrinkage_kappa %||% list(), use.names = FALSE))
    out$omega_iov_param_corr <- .fitrx_matrix_from_wire(w$iov$omega_iov_param_corr)
    out$kappa_init_as_sd <- as.logical(unlist(w$iov$kappa_init_as_sd %||% list(), use.names = FALSE))
  } else {
    out$kappa_names <- character()
    out$kappa_fixed <- logical()
    out$se_kappa <- NULL
    out$shrinkage_kappa <- numeric()
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

.fitrx_write_ebes_csv <- function(fit, path) {
  ebes <- fit$ebe_etas
  if (is.null(ebes) || nrow(ebes) == 0L) {
    # Header-only file so the entry shape is stable.
    cat("ID,ofv_contribution,n_obs\n", file = path)
    return(invisible())
  }
  # Pull ofv_contribution / n_obs out of sdtab (one row per subject) when not
  # already present on ebes — the Rust shim doesn't expose them on ebe_etas.
  if (!all(c("ofv_contribution", "n_obs") %in% names(ebes))) {
    per_subj <- .fitrx_per_subject_ofv_nobs(fit)
    if (!is.null(per_subj)) {
      # Match by string ID so character IDs work (PT001 etc.). Falls back to
      # NA on subjects whose ID isn't represented in sdtab.
      idx <- match(as.character(ebes$ID), per_subj$id)
      ebes$ofv_contribution <- per_subj$ofv[idx]
      ebes$n_obs <- per_subj$n_obs[idx]
    } else {
      ebes$ofv_contribution <- NA_real_
      ebes$n_obs <- NA_integer_
    }
  }
  utils::write.table(
    ebes, path,
    row.names = FALSE, quote = FALSE, sep = ",", na = ""
  )
}

.fitrx_write_ebes_kappa_csv <- function(fit, path) {
  k <- fit$ebe_kappas
  if (is.null(k) || nrow(k) == 0L) {
    return(invisible())
  }
  utils::write.table(
    k, path,
    row.names = FALSE, quote = FALSE, sep = ",", na = ""
  )
}

.fitrx_write_predictions_csv <- function(fit, path) {
  sdtab <- fit$sdtab
  if (is.null(sdtab) || nrow(sdtab) == 0L) {
    cat("ID,TIME,DV,PRED,IPRED,CWRES,IWRES,EBE_OFV,N_OBS\n", file = path)
    return(invisible())
  }
  # sdtab$ID coming straight from the Rust shim is a 1-based numeric subject
  # index; translate to the string IDs stored on ebe_etas so a cross-language
  # reader can join with data.csv. When sdtab$ID is already character (e.g.
  # because the fit was previously loaded from a .fitrx, where predictions.csv
  # stores string IDs), leave it untouched — the remap would be a no-op and
  # `max(sdtab$ID, ...)` would error on character input.
  if (is.numeric(sdtab$ID)) {
    string_ids <- .fitrx_subject_string_ids(fit)
    if (!is.null(string_ids) &&
        length(string_ids) >= max(sdtab$ID, 0L, na.rm = TRUE)) {
      sdtab$ID <- string_ids[as.integer(sdtab$ID)]
    }
  }
  utils::write.table(
    sdtab, path,
    row.names = FALSE, quote = FALSE, sep = ",", na = ""
  )
}

.fitrx_subject_string_ids <- function(fit) {
  ebes <- fit$ebe_etas
  if (is.null(ebes) || !"ID" %in% names(ebes)) return(NULL)
  as.character(ebes$ID)
}

.fitrx_per_subject_ofv_nobs <- function(fit) {
  sdtab <- fit$sdtab
  if (is.null(sdtab) || nrow(sdtab) == 0L) return(NULL)
  needed <- c("ID", "EBE_OFV", "N_OBS")
  if (!all(needed %in% names(sdtab))) return(NULL)
  ord <- !duplicated(sdtab$ID)
  list(
    id = as.character(sdtab$ID[ord]),
    ofv = as.numeric(sdtab$EBE_OFV[ord]),
    n_obs = as.integer(sdtab$N_OBS[ord])
  )
}

.fitrx_has_kappa <- function(fit) {
  k <- fit$ebe_kappas
  !is.null(k) && (is.data.frame(k) && nrow(k) > 0L)
}

.fitrx_method_to_token <- function(label) {
  if (is.null(label) || !nzchar(label)) return("focei")
  switch(tolower(label),
    "foce" = "foce",
    "focei" = "focei",
    "foce-gn" = "foce_gn",
    "foce_gn" = "foce_gn",
    "foce-gn-hybrid" = "foce_gn_hybrid",
    "foce_gn_hybrid" = "foce_gn_hybrid",
    "saem" = "saem",
    tolower(label)
  )
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

.fitrx_covariance_status_to_token <- function(s) {
  if (is.null(s) || !nzchar(s)) return("not_requested")
  switch(tolower(s),
    "computed" = "computed",
    "failed" = "failed",
    "not_requested" = "not_requested",
    "notrequested" = "not_requested",
    tolower(s)
  )
}

.fitrx_covariance_status_label <- function(token) {
  if (is.null(token)) return("NotRequested")
  switch(as.character(token),
    "computed" = "Computed",
    "failed" = "Failed",
    "not_requested" = "NotRequested",
    as.character(token)
  )
}

.fitrx_error_model_from_sigma_types <- function(types) {
  if (is.null(types) || length(types) == 0L) return("proportional")
  t <- tolower(as.character(types))
  if (all(t == "proportional")) return("proportional")
  if (all(t == "additive")) return("additive")
  "combined"
}

.fitrx_matrix_to_wire <- function(m) {
  if (is.null(m)) return(NULL)
  if (!is.matrix(m)) {
    # Try to coerce flat vector with an attribute; otherwise NULL.
    return(NULL)
  }
  list(
    rows = as.integer(nrow(m)),
    cols = as.integer(ncol(m)),
    data = as.numeric(as.vector(t(m))) # row-major
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

.fitrx_build_sir_wire <- function(fit) {
  has_any <- !is.null(fit$sir_ess) || !is.null(fit$sir_ci_theta) ||
    !is.null(fit$sir_ci_omega) || !is.null(fit$sir_ci_sigma)
  if (!has_any) return(NULL)
  list(
    ci_theta = .fitrx_ci_to_wire(fit$sir_ci_theta),
    ci_omega = .fitrx_ci_to_wire(fit$sir_ci_omega),
    ci_sigma = .fitrx_ci_to_wire(fit$sir_ci_sigma),
    ess = .fitrx_opt_num(fit$sir_ess),
    resamples_packed = NULL
  )
}

# CIs on the R side are either NULL, a 2-column matrix (rows = params), or a
# flat numeric of length 2k (legacy FFI layout). Normalise to a list of
# 2-element lists, which JSON encodes as an array of [lo, hi] pairs.
.fitrx_ci_to_wire <- function(ci) {
  if (is.null(ci) || length(ci) == 0L) return(NULL)
  if (is.matrix(ci) && ncol(ci) == 2L) {
    return(lapply(seq_len(nrow(ci)), function(i) c(ci[i, 1L], ci[i, 2L])))
  }
  flat <- as.numeric(ci)
  if (length(flat) %% 2L != 0L) return(NULL)
  pairs <- split(flat, ceiling(seq_along(flat) / 2))
  unname(pairs)
}

.fitrx_unwrap_ci <- function(wire_ci) {
  if (is.null(wire_ci) || length(wire_ci) == 0L) return(NULL)
  mat <- do.call(rbind, lapply(wire_ci, function(p) as.numeric(unlist(p, use.names = FALSE))))
  if (is.null(mat) || ncol(mat) != 2L) return(NULL)
  mat
}

.fitrx_build_iov_wire <- function(fit) {
  if (is.null(fit$omega_iov)) return(NULL)
  list(
    kappa_names = as.character(fit$kappa_names %||% character()),
    kappa_fixed = as.logical(fit$kappa_fixed %||% rep(FALSE, length(fit$kappa_names))),
    se_kappa = .fitrx_opt_num_vec(fit$se_kappa),
    shrinkage_kappa = as.numeric(fit$shrinkage_kappa %||% numeric()),
    omega_iov = .fitrx_matrix_to_wire(fit$omega_iov),
    omega_iov_param_corr = .fitrx_matrix_to_wire(fit$omega_iov_param_corr)
  )
}

.fitrx_build_eta_param_info <- function(fit) {
  eta_names <- as.character(fit$eta_names %||% character())
  if (length(eta_names) == 0L) return(list())
  types <- as.character(fit$eta_param_types %||% rep("custom", length(eta_names)))
  linked <- as.character(fit$eta_linked_theta %||% rep("", length(eta_names)))
  mapply(
    function(name, type, link) {
      list(
        eta_name = name,
        param_type = type,
        linked_theta = if (nzchar(link)) link else NULL,
        individual_param_name = "" # not exposed by the R shim yet
      )
    },
    eta_names, types, linked,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
}

# Fields that have no place in the cross-language schema but matter for R
# users (R-side derived labels, call settings, model file settings, etc.).
# Stored under fit.json:r_extras so the Rust loader silently ignores them.
.fitrx_collect_r_extras <- function(fit) {
  r_only_keys <- c(
    "call_settings", "model_file_settings", "model_structure",
    "data_name", "gradient", "gradient_used", "model_file_path"
  )
  out <- list()
  for (k in r_only_keys) {
    if (!is.null(fit[[k]])) out[[k]] <- fit[[k]]
  }
  out
}

.fitrx_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.numeric(x)
  if (is.na(v)) return(NULL)
  v
}

.fitrx_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) return(NULL)
  v
}

.fitrx_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.character(x)
  if (is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(x)
  if (length(v) == 0L) return(NULL)
  v
}

.fitrx_unwrap_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.character(x)
  if (length(v) == 0L || is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(unlist(x, use.names = FALSE))
  if (length(v) == 0L) return(NULL)
  v
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

.fitrx_unwrap_named_se <- function(se_wire, names_wire) {
  if (is.null(se_wire)) return(NULL)
  v <- as.numeric(unlist(se_wire, use.names = FALSE))
  if (length(v) == 0L) return(NULL)
  nm <- unlist(names_wire, use.names = FALSE)
  if (length(nm) == length(v)) stats::setNames(v, nm) else v
}
