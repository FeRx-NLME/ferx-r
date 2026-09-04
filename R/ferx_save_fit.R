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
#' @param output Path to write. When the path has no file extension, the
#'   conventional \code{.fitrx} extension is appended automatically (so
#'   \code{"results/run1"} is written as \code{"results/run1.fitrx"});
#'   \code{\link{ferx_load_fit}} accepts either the bare or the extended
#'   path. An explicit extension is honoured as given.
#' @param include_data Logical. When \code{TRUE}, the input NONMEM CSV used
#'   to fit the model is embedded verbatim inside the archive (as
#'   \code{data.csv}). Requires that the file is still accessible at the
#'   path captured at fit time (\code{fit$data_path}). Default \code{FALSE}.
#' @return Invisibly returns \code{fit}, so the function can be used inside a
#'   pipe without breaking the chain:
#'   \preformatted{
#'   fit |> ferx_save_fit("run1.fitrx") |> (\(f) f$estimates)()
#'   }
#' @family utilities
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

  # Default the archive extension to `.fitrx` when the caller passes a bare
  # path with no extension, so `ferx_save_fit(fit, "results/run1")` lands at
  # "results/run1.fitrx" - the documented convention - rather than at a
  # bare "results/run1" or (via Info-ZIP's own habit) "results/run1.zip".
  # `ferx_load_fit()` mirrors this by falling back to `<path>.fitrx`.
  if (!nzchar(tools::file_ext(output))) {
    output <- paste0(output, ".fitrx")
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

  # ebes_kappa.csv - present only when the model uses IOV
  if (.fitrx_has_kappa(fit)) {
    .fitrx_write_ebes_kappa_csv(fit, file.path(staging, "ebes_kappa.csv"))
    entries <- c(entries, "ebes_kappa.csv")
  }

  # conddist.csv - present only when the SAEM conditional-distribution pass
  # ran (`conddist = TRUE`, ferx-core #257). Long format (ID, ETA, COND_MEAN,
  # COND_SD, COND_MODE), mirroring ferx-core's own `conddist.csv` .fitrx entry
  # (ferx-core #675) so bundles are interoperable, but this reimplementation
  # doesn't call into ferx-core's Rust bundle writer at all - see #244.
  if (.fitrx_has_cond_dist(fit)) {
    .fitrx_write_conddist_csv(fit, file.path(staging, "conddist.csv"))
    entries <- c(entries, "conddist.csv")
  }

  # predictions.csv
  .fitrx_write_predictions_csv(fit, file.path(staging, "predictions.csv"))
  entries <- c(entries, "predictions.csv")

  # covtab.csv - present only when the model declares a [covariates] block
  if (!is.null(fit$covtab) && nrow(fit$covtab) > 0L) {
    .fitrx_write_covtab_csv(fit, file.path(staging, "covtab.csv"))
    entries <- c(entries, "covtab.csv")
  }

  # trace.csv - present only when optimizer_trace = TRUE was used. Bundled so
  # fit$trace survives the round-trip even though fit$trace_path (a temp
  # file) usually won't. Bundle on `is.data.frame()`, not `nrow() > 0L`, so a
  # valid header-only (zero-iteration) trace still round-trips.
  if (is.data.frame(fit[["trace"]])) {
    utils::write.csv(fit[["trace"]], file.path(staging, "trace.csv"), row.names = FALSE, na = "")
    entries <- c(entries, "trace.csv")
  }

  # model.ferx
  #
  # Prefer `file.copy(fit$model_path, ...)` so the bundled bytes match
  # the original file exactly - including line endings, BOM, and the
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
  # is available (e.g. a hand-constructed fit) - in that case the bundle
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

  # manifest.json - last so it can list every entry
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
  # Info-ZIP appends `.zip` when the archive name has no extension, so a
  # bare `output` like "results/run1" would silently land at "results/run1.zip"
  # and `ferx_load_fit("results/run1")` would then fail with "File does not
  # exist". Zip into a temp `.zip` (extension present -> no rename by zip) and
  # move it onto the exact requested path afterwards.
  zip_tmp <- tempfile("fitrx_zip_", fileext = ".zip")
  old_wd <- getwd()
  setwd(staging)
  zip_status <- tryCatch(
    utils::zip(zip_tmp, files_to_zip, flags = "-q9X"),
    finally = setwd(old_wd)
  )
  if (!is.null(zip_status) && is.numeric(zip_status) && zip_status != 0L) {
    stop(sprintf("zip command failed with status %d while writing %s", zip_status, output))
  }
  if (!file.exists(zip_tmp)) {
    stop(sprintf("zip command did not produce an archive while writing %s", output))
  }
  if (!file.copy(zip_tmp, out_abs, overwrite = TRUE)) {
    unlink(zip_tmp)
    stop(sprintf("failed to move zip archive into place at %s", out_abs))
  }
  unlink(zip_tmp)

  invisible(fit)
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
    # Free-parameter tally by Delattre class (ferx-core #1177) - what makes a
    # stored fit rankable by ferx_bic() without re-parsing the model. The
    # engine reads it back with `serde(default)`, so omitting it degrades to an
    # all-zero tally rather than a load error; write it whenever we have it.
    bic_inputs = .fitrx_bic_inputs_to_wire(fit$bic_inputs),
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
    covariance_n_evals_estimated = .fitrx_opt_num(fit$covariance_n_evals_estimated),
    trace_path = .fitrx_opt_chr(fit$trace_path),
    ebe_convergence_warnings = as.integer(fit$ebe_convergence_warnings %||% 0L),
    max_unconverged_subjects = as.integer(fit$max_unconverged_subjects %||% 0L),
    total_ebe_fallbacks = as.integer(fit$total_ebe_fallbacks %||% 0L),
    warnings = as.character(fit$warnings %||% character()),
    saem_mu_ref_m_step_evals_saved = .fitrx_opt_num(fit$saem_mu_ref_m_step_evals_saved),

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
      shrinkage = as.numeric(fit$shrinkage_eta %||% numeric()),
      init_as_sd = as.logical(fit$omega_init_as_sd %||% rep(FALSE, NROW(fit$omega)))
    ),
    sigma = list(
      names = as.character(fit$sigma_names %||% character()),
      estimates = as.numeric(fit$sigma %||% numeric()),
      se = .fitrx_opt_num_vec(fit$se_sigma),
      fixed = as.logical(fit$sigma_fixed %||% rep(FALSE, length(fit$sigma %||% numeric()))),
      types = as.character(fit$sigma_types %||% rep("proportional", length(fit$sigma %||% numeric()))),
      init_as_sd = as.logical(fit$sigma_init_as_sd %||% rep(FALSE, length(fit$sigma %||% numeric())))
    ),
    error_model = .fitrx_error_model_from_sigma_types(fit$sigma_types),
    shrinkage_eps = as.numeric(fit$shrinkage_eps %||% NA_real_),
    dw_statistic = .fitrx_opt_num(fit$dw_statistic),
    iwres_lag1_r = .fitrx_opt_num(fit$iwres_lag1_r),
    covariance_matrix = .fitrx_matrix_to_wire(fit$cov_matrix),
    cov_eigenvalues = .fitrx_opt_num_vec(fit$cov_eigenvalues),
    cov_condition_number = .fitrx_opt_num(fit$cov_condition_number),
    # The outer optimizer's own init-escape verdict (ferx-core #1177), which is
    # the reading `check_strictness(reject_init_stall = TRUE)` prefers. NULL
    # when the fit recorded none.
    left_init = .fitrx_opt_lgl(fit$left_init),

    sir = .fitrx_build_sir_wire(fit),
    bayes = .fitrx_build_bayes_wire(fit),
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

    input_columns    = as.character(fit$input_columns    %||% character()),
    covariate_names  = as.character(fit$covariate_names  %||% character()),

    exclusions = if (!is.null(fit$exclusions)) {
      ex <- fit$exclusions
      list(
        n_records_total      = as.integer(ex$n_records_total  %||% 0L),
        n_obs_excluded       = as.integer(ex$n_obs_excluded   %||% 0L),
        n_dose_excluded      = as.integer(ex$n_dose_excluded  %||% 0L),
        n_other_excluded     = as.integer(ex$n_other_excluded %||% 0L),
        excluded_subject_ids = as.character(ex$excluded_subject_ids %||% character()),
        fired_ignore         = as.character(ex$fired_ignore          %||% character()),
        fired_accept         = as.character(ex$fired_accept          %||% character())
      )
    } else {
      NULL
    },

    # Named covariate -> "continuous"/"categorical" map, as a JSON object.
    # `as.list` keeps the names; NULL when no [covariates] block was declared.
    covariate_types = if (length(fit$covariate_types)) {
      as.list(fit$covariate_types)
    } else {
      NULL
    },

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

.fitrx_write_ebes_csv <- function(fit, path) {
  ebes <- fit$ebe_etas
  if (is.null(ebes) || nrow(ebes) == 0L) {
    # Header-only file so the entry shape is stable.
    cat("ID,ofv_contribution,n_obs\n", file = path)
    return(invisible())
  }
  # Pull ofv_contribution / n_obs out of sdtab (one row per subject) when not
  # already present on ebes - the Rust shim doesn't expose them on ebe_etas.
  if (!all(c("ofv_contribution", "n_obs") %in% names(ebes))) {
    per_subj <- .fitrx_per_subject_ofv_nobs(fit)
    if (!is.null(per_subj)) {
      # ebe_etas and the per-subject sdtab summary are both one row per subject
      # in subject order, so join positionally: a dataset that reuses a subject
      # ID in a non-contiguous block has two subjects sharing a textual ID, and
      # an ID match would collapse them (both ebes rows get the first subject's
      # n_obs, so ebes.csv disagrees with predictions.csv and the ferx-core
      # loader rejects the bundle). Fall back to a string-ID match when the row
      # counts don't line up (character IDs like PT001 still work; unrepresented
      # subjects fall to NA).
      idx <- if (length(per_subj$id) == nrow(ebes)) {
        seq_len(nrow(ebes))
      } else {
        match(as.character(ebes$ID), per_subj$id)
      }
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

.fitrx_write_conddist_csv <- function(fit, path) {
  utils::write.table(
    fit$cond_dist$data, path,
    row.names = FALSE, quote = FALSE, sep = ",", na = ""
  )
}

.fitrx_has_cond_dist <- function(fit) {
  cd <- fit$cond_dist
  !is.null(cd) && is.data.frame(cd$data) && nrow(cd$data) > 0L
}

.fitrx_write_predictions_csv <- function(fit, path) {
  sdtab <- fit$sdtab
  if (is.null(sdtab) || nrow(sdtab) == 0L) {
    cat("ID,TIME,DV,PRED,IPRED,CWRES,IWRES,EBE_OFV,N_OBS\n", file = path)
    return(invisible())
  }
  # sdtab$ID is the original numeric subject ID from the NONMEM CSV.
  utils::write.table(
    sdtab, path,
    row.names = FALSE, quote = FALSE, sep = ",", na = ""
  )
}

.fitrx_write_covtab_csv <- function(fit, path) {
  covtab <- fit$covtab
  if (is.null(covtab) || nrow(covtab) == 0L) {
    return(invisible())
  }
  # Quote so a free-form character ID containing a comma/quote round-trips
  # intact; missing covariate values (NA/NaN) are written as empty cells.
  utils::write.csv(covtab, path, row.names = FALSE, na = "")
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
  # First row of each subject block, keyed on the ordinal subject index rather
  # than on the raw ID - `!duplicated(sdtab$ID)` would drop the second of two
  # subjects that reuse a textual ID in a non-contiguous block.
  subj <- .ferx_subject_index(sdtab$ID)
  first <- !duplicated(subj)
  list(
    id = as.character(sdtab$ID[first]),
    ofv = as.numeric(sdtab$EBE_OFV[first]),
    n_obs = as.integer(sdtab$N_OBS[first])
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

.fitrx_covariance_status_to_token <- function(s) {
  if (is.null(s) || !nzchar(s)) return("not_requested")
  switch(tolower(s),
    "computed" = "computed",
    "failed" = "failed",
    "not_requested" = "not_requested",
    "notrequested" = "not_requested",
    "sir_fallback" = "sir_fallback",
    "sirfallback" = "sir_fallback",
    tolower(s)
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

# Bayes posterior summary (method = "bayes"). Stored as scalar metadata plus
# parallel per-parameter vectors, mirroring the in-memory `$bayes` list so the
# posterior survives a save/load round-trip (the engine's binary .fitrx does
# not carry it). NULL for non-Bayes fits.
.fitrx_build_bayes_wire <- function(fit) {
  b <- fit$bayes
  if (is.null(b)) return(NULL)
  list(
    n_chains = .fitrx_opt_num(b$n_chains),
    n_warmup = .fitrx_opt_num(b$n_warmup),
    n_draws_per_chain = .fitrx_opt_num(b$n_draws_per_chain),
    n_divergent = .fitrx_opt_num(b$n_divergent),
    max_rhat = .fitrx_opt_num(b$max_rhat),
    param_names = as.list(as.character(b$param_names %||% character())),
    mean = .fitrx_opt_num_vec(b$mean),
    sd = .fitrx_opt_num_vec(b$sd),
    q025 = .fitrx_opt_num_vec(b$q025),
    median = .fitrx_opt_num_vec(b$median),
    q975 = .fitrx_opt_num_vec(b$q975),
    rhat = .fitrx_opt_num_vec(b$rhat),
    ess_bulk = .fitrx_opt_num_vec(b$ess_bulk),
    ess_tail = .fitrx_opt_num_vec(b$ess_tail),
    mcse = .fitrx_opt_num_vec(b$mcse)
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

.fitrx_build_iov_wire <- function(fit) {
  if (is.null(fit$omega_iov)) return(NULL)
  wire <- list(
    kappa_names = as.character(fit$kappa_names %||% character()),
    kappa_fixed = as.logical(fit$kappa_fixed %||% rep(FALSE, length(fit$kappa_names))),
    se_kappa = .fitrx_opt_num_vec(fit$se_kappa),
    shrinkage_kappa = as.numeric(fit$shrinkage_kappa %||% numeric()),
    shrinkage_kappa_by_occ = if (is.data.frame(fit$shrinkage_kappa_by_occ) &&
                                   nrow(fit$shrinkage_kappa_by_occ) > 0L) {
      kap_cols <- setdiff(colnames(fit$shrinkage_kappa_by_occ), "occ")
      lapply(seq_len(nrow(fit$shrinkage_kappa_by_occ)), function(i) {
        as.numeric(fit$shrinkage_kappa_by_occ[i, kap_cols])
      })
    } else {
      list()
    },
    omega_iov = .fitrx_matrix_to_wire(fit$omega_iov),
    omega_iov_param_corr = .fitrx_matrix_to_wire(fit$omega_iov_param_corr),
    kappa_init_as_sd = as.logical(fit$kappa_init_as_sd %||% rep(FALSE, length(fit$kappa_names)))
  )
  # Sample-size-weighted IOV (ferx-core #1031). The two keys must be *absent*
  # (not JSON `null`) for an unweighted model: the engine declares them
  # `#[serde(default)] Vec<Option<..>>`, and `serde(default)` only fires for a
  # missing key - an explicit `null` fails with "invalid type: null, expected a
  # sequence". Assigning NULL into a list() literal keeps the name, so the
  # fields are appended here only when a weight exists. That also keeps the
  # bundle byte-identical to a pre-#1031 one for every unweighted model.
  # `as.list()` keeps them JSON arrays under `auto_unbox = TRUE` (a single
  # weighted kappa is the common MBMA case), and `na = "null"` maps an
  # unweighted slot to `null`, which is what `Vec<Option<..>>` expects.
  if (!is.null(fit$kappa_weights)) {
    wire$kappa_weights <- as.list(as.character(fit$kappa_weights))
    wire$kappa_weight_typical <- as.list(as.numeric(
      fit$kappa_weight_typical %||% rep(NA_real_, length(fit$kappa_weights))
    ))
  }
  wire
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
    "data_name", "gradient", "gradient_used", "model_file_path",
    # Strictness-gate ingredients (ferx-core #1177). The engine recomputes each
    # of these from a `FitResult` it holds; R has no `FitResult`, so a loaded
    # fit has to carry them or `check_strictness()` would silently skip the
    # correlation, boundary and init-stall gates on it. `bic_inputs` and
    # `left_init` are proper wire fields and are written above.
    "max_abs_correlation", "near_boundary", "boundary_estimate_message",
    "stalled_at_init"
  )
  out <- list()
  for (k in r_only_keys) {
    if (!is.null(fit[[k]])) out[[k]] <- fit[[k]]
  }
  out
}
