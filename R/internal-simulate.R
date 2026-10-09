# Internal: pull theta/omega/sigma/omega_iov out of a ferx_fit result for FFI.
# Flattens the matrices row-major for the Rust side.
validate_fit_for_params <- function(fit) {
  if (!is.list(fit) || is.null(fit$theta) || is.null(fit$omega) || is.null(fit$sigma)) {
    stop("`fit` must be a ferx_fit result with theta, omega, and sigma components.")
  }
  theta <- as.numeric(fit$theta)
  sigma <- as.numeric(fit$sigma)
  omega <- fit$omega
  if (!is.matrix(omega) || nrow(omega) != ncol(omega)) {
    stop("`fit$omega` must be a square matrix.")
  }
  c(list(
    theta = theta,
    omega_flat = as.numeric(t(omega)),  # row-major
    omega_dim = as.integer(nrow(omega)),
    sigma = sigma
  ), .ferx_omega_iov_args(fit), list(
    # Fitted `block_sigma` residual correlations, in model declaration order.
    # A plain (non-FIX) block estimates rho, so passing them is what keeps the
    # engine from rebuilding this fit at the model file's declared correlation.
    # Empty when the fit carries none: the engine accepts that only when every
    # correlation is FIX, and refuses an estimated one (#480).
    residual_rho = .ferx_residual_rho_vec(fit)
  ), .ferx_fit_binding_args(fit))
}

# Fitted IOV (kappa) covariance, flattened row-major for the FFI. A `kappa` model
# needs it downstream: the engine draws (simulate) or conditions on (predict /
# npde) one kappa vector per occasion from it, and dropping it panicked
# ferx-core's simulate path (#1019), and the engine would rebuild the SIR and
# covariance skeletons at the model file's *initial* kappa without it, so the
# glue refuses a kappa fit that lacks it (#465, #473). NULL for a non-IOV fit,
# which maps to an empty vector + dim 0 ("no IOV") on the Rust side. The one
# source of these two arguments for every from-fit entry point.
.ferx_omega_iov_args <- function(fit) {
  omega_iov <- fit$omega_iov
  if (!is.null(omega_iov) &&
      (!is.matrix(omega_iov) || nrow(omega_iov) != ncol(omega_iov))) {
    stop("`fit$omega_iov` must be a square matrix (or NULL for a model without IOV).")
  }
  list(
    omega_iov_flat = if (is.null(omega_iov)) numeric(0) else as.numeric(t(omega_iov)),
    omega_iov_dim = if (is.null(omega_iov)) 0L else as.integer(nrow(omega_iov))
  )
}

# The data-derived bindings the fit was made with, as the one `fit_bindings`
# list the glue decodes into core's `DataBindings` (#412): the theta level-block
# layout (`fit$theta_levels`, #370) and the `[covariate_model]` statistics
# (`fit$covariate_stats`), flattened column by column. Both halves travel
# together so neither can be dropped on its way to the engine. Empty for a
# model that declares nothing data-derived - and for a fit that predates a
# field, which the engine refuses on a model that needs it. The one source of
# this argument for every from-fit entry point (predict / simulate / npde via
# `validate_fit_for_params()`, `ferx_sir()`, `ferx_covariance()`), so the
# paths cannot drift apart.
.ferx_fit_binding_args <- function(fit) {
  tl <- fit$theta_levels
  cs <- fit$covariate_stats
  list(fit_bindings = list(
    level_block = as.character(tl$block),
    level_index = as.integer(tl$index),
    level_label = as.character(tl$label),
    level_group = as.integer(tl$group),
    level_contrast = as.character(tl$contrast),
    stat_covariate = as.character(cs$covariate),
    stat_median = as.numeric(cs$median),
    stat_mean = as.numeric(cs$mean),
    stat_min = as.numeric(cs$min),
    stat_max = as.numeric(cs$max),
    stat_mode = as.numeric(cs$mode),
    stat_levels = lapply(cs$levels, as.numeric)
  ))
}

# The fit's subject IDs, verbatim and in fit order, for the skeleton FitResult
# that `ferx_sir()` / `ferx_covariance()` hand the engine, which checks them by
# position against the population it re-reads from the data (#468). Read from
# `fit$ebe_etas` - the rows the EBE matrix is built from, so the two agree by
# construction - and, on a fit without random effects (no EBE rows), from
# `fit$individual_estimates`. `fit$sdtab$ID` is numeric and not an ID source.
.ferx_fit_subject_ids <- function(fit, caller) {
  ids <- NULL
  source <- NULL
  if (!is.null(fit$ebe_etas) && nrow(fit$ebe_etas) > 0L) {
    # EBE rows are what the warm-start is built from, so their IDs are the only
    # ones known to be in the same order; never borrow another table's.
    if (is.null(fit$ebe_etas$ID)) {
      stop(
        caller, ": fit$ebe_etas has no ID column, so its rows cannot be ",
        "matched to the data's subjects.",
        call. = FALSE
      )
    }
    ids <- fit$ebe_etas$ID
    source <- "fit$ebe_etas$ID"
  } else if (!is.null(fit$individual_estimates$ID) &&
             length(fit$individual_estimates$ID) > 0L) {
    ids <- fit$individual_estimates$ID
    source <- "fit$individual_estimates$ID"
  }
  if (is.null(ids)) {
    stop(
      caller, ": the fit carries no subject IDs (neither fit$ebe_etas$ID nor ",
      "fit$individual_estimates$ID), so it cannot be matched to the data. ",
      "Re-fit via ferx_fit(model, data).",
      call. = FALSE
    )
  }
  n_na <- sum(is.na(ids))
  if (n_na > 0L) {
    stop(
      caller, ": ", n_na, " of the ", length(ids), " subject IDs in ", source,
      " are NA. The IDs are matched to the data subject by subject, so every ",
      "one must be present.",
      call. = FALSE
    )
  }
  as.character(ids)
}

# The fit's reader settings and population fingerprint (#462), as the JSON
# strings the glue reads back; "" for a fit that records neither. The one
# source of these arguments for `ferx_sir()` and `ferx_covariance()`;
# `ferx_calc_npde()` passes the settings alone (the engine cannot verify a
# fingerprint outside its own post-hoc steps).
.ferx_fit_reader_args <- function(fit) {
  list(
    reader_settings = as.character(fit$reader_settings %||% ""),
    population_fingerprint = as.character(fit$population_fingerprint %||% "")
  )
}

# How the fit was scored (#511, #472): `fit$scoring_settings` and
# `fit$sir_settings` (the record lists, NULL when the fit carries none) and
# `fit$packed_estimate` (numeric(0) when none). The glue puts them on the
# skeleton both `ferx_sir()` and `ferx_covariance()` build, so with no
# arguments those steps score what the fit scored.
.ferx_fit_record_args <- function(fit) {
  list(
    scoring_settings = fit$scoring_settings,
    sir_settings = fit$sir_settings,
    packed_estimate = as.numeric(fit$packed_estimate %||% numeric())
  )
}

# Refuse a post-hoc step on a fit that records no reader settings when it was
# read with settings its model file does not state (#462 / #416).
#
# A fit made before `fit$reader_settings` existed is re-read with the model
# file's reader settings only, so a selection passed through `ferx_fit(ignore =,
# accept =, ignore_ids =)`, or an `iov_column` passed through `settings =`,
# would be silently lost and the step would score other data than the fit.
# `fit$exclusions` names every clause that removed a row and
# `fit$call_settings` the call's `iov_column`; the glue subtracts what the
# model file states, labelled through the engine's own parse. A fit with no
# record of either shows nothing to check and runs as before.
#
# `model_path` / `model_hash` are the file the step reads and the hash it must
# have: the fit's own (`fit$model_path`, `fit$model_hash`), or a `model =`
# override with no hash to check. The glue reads that file once and refuses an
# edited, unreadable or unparseable one with core's own text, so the comparison
# is never made against the wrong file - `ferx_calc_npde()` has no hash check of
# its own to say so (#526 review 1). `data` names a `data =` override, which
# cannot be checked against the fit's selection from here. `remedy` is the entry
# point's own way to get the step from a fresh fit.
.ferx_refuse_unrecorded_selection <- function(fit, caller, remedy,
                                              model_path = fit$model_path,
                                              model_hash = fit$model_hash,
                                              data = NULL) {
  if (!is.null(fit$reader_settings)) {
    return(invisible(NULL))
  }
  fired_ignore <- as.character(fit$exclusions$fired_ignore %||% character())
  fired_accept <- as.character(fit$exclusions$fired_accept %||% character())
  call_iov <- fit$call_settings$iov_column
  call_iov <- if (is.character(call_iov) && length(call_iov) == 1L && !is.na(call_iov)) {
    call_iov
  } else {
    ""
  }
  if (length(fired_ignore) + length(fired_accept) == 0L && !nzchar(call_iov)) {
    return(invisible(NULL))
  }
  if (is.null(model_path) || is.na(model_path) || !nzchar(model_path)) {
    return(invisible(NULL))
  }
  model_path <- normalizePath(model_path, mustWork = FALSE)
  unstated <- .ferx_engine_call(
    ferx_rust_unstated_selection(
      entry_point = caller,
      model_path = model_path,
      model_hash = .ferx_hash_arg(model_hash),
      fired_ignore = fired_ignore,
      fired_accept = fired_accept,
      call_iov_column = call_iov
    ),
    model_path, data %||% fit$data_path
  )
  if (length(unstated) == 0L) {
    return(invisible(NULL))
  }
  subject_prefix <- "ignore_subjects: "
  iov_prefix <- "iov_column: "
  is_subject <- startsWith(unstated, subject_prefix)
  is_iov <- startsWith(unstated, iov_prefix)
  is_clause <- !is_subject & !is_iov
  parts <- character()
  if (any(is_clause)) {
    parts <- c(parts, paste0("`", unstated[is_clause], "`"))
  }
  if (any(is_subject)) {
    ids <- substring(unstated[is_subject], nchar(subject_prefix) + 1L)
    parts <- c(parts, paste0(
      "subject(s) ", paste(ids, collapse = ", "),
      " (from `ignore_ids` / `ignore_subjects`)"
    ))
  }
  if (any(is_iov)) {
    col <- substring(unstated[is_iov], nchar(iov_prefix) + 1L)
    parts <- c(parts, paste0("`settings = list(iov_column = \"", col, "\")`"))
  }
  data_note <- if (!is.null(data)) {
    paste0(
      " `data = \"", data, "\"` cannot be checked against the fit's selection ",
      "from here either."
    )
  } else {
    ""
  }
  stop(
    caller, ": this fit predates the record of the data selection it was ",
    "read with, and it was read with reader settings its model file does not ",
    "state: ", paste(parts, collapse = "; "), ". Re-reading the data with the ",
    "model file's settings only would not reproduce the data the fit was ",
    "scored on.", data_note, " ", remedy,
    call. = FALSE
  )
}

# The one constructor of `fit$theta_levels` (#370), shared by `ferx_fit()` and
# `ferx_load_fit()` so a fresh fit and a reloaded one are `identical()`. Takes
# the columns as plain vectors; a NULL column reads as zero rows. `value` (#430)
# is the engine's fitted value per level; NULL there on a non-empty table (a
# bundle saved before ferx recorded it) reads as NA, unknown, never 0.
.ferx_theta_levels_frame <- function(block = NULL, index = NULL, label = NULL,
                                     group = NULL, contrast = NULL,
                                     theta_name = NULL, value = NULL) {
  if (is.null(value)) value <- rep(NA_real_, length(block))
  data.frame(
    block = as.character(block),
    index = as.integer(index),
    label = as.character(label),
    group = as.integer(group),
    contrast = as.character(contrast),
    theta_name = as.character(theta_name),
    value = as.numeric(value),
    stringsAsFactors = FALSE
  )
}

# The one constructor of `fit$covariate_stats` (#412), shared by `ferx_fit()`
# and `ferx_load_fit()` so a fresh fit and a reloaded one are `identical()`:
# one row per covariate a `[covariate_model]` relation reads, with the
# statistics the fit's symbolic centres resolved against, and `levels` a list
# column of the distinct values (`levels = auto`). A NULL column reads as zero
# rows.
.ferx_covariate_stats_frame <- function(covariate = NULL, median = NULL,
                                        mean = NULL, min = NULL, max = NULL,
                                        mode = NULL, levels = NULL) {
  out <- data.frame(
    covariate = as.character(covariate),
    median = as.numeric(median),
    mean = as.numeric(mean),
    min = as.numeric(min),
    max = as.numeric(max),
    mode = as.numeric(mode),
    stringsAsFactors = FALSE
  )
  out$levels <- lapply(unname(as.list(levels)), as.numeric)
  out
}

# The fitted residual correlations as a bare numeric vector for the FFI.
# Reads the tidy frame ferx_fit() builds, and stays empty when the fit carries
# none (NULL or zero rows). The engine decides what empty means, because only it
# has the model's FIX flags: declared values when every correlation is FIX, a
# refusal naming the estimated ones otherwise (#480).
.ferx_residual_rho_vec <- function(fit) {
  rc <- fit$residual_correlations
  if (is.data.frame(rc) && nrow(rc) > 0L) return(as.numeric(rc$rho))
  numeric(0)
}

# The prior half of a fit's objective, for the FFI.
#
# `fit$ofv` is the *penalized* total once the model declares a `prior(...)`
# (ferx-core #254). Both the SIR and the standalone covariance binding rebuild a
# skeleton FitResult from the primitives flattened out of the fit, and SIR takes
# its reference objective as `ofv - ofv_prior` - so handing it a penalized `ofv`
# with a zero prior half counts the penalty twice (ferx-r #366).
#
# Returns 0 for an unpriored fit, and for one that predates the field (an older
# .fitrx bundle), which is the truth in both cases: no prior was applied, so
# `ofv` is already the data half.
.ferx_ofv_prior <- function(fit) {
  v <- suppressWarnings(as.numeric(fit$ofv_prior %||% 0))
  if (length(v) != 1L || !is.finite(v)) 0 else v
}
