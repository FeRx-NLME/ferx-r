#' Simulation-based NPDE / NPD diagnostics from a fit
#'
#' Computes Normalized Prediction Distribution Errors (\code{NPDE}, decorrelated
#' within subject) and Normalized Prediction Discrepancies (\code{NPD}) post-hoc
#' by Monte-Carlo simulation under the fitted model (Brendel et al. 2006; Comets
#' et al. 2008). Use this when a model was fitted without
#' \code{[fit_options] npde_nsim} and you want the diagnostics without re-running
#' \code{\link{ferx_fit}}. Unlike CWRES, NPDE/NPD are robust to model
#' nonlinearity and non-Gaussian random effects, and follow N(0, 1) under a
#' correctly specified model.
#'
#' The data are re-read with the selection the fit was made with
#' (\code{fit$reader_settings}): a \code{ferx_fit(ignore =, accept =,
#' ignore_ids =)} selection applies as well as the model file's
#' \code{[data_selection]}, also to a file passed as \code{data}. On a model
#' passed as \code{model}, only the fit's selection is replayed: that file's
#' own \code{[data]}, \code{[covariates]} and \code{iov_column} decide how the
#' data are read. A fit made before the selection was recorded is refused when
#' its \code{fit$exclusions} show a clause, or its \code{fit$call_settings} an
#' \code{iov_column}, that the model file does not state; that model file
#' (\code{model}, or \code{fit$model_path}, which must still have the fit's
#' hash) is the one compared.
#'
#' When the model and data are the fit's own (\code{model} and \code{data}
#' left \code{NULL}, or naming \code{fit$model_path} / \code{fit$data_path}),
#' the re-read data are checked against the fit's
#' \code{fit$population_fingerprint}, as
#' \code{\link{ferx_covariance}} and \code{\link{ferx_sir}} check theirs: a data
#' file edited since the fit is refused with code
#' \code{E_POPULATION_MISMATCH}, saying what differs.
#'
#' @param fit A \code{ferx_fit} result, carrying \code{theta}, \code{omega},
#'   \code{sigma}, and (unless overridden) the \code{model_path} / \code{data_path}
#'   captured at fit time.
#' @param nsim Number of Monte-Carlo replicates per subject (default \code{1000}).
#'   NPDE needs \code{nsim} greater than each subject's observation count for a
#'   full-rank simulated covariance; subjects that fail this get \code{NA} NPDE
#'   (NPD is still computed).
#' @param seed Optional integer RNG seed for reproducibility. \code{NULL}
#'   (default) uses the engine's built-in default seed.
#' @param model Path to the \code{.ferx} model file. Defaults to
#'   \code{fit$model_path}.
#' @param data Path to the NONMEM-format CSV, or a data.frame (see
#'   \code{\link{ferx_fit}}). Defaults to the fit's own data.
#'
#' @return A new fit with new bytes, the input \code{fit} being left as it was,
#'   with \code{NPDE} and \code{NPD} columns added to its \code{sdtab} data
#'   frame (replacing any existing ones). They are computed for the Gaussian
#'   observations only: on the discrete (e.g. binary) and time-to-event rows of
#'   \code{sdtab} they are \code{NA}. They are in the fit's bundle, so they
#'   survive \code{\link{ferx_save_fit}}. A fit with no bundle, such as one from
#'   an \code{.rds} written by an older ferx, is refused with an error of class
#'   \code{ferx_no_bundle_error}: refit it. Because the diagnostics live in
#'   \code{fit$sdtab}, downstream consumers such as
#'   \code{\link{ferx_xpose}} and goodness-of-fit plots pick them up
#'   automatically.
#'
#' @inheritSection ferx_simulate Errors raised by the engine
#' @inheritSection ferx_simulate Theta level blocks
#'
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' fit <- ferx_calc_npde(fit, nsim = 1000L, seed = 12345L)
#' head(fit$sdtab[, c("ID", "TIME", "NPDE", "NPD")])
#'
#' @family simulation
#' @export
ferx_calc_npde <- function(fit, nsim = 1000L, seed = NULL, model = NULL, data = NULL) {
  # Reject a non-fit before asking it for bytes.
  validate_fit_for_params(fit)
  .ferx_fit_bundle(fit, "ferx_calc_npde")

  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim <= 0L) {
    stop("`nsim` must be a single positive integer.")
  }
  # -1 is the FFI sentinel for "use the engine default seed". Reject negative
  # user seeds: the engine maps any negative value to the default, so they would
  # silently collide rather than seed distinct draws.
  seed_int <- if (is.null(seed)) -1L else as.integer(seed)
  if (length(seed_int) != 1L || is.na(seed_int)) {
    stop("`seed` must be a single integer or NULL.")
  }
  if (!is.null(seed) && seed_int < 0L) {
    stop("`seed` must be a non-negative integer (or NULL for the engine default).")
  }

  # Everything below reads the fit as its bytes say it is, not the caller's list:
  # diagnostics computed at an edited theta would land on a fit that does not have it.
  fit <- .ferx_fit_canonical(fit, "ferx_calc_npde")
  fit_pieces <- validate_fit_for_params(fit)
  if (is.null(fit$sdtab) || !is.data.frame(fit$sdtab) || nrow(fit$sdtab) == 0L) {
    stop("`fit$sdtab` is empty; cannot attach NPDE/NPD. Refit so the fit carries an sdtab.")
  }

  # A `model =` override keeps its own reader settings and has no recorded
  # hash; only the fit's selection is replayed on it (#526 review 3, 4).
  model_is_fits <- is.null(model) ||
    identical(normalizePath(model, mustWork = FALSE),
              normalizePath(fit$model_path %||% "", mustWork = FALSE))
  data <- .ferx_data_arg(data)
  data_override <- data
  # `data =` naming the fit's own file is still the fit's data; so is a frame
  # whose CSV text is the fit's `data_source`, whatever path it was written to
  # (a later session, a loaded fit, #419).
  data_is_fits <- is.null(data) ||
    identical(normalizePath(data, mustWork = FALSE),
              normalizePath(fit$data_path %||% "", mustWork = FALSE)) ||
    (!is.null(fit$data_source) && file.exists(data) &&
       identical(.ferx_read_data_text(data), fit$data_source))
  model <- model %||% fit$model_path
  data  <- data  %||% .ferx_fit_data_path(fit)
  if (is.null(model) || is.na(model) || !file.exists(model)) {
    stop("No usable model file: pass `model=` or refit so `fit$model_path` is set.")
  }
  if (is.null(data) || is.na(data) || !file.exists(data)) {
    stop("No usable data file: pass `data=` or refit so `fit$data_path` is set.")
  }

  # On the fit's own model and data the engine checks the re-read population
  # against the fit's fingerprint (ferx-core #1814): a data file edited since
  # the fit is refused with E_POPULATION_MISMATCH. That includes `model =` /
  # `data =` naming the fit's own files. Another model or data file is not
  # the fit's population by construction, so it is not checked.
  reader_args <- .ferx_fit_reader_args(fit)
  fingerprint <- if (model_is_fits && data_is_fits) {
    reader_args$population_fingerprint
  } else {
    ""
  }

  .ferx_refuse_unrecorded_selection(
    fit, "ferx_calc_npde",
    paste0(
      "Refit with `ferx_fit(..., settings = list(npde_nsim = ", nsim, "))`, ",
      "or refit and call ferx_calc_npde() on the new fit."
    ),
    model_path = model,
    model_hash = if (model_is_fits) fit$model_hash else "",
    data = data_override,
    data_path = data,
    fingerprint = fingerprint
  )

  # A refusal (bad params, unreadable data, ...) is an R error, classed like a
  # refused `ferx_fit()` (#385).
  npde_tbl <- .ferx_engine_call(
    ferx_rust_npde_from_fit(
      model_path = normalizePath(model),
      data_path  = normalizePath(data),
      theta      = fit_pieces$theta,
      omega_flat = fit_pieces$omega_flat,
      omega_dim  = fit_pieces$omega_dim,
      sigma      = fit_pieces$sigma,
      omega_iov_flat = fit_pieces$omega_iov_flat,
      omega_iov_dim  = fit_pieces$omega_iov_dim,
      residual_rho = fit_pieces$residual_rho,
      residual_rho_pairs = fit_pieces$residual_rho_pairs,
      mixture_overrides = fit_pieces$mixture_overrides,
      fit_bindings = fit_pieces$fit_bindings,
      # The fit's own selection, replayed on whichever file is read (#416).
      reader_settings = reader_args$reader_settings,
      model_is_fits = model_is_fits,
      population_fingerprint = fingerprint,
      # The fit's IOV occasion rule: a `dose` / `time(...)` rule derives the
      # occasion labels the kappa reference draw needs (ferx-core #1980).
      iov_occasion = as.character(fit$iov_occasion %||% ""),
      nsim       = nsim,
      seed       = seed_int
    ),
    model, data
  )
  # The glue raises rather than returning NULL, so this only guards the alignment
  # step against a table of the wrong shape.
  if (is.null(npde_tbl) || !is.data.frame(npde_tbl)) {
    stop("ferx_calc_npde: the engine returned no NPDE table.", call. = FALSE)
  }
  # The fingerprint check's notes, e.g. a fingerprint from another ferx
  # version's encoding that could not be compared.
  for (note in attr(npde_tbl, "npde_notes", exact = TRUE)) {
    warning("ferx_calc_npde: ", note, call. = FALSE)
  }

  # A new fit from new bytes, the columns in them; `fit` is left as it was.
  sdtab <- .ferx_attach_npde(fit$sdtab, npde_tbl,
                             gaussian_rows = .fitrx_sdtab_blocks(fit)$gaussian)
  .ferx_fit_with_bundle(
    fit, ferx_rust_bundle_with_sdtab(fit$bundle, lapply(sdtab, as.numeric))
  )
}

# Internal: splice NPDE/NPD onto an sdtab by position.
#
# The engine emits the NPDE table in exactly the same subject/observation order
# as the Gaussian rows of `fit$sdtab` (both iterate the same freshly read
# population) and emits ID and TIME the same way `io::output::sdtab` does, so a
# row-for-row positional copy is the correct alignment. We *assert* that
# invariant (equal row count, matching per-row ID and TIME) rather than
# silently re-joining: a mismatch means the `data`/`model` differ from the
# fit, and should be a clear error, not a quietly wrong or NA-filled column.
#
# The NPDE table covers the Gaussian observations only (ferx-core
# `compute_npde_npd` loops over `subject.observations`), while `fit$sdtab`
# holds the Gaussian rows, then the discrete (binary, ferx-core #900) rows,
# then the time-to-event residual rows (#901). `gaussian_rows` are the
# Gaussian rows as `.fitrx_sdtab_blocks()` finds them (#696): NPDE/NPD go on
# those, checked against the table by ID/TIME, and are NA on every other row
# (#708).
.ferx_attach_npde <- function(sdtab, npde_tbl,
                              gaussian_rows = seq_len(nrow(sdtab))) {
  gaussian_rows <- as.integer(gaussian_rows)
  if (!identical(gaussian_rows, seq_len(nrow(sdtab)))) {
    gauss <- .ferx_attach_npde(sdtab[gaussian_rows, , drop = FALSE], npde_tbl)
    sdtab$NPDE <- NA_real_
    sdtab$NPD  <- NA_real_
    sdtab$NPDE[gaussian_rows] <- gauss$NPDE
    sdtab$NPD[gaussian_rows]  <- gauss$NPD
    return(sdtab)
  }
  if (nrow(sdtab) != nrow(npde_tbl)) {
    stop(sprintf(paste0(
      "ferx_calc_npde: the NPDE table has %d row(s) but fit$sdtab has %d; cannot align. ",
      "This happens when `model`/`data` differ from the fit."),
      nrow(npde_tbl), nrow(sdtab)), call. = FALSE)
  }
  id_ok   <- isTRUE(all.equal(as.numeric(sdtab$ID),   as.numeric(npde_tbl$ID)))
  time_ok <- isTRUE(all.equal(as.numeric(sdtab$TIME), as.numeric(npde_tbl$TIME)))
  if (!id_ok || !time_ok) {
    stop("ferx_calc_npde: NPDE rows do not line up with fit$sdtab by ID/TIME; ",
         "refusing to attach possibly-misaligned diagnostics.", call. = FALSE)
  }
  sdtab$NPDE <- npde_tbl$NPDE
  sdtab$NPD  <- npde_tbl$NPD
  sdtab
}
