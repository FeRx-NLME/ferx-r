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
#' The same simulation also gives \code{EPRED}, the Monte-Carlo expected
#' prediction (the mean of the simulated replicates, NONMEM's \code{EPRED}), and
#' \code{PWRES}, the population-weighted residual
#' \eqn{\Sigma_{sim}^{-1/2}(y - EPRED)} (Monolix's \code{PWRES}), also written as
#' \code{EWRES} (NONMEM's name for it). PWRES is the simulation-based analogue of
#' CWRES: it needs no FOCE linearization, so it is valid after any estimator, and
#' it keeps the size of an outlier that NPDE, a rank, caps. It uses only the first
#' two moments, so it is N(0, 1) under the true model only when the marginal
#' distribution of the observations is Gaussian; NPDE stays the primary
#' simulation-based check. PWRES is decorrelated with the symmetric inverse
#' square root of the simulated covariance (NONMEM's default for \code{EWRES}),
#' NPDE with its Cholesky factor (the \pkg{npde} package's convention).
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
#' @section Time-to-event records:
#' A single-event TTE record is scored against the model's population
#' predictive distribution of the event time (ferx-core #1889): with
#' \eqn{\hat F(y)} the fraction of \code{nsim} simulated replicates of the
#' subject with an event before \eqn{y}, an observed event at \eqn{t} has
#' \eqn{pd = \hat F(t)}, and a time right-censored at \eqn{c} has \eqn{pd}
#' drawn uniformly on \eqn{[\hat F(c), 1]}. The score is \eqn{\Phi^{-1}(pd)},
#' so an event earlier than the model predicts is negative. One event per
#' subject leaves nothing to decorrelate, so \code{NPDE} and \code{NPD} hold
#' the same number. The engine is called as a fit with
#' \code{settings = list(npde_nsim = nsim, npde_seed = seed)} calls it, so the
#' scores equal that fit's.
#'
#' Repeated events, competing risks, left truncation and interval-censored
#' records are not scored: their rows stay \code{NA}, the Gaussian rows are
#' scored as usual, and a \code{W_NPDE_TTE_SKIPPED} warning says why.
#'
#' @param fit A \code{ferx_fit} result, carrying \code{theta}, \code{omega},
#'   \code{sigma}, and (unless overridden) the \code{model_path} / \code{data_path}
#'   captured at fit time.
#' @param nsim Number of Monte-Carlo replicates per subject (default \code{1000}).
#'   NPDE needs \code{nsim} greater than each subject's observation count for a
#'   full-rank simulated covariance; subjects that fail this get \code{NA} NPDE
#'   and PWRES (NPD and EPRED are still computed).
#' @param seed Optional integer RNG seed for reproducibility. \code{NULL}
#'   (default) uses the engine's built-in default seed.
#' @param model Path to the \code{.ferx} model file. Defaults to
#'   \code{fit$model_path}.
#' @param data Path to the NONMEM-format CSV, or a data.frame (see
#'   \code{\link{ferx_fit}}). Defaults to the fit's own data.
#'
#' @return A new fit with new bytes, the input \code{fit} being left as it was,
#'   with \code{NPDE}, \code{NPD}, \code{EPRED}, \code{EWRES} and \code{PWRES}
#'   columns added to its \code{sdtab} data frame (replacing any existing
#'   ones). They are computed for the Gaussian observations; on the discrete
#'   (e.g. binary) rows of \code{sdtab} they are \code{NA}. On the
#'   time-to-event rows of a single-event TTE endpoint, \code{NPDE} and
#'   \code{NPD} hold the record's simulation-based npde (see the section on
#'   time-to-event records)
#'   and \code{EPRED} / \code{EWRES} / \code{PWRES} are \code{NA}. They are in
#'   the fit's bundle, so they
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
#' head(fit$sdtab[, c("ID", "TIME", "NPDE", "NPD", "EPRED", "PWRES")])
#'
#' # A single-event time-to-event model: NPDE on its TTE rows.
#' tte <- ferx_example("tte_weibull")
#' tfit <- ferx_fit(tte$model, tte$data, covariance = FALSE)
#' tfit <- ferx_calc_npde(tfit, nsim = 500L, seed = 1L)
#' summary(tfit$sdtab$NPDE)
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
  blocks <- .fitrx_sdtab_blocks(fit)
  sdtab <- .ferx_attach_npde(fit$sdtab, npde_tbl, gaussian_rows = blocks$gaussian)
  # The single-event time-to-event records' NPDE (ferx-core #1889), the same
  # number in NPDE and NPD; NULL when the model has none or the engine could not
  # score them (its W_NPDE_TTE_SKIPPED note was raised above).
  sdtab <- .ferx_attach_tte_npde(sdtab, attr(npde_tbl, "tte_npde", exact = TRUE),
                                 blocks$tte)
  # `ID` is the subjects' text (#786), which the glue writes from the fit; every
  # other column goes over as numbers.
  columns <- lapply(names(sdtab), function(name) {
    if (identical(name, "ID")) as.character(sdtab[[name]]) else as.numeric(sdtab[[name]])
  })
  names(columns) <- names(sdtab)
  .ferx_fit_with_bundle(fit, ferx_rust_bundle_with_sdtab(fit$bundle, columns))
}

# Internal: put the engine's TTE npde table (one row per scored record: ID,
# TIME, CMT, NPDE; ferx-core #1889) on the sdtab's time-to-event rows,
# `tte_rows` as `.fitrx_sdtab_blocks()` finds them, into both NPDE and NPD. The
# engine emits the records in the sdtab's TTE-row order, so the copy is
# positional, checked by ID and TIME as `.ferx_attach_npde()` checks the
# Gaussian rows. EPRED / PWRES stay NA there: a single event has no moment to
# weight. A NULL or empty table leaves the rows as they are (NA).
# The refusal when the engine's table names other subjects than fit$sdtab on
# the same rows, naming the first pair; NULL when the IDs agree. IDs are the
# subjects' text since #786, so a `data` frame read with `read.csv()` turns
# "001" into 1 and no longer matches a fit made from the file.
.ferx_npde_id_mismatch <- function(fit_ids, table_ids, what) {
  fit_ids <- as.character(fit_ids)
  table_ids <- as.character(table_ids)
  if (identical(fit_ids, table_ids)) return(NULL)
  at <- which(is.na(fit_ids) != is.na(table_ids) | fit_ids != table_ids)[1L]
  sprintf(paste0(
    "ferx_calc_npde: %s row %d is subject \"%s\" in fit$sdtab but \"%s\" in the ",
    "engine's table, so the rows do not line up by ID; refusing to attach ",
    "possibly-misaligned diagnostics. Subject IDs are compared as text: if ",
    "`data` is a data frame, read its ID column as text (for example ",
    "`read.csv(path, colClasses = c(ID = \"character\"))`) or pass the file path."),
    what, at, fit_ids[at], table_ids[at])
}

.ferx_attach_tte_npde <- function(sdtab, tte_tbl, tte_rows) {
  if (!is.data.frame(tte_tbl) || nrow(tte_tbl) == 0L) return(sdtab)
  tte_rows <- as.integer(tte_rows)
  if (length(tte_rows) != nrow(tte_tbl)) {
    stop(sprintf(paste0(
      "ferx_calc_npde: the engine scored %d time-to-event record(s) but ",
      "fit$sdtab has %d time-to-event row(s); cannot align. This happens when ",
      "`model`/`data` differ from the fit."),
      nrow(tte_tbl), length(tte_rows)), call. = FALSE)
  }
  rows <- sdtab[tte_rows, , drop = FALSE]
  id_error <- .ferx_npde_id_mismatch(rows$ID, tte_tbl$ID, "time-to-event")
  if (!is.null(id_error)) stop(id_error, call. = FALSE)
  if (!isTRUE(all.equal(as.numeric(rows$TIME), as.numeric(tte_tbl$TIME)))) {
    stop("ferx_calc_npde: time-to-event NPDE rows do not line up with ",
         "fit$sdtab by ID/TIME; refusing to attach possibly-misaligned ",
         "diagnostics.", call. = FALSE)
  }
  for (col in c("NPDE", "NPD")) {
    if (is.null(sdtab[[col]])) sdtab[[col]] <- NA_real_
    sdtab[[col]][tte_rows] <- as.numeric(tte_tbl$NPDE)
  }
  sdtab
}

# The columns the engine's NPDE table carries beyond ID/TIME, in the order
# ferx-core's sdtab writes them. A table without some of them (an older engine)
# attaches the ones it has.
.ferx_npde_columns <- c("NPDE", "NPD", "EPRED", "EWRES", "PWRES")

# Internal: splice NPDE/NPD (and EPRED / EWRES / PWRES) onto an sdtab by position.
#
# The engine emits the NPDE table in exactly the same subject/observation order
# as the Gaussian rows of `fit$sdtab` (both iterate the same freshly read
# population) and emits ID (the subject's text, #786) and TIME as `fit$sdtab`
# holds them, so a
# row-for-row positional copy is the correct alignment. We *assert* that
# invariant (equal row count, matching per-row ID and TIME) rather than
# silently re-joining: a mismatch means the `data`/`model` differ from the
# fit, and should be a clear error, not a quietly wrong or NA-filled column.
#
# The NPDE table covers the Gaussian observations only (ferx-core
# `compute_npde_npd` loops over `subject.observations`), while `fit$sdtab`
# holds the Gaussian rows, then the discrete (binary, ferx-core #900) rows,
# then the time-to-event residual rows (#901). `gaussian_rows` are the
# Gaussian rows as `.fitrx_sdtab_blocks()` finds them (#696): the columns go
# on those, checked against the table by ID/TIME, and are NA on every other
# row (#708).
.ferx_attach_npde <- function(sdtab, npde_tbl,
                              gaussian_rows = seq_len(nrow(sdtab))) {
  cols <- intersect(.ferx_npde_columns, names(npde_tbl))
  gaussian_rows <- as.integer(gaussian_rows)
  if (!identical(gaussian_rows, seq_len(nrow(sdtab)))) {
    gauss <- .ferx_attach_npde(sdtab[gaussian_rows, , drop = FALSE], npde_tbl)
    for (col in cols) {
      sdtab[[col]] <- NA_real_
      sdtab[[col]][gaussian_rows] <- gauss[[col]]
    }
    return(sdtab)
  }
  if (nrow(sdtab) != nrow(npde_tbl)) {
    stop(sprintf(paste0(
      "ferx_calc_npde: the NPDE table has %d row(s) but fit$sdtab has %d; cannot align. ",
      "This happens when `model`/`data` differ from the fit."),
      nrow(npde_tbl), nrow(sdtab)), call. = FALSE)
  }
  # Both ID columns are the subjects' text (#786): compared as text, so "S1"
  # and "S2" cannot both read NA and pass.
  id_error <- .ferx_npde_id_mismatch(sdtab$ID, npde_tbl$ID, "observation")
  if (!is.null(id_error)) stop(id_error, call. = FALSE)
  if (!isTRUE(all.equal(as.numeric(sdtab$TIME), as.numeric(npde_tbl$TIME)))) {
    stop("ferx_calc_npde: NPDE rows do not line up with fit$sdtab by ID/TIME; ",
         "refusing to attach possibly-misaligned diagnostics.", call. = FALSE)
  }
  for (col in cols) sdtab[[col]] <- npde_tbl[[col]]
  sdtab
}
