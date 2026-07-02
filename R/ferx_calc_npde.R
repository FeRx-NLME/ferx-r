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
#' @param data Path to the NONMEM-format CSV. Defaults to \code{fit$data_path}.
#'
#' @return The input \code{fit}, with \code{NPDE} and \code{NPD} columns added
#'   to its \code{sdtab} data frame (replacing any existing ones). Because the
#'   diagnostics live in \code{fit$sdtab}, downstream consumers such as
#'   \code{\link{ferx_xpose}} and goodness-of-fit plots pick them up
#'   automatically.
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
  fit_pieces <- validate_fit_for_params(fit)
  if (is.null(fit$sdtab) || !is.data.frame(fit$sdtab) || nrow(fit$sdtab) == 0L) {
    stop("`fit$sdtab` is empty; cannot attach NPDE/NPD. Refit so the fit carries an sdtab.")
  }

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

  model <- model %||% fit$model_path
  data  <- data  %||% fit$data_path
  if (is.null(model) || is.na(model) || !file.exists(model)) {
    stop("No usable model file: pass `model=` or refit so `fit$model_path` is set.")
  }
  if (is.null(data) || is.na(data) || !file.exists(data)) {
    stop("No usable data file: pass `data=` or refit so `fit$data_path` is set.")
  }

  npde_tbl <- ferx_rust_npde_from_fit(
    model_path = normalizePath(model),
    data_path  = normalizePath(data),
    theta      = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim  = fit_pieces$omega_dim,
    sigma      = fit_pieces$sigma,
    nsim       = nsim,
    seed       = seed_int
  )
  # The engine prints its error and returns NULL on failure (bad params, unreadable
  # data, ...). Surface that as a clean R error instead of letting the alignment
  # step fail cryptically on a NULL table.
  if (is.null(npde_tbl) || !is.data.frame(npde_tbl)) {
    stop("ferx_calc_npde: the engine returned no NPDE table (see the message above).",
         call. = FALSE)
  }

  fit$sdtab <- .ferx_attach_npde(fit$sdtab, npde_tbl)
  fit
}

# Internal: splice NPDE/NPD onto an sdtab by position.
#
# The engine emits the NPDE table in exactly the same subject/observation order
# as `fit$sdtab` (both iterate the same freshly read population) and emits ID and
# TIME the same way `io::output::sdtab` does, so a row-for-row positional copy is
# the correct alignment. We *assert* that invariant (equal row count, matching
# per-row ID and TIME) rather than silently re-joining: a mismatch means the
# `data`/`model` differ from the fit, or the model has non-Gaussian (e.g. TTE)
# rows that sdtab and the NPDE table count differently - both cases should be a
# clear error, not a quietly wrong or NA-filled column.
.ferx_attach_npde <- function(sdtab, npde_tbl) {
  if (nrow(sdtab) != nrow(npde_tbl)) {
    stop(sprintf(paste0(
      "ferx_calc_npde: the NPDE table has %d row(s) but fit$sdtab has %d; cannot align. ",
      "This happens when `model`/`data` differ from the fit, or for non-Gaussian ",
      "(e.g. TTE) endpoints, which are not yet supported."),
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
