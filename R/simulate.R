#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error. When a \code{fit} is supplied, the fitted theta, omega,
#' and sigma replace the model file's initial values - which is the usual flow
#' after \code{\link{ferx_fit}} (e.g. for posterior-predictive checks or VPCs).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure: doses, obs times)
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#' @param fit Optional \code{ferx_fit} result. When provided, simulation uses
#'   \code{fit$theta}, \code{fit$omega}, and \code{fit$sigma} instead of the
#'   model file's initial values.
#' @param match Logical (default \code{FALSE}). When \code{TRUE}, each
#'   replicate's drawn etas are reassigned to subjects by \strong{propensity-score
#'   matching} against the subjects' fitted (posthoc) etas - optimal Mahalanobis
#'   matching under the model omega. This pairs each subject's observed
#'   dosing/sampling design with a similar drawn eta, correcting VPC bias from
#'   treatment adaptation in real-world data (e.g. longer dosing intervals for
#'   high-clearance patients). Requires \code{data} to be real observed data
#'   (every subject must have observations, so its posthoc eta can be computed).
#'   The posthoc etas are computed at the fitted parameters when \code{fit} is
#'   supplied, otherwise at the model file's initial values.
#'
#' @return A data.frame with columns: SIM, ID, TIME, IPRED, DV_SIM
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' sim <- ferx_simulate(ex$model, ex$data, n_sim = 10L, seed = 1L, fit = fit)
#' head(sim)
#'
#' # Propensity-score-matched simulation for a real-world-data VPC:
#' sim_pm <- ferx_simulate(ex$model, ex$data, n_sim = 10L, seed = 1L,
#'                         fit = fit, match = TRUE)
#'
#' @family simulation
#' @export
ferx_simulate <- function(model, data, n_sim = 1L, seed = 42L, fit = NULL,
                          match = FALSE) {
  stopifnot(file.exists(model), file.exists(data))
  propensity_match <- isTRUE(as.logical(match))

  if (is.null(fit)) {
    return(ferx_rust_simulate(
      model_path = normalizePath(model),
      data_path = normalizePath(data),
      n_sim = as.integer(n_sim),
      seed = as.integer(seed),
      propensity_match = propensity_match
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_simulate_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma,
    n_sim = as.integer(n_sim),
    seed = as.integer(seed),
    propensity_match = propensity_match
  )
}

#' Population predictions from a NLME model
#'
#' Computes population-level predictions (eta = 0) for all subjects.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV
#' @param fit Optional \code{ferx_fit} result. When provided, predictions use
#'   \code{fit$theta} instead of the model file's initial estimate for theta.
#'
#' @return A data.frame with columns: ID, TIME, PRED
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' preds <- ferx_predict(ex$model, ex$data, fit = fit)
#' head(preds)
#'
#' @family simulation
#' @export
ferx_predict <- function(model, data, fit = NULL) {
  stopifnot(file.exists(model), file.exists(data))

  if (is.null(fit)) {
    return(ferx_rust_predict(
      model_path = normalizePath(model),
      data_path = normalizePath(data)
    ))
  }

  fit_pieces <- validate_fit_for_params(fit)
  ferx_rust_predict_from_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    theta = fit_pieces$theta,
    omega_flat = fit_pieces$omega_flat,
    omega_dim = fit_pieces$omega_dim,
    sigma = fit_pieces$sigma
  )
}

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
#' fit <- ferx_npde(fit, nsim = 1000L, seed = 12345L)
#' head(fit$sdtab[, c("ID", "TIME", "NPDE", "NPD")])
#'
#' @family simulation
#' @export
ferx_npde <- function(fit, nsim = 1000L, seed = NULL, model = NULL, data = NULL) {
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
    stop("ferx_npde: the engine returned no NPDE table (see the message above).",
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
      "ferx_npde: the NPDE table has %d row(s) but fit$sdtab has %d; cannot align. ",
      "This happens when `model`/`data` differ from the fit, or for non-Gaussian ",
      "(e.g. TTE) endpoints, which are not yet supported."),
      nrow(npde_tbl), nrow(sdtab)), call. = FALSE)
  }
  id_ok   <- isTRUE(all.equal(as.numeric(sdtab$ID),   as.numeric(npde_tbl$ID)))
  time_ok <- isTRUE(all.equal(as.numeric(sdtab$TIME), as.numeric(npde_tbl$TIME)))
  if (!id_ok || !time_ok) {
    stop("ferx_npde: NPDE rows do not line up with fit$sdtab by ID/TIME; ",
         "refusing to attach possibly-misaligned diagnostics.", call. = FALSE)
  }
  sdtab$NPDE <- npde_tbl$NPDE
  sdtab$NPD  <- npde_tbl$NPD
  sdtab
}

# Internal: pull theta/omega/sigma out of a ferx_fit result for FFI.
# Flattens omega row-major for the Rust side.
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
  list(
    theta = theta,
    omega_flat = as.numeric(t(omega)),  # row-major
    omega_dim = as.integer(nrow(omega)),
    sigma = sigma
  )
}

#' Simulate with parameter-uncertainty propagation
#'
#' For each parameter set drawn from the uncertainty distribution, the
#' per-subject random-effect / residual-error simulator runs
#' \code{n_sim_per_draw} times. The result includes both individual
#' variability (etas, epsilons) and parameter uncertainty - useful for
#' uncertainty-aware VPCs, dose-recommendation intervals, and any analysis
#' where treating the ML estimates as fixed would understate variability.
#'
#' Two uncertainty sources are supported:
#' \describe{
#'   \item{\code{method = "asymptotic"}}{Multivariate normal around the ML
#'     estimate in the engine's packed (log-theta, Cholesky-omega, log-sigma)
#'     parameter space, using \code{fit$cov_matrix}. Requires \code{fit}
#'     to come from a \code{ferx_fit()} call with \code{covariance = TRUE}.}
#'   \item{\code{method = "sir"}}{Sample with replacement from
#'     \code{fit$sir_resamples}. Requires the fit to have been run with
#'     \code{sir = TRUE} and \code{sir_keep_samples = TRUE} (passed via
#'     \code{settings}).}
#' }
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure)
#' @param fit A \code{ferx_fit} result. Must carry either \code{cov_matrix}
#'   (asymptotic) or \code{sir_resamples} (SIR) depending on \code{method}.
#' @param n_uncertainty_draws Number of parameter sets to draw from the
#'   uncertainty distribution
#' @param n_sim_per_draw Number of eta/eps replicates per parameter draw
#' @param method Either \code{"asymptotic"} (default) or \code{"sir"}
#' @param seed Random seed for reproducibility
#'
#' @return A data.frame with columns: DRAW, SIM, ID, TIME, IPRED, DV_SIM.
#'   Row count: \code{n_uncertainty_draws * n_sim_per_draw * n_obs}.
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE)
#'
#' # Asymptotic (default): fast, MVN draws around the ML estimate
#' sims <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10
#' )
#' head(sims)         # SIM, ID, TIME, IPRED, DV_SIM
#' range(sims$DV_SIM) # sanity check on simulated range
#'
#' # SIR method: requires sir = TRUE + sir_keep_samples = TRUE at fit time
#' fit_sir <- ferx_fit(ex$model, ex$data,
#'                     covariance = TRUE, sir = TRUE,
#'                     settings = list(sir_keep_samples = TRUE))
#' sims_sir <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit_sir,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10,
#'   method = "sir"
#' )
#'
#' # Summarise: 90% prediction interval per time point
#' pi90 <- aggregate(DV_SIM ~ TIME, data = sims,
#'                   FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
#' }
#'
#' @family simulation
#' @export
ferx_simulate_with_uncertainty <- function(model, data, fit,
                                           n_uncertainty_draws = 100L,
                                           n_sim_per_draw = 1L,
                                           method = c("asymptotic", "sir"),
                                           seed = 42L) {
  stopifnot(file.exists(model), file.exists(data))
  method <- match.arg(method)

  # Check the cheap, user-facing arg constraints before the fit-object
  # validation so callers passing `n_uncertainty_draws = 0` see the obvious
  # error rather than a follow-on complaint about `cov_matrix` / `sir_resamples`.
  if (n_uncertainty_draws < 1L) {
    stop("`n_uncertainty_draws` must be >= 1.")
  }
  if (n_sim_per_draw < 1L) {
    stop("`n_sim_per_draw` must be >= 1.")
  }

  fit_pieces <- validate_fit_for_params(fit)
  unc_pieces <- validate_fit_for_uncertainty(fit, method)

  ferx_rust_simulate_with_uncertainty(
    model_path           = normalizePath(model),
    data_path            = normalizePath(data),
    theta                = fit_pieces$theta,
    omega_flat           = fit_pieces$omega_flat,
    omega_dim            = fit_pieces$omega_dim,
    sigma                = fit_pieces$sigma,
    method               = method,
    cov_matrix_flat      = unc_pieces$cov_matrix_flat,
    cov_matrix_dim       = unc_pieces$cov_matrix_dim,
    sir_resamples_flat   = unc_pieces$sir_resamples_flat,
    sir_resamples_n      = unc_pieces$sir_resamples_n,
    sir_resamples_dim    = unc_pieces$sir_resamples_dim,
    n_uncertainty_draws  = as.integer(n_uncertainty_draws),
    n_sim_per_draw       = as.integer(n_sim_per_draw),
    seed                 = as.integer(seed)
  )
}

# Internal: pull the uncertainty payload out of a ferx_fit result, validate it
# matches the requested method, and return flat representations for FFI. The
# empty branches give the Rust side something to ignore - it switches on
# `method` and only inspects the relevant arrays.
validate_fit_for_uncertainty <- function(fit, method) {
  if (method == "asymptotic") {
    cov <- fit$cov_matrix
    if (is.null(cov) || length(cov) == 0L) {
      stop("`fit$cov_matrix` is empty - re-fit with `covariance = TRUE` for ",
           "asymptotic uncertainty.")
    }
    # `cov` is a square R matrix after `process_fit_result()`. Flatten
    # row-major to match the engine's `DMatrix::from_row_slice` reader on
    # the Rust side - same convention used for `omega` above. The
    # transpose is mathematically a no-op for any well-formed covariance
    # (which must be symmetric), but the explicit `t()` keeps the FFI
    # contract uniform with `omega_flat` and prevents a future footgun
    # if the Rust reader ever changes.
    if (!is.matrix(cov) || nrow(cov) != ncol(cov)) {
      stop("`fit$cov_matrix` must be a square matrix.")
    }
    list(
      cov_matrix_flat    = as.numeric(t(cov)),
      cov_matrix_dim     = as.integer(nrow(cov)),
      sir_resamples_flat = numeric(0),
      sir_resamples_n    = 0L,
      sir_resamples_dim  = 0L
    )
  } else {
    resamples <- fit$sir_resamples
    n <- fit$sir_resamples_n
    d <- fit$sir_resamples_dim
    if (is.null(resamples) || length(resamples) == 0L ||
        is.null(n) || n == 0L || is.null(d) || d == 0L) {
      stop("`fit$sir_resamples` is empty - re-fit with `sir = TRUE` and ",
           "`sir_keep_samples = TRUE` in `settings` for SIR uncertainty.")
    }
    list(
      cov_matrix_flat    = numeric(0),
      cov_matrix_dim     = 0L,
      sir_resamples_flat = as.numeric(resamples),
      sir_resamples_n    = as.integer(n),
      sir_resamples_dim  = as.integer(d)
    )
  }
}
