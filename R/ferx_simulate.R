#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error. When a \code{fit} is supplied, the fitted theta, omega,
#' and sigma replace the model file's initial values - which is the usual flow
#' after \code{\link{ferx_fit}} (e.g. for posterior-predictive checks or VPCs).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure: doses,
#'   obs times). When omitted, the model file's \code{[data]} block (\code{path
#'   = ...}) is used.
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#' @param fit Optional \code{ferx_fit} result. When provided, simulation uses
#'   \code{fit$theta}, \code{fit$omega}, and \code{fit$sigma} instead of the
#'   model file's initial values.
#' @param match Propensity-score matching method (default \code{FALSE}). When
#'   enabled, each replicate's drawn etas are reassigned to subjects by
#'   \strong{propensity-score matching} against the subjects' fitted (posthoc)
#'   etas - Mahalanobis matching under the model omega. This pairs each subject's
#'   observed dosing/sampling design with a similar drawn eta, correcting VPC
#'   bias from treatment adaptation in real-world data (e.g. longer dosing
#'   intervals for high-clearance patients). Accepts:
#'   \describe{
#'     \item{\code{FALSE} / \code{"none"}}{No matching (default).}
#'     \item{\code{TRUE} / \code{"optimal"}}{Global linear-assignment minimum
#'       Mahalanobis distance (\code{MatchIt(method = "optimal")}); best on
#'       average in simulation, the recommended method.}
#'     \item{\code{"nearest"}}{Greedy nearest-neighbour
#'       (\code{MatchIt(method = "nearest", distance = "mahalanobis")}).}
#'     \item{\code{"rank"}}{Pair by the rank of each eta's Mahalanobis norm.}
#'   }
#'   Requires \code{data} to be real observed data (every subject must have
#'   observations, so its posthoc eta can be computed). The posthoc etas are
#'   computed at the fitted parameters when \code{fit} is supplied, otherwise at
#'   the model file's initial values.
#'
#' @param horizon Optional administrative censoring time for time-to-event (TTE)
#'   endpoints. A finite, positive \code{horizon} is \strong{required} to simulate
#'   a drug-driven (joint PK-TTE) model: the augmented hazard ODE is integrated
#'   until the cumulative hazard reaches \code{-log U}, censoring at
#'   \code{horizon} if no event fires (ferx-core #564). Purely-Gaussian models
#'   ignore it. \code{NULL} (default) leaves it unset.
#' @return A data.frame. Gaussian rows carry DRAW, SIM, ID, TIME, CMT, IPRED,
#'   DV_SIM (with \code{OBSERVED = NA}). For a joint PK-TTE model each subject
#'   also yields a TTE row on the event CMT, where TIME is the sampled
#'   event/censor time and \code{OBSERVED} is 1 (event before \code{horizon}) or
#'   0 (right-censored at it); its IPRED and DV_SIM are \code{NA}. Use
#'   \code{is.na(OBSERVED)} to separate continuous rows from event rows.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' sim <- ferx_simulate(ex$model, ex$data, n_sim = 10L, seed = 1L, fit = fit)
#' head(sim)
#'
#' # Propensity-score-matched simulation for a real-world-data VPC:
#' sim_pm <- ferx_simulate(ex$model, ex$data, n_sim = 10L, seed = 1L,
#'                         fit = fit, match = "optimal")
#'
#' @family simulation
#' @export
ferx_simulate <- function(model, data = NULL, n_sim = 1L, seed = 42L, fit = NULL,
                          match = FALSE, horizon = NULL) {
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  stopifnot(file.exists(model), file.exists(data))
  match_method <- normalize_match_method(match)
  # A finite, positive `horizon` is required to sample drug-driven (joint PK-TTE)
  # event times (ferx-core #564); it is ignored by purely-Gaussian models.
  # `NULL` maps to the sentinel -1, which the Rust side reads as "unset".
  if (!is.null(horizon)) {
    if (length(horizon) != 1L || !is.numeric(horizon) || !is.finite(horizon) ||
        horizon <= 0) {
      stop("`horizon` must be a single finite positive number (or NULL).",
           call. = FALSE)
    }
  }
  horizon_arg <- if (is.null(horizon)) -1 else as.numeric(horizon)

  if (is.null(fit)) {
    return(ferx_rust_simulate(
      model_path = normalizePath(model),
      data_path = normalizePath(data),
      n_sim = as.integer(n_sim),
      seed = as.integer(seed),
      match_method = match_method,
      horizon = horizon_arg
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
    match_method = match_method,
    horizon = horizon_arg
  )
}

# Internal: normalize the user-facing `match` argument to the string token the
# Rust side expects ("none" | "optimal" | "nearest" | "rank"). Accepts a logical
# (FALSE -> "none", TRUE -> "optimal" for backward compatibility) or one of the
# method strings (case-insensitive).
normalize_match_method <- function(match) {
  if (length(match) != 1L || is.na(match)) {
    stop("`match` must be a single value: FALSE/TRUE or one of ",
         "\"none\", \"optimal\", \"nearest\", \"rank\".", call. = FALSE)
  }
  if (is.logical(match)) {
    return(if (isTRUE(match)) "optimal" else "none")
  }
  method <- tolower(as.character(match))
  valid <- c("none", "optimal", "nearest", "rank")
  if (!method %in% valid) {
    stop("`match` must be FALSE/TRUE or one of ",
         paste0("\"", valid, "\"", collapse = ", "), ".", call. = FALSE)
  }
  method
}
