#' Simulate from a NLME model
#'
#' Simulates observations from a parsed model with between-subject variability
#' and residual error. When a \code{fit} is supplied, the fitted theta, omega,
#' and sigma replace the model file's initial values - which is the usual flow
#' after \code{\link{ferx_fit}} (e.g. for posterior-predictive checks or VPCs).
#'
#' @section Which predictive distribution you get:
#' \code{ferx_simulate()} draws a \strong{fresh set of random effects for every
#' ID in \code{data}, in every replicate}, and adds residual error on top. At
#' the default \code{match = FALSE} it does not read the observed \code{DV} at
#' all (see the \code{data} argument), so each subject's draw is unconditional:
#' the \emph{prior} predictive under the estimated parameters, not a posterior
#' one.
#'
#' \code{match} is the exception, on the assignment rather than on the draw.
#' The pool of etas per replicate is the same unconditional sample with the same
#' marginal distribution, but \emph{which} member of it lands on a given
#' subject is chosen by matching against that subject's posthoc eta - computed
#' from that subject's observed \code{DV}, which is why matching requires every
#' subject to carry observations. So under \code{match} the eta paired with a
#' particular observed design \strong{is} informed by that design's data; that
#' is the point of it (it restores the design-eta association that adaptive
#' dosing puts in real-world data), and it means the per-subject draw is no
#' longer the unconditional one described above.
#'
#' What that distribution is \emph{of} follows from what one ID means in your
#' data and what level the etas were estimated at. The two are the same thing in
#' individual-level PK, where an ID is a patient - but not in a model-based
#' meta-analysis, where a row is a trial-arm summary, an ID is a study or an arm
#' and the etas are \strong{between-study} random effects. There, each replicate
#' is a set of \strong{new studies}, and \code{DV_SIM} is the predictive
#' distribution of the \emph{next trial's readout} - not of the next patient. A
#' patient-level interval needs a patient-level random effect in the model and
#' rows at patient level in the data; ferx will not manufacture one, and quoting
#' a between-study interval as a between-patient one is the error this note
#' exists to prevent.
#'
#' Within a run, pick the column that matches the question:
#' \describe{
#'   \item{\code{IPRED}}{the mean response of the newly drawn study/subject -
#'     random effects, no residual error.}
#'   \item{\code{DV_SIM}}{that study's/subject's \emph{observed} readout -
#'     random effects plus residual error. This is the column a VPC and a
#'     predictive interval are built from.}
#' }
#' Two neighbouring functions cover the other levels: \code{\link{ferx_predict}}
#' gives the typical-value curve (all etas at zero, no residual error), and
#' \code{\link{ferx_simulate_with_uncertainty}} adds \strong{parameter}
#' uncertainty (a theta/omega/sigma draw per replicate) on top of the random
#' effects, which is what a decision interval for a future trial usually wants.
#' For the individual etas conditioned on observed data, see
#' \code{fit$ebe_etas} and \code{\link{ferx_conddist}}.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure: doses,
#'   obs times). When omitted, the model file's \code{[data]} block (\code{path
#'   = ...}) is used.
#'   The \code{DV} column may be left empty (\code{.} / \code{NA}) on the
#'   sampling rows - the DV is what the simulation produces, so an empty cell
#'   means "simulate here" (a placeholder value is not needed). Rows marked
#'   \code{MDV = 1} are excluded, as always. Kept empty-DV records are reported
#'   in the \code{simulation_warnings} attribute and re-emitted as an R warning:
#'   \code{ferx_fit()} skips those same records, so simulated rows at those
#'   times have no counterpart in a fit's \code{sdtab} (do not overlay the two,
#'   e.g. in a VPC). Every other data-reader diagnostic the engine raises for
#'   \code{data} - a dose that never landed, a covariate with no value for some
#'   subjects, and so on - travels the same channel.
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
#' @return A data.frame, one row per simulated record, with these columns:
#'   \describe{
#'     \item{\code{DRAW}}{Parameter draw. Always 1 here; it is the grouping
#'       variable in \code{\link{ferx_simulate_with_uncertainty}}.}
#'     \item{\code{SIM}}{Replicate index, 1 to \code{n_sim}.}
#'     \item{\code{ID}}{Subject identifier, as a character string, exactly as
#'       the \code{ID} column of \code{data} spelled it.}
#'     \item{\code{TIME}}{Sampling time from \code{data} - or, on a TTE row,
#'       the \emph{sampled} event or censoring time (see \code{OBSERVED}).}
#'     \item{\code{CMT}}{Observation compartment: the data file's \code{CMT}
#'       for a continuous row, the \code{[event_model]} / \code{[binary_model]}
#'       compartment for an event or categorical row.}
#'     \item{\code{IPRED}}{Individual prediction - the drawn subject's mean
#'       response, random effects but no residual error.}
#'     \item{\code{DV_SIM}}{The simulated observation: \code{IPRED} plus
#'       residual error, and the column a VPC or a predictive interval is built
#'       from.}
#'     \item{\code{OBSERVED}}{Event indicator for a time-to-event row only -
#'       1 = the event fired before \code{horizon}, 0 = administratively
#'       right-censored at it - and \code{NA} on every other row. It is
#'       \strong{not} the input \code{DV}: nothing in this frame echoes the
#'       data file's observations, since simulation is what produces that
#'       column. Compare against the observed data by joining \code{data} on
#'       ID/TIME yourself.}
#'   }
#'   Gaussian rows carry DRAW, SIM, ID, TIME, CMT, IPRED and DV_SIM, with
#'   \code{OBSERVED = NA}. For a joint PK-TTE model each subject also yields a
#'   TTE row on the event CMT, whose IPRED and DV_SIM are \code{NA}; use
#'   \code{is.na(OBSERVED)} to separate continuous rows from event rows. For a
#'   \code{[binary_model]} endpoint the categorical row on the binary CMT carries
#'   the simulated 0/1 outcome in \code{DV_SIM} (coded as the input CSV codes DV),
#'   with \code{IPRED} and \code{OBSERVED} both \code{NA}; select it by its CMT.
#'
#'   The returned frame carries a \code{simulation_warnings} attribute (a
#'   character vector, empty for a clean run) listing any per-subject simulation
#'   diagnostics from ferx-core (e.g. a degenerate or pathological hazard that
#'   censored a subject with no event) together with the engine's data-reader
#'   diagnostics for \code{data} - the same ones \code{ferx_fit()} returns in
#'   \code{fit$warnings}; these are also raised as an R warning.
#'
#' @section Errors raised by the engine:
#' A model or dataset the engine refuses is an R error, never a \code{NULL}
#' result. Where the engine's validation pass names the finding, the error has
#' class \code{ferx_engine_error} and carries \code{code}, \code{block},
#' \code{line} and \code{suggestion} - the same condition a refused
#' \code{\link{ferx_fit}} raises, so one
#' \code{tryCatch(..., ferx_engine_error = )} handler covers every entry
#' point. A failure the validation pass cannot attribute to one diagnostic is
#' raised as an ordinary error carrying the engine's message.
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

  # A refusal is an R error, classed like a refused `ferx_fit()` (#385).
  res <- if (is.null(fit)) {
    .ferx_engine_call(
      ferx_rust_simulate(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        n_sim = as.integer(n_sim),
        seed = as.integer(seed),
        match_method = match_method,
        horizon = horizon_arg
      ),
      model, data
    )
  } else {
    fit_pieces <- validate_fit_for_params(fit)
    .ferx_engine_call(
      ferx_rust_simulate_from_fit(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        theta = fit_pieces$theta,
        omega_flat = fit_pieces$omega_flat,
        omega_dim = fit_pieces$omega_dim,
        sigma = fit_pieces$sigma,
        omega_iov_flat = fit_pieces$omega_iov_flat,
        omega_iov_dim = fit_pieces$omega_iov_dim,
        residual_rho = fit_pieces$residual_rho,
        n_sim = as.integer(n_sim),
        seed = as.integer(seed),
        match_method = match_method,
        horizon = horizon_arg
      ),
      model, data
    )
  }

  # Surface per-subject simulation diagnostics from ferx-core (#762/#763): a
  # degenerate or pathological hazard that would otherwise censor a subject with
  # no signal is attached by the Rust glue as a `simulation_warnings` attribute.
  .ferx_surface_sim_warnings(res)
}

# Emit the diagnostics attached by the Rust glue as a `simulation_warnings`
# character-vector attribute - ferx-core's per-subject simulation diagnostics
# (#762/#763) plus every data-reader diagnostic the engine raised for the dataset
# (#283) - as a single R warning, so they are not silently lost; return `res`
# unchanged (the attribute is left in place for programmatic access).
.ferx_surface_sim_warnings <- function(res, fn = "ferx_simulate") {
  w <- attr(res, "simulation_warnings", exact = TRUE)
  if (length(w) > 0L) {
    warning(
      fn, " produced ", length(w), " diagnostic",
      if (length(w) > 1L) "s" else "", ":\n  ",
      paste(w, collapse = "\n  "),
      call. = FALSE
    )
  }
  res
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
