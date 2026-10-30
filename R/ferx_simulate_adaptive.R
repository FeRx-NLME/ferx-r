#' Simulate state-reactive (adaptive / feedback) dosing
#'
#' Runs a forward simulation in which the dosing regimen is decided at run time
#' by the model file's \code{[adaptive_dosing]} block - a declarative
#' first-matching-rule controller that reads the simulated (optionally
#' assay-noised) state at each decision time and titrates the next dose. This is
#' the file-driven therapeutic-drug-monitoring / feedback-dosing path; the plain
#' \code{\link{ferx_simulate}} replays a fixed regimen from the data instead.
#'
#' The data may carry a \strong{pre-scheduled base regimen} - ordinary dose rows
#' (\code{EVID = 1} or \code{4}, with \code{AMT}, and optionally \code{RATE},
#' \code{SS}, \code{II}) such as a loading or maintenance dose. The controller
#' then \emph{augments} that standing prescription: at each decision time it
#' reads the simulated (optionally assay-noised) state and appends its own dose.
#' A dose-free base subject (observation grid only) is equally valid - the
#' controller simply supplies the entire regimen. A model without an
#' \code{[adaptive_dosing]} block is an error.
#'
#' \strong{The dose ledger (\code{$doses}) holds controller-issued doses only.}
#' Pre-scheduled base doses never appear in it, so \code{$metrics$CUM_DOSE} and
#' \code{$metrics$N_DOSES} \emph{exclude} the base regimen (e.g. a 1500 mg
#' loading dose is not counted). The full realized regimen is the union of
#' \code{$doses} with the input CSV's \code{EVID = 1/4} rows. In contrast
#' \code{$trajectories}, \code{$decisions}, \code{PCT_TIME_IN_WINDOW}, and the
#' \code{auc_target} attainment metric \emph{do} reflect the base regimen. A base
#' dose scheduled exactly at a decision time is read \strong{pre-dose} (the
#' trough), and a base dose scheduled after a \code{stop} still lands - a
#' \code{stop} ends controller augmentation, not the standing prescription.
#'
#' \strong{System resets.} An \code{EVID = 3} row zeroes the compartments (and
#' re-seeds any \code{init(state) = expr}) at its time; \code{EVID = 4} does the
#' same and then doses. Ordering at a shared time is reset, then dose, then
#' observation. Some combinations remain rejected with a typed error rather than
#' silently mis-integrated: a base regimen combined with a reset \emph{under} a
#' time-varying covariate or IOV; a steady-state, lag-timed, or modeled-rate
#' (\code{RATE}) base dose \emph{under} a time-varying covariate or IOV (these
#' are supported on the constant-covariate path); and an \code{auc_target} metric
#' combined with a time-varying covariate, a \code{TIME}-dependent PK parameter,
#' IOV, or a reset.
#'
#' @param model Path to a .ferx model file containing an \code{[adaptive_dosing]}
#'   block
#' @param data Path to a NONMEM-format CSV, or a data.frame (see \code{\link{ferx_fit}}): an observation grid, any covariates,
#'   and optionally a pre-scheduled base regimen (dose rows). When omitted, the
#'   data \code{fit} was fitted to (when \code{fit} is given and its data can
#'   still be found), otherwise the model file's \code{[data]} block
#'   (\code{path = ...}).
#'   The \code{DV} column may be left empty (\code{.} / \code{NA}) on the
#'   sampling rows - the DV is what the simulation produces, so an empty cell
#'   means "simulate here" (a placeholder value is not needed). Rows marked
#'   \code{MDV = 1} are excluded, as always. Kept empty-DV records are reported
#'   in the \code{simulation_warnings} attribute and re-emitted as an R warning:
#'   \code{ferx_fit()} skips those same records, so simulated rows at those
#'   times have no counterpart in a fit's \code{sdtab} (do not overlay the two,
#'   e.g. in a VPC). Every other data-reader diagnostic the engine raises for
#'   \code{data} - a dose that never landed, a covariate with no value for some
#'   subjects, and so on - travels the same channel, as do the model/data and
#'   ODE-solver findings \code{ferx_simulate()} reports. A dataset with no dose
#'   rows raises no \code{W_NO_DOSES} here: the controller supplies the regimen.
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#' @param verify Run the frozen-schedule replay verifier after every
#'   (subject, replicate) (default \code{TRUE}). A divergence is a hard error,
#'   not a warning: it means the realized dose ledger does not reproduce the
#'   trajectory, so the result is not trustworthy.
#' @param max_decisions Per-run cap on the number of decision points (the
#'   runaway / closed-loop guard). \code{0} (default) keeps the engine default.
#' @param fit Optional \code{ferx_fit} result. When provided, the simulation
#'   runs at the fitted estimates (\code{fit$theta}, \code{fit$omega},
#'   \code{fit$sigma}, the IOV omega when the model has one, and the fitted
#'   \code{block_sigma} correlations) instead of the model file's initial
#'   values, as in \code{\link{ferx_simulate}}. A theta
#'   level block is laid out as the fit laid it out (\code{fit$theta_levels}),
#'   so a design that observes only some of the fit's levels reads each level's
#'   own fitted value; a level the fit never observed is an error. Random
#'   effects are still drawn from the fitted omega - a subject's own
#'   individual estimates are not carried over, and the per-occasion kappas
#'   are drawn per decision window, so the fit's IOV occasion rule plays no
#'   part.
#'
#' @return A list of four data frames:
#'   \describe{
#'     \item{\code{trajectories}}{DRAW, SIM, ID, TIME, CMT, IPRED, DV_SIM,
#'       OBSERVED - the per-observation predictions, as in
#'       \code{\link{ferx_simulate}} (\code{OBSERVED} is \code{NA} here: the
#'       adaptive path has no time-to-event endpoint). A \code{[mixture]}
#'       model adds an integer \code{MIXNUM} column, the class each subject
#'       was drawn into for that replicate; it is absent otherwise.}
#'     \item{\code{doses}}{The realized dose ledger, \strong{controller-issued
#'       doses only} (pre-scheduled base doses are excluded): DRAW, SIM, ID, TIME,
#'       AMT, CMT, RATE, DECISION, SIGNAL (the value the controller titrated on -
#'       assay-noised when \code{with_assay_error}), RULE (the rule that fired).}
#'     \item{\code{decisions}}{One row per decision including holds: DRAW, SIM,
#'       ID, DECISION, TIME, SIGNAL, OUTCOME ("dosed" / "hold" / "stop"),
#'       N_DOSED.}
#'     \item{\code{metrics}}{One row per realized run (subject x draw x
#'       replicate) - the per-subject outcome summary: DRAW, SIM, ID, CUM_DOSE
#'       (sum of controller-issued doses; excludes the base regimen), N_DOSES
#'       (controller doses only), N_INCREASES / N_DECREASES (realized
#'       dose-change counts, by dose delta), N_HOLDS, DISCONTINUED,
#'       TIME_TO_DISCONT, SIGNAL_MIN / SIGNAL_MAX / SIGNAL_MEAN (over the
#'       decisions that recorded a signal), and PCT_TIME_IN_WINDOW - the fraction
#'       of signal-bearing decisions whose observed value fell inside the model's
#'       \code{[adaptive_dosing]} \code{target_window} (NA when no
#'       \code{target_window} is declared).}
#'   }
#'
#' @inheritSection ferx_simulate Errors raised by the engine
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("adaptive_tdm")
#' res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 20L, seed = 1L)
#' head(res$doses)         # how the dose was titrated per decision
#' head(res$trajectories)  # the resulting concentration-time profiles
#' }
#'
#' @family simulation
#' @export
ferx_simulate_adaptive <- function(model, data = NULL, n_sim = 1L, seed = 42L,
                                   verify = TRUE, max_decisions = 0L,
                                   fit = NULL) {
  if (is.null(data)) data <- .ferx_default_data(model, fit)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  data <- .ferx_data_arg(data)
  stopifnot(file.exists(model), file.exists(data))
  if (length(verify) != 1L || is.na(verify) || !is.logical(verify)) {
    stop("`verify` must be a single TRUE or FALSE.", call. = FALSE)
  }
  # Validate `n_sim` here: the Rust side casts it to an unsigned count, so a
  # negative value would wrap to an enormous number of replicates.
  n_sim <- as.integer(n_sim)
  if (length(n_sim) != 1L || is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a single positive integer.", call. = FALSE)
  }
  # A refusal is an R error, classed like a refused `ferx_fit()` (#385, #390).
  # This glue prefixes its parse / data-read messages with the function's own
  # name ("ferx_simulate_adaptive: error parsing model: ..."), which the
  # single-error fallback in `.ferx_engine_call()` allows for; a controller or
  # regimen refusal the validation pass cannot name still gets the engine's
  # prose and no code, never another finding's code.
  verify_arg <- if (isTRUE(verify)) "true" else "false"
  res <- if (is.null(fit)) {
    .ferx_engine_call(
      ferx_rust_simulate_adaptive(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        n_sim = n_sim,
        seed = as.integer(seed),
        verify = verify_arg,
        max_decisions = as.integer(max_decisions)
      ),
      model, data
    )
  } else {
    # The fitted estimates, placed on the fit's theta layout exactly as
    # `ferx_simulate(fit = )` places them (ferx-core #1594). No occasion
    # rule is passed: the adaptive engine draws kappa per decision window.
    fit_pieces <- validate_fit_for_params(fit)
    .ferx_engine_call(
      ferx_rust_simulate_adaptive_from_fit(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        theta = fit_pieces$theta,
        omega_flat = fit_pieces$omega_flat,
        omega_dim = fit_pieces$omega_dim,
        sigma = fit_pieces$sigma,
        omega_iov_flat = fit_pieces$omega_iov_flat,
        omega_iov_dim = fit_pieces$omega_iov_dim,
        residual_rho = fit_pieces$residual_rho,
        residual_rho_pairs = fit_pieces$residual_rho_pairs,
        mixture_overrides = fit_pieces$mixture_overrides,
        fit_bindings = fit_pieces$fit_bindings,
        n_sim = n_sim,
        seed = as.integer(seed),
        verify = verify_arg,
        max_decisions = as.integer(max_decisions)
      ),
      model, data
    )
  }

  # Same `simulation_warnings` channel `ferx_simulate()` uses - here it carries
  # the engine's warning bundle for the run (ferx-r #426).
  .ferx_surface_sim_warnings(res, "ferx_simulate_adaptive")
}
