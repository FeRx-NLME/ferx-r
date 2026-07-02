#' Simulate state-reactive (adaptive / feedback) dosing
#'
#' Runs a forward simulation in which the dosing regimen is decided at run time
#' by the model file's \code{[adaptive_dosing]} block - a declarative
#' first-matching-rule controller that reads the simulated (optionally
#' assay-noised) state at each decision time and titrates the next dose. This is
#' the file-driven therapeutic-drug-monitoring / feedback-dosing path; the plain
#' \code{\link{ferx_simulate}} replays a fixed regimen from the data instead.
#'
#' The base subjects must be \strong{dose-free}: the controller supplies every
#' dose, so the data provides only the observation grid and any covariates. A
#' model without an \code{[adaptive_dosing]} block is an error.
#'
#' @param model Path to a .ferx model file containing an \code{[adaptive_dosing]}
#'   block
#' @param data Path to a NONMEM-format CSV (dose-free base subjects + observation
#'   times)
#' @param n_sim Number of simulation replicates
#' @param seed Random seed for reproducibility
#' @param verify Run the frozen-schedule replay verifier after every
#'   (subject, replicate) (default \code{TRUE}). A divergence is a hard error,
#'   not a warning: it means the realized dose ledger does not reproduce the
#'   trajectory, so the result is not trustworthy.
#' @param max_decisions Per-run cap on the number of decision points (the
#'   runaway / closed-loop guard). \code{0} (default) keeps the engine default.
#'
#' @return A list of four data frames:
#'   \describe{
#'     \item{\code{trajectories}}{DRAW, SIM, ID, TIME, CMT, IPRED, DV_SIM,
#'       OBSERVED - the per-observation predictions, as in
#'       \code{\link{ferx_simulate}} (\code{OBSERVED} is \code{NA} here: the
#'       adaptive path has no time-to-event endpoint).}
#'     \item{\code{doses}}{The realized dose ledger: DRAW, SIM, ID, TIME, AMT,
#'       CMT, RATE, DECISION, SIGNAL (the value the controller titrated on -
#'       assay-noised when \code{with_assay_error}), RULE (the rule that fired).}
#'     \item{\code{decisions}}{One row per decision including holds: DRAW, SIM,
#'       ID, DECISION, TIME, SIGNAL, OUTCOME ("dosed" / "hold" / "stop"),
#'       N_DOSED.}
#'     \item{\code{metrics}}{One row per realized run (subject x draw x
#'       replicate) - the per-subject outcome summary: DRAW, SIM, ID, CUM_DOSE
#'       (sum of realized doses), N_DOSES, N_INCREASES / N_DECREASES (realized
#'       dose-change counts, by dose delta), N_HOLDS, DISCONTINUED,
#'       TIME_TO_DISCONT, SIGNAL_MIN / SIGNAL_MAX / SIGNAL_MEAN (over the
#'       decisions that recorded a signal), and PCT_TIME_IN_WINDOW - the fraction
#'       of signal-bearing decisions whose observed value fell inside the model's
#'       \code{[adaptive_dosing]} \code{target_window} (NA when no
#'       \code{target_window} is declared).}
#'   }
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
ferx_simulate_adaptive <- function(model, data, n_sim = 1L, seed = 42L,
                                   verify = TRUE, max_decisions = 0L) {
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
  ferx_rust_simulate_adaptive(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    n_sim = n_sim,
    seed = as.integer(seed),
    verify = if (isTRUE(verify)) "true" else "false",
    max_decisions = as.integer(max_decisions)
  )
}
