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
#'   part. To simulate each subject at its own estimates instead, pass
#'   \code{etas = fit$ebe_etas}.
#' @param etas Optional data frame of each known subject's own random effects:
#'   an \code{ID} column plus one column per eta, named and ordered as the
#'   model declares them (\code{fit$ebe_etas} has exactly this shape, so
#'   \code{ferx_simulate_adaptive(model, fit = fit, etas = fit$ebe_etas)}
#'   simulates every fitted subject at its EBEs). A subject with a row is
#'   simulated at those etas in every replicate instead of a draw from omega;
#'   its kappas (IOV, drawn per decision window) and the controller's assay
#'   error are still drawn, so without either its replicates are identical.
#'   Rows are matched to the data's subjects by \code{ID}; a row whose
#'   subject is not in the data is ignored. A row gives every eta of the
#'   model - there is no partial row, since under a block omega the missing
#'   eta would need a conditional draw. Each supplied subject still consumes
#'   its eta draw, so the drawn subjects get the etas the same \code{seed}
#'   gives without \code{etas}. Works with and without \code{fit}. Refused,
#'   coded \code{E_SUPPLIED_ETAS}: eta columns that are not the model's
#'   names in its order, a repeated \code{ID}, a non-finite value, and a
#'   \code{[mixture]} model (each subject's eta was estimated under its own
#'   class, but the run draws each subject's class). A repeated column name
#'   is refused before the engine is called. Supplied etas are point
#'   estimates treated as exact, which is optimistic for a subject with few
#'   observations (high shrinkage).
#' @param unmatched What to do with a data subject that has no row in
#'   \code{etas}: \code{"refuse"} (the default) errors with
#'   \code{E_SUPPLIED_ETAS}, listing the subjects; \code{"draw"} draws their
#'   etas from omega, so known patients and virtual ones run in one
#'   population. Under \code{"draw"}, an \code{etas} whose IDs match no
#'   data subject (\code{"01"} against \code{"1"}, say) is warned about,
#'   since every subject was then drawn. IDs are matched as text; a numeric
#'   \code{ID} is spelled as a data.frame \code{data} writes it. Passing
#'   \code{unmatched} without \code{etas} has no effect and warns.
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
#'       \code{target_window} is declared), and ETA_SUPPLIED (\code{TRUE}
#'       when the run used the subject's row of \code{etas}, \code{FALSE}
#'       when its eta was drawn).}
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
#'
#' # Subjects 1 and 2 at known etas; the other subjects draw theirs.
#' known <- data.frame(ID = c("1", "2"), ETA_CL = c(0.3, -0.2))
#' res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 20L, seed = 1L,
#'                               etas = known, unmatched = "draw")
#' table(res$metrics$ID, res$metrics$ETA_SUPPLIED)
#'
#' # With a fit to observed data, each fitted subject at its own EBEs:
#' # ferx_simulate_adaptive(model, data, fit = fit, etas = fit$ebe_etas)
#' }
#'
#' @family simulation
#' @export
ferx_simulate_adaptive <- function(model, data = NULL, n_sim = 1L, seed = 42L,
                                   verify = TRUE, max_decisions = 0L,
                                   fit = NULL, etas = NULL,
                                   unmatched = c("refuse", "draw")) {
  if (!missing(unmatched) && is.null(etas)) {
    warning("ferx_simulate_adaptive: `unmatched` has no effect without ",
            "`etas`; every subject's eta is drawn.", call. = FALSE)
  }
  unmatched <- match.arg(unmatched)
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
  supplied <- .ferx_adaptive_etas(etas)
  # The engine's unmatched-subject refusal names `UnmatchedSubjects::Draw`,
  # the Rust spelling of `unmatched = "draw"`.
  res <- .ferx_engine_spelling(if (is.null(fit)) {
    .ferx_engine_call(
      ferx_rust_simulate_adaptive(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        n_sim = n_sim,
        seed = as.integer(seed),
        verify = verify_arg,
        max_decisions = as.integer(max_decisions),
        supply_etas = !is.null(etas),
        eta_ids = supplied$ids,
        eta_names = supplied$names,
        eta_flat = supplied$flat,
        unmatched = unmatched
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
        max_decisions = as.integer(max_decisions),
        supply_etas = !is.null(etas),
        eta_ids = supplied$ids,
        eta_names = supplied$names,
        eta_flat = supplied$flat,
        unmatched = unmatched
      ),
      model, data
    )
  }, "E_SUPPLIED_ETAS", "`UnmatchedSubjects::Draw`", "`unmatched = \"draw\"`")

  # Under "draw" an `etas` keyed differently from the data (`01` against `1`)
  # matches nobody and every subject is drawn - the mistake "refuse" exists to
  # catch. Say so rather than return an all-drawn run as if etas were used.
  if (length(supplied$ids) > 0L && nrow(res$metrics) > 0L &&
      !any(res$metrics$ETA_SUPPLIED)) {
    some <- function(x) {
      x <- unique(x)
      paste0(paste(utils::head(x, 3L), collapse = ", "),
             if (length(x) > 3L) ", ..." else "")
    }
    warning(
      "ferx_simulate_adaptive: no row of `etas` matched a subject in `data`, ",
      "so every subject's eta was drawn. IDs are matched as text: `etas` has ",
      some(supplied$ids), "; `data` has ", some(res$metrics$ID), ".",
      call. = FALSE
    )
  }

  # Same `simulation_warnings` channel `ferx_simulate()` uses - here it carries
  # the engine's warning bundle for the run (ferx-r #426).
  .ferx_surface_sim_warnings(res, "ferx_simulate_adaptive")
}

# `etas =` flattened for the engine (ferx-core #1594): IDs as text, matched to
# the data's verbatim subject IDs; eta columns in the order given, which the
# engine checks against the model's names and order; values row-major. The
# table's own checks (finite, one row per ID) are the engine's too, so every
# refusal of the values carries `E_SUPPLIED_ETAS`.
.ferx_adaptive_etas <- function(etas) {
  if (is.null(etas)) {
    return(list(ids = character(0), names = character(0), flat = numeric(0)))
  }
  if (!is.data.frame(etas) || !"ID" %in% names(etas)) {
    stop("ferx_simulate_adaptive: `etas` must be a data frame with an ID ",
         "column and one column per eta (e.g. `fit$ebe_etas`).", call. = FALSE)
  }
  # A repeated column name would pass the engine's name check on its first
  # copy and drop the second, as a repeated `data` column would.
  dup <- unique(names(etas)[duplicated(names(etas))])
  if (length(dup) > 0L) {
    stop("ferx_simulate_adaptive: `etas` has more than one column named ",
         paste0("`", dup, "`", collapse = ", "), ".", call. = FALSE)
  }
  etas <- as.data.frame(etas)
  ids <- etas$ID
  if (anyNA(ids)) {
    stop("ferx_simulate_adaptive: `etas$ID` has ", sum(is.na(ids)), " NA ",
         "value(s); every row needs its subject's ID.", call. = FALSE)
  }
  # A double ID is spelled as a data.frame `data` writes it
  # (`.ferx_data_column_text()`), so the two match as text.
  ids <- if (is.double(ids)) .fitrx_double_text(ids) else as.character(ids)
  eta_cols <- setdiff(names(etas), "ID")
  not_num <- eta_cols[!vapply(etas[eta_cols], is.numeric, logical(1))]
  if (length(not_num) > 0L) {
    stop("ferx_simulate_adaptive: `etas` column(s) ",
         paste(not_num, collapse = ", "), " are not numeric.", call. = FALSE)
  }
  eta_mat <- as.matrix(etas[, eta_cols, drop = FALSE])
  storage.mode(eta_mat) <- "double"
  list(ids = ids, names = eta_cols, flat = as.numeric(t(eta_mat)))
}
