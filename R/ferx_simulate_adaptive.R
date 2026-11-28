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
#' \code{$trajectories}, \code{$decisions}, \code{PCT_READINGS_IN_WINDOW}, and the
#' \code{auc_target} attainment metric \emph{do} reflect the base regimen. A base
#' dose scheduled exactly at a decision time is read \strong{pre-dose} (the
#' trough), and a base dose scheduled after a \code{stop} still lands - a
#' \code{stop} ends controller augmentation, not the standing prescription.
#'
#' \strong{System resets.} An \code{EVID = 3} row zeroes the compartments (and
#' re-seeds any \code{init(state) = expr}) at its time; \code{EVID = 4} does the
#' same and then doses. Ordering at a shared time is reset, then dose, then
#' observation. A base regimen combined with resets runs on the
#' constant-covariate path and under a time-varying covariate or IOV. One
#' combination remains rejected with a typed error rather than silently
#' mis-integrated: a steady-state, lag-timed, or modeled-rate (\code{RATE})
#' base dose \emph{under} a time-varying covariate or IOV (supported on the
#' constant-covariate path). An \code{auc_target} metric is accepted on a
#' subject with resets, under a time-varying covariate, a \code{TIME}-dependent
#' PK parameter, or IOV: the window AUC is integrated on the run's own
#' trajectory, so a dose or reset inside a decision window is integrated
#' exactly.
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
#'   model declares them, leaving out any \code{omega NAME ~ 0 FIX}
#'   placeholder, which is removed from the model (ferx-core #3150; a table
#'   with its column is refused with \code{E_SUPPLIED_ETAS}).
#'   \code{fit$ebe_etas} has exactly this shape (so
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
#' @param uncertainty Optional parameter uncertainty to propagate from
#'   \code{fit}: \code{NULL} (the default) simulates at the fit's point
#'   estimates; otherwise \code{list(method = , n_draws = )}, where
#'   \code{method} is \code{"asymptotic"} (the default; a multivariate normal
#'   on \code{fit$cov_matrix}, so the fit needs its covariance step) or
#'   \code{"sir"} (resampled from the fit's kept SIR draws, so the fit needs
#'   \code{sir_keep_samples = TRUE}), and \code{n_draws} (a positive whole
#'   number) is the number of parameter sets. Every subject and each of the
#'   \code{n_sim} replicates is simulated under every draw, so a run has
#'   \code{n_draws * n_sim} replicates per subject. The draws are the
#'   sampler \code{\link{ferx_simulate_with_uncertainty}} uses, with a FIX or
#'   held-at-bound parameter kept at its estimate. They come from a stream of
#'   their own: draw 1 uses the random effects the same \code{seed} gives
#'   without uncertainty (under draw 1's omega), so draw 1 does not depend on
#'   \code{n_draws}, and each further draw simulates new virtual patients.
#'   Requires \code{fit}, and cannot be combined with \code{etas}: each
#'   supplied eta was estimated at the fit's theta, so pairing it with a drawn
#'   theta would not describe the same patient (refused, coded
#'   \code{E_SUPPLIED_ETAS}). Dose expressions that read only theta are
#'   evaluated at the fit's estimate in every draw, while a
#'   \code{dose_knowledge = oracle} subject is dosed from its own draw's
#'   parameters. Warnings and refusals name the replicate as
#'   \code{"(draw d, sim s)"} only in a run with uncertainty.
#' @param n_reps Assay replicates per subject run (ferx-core #946), a whole
#'   number >= 1 (default \code{1L}). Each subject's eta (and \code{[mixture]}
#'   class) is drawn once per \code{(DRAW, SIM)} and held fixed across
#'   replicates \code{1..n_reps}, while the IOV kappas and the controller's
#'   assay noise are drawn afresh for each, so the controller re-runs on the
#'   same patient under new within-subject noise. \code{doses},
#'   \code{decisions} and \code{metrics} carry each row's \code{REP};
#'   \code{trajectories} holds replicate 1 only (join it on
#'   \code{REP == 1}). \code{n_reps = 1} is the run without replicates, bit
#'   for bit. \code{decompose = TRUE} needs \code{n_reps >= 2}.
#' @param summary_probs The three probabilities of every band and prediction
#'   interval in the population summary (the \code{LO}, \code{MED} and
#'   \code{HI} columns), strictly inside (0, 1) and strictly increasing
#'   (default \code{c(0.05, 0.5, 0.95)}; otherwise refused before the run, coded
#'   \code{E_ADAPTIVE_OPTIONS}).
#' @param summary_thresholds Optional thresholds for the population summary:
#'   each adds a \code{"fraction_at_least"} statistic, the fraction of runs
#'   whose value is at least the threshold. Either a named list (or named
#'   numeric vector) from a \code{metrics} column name to one or more
#'   thresholds, e.g. \code{list(PCT_READINGS_IN_WINDOW = 0.8)}, or a data
#'   frame with \code{METRIC} and \code{THRESHOLD} columns, plus
#'   \code{ANALYTE} and \code{DECISION} for a \code{"SIGNAL"} threshold on one
#'   monitor's latent signal at one decision (\code{DECISION} numbered as in
#'   \code{decisions$DECISION}). The metrics are \code{CUM_DOSE},
#'   \code{N_DOSES}, \code{N_INCREASES}, \code{N_DECREASES}, \code{N_HOLDS},
#'   \code{N_NOT_DUE}, \code{DISCONTINUED} (as 1/0), \code{SIGNAL_MIN},
#'   \code{SIGNAL_MAX}, \code{SIGNAL_MEAN}, \code{PCT_READINGS_IN_WINDOW},
#'   \code{AUC_TARGET_ATTAINMENT} and \code{SIGNAL}. A non-finite threshold is
#'   refused before the run, coded \code{E_ADAPTIVE_OPTIONS}. A
#'   \code{"SIGNAL"} threshold whose analyte and decision no run recorded (a
#'   mistyped \code{ANALYTE}, a \code{DECISION} past the schedule) gives
#'   \code{NA} rows and a warning naming what the runs did record.
#' @param decompose If \code{TRUE}, also split the variance of every summary
#'   metric into its parameter-uncertainty, between-subject and within-subject
#'   layers (\code{variance_components}, \code{variance_skipped}). The design
#'   has to separate the layers: \code{n_reps >= 2}, \code{n_sim >= 2}, and
#'   with \code{uncertainty} at least two draws; otherwise the call is refused
#'   before the simulation runs, coded \code{E_ADAPTIVE_OPTIONS}. Default
#'   \code{FALSE}.
#' @param reference \code{NULL} (default) for one run under the controller,
#'   or \code{"fixed"} for the controller's effect as a contrast (ferx-core
#'   #946): the same model, data, seed and options are run twice, once under
#'   the \code{[adaptive_dosing]} rules (the adaptive arm) and once under the
#'   fixed reference regimen, the block with its rules cleared, which gives
#'   each subject the block's start dose at every due point (its own value of
#'   a \code{start_dose} expression, on the \code{round_to} grid, within
#'   \code{dose_bounds}, on \code{start_interval}). The two arms share every
#'   random number - parameter draws, etas, mixture classes, IOV kappas and
#'   assay noise - so their difference is the controller's effect, not noise.
#'   It is not a fair reference for a controller with a loading dose: say so
#'   in \code{arm_labels}.
#' @param arm_labels With \code{reference}: what each arm is, as
#'   \code{c(adaptive = ..., reference = ...)}. Defaults to
#'   \code{c(adaptive = "[adaptive_dosing] rules",
#'   reference = "start dose every due point")}. An empty label is refused,
#'   coded \code{E_ADAPTIVE_OPTIONS}: the reference is never implicit.
#'
#' @return A list of data frames (five run tables, three summary tables and
#'   two decomposition tables):
#'   \describe{
#'     \item{\code{trajectories}}{DRAW, SIM, ID, TIME, CMT, IPRED, DV_SIM,
#'       OBSERVED - the per-observation predictions, as in
#'       \code{\link{ferx_simulate}} (\code{OBSERVED} is \code{NA} here: the
#'       adaptive path has no time-to-event endpoint). A \code{[mixture]}
#'       model adds an integer \code{MIXNUM} column, the class each subject
#'       was drawn into for that replicate; it is absent otherwise.}
#'     \item{\code{doses}}{The realized dose ledger, \strong{controller-issued
#'       doses only} (pre-scheduled base doses are excluded): DRAW, SIM, REP, ID, TIME,
#'       AMT, CMT, RATE, DECISION, SIGNAL (the value the controller titrated on -
#'       assay-noised when \code{with_assay_error}), LATENT (the noise-free
#'       monitor value that reading was drawn from; equal to SIGNAL when the
#'       monitor has no assay error), RULE (the rule that fired).}
#'     \item{\code{decisions}}{One row per review-grid point including holds:
#'       DRAW, SIM, REP, ID, DECISION, TIME, SIGNAL, LATENT, OUTCOME ("dosed" / "hold" /
#'       "not_due" / "stop"), N_DOSED, INTERVAL. "not_due" marks a grid point
#'       at which no dose was due under the dosing interval in force (the
#'       block declares \code{intervals} and \code{start_interval}): the
#'       controller was not consulted and gave no dose, but SIGNAL still holds
#'       what the monitor read. INTERVAL is the dosing interval in force after
#'       the decision, NA when the block declares no \code{intervals}.}
#'     \item{\code{metrics}}{One row per realized run (subject x draw x
#'       replicate) - the per-subject outcome summary: DRAW, SIM, REP, ID, CUM_DOSE
#'       (sum of controller-issued doses; excludes the base regimen), N_DOSES
#'       (controller doses only), N_INCREASES / N_DECREASES (realized
#'       dose-change counts, by dose delta), N_HOLDS (consulted decisions that
#'       gave no dose; not-due points are not counted), N_NOT_DUE (grid points
#'       that were not due), DISCONTINUED,
#'       TIME_TO_DISCONT, SIGNAL_MIN / SIGNAL_MAX / SIGNAL_MEAN (of LATENT over
#'       the decisions that recorded a signal), PCT_READINGS_IN_WINDOW - the
#'       fraction of signal-bearing decision readings whose LATENT value fell
#'       inside the model's \code{[adaptive_dosing]} \code{target_window} (NA
#'       when no \code{target_window} is declared; it counts readings, not
#'       elapsed time), PCT_TIME_IN_WINDOW (deprecated: the same values under
#'       the old name, to be removed in a later release),
#'       AUC_TARGET_ATTAINMENT - the
#'       fraction of windows between consecutive realized decisions whose area
#'       under the monitored signal fell inside the \code{auc_target} band
#'       (NA when no \code{auc_target} is declared, or with fewer than two
#'       decisions), ETA_SUPPLIED (\code{TRUE} when the run used the
#'       subject's row of \code{etas}, \code{FALSE} when its eta was drawn),
#'       and DOSE_ORACLE (\code{TRUE} for a subject whose eta was drawn and
#'       whose dose expressions read its drawn individual parameters, under
#'       the block's \code{dose_knowledge = oracle}: a perfect-information
#'       comparison arm, never a dosing rule).}
#'     In \code{doses}, \code{decisions} and \code{metrics}, \code{REP} is the
#'       1-based assay replicate within a \code{(DRAW, SIM)} (ferx-core #946):
#'       \code{1} on every row of a run without replicates.
#'     \item{\code{parameters}}{The parameter set each draw was simulated
#'       under: one row per draw, \code{DRAW} and then one column per
#'       estimated quantity - theta, the omega elements, sigma, the IOV
#'       (kappa) omega elements, a \code{[mixture]} model's per-class
#'       overrides and the \code{block_sigma} correlations - named as the
#'       rows of \code{fit$cov_matrix} (\code{TVCL}, \code{ETA_CL},
#'       \code{ETA_V ~ ETA_CL}, \code{EPS_PROP}, \code{ETA_CL (class 2)},
#'       ...) and on the scale \code{fit$theta}, \code{fit$omega} and
#'       \code{fit$sigma} report them. Join it to the other tables by
#'       \code{DRAW}; a parameter itself named \code{DRAW} is renamed
#'       \code{DRAW.1}, as \code{make.unique()} would. Without \code{uncertainty} it has these columns and
#'       no rows.}
#'     \item{\code{summary_rows}}{The population summary (ferx-core #946),
#'       one row per metric, statistic and draw: DRAW, METRIC (a
#'       \code{metrics} column name, or \code{"SIGNAL"} for one monitor's
#'       latent signal at one decision), ANALYTE and DECISION (the monitor's
#'       name and the decision index, as in \code{decisions$DECISION}, for
#'       \code{"SIGNAL"}; \code{NA} otherwise), STATISTIC (\code{"mean"},
#'       \code{"median"} or \code{"fraction_at_least"}), THRESHOLD (the
#'       \code{summary_thresholds} value of a \code{"fraction_at_least"} row,
#'       \code{NA} otherwise), VALUE (over the draw's runs that carry the
#'       metric; \code{NA} when none does), MC_SE (the Monte Carlo standard
#'       error of VALUE, ferx-core #3215: the delete-one jackknife over the
#'       draw's sims that carry the metric, a sim being the independent
#'       replicate since its runs share each subject's eta; read it with a t
#'       on S - 1 degrees of freedom for S such sims. \code{NA} on a
#'       \code{"median"} row, where the jackknife is inconsistent, and with
#'       fewer than two such sims), N_RUNS (runs in the draw, every
#'       subject, SIM and REP counted equally) and N_MISSING (runs carrying no
#'       value, which are skipped, never imputed). \code{DISCONTINUED} is
#'       summarised as 1/0, so its mean is the fraction discontinued;
#'       \code{TIME_TO_DISCONT} is censored and not summarised.}
#'     \item{\code{summary_bands}}{The parameter-uncertainty band of each
#'       per-draw statistic, its quantiles across draws: METRIC, ANALYTE,
#'       DECISION, STATISTIC, THRESHOLD, N_DRAWS_PRESENT, N_DRAWS_MISSING,
#'       STATUS (\code{"not_drawn"} for a run with one draw, which has no
#'       band and is never reported as a zero-width one, else
#'       \code{"quantiles"}), LO, MED, HI at \code{summary_probs} (type 7;
#'       \code{NA} until there are enough values to resolve a quantile, e.g.
#'       19 for 5\% / 95\%). The probabilities are its \code{probs}
#'       attribute.}
#'     \item{\code{prediction_intervals}}{Each metric's individual values
#'       pooled over every run of every draw - where one simulated patient
#'       falls, parameter uncertainty included: METRIC, ANALYTE, DECISION, LO,
#'       MED, HI, N_RUNS, N_MISSING, with the same quantile rule and
#'       \code{probs} attribute as \code{summary_bands}.}
#'     \item{\code{variance_components}}{With \code{decompose = TRUE}
#'       (\code{NULL} otherwise): the nested variance decomposition of each
#'       metric, long form. METRIC, ANALYTE, DECISION, SCOPE
#'       (\code{"subject"} or \code{"population"}), ID (the subject;
#'       \code{NA} on a population row), LAYER (\code{"parameter"},
#'       \code{"bsv"}, \code{"within"} - IOV and assay noise -,
#'       \code{"between_subject"} on population rows, and \code{"total"}),
#'       STATUS (\code{"estimated"}, \code{"not_drawn"} for the parameter
#'       layer of a run without \code{uncertainty}, \code{"supplied"} for the
#'       BSV of a subject run at its own \code{etas}), VALUE (as estimated: a
#'       negative estimate is not truncated), MC_SE (its jackknife Monte Carlo
#'       standard error; \code{NA} when the design cannot form one) and
#'       N_SUBJECTS. A subject's total is the sum of its estimated layers; the
#'       population total is the mean subject total plus the between-subject
#'       line, which is the sum of the population lines when no subject's eta
#'       was supplied. Attributes \code{n_draws}, \code{n_sims} and
#'       \code{n_reps} give the design.}
#'     \item{\code{variance_skipped}}{With \code{decompose = TRUE}: the
#'       (metric, subject) pairs left out because the metric is missing in
#'       some of the subject's runs: METRIC, ANALYTE, DECISION, ID,
#'       N_MISSING.}
#'   }
#'
#'   \code{DRAW} on every table is the parameter draw a row was simulated
#'   under: \code{1} throughout without \code{uncertainty}, \code{1} to
#'   \code{n_draws} with it.
#'
#'   With \code{reference}, every table gains a first column ARM:
#'   \code{"adaptive"} and \code{"reference"} rows on all of them (each
#'   arm's rows are exactly what that arm alone gives), and on the summary
#'   and decomposition tables also \code{"difference"} rows, adaptive minus
#'   reference. The difference is paired: \code{summary_rows} holds the
#'   summary of the per-run differences, so its \code{"mean"} row is the mean
#'   effect with a paired MC_SE (one sim left out of both arms at once), its
#'   band the parameter uncertainty of that effect and its prediction
#'   interval the spread of one patient's effect; a \code{"fraction_at_least"}
#'   difference row is P(A - R >= threshold), while the change in attainment
#'   is the adaptive row minus the reference row. On
#'   \code{variance_components} a negative \code{"bsv"} difference at
#'   population scope is between-patient variance the controller removed. A
#'   run (or subject) whose metric is missing in either arm is missing in the
#'   difference. The labels are the list's \code{arm_labels} element.
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
#' # Three assay replicates per patient, the fraction of runs with at least
#' # 80% of readings in the target window, and the variance decomposition.
#' res <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 20L, seed = 1L,
#'                               n_reps = 3L, decompose = TRUE,
#'                               summary_thresholds = list(PCT_READINGS_IN_WINDOW = 0.8))
#' subset(res$summary_rows, METRIC == "PCT_READINGS_IN_WINDOW")
#' subset(res$variance_components, SCOPE == "population" & METRIC == "CUM_DOSE")
#'
#' # The controller's effect: the same patients under the fixed start dose.
#' con <- ferx_simulate_adaptive(ex$model, ex$data, n_sim = 50L, seed = 1L,
#'                               reference = "fixed")
#' subset(con$summary_rows, METRIC == "PCT_READINGS_IN_WINDOW" &
#'                          STATISTIC == "mean")
#'
#' # With a fit to observed data, each fitted subject at its own EBEs:
#' # ferx_simulate_adaptive(model, data, fit = fit, etas = fit$ebe_etas)
#'
#' # Or propagate the fit's parameter uncertainty: 50 parameter draws.
#' # res <- ferx_simulate_adaptive(model, data, fit = fit, seed = 1L,
#' #                               uncertainty = list(n_draws = 50L))
#' # merge(res$metrics, res$parameters, by = "DRAW")
#' }
#'
#' @family simulation
#' @export
ferx_simulate_adaptive <- function(model, data = NULL, n_sim = 1L, seed = 42L,
                                   verify = TRUE, max_decisions = 0L,
                                   fit = NULL, etas = NULL,
                                   unmatched = c("refuse", "draw"),
                                   uncertainty = NULL, n_reps = 1L,
                                   summary_probs = c(0.05, 0.5, 0.95),
                                   summary_thresholds = NULL,
                                   decompose = FALSE, reference = NULL,
                                   arm_labels = NULL) {
  uncertainty <- .ferx_adaptive_uncertainty(uncertainty, fit)
  arm_labels <- .ferx_adaptive_arm_labels(reference, arm_labels)
  run_options <- .ferx_adaptive_run_options(n_reps, summary_probs,
                                            summary_thresholds, decompose,
                                            arm_labels)
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
  # Engine refusals name Rust spellings; here they are this function's
  # arguments. E_SUPPLIED_ETAS: the refusal of supplied etas under parameter
  # uncertainty (ferx-core #946) names the two `AdaptiveInputs` builders, and
  # the unmatched-subject refusal names `UnmatchedSubjects::Draw`.
  # E_ADAPTIVE_OPTIONS: the decomposition's refusal of a design that cannot
  # separate a layer names the Rust options (the glue now refuses those
  # designs before the run, but core's check stays behind it).
  supplied_etas_spelling <- function(call) {
    .ferx_engine_spelling(
      call, "E_SUPPLIED_ETAS",
      c("`with_uncertainty`", "`with_etas`", "`UnmatchedSubjects::Draw`"),
      c("`uncertainty`", "`etas`", "`unmatched = \"draw\"`")
    )
  }
  options_spelling <- function(call) {
    .ferx_engine_spelling(
      call, "E_ADAPTIVE_OPTIONS",
      c("`AdaptiveSimulateOptions::n_reps`", "in `UncertaintySpec`"),
      c("`n_reps`", "in `uncertainty`")
    )
  }
  res <- options_spelling(supplied_etas_spelling(if (is.null(fit)) {
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
        unmatched = unmatched,
        run_options = run_options
      ),
      model, data
    )
  } else {
    # The fitted estimates, placed on the fit's theta layout exactly as
    # `ferx_simulate(fit = )` places them (ferx-core #1594). No occasion
    # rule is passed: the adaptive engine draws kappa per decision window.
    fit_pieces <- validate_fit_for_params(fit)
    # The fit's covariance matrix or SIR pool, refused here with the R remedy
    # when the fit lacks the one `method` needs, as `ferx_simulate_with_
    # uncertainty()` refuses it.
    unc_pieces <- if (is.null(uncertainty)) {
      validate_fit_for_uncertainty(NULL, "none")
    } else {
      # Re-raised under this function's name: the helper's own call is
      # internal and means nothing to the caller.
      tryCatch(
        validate_fit_for_uncertainty(fit, uncertainty$method),
        error = function(e) {
          stop("ferx_simulate_adaptive: ", conditionMessage(e), call. = FALSE)
        }
      )
    }
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
        unmatched = unmatched,
        uncertainty_method = if (is.null(uncertainty)) "" else uncertainty$method,
        cov_matrix_flat = unc_pieces$cov_matrix_flat,
        cov_matrix_dim = unc_pieces$cov_matrix_dim,
        sir_resamples_flat = unc_pieces$sir_resamples_flat,
        sir_resamples_n = unc_pieces$sir_resamples_n,
        sir_resamples_dim = unc_pieces$sir_resamples_dim,
        n_draws = if (is.null(uncertainty)) 0L else uncertainty$n_draws,
        run_options = run_options
      ),
      model, data
    )
  }))

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

  if (length(arm_labels) > 0L) res <- .ferx_adaptive_stack_arms(res, arm_labels)

  .ferx_adaptive_warn_unrecorded_signals(res, run_options)

  # Same `simulation_warnings` channel `ferx_simulate()` uses - here it carries
  # the engine's warning bundle for the run (ferx-r #426).
  .ferx_surface_sim_warnings(res, "ferx_simulate_adaptive")
}

# `n_reps`, `summary_probs`, `summary_thresholds` and `decompose` checked for
# shape and flattened into the `run_options` list the glue reads (ferx-core
# #946 PR2-4). Only what the glue cannot pass on is refused here: `n_reps` that
# is not a whole number >= 1 (the glue casts it to an unsigned count), a
# `summary_probs` that is not three numbers, a malformed threshold table. The
# values themselves (probabilities strictly inside (0, 1) and increasing, a
# finite threshold, a known metric name) are the engine's and the glue's to
# refuse, so the rules live in one place.
.ferx_adaptive_run_options <- function(n_reps = 1L,
                                       summary_probs = c(0.05, 0.5, 0.95),
                                       summary_thresholds = NULL,
                                       decompose = FALSE,
                                       arm_labels = character(0)) {
  if (!is.numeric(n_reps) || length(n_reps) != 1L || !is.finite(n_reps) ||
      n_reps != round(n_reps) || n_reps < 1 || n_reps > .Machine$integer.max) {
    stop("ferx_simulate_adaptive: `n_reps` must be a single whole number ",
         ">= 1 (assay replicates per subject run; 1 is a run without ",
         "replicates).", call. = FALSE)
  }
  if (!is.numeric(summary_probs) || length(summary_probs) != 3L) {
    stop("ferx_simulate_adaptive: `summary_probs` must be three ",
         "probabilities, e.g. c(0.05, 0.5, 0.95).", call. = FALSE)
  }
  if (length(decompose) != 1L || is.na(decompose) || !is.logical(decompose)) {
    stop("ferx_simulate_adaptive: `decompose` must be a single TRUE or FALSE.",
         call. = FALSE)
  }
  thr <- .ferx_adaptive_thresholds(summary_thresholds)
  list(
    n_reps = as.integer(n_reps),
    probs = as.double(summary_probs),
    threshold_metric = thr$metric,
    threshold_value = thr$value,
    threshold_analyte = thr$analyte,
    threshold_decision = thr$decision,
    decompose = decompose,
    arm_labels = unname(arm_labels)
  )
}

# `reference` and `arm_labels` as the glue's `arm_labels`: `character(0)`
# for one run, `c(adaptive, reference)` for a contrast. Only the shape is
# checked here; an empty label is core's to refuse (E_ADAPTIVE_OPTIONS).
.ferx_adaptive_arm_labels <- function(reference, arm_labels) {
  if (is.null(reference)) {
    if (!is.null(arm_labels)) {
      warning("ferx_simulate_adaptive: `arm_labels` has no effect without ",
              "`reference`.", call. = FALSE)
    }
    return(character(0))
  }
  if (!identical(reference, "fixed")) {
    stop("ferx_simulate_adaptive: `reference` must be NULL or \"fixed\" ",
         "(the start dose at every due point).", call. = FALSE)
  }
  if (is.null(arm_labels)) {
    return(c(adaptive = "[adaptive_dosing] rules",
             reference = "start dose every due point"))
  }
  if (!is.character(arm_labels) || length(arm_labels) != 2L ||
      !setequal(names(arm_labels), c("adaptive", "reference")) ||
      anyNA(arm_labels)) {
    stop("ferx_simulate_adaptive: `arm_labels` must be a character vector ",
         "c(adaptive = ..., reference = ...).", call. = FALSE)
  }
  arm_labels[c("adaptive", "reference")]
}

# The glue's contrast list (`adaptive`, `reference` - each one arm's full
# result - and `difference`, the difference's outcome tables) as one result
# shaped like a single run's: every table stacked by arm under a leading ARM
# column. Attributes (`probs`, the design of `variance_components`) are the
# same in every arm and kept from the first.
.ferx_adaptive_stack_arms <- function(res, arm_labels) {
  arms <- list(adaptive = res$adaptive, reference = res$reference,
               difference = res$difference)
  out <- lapply(stats::setNames(nm = names(res$adaptive)), function(tab) {
    parts <- Filter(Negate(is.null), lapply(names(arms), function(a) {
      df <- arms[[a]][[tab]]
      if (is.null(df)) return(NULL)
      cbind(data.frame(ARM = rep(a, nrow(df))), df)
    }))
    if (length(parts) == 0L) return(NULL)
    stacked <- do.call(rbind, parts)
    rownames(stacked) <- NULL
    first <- res$adaptive[[tab]]
    for (a in setdiff(names(attributes(first)), c("names", "row.names", "class"))) {
      attr(stacked, a) <- attr(first, a)
    }
    stacked
  })
  out$arm_labels <- arm_labels
  warnings <- unique(c(attr(res$adaptive, "simulation_warnings", exact = TRUE),
                       attr(res$reference, "simulation_warnings", exact = TRUE)))
  attr(out, "simulation_warnings") <- warnings
  out
}

# A `SIGNAL` threshold whose analyte and decision no run recorded - a
# mistyped ANALYTE, a DECISION past the schedule - gives summary rows that are
# NA throughout, as core reports a metric no run carries. Not an error: a late
# decision every run stopped before is a real, empty answer. But a typo should
# not pass as one, so warn, naming what was recorded (#779 review, finding 2).
.ferx_adaptive_warn_unrecorded_signals <- function(res, run_options) {
  is_sig <- run_options$threshold_metric == "SIGNAL"
  if (!any(is_sig)) return(invisible())
  rows <- res$summary_rows
  sig <- rows[rows$METRIC == "SIGNAL" & rows$STATISTIC == "mean", , drop = FALSE]
  recorded <- sig[sig$N_MISSING < sig$N_RUNS, , drop = FALSE]
  key <- function(a, d) paste0("ANALYTE \"", a, "\" at DECISION ", d)
  asked <- unique(key(run_options$threshold_analyte[is_sig],
                      run_options$threshold_decision[is_sig]))
  missing <- setdiff(asked, unique(key(recorded$ANALYTE, recorded$DECISION)))
  if (length(missing) == 0L) return(invisible())
  seen <- if (nrow(recorded) == 0L) {
    "no signal at all"
  } else {
    ranges <- vapply(split(recorded$DECISION, recorded$ANALYTE), function(d) {
      paste0(range(d), collapse = "-")
    }, character(1))
    paste0("ANALYTE \"", names(ranges), "\" at DECISION ", ranges,
           collapse = ", ")
  }
  warning(
    "ferx_simulate_adaptive: no run recorded a signal for the ",
    "`summary_thresholds` row(s) ", paste(missing, collapse = "; "),
    ", so their summary rows are NA. The runs recorded ", seen, ".",
    call. = FALSE
  )
  invisible()
}

# `summary_thresholds =` as four parallel vectors. A named list or named
# numeric vector maps a `metrics` column name to one or more thresholds
# (`list(PCT_READINGS_IN_WINDOW = 0.8, CUM_DOSE = c(5000, 8000))`); a data
# frame has METRIC and THRESHOLD columns, plus ANALYTE and DECISION for a
# `SIGNAL` threshold (one analyte's latent signal at one decision).
.ferx_adaptive_thresholds <- function(x) {
  empty <- list(metric = character(0), value = numeric(0),
                analyte = character(0), decision = integer(0))
  if (is.null(x) || length(x) == 0L) return(empty)
  bad <- function() {
    stop("ferx_simulate_adaptive: `summary_thresholds` must be NULL, a ",
         "named list (or named numeric vector) of thresholds per metrics ",
         "column, e.g. `list(PCT_READINGS_IN_WINDOW = 0.8)`, or a data frame ",
         "with METRIC and THRESHOLD columns (plus ANALYTE and DECISION for a ",
         "SIGNAL threshold).", call. = FALSE)
  }
  if (is.data.frame(x)) {
    if (!all(c("METRIC", "THRESHOLD") %in% names(x)) ||
        !is.numeric(x$THRESHOLD)) bad()
    n <- nrow(x)
    col_or_na <- function(name, na) {
      if (name %in% names(x)) x[[name]] else rep(na, n)
    }
    decision <- col_or_na("DECISION", NA_integer_)
    if (!is.numeric(decision) && !all(is.na(decision))) bad()
    if (any(!is.na(decision) & decision != round(decision))) bad()
    # A table read with `read.csv()` has "" where a scalar row has no
    # analyte: the same as NA (#779 review, finding 7).
    analyte <- as.character(col_or_na("ANALYTE", NA_character_))
    analyte[!is.na(analyte) & analyte == ""] <- NA_character_
    return(list(
      metric = as.character(x$METRIC),
      value = as.double(x$THRESHOLD),
      analyte = analyte,
      decision = as.integer(decision)
    ))
  }
  nms <- names(x)
  if (is.null(nms) || any(is.na(nms) | nms == "")) bad()
  if (is.list(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) bad()
  } else if (!is.numeric(x)) {
    bad()
  }
  lens <- lengths(x)
  list(
    metric = rep(nms, lens),
    value = as.double(unlist(x, use.names = FALSE)),
    analyte = rep(NA_character_, sum(lens)),
    decision = rep(NA_integer_, sum(lens))
  )
}

# `uncertainty =` checked and normalised to list(method, n_draws), or NULL
# (ferx-core #946). Refused here: an uncertainty without `fit` (there is no
# estimate to draw around), a malformed list, and `n_draws` that is not a
# positive whole number (the glue casts it to an unsigned count). Combined with
# `etas` it reaches the engine, which refuses it coded `E_SUPPLIED_ETAS`.
.ferx_adaptive_uncertainty <- function(uncertainty, fit) {
  if (is.null(uncertainty)) return(NULL)
  if (is.null(fit)) {
    stop("ferx_simulate_adaptive: `uncertainty` needs `fit`: the parameter ",
         "draws come from a fit's covariance matrix or SIR samples. Pass ",
         "`fit = `, or drop `uncertainty`.", call. = FALSE)
  }
  if (!is.list(uncertainty) || is.null(names(uncertainty)) ||
      any(names(uncertainty) == "")) {
    stop("ferx_simulate_adaptive: `uncertainty` must be NULL or a named ",
         "list, e.g. `list(method = \"asymptotic\", n_draws = 100)`.",
         call. = FALSE)
  }
  unknown <- setdiff(names(uncertainty), c("method", "n_draws"))
  if (length(unknown) > 0L) {
    stop("ferx_simulate_adaptive: unknown `uncertainty` element(s) ",
         paste0("`", unknown, "`", collapse = ", "),
         "; expected `method` and `n_draws`.", call. = FALSE)
  }
  method <- uncertainty$method
  if (is.null(method)) method <- "asymptotic"
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("asymptotic", "sir")) {
    stop("ferx_simulate_adaptive: `uncertainty$method` must be ",
         "\"asymptotic\" or \"sir\".", call. = FALSE)
  }
  n_draws <- uncertainty$n_draws
  if (is.null(n_draws)) {
    stop("ferx_simulate_adaptive: `uncertainty$n_draws` is required: the ",
         "number of parameter draws, e.g. `list(n_draws = 100)`.",
         call. = FALSE)
  }
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws != round(n_draws) || n_draws < 1 ||
      n_draws > .Machine$integer.max) {
    stop("ferx_simulate_adaptive: `uncertainty$n_draws` must be a single ",
         "positive whole number (the number of parameter draws).",
         call. = FALSE)
  }
  list(method = method, n_draws = as.integer(n_draws))
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
