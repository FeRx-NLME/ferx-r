#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend.
#'
#' @param model Path to a \code{.ferx} model file, or a \code{\link{ferx_model}} object.
#' @param data Path to a NONMEM-format CSV file. When omitted, it defaults to
#'   the dataset declared in the model file's \code{[data]} block
#'   (\code{path = ...}); passing \code{data} here overrides that. Required
#'   columns: ID, TIME,
#'   DV, EVID, AMT, CMT. Optional columns recognised by the engine: RATE
#'   (infusion rate; `RATE = -1` infuses `AMT` at a modeled rate given by a
#'   per-subject parameter `R{n}`, and `RATE = -2` over a modeled duration given
#'   by `D{n}`, on dose compartment `n`),
#'   MDV (missing-DV flag), II (dosing interval, required when
#'   SS > 0), SS (steady-state flag: 1 = pre-dose at steady state, 2 = add SS
#'   concentration to current state), CENS (LOQ censoring flag for M3 method:
#'   \code{1} below LLOQ, \code{-1} above ULOQ), OCC (occasion index for IOV),
#'   and any covariate columns referenced in the model. See the steady-state
#'   section below for SS/II details.
#' @param method Estimation method(s). \code{NULL} (the default) uses whatever
#'   the model file's \code{[fit_options] method} specifies, falling back to
#'   \code{"focei"} if the model file sets none (the engine warns when it
#'   defaults). Pass a value here to override the model file: either a single
#'   string or a character
#'   vector of methods to run in sequence. Each stage is seeded with the
#'   previous stage's converged parameters, and only the final stage produces
#'   the reported covariance/diagnostics. Supported methods: \code{"foce"},
#'   \code{"focei"},
#'   \code{"laplace"} (alias \code{"laplacian"}; the Laplace approximation with the
#'   \emph{exact} Hessian - NONMEM \code{$EST METHOD=1 LAPLACIAN INTER}. This is not the
#'   same estimator as \code{"focei"}, which builds its Gaussian from the
#'   Gauss-Newton Hessian and reports a different OFV. At the default
#'   \code{n_agq = 1} it is a single node; \code{settings = list(n_agq = N)} with
#'   \code{N > 1} turns it into adaptive Gauss-Hermite quadrature - the exact
#'   conditional likelihood evaluated on a Gauss-Hermite grid around each subject's
#'   empirical-Bayes mode, which makes no Gaussian-residual assumption and so
#'   handles non-Gaussian endpoints (time-to-event, categorical). Cost is
#'   \code{n_agq^n_eta} per subject per iteration, so higher node counts suit
#'   models with few random effects; \code{n_agq} is capped at 21, and a fit whose
#'   tensor grid would exceed 100000 nodes is rejected at check time. Odd values
#'   are conventional. Supports IOV at any occasion count. There is
#'   no separate \code{"agq"} method - adaptive quadrature is this method with
#'   \code{n_agq > 1}; \code{"focei"} likewise accepts \code{n_agq > 1} for the
#'   Gauss-Newton-anchored quadrature),
#'   \code{"saem"}, \code{"gn"} (Gauss-Newton / BHHH),
#'   \code{"gn_hybrid"} (Gauss-Newton followed by a FOCEI polish step),
#'   \code{"imp"} (also accepted as \code{"importance_sampling"} or
#'   \code{"importance-sampling"}; the NONMEM \code{METHOD=IMP}
#'   importance-sampling Monte-Carlo EM estimator), or \code{"impmap"} (also
#'   accepted as \code{"importance_sampling_map"}; Importance Sampling assisted
#'   by Mode A Posteriori, the NONMEM \code{METHOD=IMPMAP} Monte-Carlo EM
#'   estimator).
#'   \code{"imp"} is an estimator by default (it updates parameters) and may run
#'   standalone (\code{method = "imp"}), lead, or sit mid-chain. Set
#'   \code{settings = list(imp_eval_only = TRUE)} (NONMEM \code{EONLY=1}) to make
#'   it instead *evaluate* the marginal \code{-2 log L} at the fixed input
#'   parameters; in that mode it must be the *last* entry of the chain, e.g.
#'   \code{c("focei", "imp")}. Plain \code{"imp"} re-centers its proposal from
#'   the previous iteration's samples and so is fragile on rich data; prefer
#'   \code{"impmap"} or warm-start with \code{c("focei", "imp")} there.
#'   \code{"impmap"} is a full estimator: it may run standalone
#'   (\code{method = "impmap"}) or as a chain stage (\code{c("focei", "impmap")}),
#'   requires a mu-referenced parameterization, and does not yet support IOV.
#'   Example chain: \code{c("saem", "focei")}.
#'   \code{"bayes"} (also accepted as \code{"bayesian"} or \code{"mcmc"}) runs
#'   full MCMC Bayesian estimation (Gibbs-within-HMC, NONMEM \code{METHOD=BAYES}
#'   parity): it returns posterior means with credible intervals and
#'   convergence diagnostics on \code{$bayes} rather than a point estimate, and
#'   runs standalone. Supports BSV and zero-mean inter-occasion variability
#'   (per-occasion \code{kappa}); the IOV variance posterior appears as
#'   \code{OMEGA_IOV(...)} in \code{$bayes}.
#'   SAEM fully supports inter-occasion variability (IOV / kappa) models.
#' @param covariance Logical, or \code{NULL} (the default). Whether to compute
#'   the covariance step for standard errors. \code{NULL} uses the model file's
#'   \code{[fit_options] covariance} (engine default \code{TRUE} when unset); a
#'   logical overrides the model file. Previously defaulted to \code{TRUE} and
#'   silently overrode a model file that set \code{covariance = false} (#558).
#' @param verbose Logical, or \code{NULL} (the default). Print progress during
#'   estimation. \code{NULL} uses the model file's \code{[fit_options] verbose}
#'   (engine default \code{TRUE} when unset); a logical overrides it.
#' @param bloq_method Handling of observations outside quantification limits.
#'   \code{NULL} (default) keeps whatever the model file specified;
#'   \code{"m3"} enables Beal's M3 likelihood (requires a \code{CENS} column
#'   in the data, with \code{DV} carrying the limit value: LLOQ on
#'   \code{CENS=1} rows and ULOQ on \code{CENS=-1} rows); \code{"drop"}
#'   disables M3 and treats censored rows as ordinary observations -- it does
#'   \emph{not} remove them. Each censored row is fitted at the limit value in
#'   \code{DV} as though it had been measured there, which biases the fit. To
#'   genuinely exclude them, use \code{ignore = "CENS == 1"}.
#' @param threads Number of worker threads for the per-subject parallel loops
#'   in the Rust backend (inner EBE search, SAEM, SIR). \code{NULL} (default)
#'   uses the engine's default: available cores minus one (floored at 1), capped
#'   at 8 -- most fits gain little from spreading across every core, and not all
#'   cores are equal on asymmetric platforms (e.g. Apple Silicon E-cores). Pass
#'   a positive integer to pin the count. The setting is per-call,
#'   so successive fits in the same R session can use different values.
#' @param mu_referencing Logical, or \code{NULL} (the default). When \code{TRUE},
#'   automatically
#'   detects mu-referencing from the model structure for faster and more
#'   accurate convergence. \code{NULL} uses the model file's
#'   \code{[fit_options] mu_referencing} (engine default \code{TRUE} when unset).
#'   Applies to all estimation methods. Set to
#'   \code{FALSE} to disable for comparison purposes. Detection works
#'   automatically for standard parameterizations such as
#'   \code{PARAM = THETA * exp(ETA)}; unusual parameterizations fall back
#'   silently to zero-centred ETA initialisation with no error. No changes
#'   to the \code{.ferx} model file are needed. Check \code{fit$warnings}
#'   to see which ETAs were detected.
#' @param gradient Inner-loop (per-subject EBE) gradient method.
#'   One of \code{"auto"} or \code{"fd"}, or \code{NULL} (the default) to use the
#'   model file's \code{[fit_options] gradient} (engine default \code{"auto"}
#'   when unset). \code{"ad"} is \strong{no longer supported}: the Enzyme
#'   automatic-differentiation path was retired in favour of the analytic
#'   \code{Dual2} sensitivities, so a model file (or call) carrying
#'   \code{gradient = ad} now fails validation with \code{E_AD_RETIRED} rather
#'   than being tolerated.
#'
#'   The inner optimizer is BFGS; what differs is how the gradient of the
#'   individual NLL w.r.t. ETA is computed:
#'   \itemize{
#'     \item \code{"fd"}: central finite differences, \code{2 * n_eta}
#'       forward evaluations per gradient call. Always available.
#'     \item \code{"auto"} (default): let the engine choose the best
#'       available gradient path for the model. This is the right choice
#'       for almost all uses.
#'   }
#'
#'   \strong{When to set \code{"fd"} explicitly.} Mainly for validation
#'   or for reproducibility against a known finite-difference baseline.
#'
#'   Set \code{FERX_TIME_GRADIENTS=1} in the environment to print
#'   per-gradient-call timings at the end of a fit, which is the easiest
#'   way to check which method is faster on your specific model/data.
#' @param sir Logical, or \code{NULL} (the default); run Sampling Importance
#'   Resampling after the fit to
#'   produce non-parametric parameter uncertainty intervals. Requires
#'   \code{covariance = TRUE}. \code{NULL} uses the model file's
#'   \code{[fit_options] sir} (engine default \code{FALSE} when unset). Tuning
#'   knobs (\code{sir_samples},
#'   \code{sir_resamples}, \code{sir_seed}) still flow through
#'   \code{settings}.
#' @param optimizer_trace Logical. If \code{TRUE}, write a per-iteration CSV
#'   trace to a temporary file, store its path in \code{fit$trace_path}, and
#'   read it into \code{fit$trace} as a data frame. Pass the result to
#'   \code{\link{ferx_trace}} or \code{\link{plot.ferx_fit}} to inspect
#'   optimizer progress. Default \code{FALSE}.
#' @param scale_params Logical. If \code{TRUE}, apply a per-coordinate scaling
#'   layer on top of the existing log/Cholesky parameterization, dividing each
#'   transformed coordinate by \code{|x0[i]|} (when \code{|x0[i]| > 0.1},
#'   otherwise by \code{1.0}) so the outer optimizer works in a
#'   near-unit-magnitude space. Default \code{FALSE}.
#'
#'   \strong{Why it defaults to \code{FALSE}.} The scaling is \emph{not}
#'   trajectory-transparent. Although it leaves the OFV value unchanged at any
#'   fixed point, it rescales the gradient the optimizer sees, and that
#'   gradient feeds the SLSQP overshoot cap, the quasi-Newton Hessian estimate,
#'   and the xtol/ftol termination - all of which act in the scaled coordinate
#'   system. The scaling layer only ever runs on log/Cholesky-packed
#'   coordinates (it auto-disables when any identity-packed theta is present),
#'   and for those, dividing by \code{|log value|} is counterproductive: a
#'   coordinate like \code{log(V) = log(20) ~ 3} gets scale 3, so the
#'   optimizer's unit step becomes a 3-unit move in log space - an
#'   \code{e^3 ~ 20x} multiplicative jump in V. That large step both overshoots
#'   and, through the uniform gradient cap, starves the step in every other
#'   dimension (notably OMEGA), so the fit can halt well short of the minimum.
#'   See ferx-core issue #99. The \code{FALSE} default reproduces the
#'   well-tested pre-scaling-layer behaviour.
#'
#'   \strong{When to set \code{TRUE}.} Left as an opt-in for experimentation -
#'   e.g. to A/B the scaled vs unscaled trajectory on a specific model. When
#'   enabled it applies to all outer optimizers: NLopt FOCE/FOCEI (BOBYQA,
#'   SLSQP, L-BFGS, MMA), the hand-rolled BFGS, Gauss-Newton / BHHH, and the
#'   SAEM M-step.
#' @param inits_from_nca Derive NCA-based starting values from the data before
#'   the optimizer runs, overriding the model file's defaults. Either a logical
#'   (\code{FALSE}, the default, disables it; \code{TRUE} is an alias for
#'   \code{"nca_sweep"}) or one of \code{"nca"}, \code{"nca_sweep"},
#'   \code{"nca_ebe"} to pick a strategy explicitly. Most useful with
#'   \code{settings = list(optimizer = "trust_region")} or \code{method = "gn"},
#'   where bad starting values can stall the optimizer. See
#'   \code{\link{ferx_inits_from_nca}} to inspect the values without fitting.
#' @param fd_hessian_step Positive finite numeric. Relative step size for the
#'   finite-difference Hessian used in the covariance step (default \code{0.01}).
#'   The actual perturbation for parameter \emph{i} is
#'   \code{fd_hessian_step * (1 + |x_hat[i]|)}.
#'   Increase (e.g. \code{0.1}) when the fit warns about ill-conditioned Hessian
#'   entries; decrease (e.g. \code{1e-3}) on smooth OFV surfaces where FD noise
#'   is the primary concern. Has no effect when \code{covariance = FALSE}, nor
#'   when the exact analytic R-matrix serves the fit -- which is now the default
#'   for in-scope models, so on those this argument is inert. See
#'   \code{settings = list(analytic_cov_hessian = ...)} below for the scope and
#'   for how to force finite differences back on.
#' @param settings Optional named list of fine-grained options forwarded to the
#'   Rust \code{FitOptions}. Use this to tune knobs that do not have a
#'   dedicated \code{ferx_fit()} argument. Keys are validated: values that
#'   duplicate a dedicated argument (\code{method}, \code{covariance},
#'   \code{verbose}, \code{bloq_method}, \code{bloq}, \code{threads},
#'   \code{sir}, \code{gradient}, \code{gradient_method}) are
#'   rejected -- pass them via the dedicated argument. Unknown keys and
#'   malformed values also raise an error.
#'
#'   Every key the engine accepts is listed below, with two deliberate
#'   exceptions: \code{frem_predictions} and \code{frem_sigma} are structural
#'   maps written into the generated model file by
#'   \code{\link{ferx_model_to_frem}} and cannot sensibly be hand-authored.
#'
#'   \strong{Precedence}: dedicated \code{ferx_fit()} arguments win over
#'   \code{settings}, which in turn win over the model file's
#'   \code{[fit_options]} block. A \code{warning()} is issued whenever a
#'   call-time value overrides a different value from \code{[fit_options]}.
#'   Inspect \code{fit$model_file_settings} and \code{fit$call_settings} to
#'   audit the full picture.
#'
#'   \strong{Shared options (FOCE / FOCEI / Laplace / GN / GN-hybrid / SAEM /
#'   IMP / IMPMAP / Bayes)}
#'
#'   The engine accepts the two inner-loop keys for every method that runs an
#'   inner EBE solve, which includes \code{"laplace"}, \code{"imp"},
#'   \code{"impmap"} and \code{"bayes"} as well as the five above.
#'   \describe{
#'     \item{\code{inner_maxiter}}{Per-subject EBE iteration cap (default 200).}
#'     \item{\code{inner_restarts}}{(default 1) Guarded multi-start count for the
#'       inner (per-subject EBE) optimizer, to escape a multimodal individual
#'       objective -- saturable protein binding is the classic case. A subject
#'       re-solves its EBE from this many Omega-scaled alternate seeds and keeps
#'       the lowest-objective mode; an alternate seed that reconverges to the same
#'       mode is not accepted, so unimodal subjects are bit-identical to
#'       \code{inner_restarts = 0}. The covariance step reconverges with the same
#'       setting. FOCE / FOCEI / Laplace / GN / GN-hybrid.}
#'     \item{\code{inner_optimizer}}{\code{"auto"} (default), \code{"bfgs"},
#'       \code{"lbfgs"} or \code{"nelder_mead"}. Inner-loop algorithm.
#'       \code{"auto"} uses dense BFGS, switching to limited-memory L-BFGS above
#'       32 random effects; an explicit value pins one algorithm with no
#'       switching. All gradient-based variants converge to the same EBE, so this
#'       trades per-step cost against memory rather than accuracy.}
#'     \item{\code{cov_inner_tol}}{Inner EBE-reconvergence tolerance used
#'       \emph{only} by the covariance step, decoupled from \code{inner_tol}
#'       (default: whatever \code{inner_tol} is; for LTBS models
#'       \code{min(inner_tol, 1e-8)}). The covariance R-matrix is a
#'       second-difference of the reconverged objective and so is more sensitive to
#'       EBE precision than the fit itself -- on a flat surface an EBE converged
#'       only to \code{inner_tol} can visibly perturb the standard errors. Worth
#'       reaching for on heavily-censored M3 + IOV models (try \code{1e-11}).
#'       Note the engine may emit a spurious "not used by method ... will be
#'       ignored" warning for this key; the value \emph{is} applied.}
#'     \item{\code{parameter_scaling}}{\code{"auto"} (default), \code{"none"},
#'       \code{"abs"} or \code{"rescale2"}. Parameter-scaling strategy for the
#'       outer optimizer; supersedes \code{scale_params} when not \code{"none"}.
#'       \code{"auto"} applies \code{"rescale2"} to the gradient-based optimizers
#'       that benefit (\code{nlopt_lbfgs}, \code{slsqp}) and leaves the
#'       derivative-free \code{bobyqa} unscaled, where it would distort the trust
#'       region. \code{"rescale2"} is the nlmixr2-style normalisation -- each
#'       packed coordinate divided by its bound half-range -- which markedly
#'       improves cold-start convergence. \code{"abs"} is the legacy
#'       \code{|packed value|} normalisation, equivalent to
#'       \code{scale_params = TRUE}.}
#'     \item{\code{ebe_warm_start}}{(default \code{FALSE}) When an inner EBE
#'       solve fails its BFGS step and falls back to Nelder-Mead, warm-start the
#'       simplex from the BFGS partial eta-hat rather than cold-starting from
#'       \code{eta = 0}. Substantially fewer prediction walks on fallback-heavy
#'       fits (~1.7x on a 2-cpt unidentifiable-V2 benchmark). Opt-in because it
#'       moves the fallback subjects' EBEs, which perturbs the outer optimizer's
#'       trajectory -- harmless for the derivative-free default, but it can derail
#'       a gradient-based outer optimizer into a worse basin. Validate OFV and
#'       estimates on your own model before enabling.}
#'     \item{\code{checkpoint}}{(default \code{TRUE}) and
#'       \code{checkpoint_interval_secs} (default 300). Periodically save a
#'       resume point so an interrupted run restarts from where it stopped.
#'       Restart is coarse: the population estimates at the last save become the
#'       new starting point, so a resumed run may need a few extra iterations to
#'       re-converge. A change to the model or data invalidates the checkpoint via
#'       a hash check and the run starts fresh; a run shorter than the interval
#'       leaves nothing behind.}
#'     \item{\code{iov_column}}{Name of the occasion column in the dataset (e.g.
#'       \code{"OCC"}) for inter-occasion variability.}
#'     \item{\code{iov_occasion}}{Derive the occasion partition from the model
#'       instead of a data column: \code{"column"} (default, use
#'       \code{iov_column}), \code{"dose"} to start a new occasion at each
#'       administration, or \code{"time(24, 48)"} for time-window breakpoints.
#'       \code{"dose"} and \code{"time(...)"} override \code{iov_column} (with a
#'       warning). Supported by FOCE, FOCEI, Laplace and SAEM -- not IMP/IMPMAP.}
#'     \item{\code{npde_nsim}}{(default 0, i.e. disabled) Monte-Carlo replicates
#'       per subject for the simulation-based NPDE/NPD diagnostics computed after
#'       the fit; with 0 no NPDE/NPD columns are produced. A typical value is
#'       1000, and cost scales linearly. \code{npde_seed} makes the simulation
#'       reproducible. See also \code{\link{ferx_calc_npde}} to add them
#'       post hoc.}
#'     \item{\code{inner_tol}}{Gradient-norm convergence tolerance for the inner
#'       (per-subject EBE) loop (default \code{1e-5}). A looser tolerance leaves
#'       residual noise in each subject's EBE solution, which propagates into the
#'       marginal objective the outer optimizer sees; the engine moved from
#'       \code{1e-4} to \code{1e-5} to reach NONMEM's minimum, at roughly 1.5x
#'       the per-fit cost. Note this changes the converged point, so it is not a
#'       pure cost knob. Tighter is not uniformly better -- at \code{1e-6} some
#'       ill-conditioned fits over-converge the inner Hessian and the outer
#'       optimizer can land in a worse basin.}
#'     \item{\code{fd_hessian_step}}{Also available as the dedicated
#'       \code{fd_hessian_step} argument above.}
#'     \item{\code{analytic_cov_hessian}}{\code{TRUE} (default) or
#'       \code{FALSE}. When \code{TRUE} and the model is in scope, the
#'       covariance step assembles the \strong{exact analytic R-matrix} from
#'       third-order sensitivities of the closed-form prediction instead of
#'       differencing the objective. That removes both the \code{fd_hessian_step}
#'       tuning knob and the \code{eps/h^2} differencing noise, and costs
#'       \code{2 * (n_theta + n_eta) + 1} sensitivity evaluations per subject
#'       rather than roughly \code{2 * n_free^2} objective evaluations that each
#'       re-solve every subject's inner loop. Serves \code{method = "focei"} and
#'       \code{method = "foce"} alike, from two separate assemblies rather than
#'       one shared formula -- the non-interaction case is built on the
#'       Sheiner-Beal gradient and carries no \code{log|H~|} term -- so both are
#'       exact, and \code{interaction} does not silently drop you back onto
#'       finite differences. Scope is otherwise deliberately narrow --
#'       plain analytical (closed-form) Gaussian models only: no IOV, LTBS,
#'       \code{[scaling]}, Form-C readout, \code{[initial_conditions]},
#'       \code{iiv_on_ruv}, \code{block_sigma} or custom-magnitude residual, M3
#'       censoring, FREM, non-Gaussian (TTE / categorical / Markov) endpoint, and
#'       not \code{method = "laplace"} at any node count, \code{method = "focei"}
#'       with \code{n_agq > 1}, a covariate-\code{Selected} error spec, or
#'       \code{gradient = "fd"}. Anything
#'       outside it -- including a single out-of-scope subject -- silently keeps
#'       the finite-difference stencil, all-or-nothing across the population.
#'       Set \code{FALSE} to force finite differences for an in-scope model,
#'       e.g. to reproduce standard errors from an earlier run.}
#'     \item{\code{covariance_method}}{Covariance estimator, mirroring NONMEM
#'       \code{$COV MATRIX=}: \code{"r"} (inverse Hessian, the default),
#'       \code{"s"} (score cross-product), or \code{"rsr"} (the Huber-White
#'       sandwich, robust to model misspecification). \code{"s"} and \code{"rsr"} are
#'       supported for FOCEI, FOCE and IOV fits alike, all three anchored against
#'       NONMEM \code{$COV MATRIX=S}/\code{RSR} within ~10\%. Note ferx defaults
#'       to \code{"r"} while NONMEM's \code{$COVARIANCE} default is
#'       \code{"rsr"} -- set \code{"rsr"} when reconciling against a NONMEM run
#'       that used the default \code{$COV}. No effect when
#'       \code{covariance = FALSE}.}
#'     \item{\code{covariance_fallback}}{\code{"none"} (default) or \code{"sir"}.
#'       When the finite-difference Hessian is not positive definite, \code{"sir"}
#'       runs SIR with an absolute-eigenvalue-rectified proposal instead of
#'       leaving the covariance step failed; \code{covariance_status} is then
#'       \code{"sir_fallback"} and SIR-based credible intervals are reported.}
#'   }
#'
#'   \strong{ODE models: solver method and tolerance}
#'   \describe{
#'     \item{\code{ode_method}}{Which stepper integrates the \code{[odes]}
#'       block: \code{"rk45"} (default), \code{"vern7"},
#'       \code{"rosenbrock23"}, \code{"rodas4"} or \code{"rodas5p"}. There are
#'       two independent reasons a fit's step size is small, and they want
#'       opposite fixes, so diagnose before switching. If steps stay tiny
#'       whatever \code{ode_reltol} you ask for, the model is
#'       \emph{stability}-limited (stiff -- fast reversible binding / TMDD,
#'       Michaelis-Menten with \code{KM} far below observed concentrations, long
#'       transit chains, QSP cascades): use one of the linearly implicit
#'       Rosenbrock methods, \code{"rosenbrock23"} at crude tolerances,
#'       \code{"rodas4"} at a typical \code{1e-6}--\code{1e-9}, or
#'       \code{"rodas5p"} at \code{1e-9} and tighter. If instead nearly every
#'       step is accepted and it is tightening \code{ode_reltol} that drives the
#'       step count up, the fit is \emph{accuracy}-limited and a stiff method
#'       buys nothing -- \code{"vern7"} (higher order) is the lever, worth about
#'       2.3x at \code{1e-9} on ferx-core's transit benchmark but about 1.4x
#'       \emph{slower} at default tolerances. A stiff step is not free: it costs
#'       \code{n + 1} extra right-hand-side evaluations for the
#'       finite-difference Jacobian plus an \code{n x n} factorization, where
#'       \code{n} is the size of the system actually integrated (for a
#'       continuous-time Markov endpoint that is the \code{s^2} occupancy
#'       system, not the \code{[odes]} state count). Every method is a full peer
#'       -- analytic sensitivities, time-to-event and categorical endpoints,
#'       simulation and adaptive dosing all work with all of them. Can also be
#'       set in the model file's \code{[fit_options]} block.}
#'     \item{\code{ode_reltol}}{Relative tolerance for ODE models
#'       (default \code{1e-4}; ignored for analytical PK). The default
#'       reproduces analytical closed forms in PRED to about \code{1e-4}, but
#'       the FOCE objective amplifies solver error, so an ODE-form model's OFV
#'       can differ from its analytical equivalent by several units. Set tighter
#'       (e.g. \code{1e-10}) when the ODE-form OFV must match an analytical
#'       reference; expect slower fits. Can also be set in the model file's
#'       \code{[fit_options]} block.}
#'     \item{\code{ode_abstol}}{Absolute tolerance for ODE models
#'       (default \code{1e-6}).}
#'     \item{\code{ode_max_steps}}{Maximum solver steps per integration segment
#'       (default \code{10000}). Raise if a tight \code{ode_reltol} exhausts the
#'       step budget on stiff multi-compartment systems -- or switch
#'       \code{ode_method} to a stiff stepper, which is the actual fix when
#'       stiffness rather than accuracy is capping the step.}
#'   }
#'
#'   \strong{FOCE / FOCEI / Laplace / GN / GN-hybrid: iteration cap}
#'   \describe{
#'     \item{\code{maxiter}}{Maximum outer-optimizer iterations (default 500).
#'       Not applicable to SAEM, which controls iterations via
#'       \code{n_exploration} and \code{n_convergence}.}
#'   }
#'
#'   \strong{FOCE / FOCEI / Laplace / GN-hybrid: outer optimizer}
#'   \describe{
#'     \item{\code{optimizer}}{Population-level optimizer. One of
#'       \code{"auto"} (default; picks \code{"nlopt_lbfgs"} when the exact
#'       analytic FOCE/FOCEI gradient is available and \code{"bobyqa"} when
#'       only finite differences are -- the fit reports the resolved pick as
#'       \code{"auto (<resolved>)"}),
#'       \code{"bobyqa"} (derivative-free, robust),
#'       \code{"slsqp"} (sequential quadratic programming, gradient-based),
#'       \code{"lbfgs"} / \code{"nlopt_lbfgs"} (limited-memory BFGS via
#'       NLopt, gradient-based),
#'       \code{"mma"} (method of moving asymptotes, gradient-based),
#'       \code{"bfgs"} and \code{"lbfgs"} (deprecated aliases that now select
#'       \code{"nlopt_lbfgs"} -- the hand-rolled built-in BFGS is no longer
#'       reachable and is slated for removal), or
#'       \code{"trust_region"} (trust-region Newton, gradient-based).
#'       Gradient-based optimizers use the inner gradient method set by the
#'       \code{gradient} argument; \code{"bobyqa"} does not.
#'       Not accepted by pure \code{"gn"}.}
#'     \item{\code{outer_xtol}}{Relative step tolerance for the derivative-free
#'       \code{"bobyqa"} outer optimizer (NLopt \code{xtol_rel}; default
#'       \code{1e-4}). Gradient-based optimizers use their own fixed internal
#'       tolerances. Can also be set in the model file's \code{[fit_options]}.}
#'     \item{\code{outer_ftol}}{Relative objective tolerance for the derivative-free
#'       \code{"bobyqa"} outer optimizer (NLopt \code{ftol_rel}). Unset = auto:
#'       \code{1e-8} for a pure time-to-event model (its hazard objective is
#'       evaluated exactly) and \code{1e-6} otherwise. The TTE tightening lands the
#'       frailty variance on the NONMEM/nlmixr2 minimum across a near-flat
#'       omega-squared ridge (ferx-core #469); the \code{1e-6} floor elsewhere
#'       avoids grinding on noisy ODE / FD-inner objectives, where \code{1e-8} is
#'       unreachable. Set an explicit value to pin it for every model.}
#'     \item{\code{global_search}}{Logical. When \code{TRUE}, run a global
#'       search phase before local refinement (default \code{FALSE}).
#'       Not accepted by pure \code{"gn"}.}
#'     \item{\code{global_maxeval}}{Function evaluations budget for the global
#'       search phase (default \code{0}, i.e. disabled when
#'       \code{global_search = FALSE}). Not accepted by pure \code{"gn"}.}
#'     \item{\code{stagnation_guard}}{Logical (default \code{TRUE}). Terminates
#'       the NLopt outer loop early when the OFV plateau is numerically flat.
#'       Set \code{FALSE} to let SLSQP / L-BFGS run to their own xtol/ftol or
#'       \code{maxiter}. Not accepted by pure \code{"gn"}.}
#'     \item{\code{reconverge_gradient_interval}}{Integer (default \code{0}).
#'       How often to re-solve each subject's inner EBE during the population
#'       gradient, instead of holding it fixed. \code{0} = never (cheap
#'       fixed-EBE gradient); \code{1} = every gradient evaluation; \code{N}
#'       = every \code{N}-th. Reconverging recovers the full gradient surface
#'       at roughly 5-6x the per-call cost and can help gradient-based
#'       optimizers (e.g. \code{slsqp}) past ill-conditioned non-IOV plateaus.
#'       IOV models always reconverge and ignore this setting. Not accepted by
#'       pure \code{"gn"}.}
#'   }
#'
#'   \strong{Trust-region optimizer (\code{optimizer = "trust_region"})}
#'   \describe{
#'     \item{\code{steihaug_max_iters}}{Conjugate-gradient iteration budget for
#'       the trust-region subproblem (default chosen by the engine). Only
#'       consumed when \code{optimizer = "trust_region"} is set under FOCE /
#'       FOCEI / GN-hybrid. Increase if the subproblem solver exits too early
#'       on ill-conditioned problems.}
#'   }
#'
#'   \strong{SAEM}
#'   \describe{
#'     \item{\code{n_exploration}}{Stochastic exploration phase iterations
#'       (default 150).}
#'     \item{\code{conddist}}{(default \code{FALSE}; alias
#'       \code{saem_conddist}) Run a post-fit conditional-distribution pass that
#'       estimates each subject's \code{p(eta_i | y_i)} by MCMC, reporting
#'       per-subject conditional mean and SD plus distribution-based
#'       eta-shrinkage -- the shrinkage-unbiased basis for eta diagnostics. Read
#'       it with \code{\link{ferx_conddist}}. The three keys below apply only
#'       when this is \code{TRUE}.}
#'     \item{\code{conddist_nsamp}}{Retained draws per subject (default 200;
#'       production diagnostics usually want thousands).}
#'     \item{\code{conddist_burnin}}{Draws discarded before accumulation
#'       (default 20), to forget the EBE-mode warm start.}
#'     \item{\code{conddist_keep_samples}}{(default \code{FALSE}) Retain the raw
#'       draws rather than just the summaries.}
#'     \item{\code{n_convergence}}{Averaging / convergence phase iterations
#'       (default 250).}
#'     \item{\code{n_mh_steps}}{(default 20) Metropolis-Hastings steps per subject
#'       per iteration. Also sizes the componentwise decorrelating kernel that
#'       prevents block-Omega collapse. When \code{n_leapfrog > 0} it applies to
#'       the subjects that fall back to MH. Consumed by \code{"saem"} and
#'       \code{"bayes"}.}
#'     \item{\code{adapt_interval}}{How often (in iterations) the MH proposal
#'       covariance is adapted (default 50).}
#'     \item{\code{omega_burnin}}{Initial iterations during which Omega is held
#'       fixed while the sampler warms up (default 20). Guards against Omega
#'       collapse on sparse data. Set to \code{0} to disable.}
#'     \item{\code{seed} / \code{saem_seed}}{RNG seed for the SAEM
#'       Metropolis-Hastings sampler (default 12345). Independent of
#'       \code{multi_start_seed}.}
#'     \item{\code{n_leapfrog} / \code{saem_n_leapfrog}}{Leapfrog steps for
#'       HMC proposals (default \code{0} = Metropolis-Hastings). A positive
#'       value replaces MH with one HMC proposal per subject per iteration.}
#'   }
#'
#'   \strong{Bayes (\code{"bayes"})}
#'   \describe{
#'     \item{\code{bayes_warmup}}{Warmup (burn-in + adaptation) sweeps per chain,
#'       discarded from the posterior (default \code{1000}).}
#'     \item{\code{bayes_iters}}{Retained sampling sweeps per chain, before
#'       thinning (default \code{1000}).}
#'     \item{\code{bayes_chains}}{Number of independent chains (default
#'       \code{4}); used for split-R-hat.}
#'     \item{\code{bayes_thin}}{Keep every \code{bayes_thin}-th sampling draw
#'       (default \code{1}).}
#'     \item{\code{bayes_seed}}{Base RNG seed for the Bayes sampler. Independent
#'       of \code{seed} / \code{saem_seed}.}
#'   }
#'
#'   \strong{Gauss-Newton (\code{"gn"} / \code{"gn_hybrid"})}
#'   \describe{
#'     \item{\code{gn_lambda}}{Levenberg-Marquardt damping factor (default
#'       \code{0.01}). Larger values make steps more conservative. Accepted
#'       by both \code{"gn"} and \code{"gn_hybrid"}.}
#'   }
#'
#'   \strong{Importance Sampling (\code{"imp"})}
#'
#'   By default \code{"imp"} is a Monte-Carlo EM estimator (NONMEM
#'   \code{METHOD=IMP}); set \code{imp_eval_only = TRUE} to evaluate the marginal
#'   \code{-2 log L} at fixed parameters instead (NONMEM \code{EONLY=1}).
#'   \describe{
#'     \item{\code{imp_eval_only}}{Logical; \code{TRUE} evaluates \code{-2 log L}
#'       at the fixed input parameters without estimating (NONMEM \code{EONLY=1};
#'       must be the terminal chain stage). \code{FALSE} (default) estimates.}
#'     \item{\code{imp_iterations}}{Number of Monte-Carlo EM iterations, ignored
#'       when \code{imp_eval_only} (default 200).}
#'     \item{\code{imp_averaging}}{Number of final iterations whose parameters are
#'       averaged into the reported estimate, ignored when \code{imp_eval_only}
#'       (default 50).}
#'     \item{\code{imp_samples}}{Importance samples drawn per subject (default
#'       1000). Halving the Monte-Carlo SE requires quadrupling this value.}
#'     \item{\code{imp_proposal_df}}{Degrees of freedom for the Student-t
#'       proposal distribution (default \code{5}); the string \code{"normal"}
#'       (or \code{"mvn"}) selects a multivariate-normal proposal.}
#'     \item{\code{imp_seed}}{RNG seed for the IS step (default chosen by the
#'       engine).}
#'     \item{\code{imp_auto}}{(default \code{TRUE}) Adaptive sample count (NONMEM
#'       \code{AUTO}). With this on, \code{imp_samples} is the \emph{starting}
#'       count and is ramped up (doubling per iteration, capped at 10000) while the
#'       objective's Monte-Carlo SE exceeds 1.0. Strongly recommended for FREM and
#'       other high-dimensional models, where a fixed count leaves a
#'       sample-count-dependent bias in the typical-value and Omega estimates.
#'       Note this makes the "quadruple the samples to halve the MC SE" rule of
#'       thumb apply only when it is \code{FALSE}.}
#'     \item{\code{imp_defensive_alpha}}{(default \code{0}, i.e. off; alias
#'       \code{impmap_defensive_alpha}) Defensive-mixture weight in
#'       \code{[0, 1)}. Each subject draws this fraction of its samples from the
#'       prior \code{N(0, Omega)} rather than the mode-centred proposal, and every
#'       sample is scored under the mixture density. That bounds the importance
#'       weights, so a weakly-identified subject cannot hijack the weighted M-step
#'       and walk theta to its bounds. Try \code{0.1}. Enabling it disables Sobol
#'       QMC and raises the per-subject ESS floor.}
#'     \item{\code{iscale_min}}{(default \code{0.1}) and \code{iscale_max}
#'       (default \code{10}). Bounds on the proposal scaling factor for adaptive
#'       importance sampling (NONMEM \code{ISCALE_MIN} / \code{ISCALE_MAX}): the
#'       proposal covariance is multiplied by \code{s^2} with \code{s} chosen from
#'       this interval to maximise per-subject ESS. Set both to \code{1} to
#'       disable.}
#'     \item{\code{frem_rao_blackwell}}{(default \code{TRUE}) FREM only:
#'       Rao-Blackwellise the covariate ETAs -- integrate them analytically and
#'       importance-sample only the PK ETAs. Strongly recommended, since
#'       brute-force sampling of the near-singular covariate dimensions has very
#'       poor ESS. Set \code{FALSE} only to diagnose the Rao-Blackwell path
#'       against the full-dimensional sampler.}
#'     \item{\code{imp_low_ess_threshold}}{ESS fraction below which a subject is
#'       flagged in \code{fit$importance_sampling$low_ess_subject_ids} (default
#'       \code{0.1}, i.e. \code{10\%} of \code{imp_samples}).}
#'   }
#'
#'   \strong{IMPMAP (\code{"impmap"} estimator)}
#'   \describe{
#'     \item{\code{impmap_iterations}}{Number of Monte-Carlo EM iterations
#'       (default 200).}
#'     \item{\code{impmap_samples}}{Importance samples drawn per subject per
#'       iteration (default 300).}
#'     \item{\code{impmap_proposal_df}}{Proposal degrees of freedom (default
#'       \code{4}). A finite value \code{>= 1} gives a heavier-tailed Student-t,
#'       which is the ferx default; the string \code{"normal"} (or \code{"mvn"})
#'       gives a multivariate-normal proposal, which is \emph{NONMEM's} default.
#'       ferx diverges deliberately: a Gaussian's lighter tails under-cover the
#'       posterior of weakly-identified parameters, so the importance weights blow
#'       up in the tail and bias the M-step moments.}
#'     \item{\code{impmap_averaging}}{Number of final iterations whose parameters
#'       are averaged into the reported estimate (default 50).}
#'     \item{\code{impmap_seed}}{RNG seed for the IMPMAP sampling (default chosen
#'       by the engine).}
#'     \item{\code{impmap_auto}}{(default \code{TRUE}) Adaptive sample count, as
#'       \code{imp_auto} above -- \code{impmap_samples} is the \emph{starting}
#'       count.}
#'     \item{\code{impmap_mceta}}{(default 0) Additional random starting points
#'       for the per-subject MAP optimization (NONMEM \code{MCETA}). Each start
#'       draws eta from \code{N(0, Omega)} and the lowest individual NLL wins; 0
#'       means a single warm start. 3 is a good choice for high-dimensional models
#'       (e.g. FREM with 5 or more ETAs).}
#'     \item{\code{impmap_sobol}}{(default \code{FALSE}) Use Sobol quasi-random
#'       sequences (with Cranley-Patterson randomization) for the importance draws
#'       instead of pseudo-random, giving more uniform posterior coverage per
#'       sample. Applies only to multivariate-normal proposals
#'       (\code{impmap_proposal_df = "normal"}); under the Student-t default it is
#'       inert.}
#'     \item{\code{impmap_low_ess_threshold}}{ESS fraction below which a subject
#'       is flagged as poorly sampled (default \code{0.1}).}
#'     \item{\code{impmap_trace}}{Logical; when \code{TRUE}, collect
#'       per-iteration parameter values into \code{fit$impmap_trace}
#'       (analogous to NONMEM \code{.ext} output). Default \code{FALSE}.}
#'   }
#'
#'   \strong{SIR (Sampling Importance Resampling)}
#'   \describe{
#'     \item{\code{sir_samples}}{Candidate draws for SIR (default 1000).}
#'     \item{\code{sir_resamples}}{Resampled draws retained for CI computation
#'       (default 250).}
#'     \item{\code{sir_seed}}{RNG seed for the SIR step (default chosen by the
#'       engine). Independent of \code{seed} / \code{saem_seed}.}
#'     \item{\code{sir_df}}{(default 5) Degrees of freedom for the SIR
#'       multivariate Student-t proposal; higher values approach a normal
#'       proposal.}
#'     \item{\code{sir_keep_samples}}{(default \code{FALSE}) Retain the resampled
#'       parameter vectors, which \code{\link{ferx_simulate_with_uncertainty}}
#'       requires.}
#'   }
#'
#'   \strong{Multi-start optimization}
#'   \describe{
#'     \item{\code{n_starts}}{Number of optimizer starts (default \code{1},
#'       i.e. single start). When \code{> 1}, starts 2..n are initialized from
#'       log-space perturbations of the model's initial values; the best OFV
#'       wins.}
#'     \item{\code{start_sigma}}{Perturbation spread in log-space (default
#'       \code{0.3}, approx. \code{30\%} CV). Log-packed thetas are multiplied by
#'       \code{exp(N(0, start_sigma))}; identity-packed thetas are shifted by
#'       \code{start_sigma * N(0,1)}.}
#'     \item{\code{multi_start_seed}}{RNG seed for the start-point perturbation
#'       (default \code{42}). Independent of \code{seed} / \code{saem_seed}.}
#'   }
#' @param ignore Character vector of filter expressions that exclude a record
#'   when the expression evaluates to true (NONMEM \code{$DATA IGNORE=}
#'   equivalent). Expressions use standard comparison operators
#'   (\code{==}, \code{!=}, \code{<}, \code{<=}, \code{>}, \code{>=}) on
#'   column names. Use \code{&&} within one expression for AND; for OR use
#'   multiple elements. These conditions are merged with any
#'   \code{[data_selection] ignore} rules in the model file.
#'   Example: \code{ignore = c("DV < 0.001", "EVID != 0 && MDV == 1")}.
#' @param accept Character vector of filter expressions. A record is kept
#'   only when all conditions pass; excluded if any fail (NONMEM
#'   \code{$DATA ACCEPT=} equivalent). Merged with model-file \code{accept}
#'   rules.
#' @param ignore_ids Numeric or character vector of subject IDs to exclude
#'   entirely. Sugar for \code{ignore = "ID == <id>"} applied per-subject.
#'   Merged with \code{[data_selection] ignore_subjects} from the model file.
#'
#' @param output Optional path to a \code{.fitrx} file. When non-\code{NULL},
#'   \code{\link{ferx_save_fit}} is invoked on the result so the fit is
#'   persisted to disk in a portable, cross-language bundle (zip of JSON +
#'   CSV) before the function returns. Equivalent to calling
#'   \code{ferx_save_fit(fit, output)} after the fit. See the format
#'   reference in the ferx-core docs for the on-disk schema.
#' @param include_data Logical. Only meaningful with \code{output}. When
#'   \code{TRUE}, embeds the input \code{data} CSV verbatim inside the
#'   \code{.fitrx} bundle so the file is self-contained. Default
#'   \code{FALSE}.
#'
#' @return A named list. Where to find commonly-needed quantities:
#'
#'   | What you want | Accessor | Shape | When present |
#'   |---|---|---|---|
#'   | Residuals, PRED, IPRED, diagnostics | `fit$sdtab` | one row per observation | always |
#'   | Covariates echoed per observation (via `[output]`) | `fit$sdtab` | one row per observation, LOCF | declared in `[output]` |
#'   | Raw covariate values for all dataset records | `fit$covtab` | one row per dataset record (doses + obs) | model has `[covariates]` block |
#'   | ETA / EBE values per observation row | `fit$sdtab` (`ETA_CL`, `ETA_V`, ...) | one row per observation, value repeated | always |
#'   | ETA / EBE values per subject | `fit$ebe_etas` | one row per subject | model declares etas |
#'   | Individual PK parameters per subject | `fit$individual_estimates` | one row per subject | always |
#'
#'   Full slot descriptions:
#'   \item{converged}{Logical; did the optimizer converge}
#'   \item{method}{Estimation method used}
#'   \item{n_iterations}{Number of outer iterations run}
#'   \item{ofv}{Objective function value (-2 log-likelihood)}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{theta}{Named numeric vector of fixed effect estimates}
#'   \item{omega}{Between-subject variability covariance matrix. Row and column
#'     names are the declared ETA names (e.g. \code{"ETA_CL"}); fallback is
#'     \code{"OMEGA(1,1)"} when names are absent.}
#'   \item{sigma}{Residual error parameter estimates on the
#'     \emph{standard-deviation} scale. For a proportional component
#'     \code{CV\% = sigma * 100}; for an additive component the value is in
#'     observation units. Variance is \code{sigma^2} for either type.}
#'   \item{sigma_names}{Names of the sigma parameters as declared in the
#'     \code{[parameters]} block of the \code{.ferx} file. Parallel to
#'     \code{sigma}.}
#'   \item{sigma_types}{Per-component error type label, parallel to
#'     \code{sigma}: \code{"proportional"} or \code{"additive"}. Combined
#'     error models report \code{c("proportional", "additive")} in that order.}
#'   \item{se_theta}{Standard errors for theta (NULL if covariance step failed)}
#'   \item{se_omega}{Standard errors for omega elements. For diagonal omega,
#'     a vector of length n_eta (one per variance). For block omega, a vector
#'     of length n_eta*(n_eta+1)/2 (full lower triangle, column-major) including
#'     off-diagonal covariance SEs.}
#'   \item{se_sigma}{Standard errors for sigma (on the SD scale, like
#'     \code{sigma} itself)}
#'   \item{sdtab}{Data frame with ID, TIME, DV, PRED, IPRED, CWRES, IWRES,
#'     EBE_OFV, N_OBS; OCC if any subject carries an occasion column; CENS
#'     if any rows are LOQ-censored; CMT if the dataset has more than one
#'     endpoint (any observation with \code{CMT != 1}).  TAFD (time after
#'     first dose) and TAD (time after last dose, SS-aware) are included
#'     automatically.  Columns declared in \code{[derived]} (per-row
#'     expressions, aggregates, integrals) and \code{[output]} (individual
#'     PK parameters or covariates echoed by name) appear at the end of the
#'     data frame in declaration order.
#'     ETA columns (\code{ETA_CL}, \code{ETA_V}, etc.) are included
#'     automatically - the same EBE value is repeated for every observation
#'     row of that subject. For a compact per-subject view use
#'     \code{fit$ebe_etas}; for individual PK parameter values use
#'     \code{fit$individual_estimates}.}
#'   \item{covtab}{Data frame echoing the declared covariate columns, present
#'     only when the model has a \code{[covariates]} block. Columns: \code{ID},
#'     \code{TIME}, \code{EVID}, then one column per declared covariate (in
#'     declaration order). Unlike \code{sdtab} (observation rows only), it has
#'     one row per input dataset record, including dose and other-event rows;
#'     missing covariate values are \code{NA}. \code{NULL} when no
#'     \code{[covariates]} block is declared.}
#'   \item{covariate_types}{Named character vector mapping each declared
#'     covariate to its type, \code{"continuous"} or \code{"categorical"} (the
#'     type from the \code{[covariates]} block). Parallel to the covariate
#'     columns of \code{covtab}. \code{NULL} when no \code{[covariates]} block
#'     is declared.}
#'   \item{ebe_etas}{Data frame with one row per subject containing the BSV
#'     empirical Bayes estimates: \code{ID} plus one column per eta named
#'     after the model's eta declarations (e.g. \code{ETA_CL}, \code{ETA_V}).
#'     \code{NULL} when the model declares no etas.}
#'   \item{individual_estimates}{Data frame with one row per subject containing
#'     each subject's individual parameter values: \code{ID} plus one column
#'     per parameter declared in the \code{[individual_parameters]} block
#'     (e.g. \code{CL}, \code{V}, \code{KA}). Computed by evaluating the
#'     individual-parameter expressions at the subject's EBE eta with kappa
#'     fixed at zero, so the result reflects each subject's typical value
#'     under the model's covariate effects. For inter-occasion variation see
#'     \code{ebe_kappas}.}
#'   \item{warnings}{Character vector of warnings}
#'   \item{warnings_structured}{Data frame with columns \code{severity}
#'     (\code{"critical"}, \code{"warning"}, or \code{"info"}),
#'     \code{category} (a fixed vocabulary such as \code{"convergence"},
#'     \code{"covariance_step"}, \code{"dw_autocorrelation"}),
#'     \code{message}, and \code{source_method} (the estimation stage that
#'     emitted the warning, e.g. \code{"FOCEI"}; empty string when not
#'     applicable). Severity/category are assigned by the engine for its
#'     own warnings and by the R diagnostics layer for the condition-number
#'     and ETA-normality checks. Inspect with \code{\link{ferx_get_warnings}}.}
#'   \item{sir_ess}{SIR effective sample size (NULL if SIR not run)}
#'   \item{sir_ci_theta, sir_ci_omega, sir_ci_sigma}{SIR \code{95\%} CI matrices
#'     with columns \code{lower} and \code{upper} (NULL if SIR not run)}
#'   \item{importance_sampling}{Importance-sampling marginal log-likelihood
#'     diagnostics, populated only when the method chain ends with an
#'     \code{"imp"} stage (\code{NULL} otherwise). A list with
#'     \code{minus2_log_likelihood} (the IS estimate of \eqn{-2 \log L},
#'     comparable across models in the same nesting family),
#'     \code{mc_standard_error} (Monte-Carlo SE; scales as
#'     \eqn{1/\sqrt{imp\_samples}}), \code{n_samples},
#'     \code{proposal_df}, \code{ess_min} and \code{ess_median}
#'     (per-subject effective sample size as a fraction of
#'     \code{imp_samples}), \code{kappa_treatment}
#'     (\code{"not_applicable"}, \code{"fixed_at_mode"} for partial-marginal
#'     IOV fits, or \code{"marginalized"}), and parallel vectors
#'     \code{low_ess_subject_ids} / \code{low_ess_subject_frac} flagging
#'     subjects whose ESS fraction fell below \code{imp_low_ess_threshold}.}
#'   \item{trace_path}{Path to the optimizer trace CSV, or \code{NULL} when
#'     \code{optimizer_trace = FALSE}. Pass to \code{\link{ferx_trace}}
#'     or \code{\link{plot.ferx_fit}}.}
#'   \item{trace}{The optimizer trace as a data frame (same columns as
#'     \code{\link{ferx_trace}}), or \code{NULL} when
#'     \code{optimizer_trace = FALSE}. Survives \code{\link{ferx_save_fit}} /
#'     \code{\link{ferx_load_fit}}, unlike \code{trace_path}.}
#'   \item{ebe_convergence_warnings}{Number of outer iterations in which too
#'     many EBEs were unconverged (step was rejected by the guard).}
#'   \item{max_unconverged_subjects}{Worst-case number of unconverged subjects
#'     in a single outer iteration.}
#'   \item{total_ebe_fallbacks}{Total Nelder-Mead fallback invocations across
#'     all subjects and iterations.}
#'   \item{covariance_status}{String: \code{"computed"}, \code{"failed"},
#'     \code{"not_requested"}, or \code{"sir_fallback"} (the finite-difference
#'     Hessian was not positive definite and \code{covariance_fallback = "sir"}
#'     produced SIR-based intervals).}
#'   \item{bic_inputs}{Named list of the free-parameter tally by Delattre
#'     class - \code{n_obs}, \code{theta_random}, \code{theta_fixed},
#'     \code{omega}, \code{kappa}, \code{sigma}, \code{sigma_random} - which
#'     \code{\link{ferx_bic}} turns into the mixed / IIV / random BIC variants.
#'     The five counts sum to \code{n_parameters}.}
#'   \item{left_init}{The outer optimizer's own verdict on whether it left the
#'     initial estimates: \code{TRUE}, \code{FALSE}, or \code{NA} when the
#'     method recorded none.}
#'   \item{stalled_at_init}{\code{TRUE} when no free parameter moved more than
#'     1\% of its initial value, \code{NA} when the fit carries nothing to
#'     judge by. Prefers \code{left_init} to comparing estimates.}
#'   \item{estimate_near_boundary}{\code{TRUE} when a theta is pinned to a
#'     declared bound - the predicate \code{\link{ferx_bootstrap}} applies as
#'     \code{skip_estimate_near_boundary}.}
#'   \item{max_abs_correlation}{Largest absolute off-diagonal of the covariance
#'     matrix in correlation form, on the natural theta / OMEGA / SIGMA scale.
#'     \code{NA} without a covariance matrix carrying two or more free
#'     parameters. See \code{\link{check_strictness}}.}
#'   \item{shrinkage_eta}{Numeric vector of ETA shrinkage per random effect
#'     (\code{1 - SD(eta_hat_k) / sqrt(omega_kk)}). \code{NA} when
#'     \code{omega_kk = 0} or fewer than 2 subjects.}
#'   \item{shrinkage_eps}{EPS shrinkage: \code{1 - SD(IWRES)}. \code{NA} when
#'     fewer than 2 valid residuals.}
#'   \item{dw_statistic}{Pooled Durbin-Watson statistic for IWRES within
#'     subjects. 2.0 = no autocorrelation; < 1.5 = positive autocorrelation
#'     (possible missing dynamics); > 2.5 = negative autocorrelation (possible
#'     over-parameterisation). \code{NA} when no subject has \eqn{\ge 2} valid
#'     IWRES values.}
#'   \item{iwres_lag1_r}{Pooled lag-1 Pearson correlation of IWRES within
#'     subjects. Positive values indicate under-fitting, negative values
#'     indicate over-fitting. \code{NA} when insufficient observations.}
#'   \item{wall_time_secs}{Total wall-clock time for the fit in seconds.}
#'   \item{model_name}{Model name from the \code{.ferx} file. Falls back to
#'     the model file's basename (without extension) when the file declares
#'     no name.}
#'   \item{data_name}{Dataset name, derived as the basename of the \code{data}
#'     path with its extension stripped.}
#'   \item{gradient}{The inner-loop gradient method as requested by the
#'     caller (usually \code{"auto"} or \code{"fd"}). The
#'     \emph{resolved} method may differ - see \code{gradient_used} and the
#'     \code{gradient} argument.}
#'   \item{gradient_used}{The inner-loop gradient method the engine actually
#'     used: \code{"analytic"}, \code{"fd"}, or \code{"N/A"} (derivative-free /
#'     sampling step). When \code{gradient = "auto"} this records which
#'     branch resolved at fit time. \code{NA} on older engine binaries that
#'     did not populate this field.}
#'   \item{gradient_method_inner}{Engine label for the inner-loop gradient
#'     method. Prefer \code{gradient_used} for display.}
#'   \item{gradient_method_outer}{Engine label for the population-level
#'     optimizer's gradient method.}
#'   \item{ferx_version}{ferx-core library version string.}
#'   \item{model_text}{Verbatim content of the \code{.ferx} model file as a
#'     character string. \code{NULL} for in-memory fits that were never
#'     associated with a file, or for fits produced by older engine versions.
#'     Used by \code{\link{ferx_runlog}}.}
#'   \item{theta_init}{Numeric vector of initial theta values as supplied to
#'     the optimizer, parallel to \code{theta} and \code{theta_names}.
#'     Empty for in-memory fits and older \code{.fitrx} bundles.}
#'   \item{omega_init}{Initial omega variance matrix, same layout as
#'     \code{omega}. A zero matrix for in-memory fits and older bundles.}
#'   \item{omega_init_dim}{Number of rows/columns of \code{omega_init}.}
#'   \item{sigma_init}{Numeric vector of initial sigma values, parallel to
#'     \code{sigma} and \code{sigma_names}.}
#'   \item{obs_time_range}{Length-2 numeric vector \code{c(min_time, max_time)}
#'     across all observation records, or \code{NULL} when there are no
#'     observations. Used by \code{\link{ferx_runlog}} for the data-summary
#'     header.}
#'   \item{final_gradient}{Numeric vector containing the gradient of the
#'     objective function at the best-OFV parameter point (packed space:
#'     log-theta, Cholesky-omega, log-sigma). \code{NULL} for derivative-free
#'     optimizers (BOBYQA), built-in BFGS, Gauss-Newton, and SAEM. Use
#'     \code{\link{ferx_runlog}} to interpret with a convergence threshold.}
#'   \item{optimizer_label}{Human-readable label for the outer optimizer used
#'     (e.g. \code{"LBFGS"}, \code{"BOBYQA"}). Mirrors the \code{optimizer}
#'     enum value from ferx-core. Used by \code{\link{ferx_runlog}}.}
#'   \item{n_starts}{Integer. Number of multi-start runs attempted. 1 for
#'     single-start fits.}
#'   \item{multi_start_seed}{Numeric scalar (or \code{NULL}) giving the random
#'     seed used for multi-start parameter perturbation.}
#'   \item{saem_seed}{Numeric scalar (or \code{NULL}) giving the SAEM
#'     stochastic seed. \code{NULL} for non-SAEM methods.}
#'   \item{sir_seed_used}{Numeric scalar (or \code{NULL}) giving the SIR
#'     resampling seed. \code{NULL} when SIR was not run.}
#'   \item{imp_seed}{Numeric scalar (or \code{NULL}) giving the importance
#'     sampling seed. \code{NULL} when IS was not run.}
#'   \item{bloq_method_label}{Character string describing the LOQ-censoring
#'     handling method used (\code{"m3"}, \code{"drop"}, etc.).}
#'   \item{outer_maxiter}{Integer. Maximum number of outer-loop iterations
#'     passed to the optimizer.}
#'   \item{outer_gtol}{Numeric. Gradient tolerance passed to the outer
#'     optimizer. \code{NA} for derivative-free optimizers.}
#'   \item{inits_from_nca}{Character string giving the NCA initialisation
#'     method used (\code{"nca"}, \code{"nca_sweep"}, or \code{"nca_ebe"}),
#'     or \code{NULL} when model-file initial values were used directly.}
#'   \item{covariate_names}{Character vector of non-standard column names
#'     present in the dataset (beyond \code{ID}, \code{TIME}, \code{DV},
#'     \code{EVID}, \code{AMT}, \code{CMT}, \code{RATE}, \code{MDV},
#'     \code{II}, \code{SS}, \code{CENS}, \code{OCC}). Empty when no
#'     covariates are present. Contains only non-standard covariate columns;
#'     see \code{input_columns} for the complete ordered column list.}
#'   \item{input_columns}{Character vector of all column headers from the
#'     data CSV in file order, analogous to NONMEM \code{$INPUT}. Includes
#'     both standard columns (\code{ID}, \code{TIME}, \code{DV}, etc.) and
#'     any covariate columns. Empty for in-memory fits that never read a file
#'     and absent after a \code{\link{ferx_save_fit}}/\code{\link{ferx_load_fit}}
#'     round-trip on older bundles. Use this field (not \code{covariate_names})
#'     when you need the complete ordered column list.}
#'   \item{cov_matrix}{Full parameter covariance matrix as a named numeric
#'     matrix (params ? params). Row/column names use declared variable names
#'     (\code{"TVCL"}, \code{"ETA_CL"}, \code{"EPS_PROP"}); fallback is
#'     \code{"OMEGA(1,1)"} / \code{"SIGMA(1)"}. \code{NULL} when covariance
#'     step was not run or failed. See \code{cor_matrix} for the derived
#'     correlation matrix.}
#'   \item{cor_matrix}{Correlation matrix derived from \code{cov_matrix}, same
#'     dimnames. A correlation close to \eqn{\pm 1} between two parameters
#'     flags a structural identifiability problem in the model. \code{NULL}
#'     when \code{cov_matrix} is \code{NULL}.}
#'   \item{estimates}{Tidy data frame of all estimated parameters (theta,
#'     omega diagonal, sigma, and - for IOV models - kappa diagonal), with
#'     columns \code{param}, \code{transform}, \code{estimate}, \code{se},
#'     \code{rse_pct}, \code{lower_95}, \code{upper_95},
#'     \code{estimate_natural}, \code{lower_95_natural},
#'     \code{upper_95_natural}, \code{init_as_sd}, \code{weight}. SE-derived
#'     and natural-scale columns are \code{NA} when not applicable or when the
#'     covariance step was not run. \code{weight} carries the sample-size
#'     weight expression of a weighted kappa row (see \code{kappa_weights})
#'     and is \code{NA} on every other row; on such a row \code{estimate} is
#'     the \emph{unweighted} variance, so the between-occasion SD at a weight
#'     \code{W} is \code{sqrt(estimate / W)}, not \code{sqrt(estimate)}.}
#'   \item{eta_cov}{Data frame of Pearson correlations between each ETA and
#'     each constant-per-subject numeric covariate in the dataset, with
#'     columns \code{eta}, \code{covariate}, \code{r}, \code{p_val},
#'     \code{flag}, sorted by descending \code{|r|}. \code{NULL} when the
#'     model declares no etas, the dataset has no usable numeric covariates,
#'     or the data file can no longer be read from \code{data_path}.}
#'   \item{eigenvalues}{Numeric vector of eigenvalues of the correlation matrix
#'     of estimated (non-fixed) parameters, sorted descending. Computed by the
#'     Rust backend. \code{NULL} when the covariance step was not run, failed,
#'     or fewer than two free parameters exist.}
#'   \item{condition_number}{Ratio of the largest to smallest eigenvalue of the
#'     correlation matrix of estimated parameters. Values above 1000 are flagged
#'     as potentially ill-conditioned and also appear in \code{warnings}.
#'     \code{Inf} when the smallest eigenvalue is non-positive. \code{NULL}
#'     when \code{eigenvalues} is \code{NULL}.}
#'   \item{eta_normality}{Data frame with Shapiro-Wilk normality test for
#'     each ETA: columns \code{eta}, \code{W}, \code{p_val}, \code{flag}.
#'     A \code{[!]} flag appears when \code{p < 0.05}. \code{NULL} only
#'     when \code{ebe_etas} is missing or has no eta columns. When an ETA
#'     has fewer than 3 finite values (e.g. tiny populations or fixed
#'     etas), its row is still returned but \code{W} and \code{p_val} are
#'     \code{NA}.}
#'   \item{omega_iov}{IOV variance matrix for kappa parameters (\code{NULL} if no IOV).}
#'   \item{kappa_names}{Names of kappa (IOV) parameters (\code{NULL} if no IOV).}
#'   \item{se_kappa}{Standard errors for kappa parameters: length \code{d}
#'     (diagonal kappa) or \code{d*(d+1)/2} (block_kappa, lower-triangle order).
#'     \code{NULL} if covariance step not run or no IOV.}
#'   \item{model_structure}{Named list returned by the Rust engine, derived
#'     from the parsed \code{CompiledModel} so it reflects exactly what
#'     ferx-core ran. Fields: \code{theta_names} (character vector of
#'     population parameter names), \code{model_type} (short label such as
#'     \code{"1-cpt oral"}, \code{"ODE"}, or \code{"compartment-free"}, or
#'     \code{NULL} when the structural form is not one of the known
#'     analytical PK families), \code{iiv}
#'     (omega names), \code{iov} (kappa names), \code{iov_weights} (per-kappa
#'     sample-size weight expressions, present only when some kappa declares
#'     \code{weight = }), \code{residual} (error type string). Use
#'     \code{\link{ferx_model_inspect}} to view this before or after
#'     fitting.}
#'   \item{shrinkage_kappa}{Named numeric vector of pooled shrinkage values for
#'     kappa EBEs (one entry per kappa parameter). \code{NULL} when the model
#'     has no IOV.}
#'   \item{shrinkage_kappa_by_occ}{Data frame with column \code{occ} (1-based
#'     occasion index) and one numeric column per kappa parameter giving the
#'     shrinkage at each occasion. Use this for per-occasion IOV diagnostics.
#'     \code{NULL} when the model has no IOV.}
#'   \item{ebe_kappas}{Data frame with columns \code{ID}, \code{OCC}, and one
#'     column per kappa parameter containing per-subject per-occasion kappa
#'     EBEs. \code{ID} carries the original subject identifier (matching
#'     \code{sdtab}); \code{OCC} carries the labeled occasion in first-seen
#'     order. \code{NULL} when no IOV.}
#'   \item{omega_param_corr}{Parameter-level correlation matrix for block omega.
#'     Uses the bivariate lognormal formula for lognormal pairs and falls back
#'     to the eta-level formula for additive or unknown parameterisations.
#'     \code{NULL} when omega is diagonal (no off-diagonal elements to report).}
#'   \item{omega_iov_param_corr}{Same as \code{omega_param_corr} but for the
#'     IOV kappa block. \code{NULL} when kappa is diagonal or absent.}
#'   \item{eta_names}{Character vector of ETA parameter names as declared in the
#'     model file (e.g. \code{"ETA_CL"}, \code{"ETA_V"}). Used to label OMEGA
#'     rows in \code{fit$estimates} and \code{print.ferx_fit()}.}
#'   \item{eta_log_transformed}{Logical vector of length \code{n_eta}; \code{TRUE}
#'     when the eta is lognormally parameterised (\code{THETA * exp(ETA)}),
#'     \code{FALSE} for additive or unknown parameterisations. \code{NULL} when
#'     the engine does not supply the information.}
#'   \item{call_settings}{Named list of the effective settings passed to Rust,
#'     with values typed back to their natural R types (logical, numeric, or
#'     character). Includes merged defaults such as \code{optimizer_trace} and
#'     \code{scale_params}. Empty list when no settings were supplied.}
#'   \item{model_file_settings}{Named list of raw \code{key = value} pairs from
#'     the model file's \code{[fit_options]} block, as character strings
#'     (before Rust type coercion). Empty list when no \code{[fit_options]}
#'     block is present. Use this alongside \code{call_settings} to audit which
#'     values the model file requested and which were overridden at the R call
#'     site. A \code{warning()} is issued automatically for each key that
#'     appears in both sources with a different value.}
#'   \item{uses_sde}{Logical; \code{TRUE} when the model file contained a
#'     \code{[diffusion]} block and the fit used the Extended Kalman Filter
#'     (EKF) likelihood. \code{FALSE} for standard ODE/analytical models.
#'     When \code{TRUE}, one or more \code{DIFF_*} entries appear in
#'     \code{theta} and \code{theta_names} - these are diffusion
#'     \emph{variances} (not standard deviations) for the named state
#'     compartments.}
#'   \item{omega_init_as_sd}{Logical vector, one entry per BSV eta. \code{TRUE}
#'     when the corresponding \code{omega} line in the model file was annotated
#'     with \code{(sd)}, meaning the initial value was specified as a standard
#'     deviation rather than a variance. \code{print.ferx_fit()} appends
#'     \code{[initial specified as SD]} to those rows.}
#'   \item{sigma_init_as_sd}{Logical vector, one entry per sigma component.
#'     Same semantics as \code{omega_init_as_sd} for residual-error parameters.}
#'   \item{kappa_init_as_sd}{Logical vector, one entry per IOV kappa parameter.
#'     Same semantics as \code{omega_init_as_sd} for inter-occasion variability.}
#'   \item{kappa_weights}{Named character vector, one entry per IOV kappa,
#'     holding the sample-size weight expression written as
#'     \code{kappa K ~ <var> weight = <expr>} - the between-treatment-arm
#'     variability (BTAV) form of every longitudinal MBMA, in which
#'     \code{kappa_ik ~ N(0, Omega_IOV / W_ik)}. \code{NA} for a kappa that
#'     carries no weight, and \code{NULL} for a model where no kappa does.
#'     The estimate in \code{omega_iov} stays the \emph{unweighted} variance
#'     (the quantity a published analysis reports); the arm's effective
#'     standard deviation is \code{sqrt(omega_iov[i, i] / W)}.}
#'   \item{kappa_weight_typical}{Named numeric vector parallel to
#'     \code{kappa_weights}: the median weight over this dataset's
#'     subject-occasions, i.e. the arm size the reported effective SD refers
#'     to. \code{NA} where the kappa carries no weight, \code{NULL} when
#'     \code{kappa_weights} is.}
#'   \item{saem_n_subjects_hmc}{Integer count of subjects whose SAEM E-step
#'     proposals were generated by HMC (i.e. \code{n_leapfrog > 0} and the AD
#'     build is active). Subjects that fell back to MH (e.g. ODE models) are
#'     not counted. \code{NA_integer_} when the method is not SAEM, when
#'     \code{n_leapfrog = 0} (MH only), or when AD is unavailable.}
#'   \item{saem_mu_ref_m_step_evals_saved}{Number of OFV evaluations saved by
#'     mu-referencing in the SAEM M-step. Only populated for SAEM fits with
#'     \code{mu_referencing = TRUE}; \code{NULL} otherwise. A diagnostic for
#'     understanding SAEM efficiency.}
#'   \item{covariance_n_evals_estimated}{Estimated number of OFV calls required
#'     for the covariance step, computed before the step runs on large models.
#'     \code{NULL} when the model is small (exact count used instead) or when
#'     \code{covariance = FALSE}.}
#'   \item{n_threads_used}{Integer. Actual number of parallel threads used by
#'     the engine during fitting. May be lower than the requested \code{threads}
#'     argument when fewer subjects are available.}
#'   \item{nlopt_missing_algorithms}{Character vector of NLopt algorithm names
#'     that are not available in the current build (empty on most platforms).
#'     Informational; the engine falls back automatically.}
#'   \item{neural_networks}{Named list of sub-lists, one per \code{[covariate_nn]}
#'     block declared in the model file. Empty list when the \code{nn} cargo
#'     feature is off or no NN blocks are present. Each sub-list contains:
#'     \code{name} (block name), \code{shape} (integer vector of layer widths),
#'     \code{hidden_activation}, \code{output_activation}, \code{n_weights}
#'     (total weight + bias count), \code{weights_offset} (0-based start index
#'     into \code{fit$theta}), \code{input_names}, \code{output_names}, and
#'     \code{weights} (numeric vector of trained values).}
#'   \item{exclusions}{Named list summarising records excluded by
#'     \code{[data_selection]} rules (model file and/or \code{ignore}/
#'     \code{accept}/\code{ignore_ids} arguments). Fields:
#'     \code{n_records_total}, \code{n_obs_excluded}, \code{n_dose_excluded},
#'     \code{n_other_excluded}, \code{excluded_subject_ids} (character),
#'     \code{fired_ignore} (character), \code{fired_accept} (character).
#'     \code{NULL} when no data-selection rules were active.}
#'
#' @section Process noise (SDE / diffusion):
#'
#' ODE-based PK/PD models occasionally produce autocorrelated IWRES when the
#' structural model is misspecified (missing compartment, unmodelled feedback).
#' Adding continuous within-subject process noise via an SDE framework - solved
#' by an Extended Kalman Filter (EKF) - absorbs this drift and yields a
#' better-calibrated likelihood.
#'
#' Add a \code{[diffusion]} block to the \code{.ferx} model file to enable SDE
#' mode. Each line declares the diffusion variance for one ODE state:
#' \preformatted{
#' [diffusion]
#'   central ~ 0.5   # initial estimate for DIFF_CENTRAL (variance units)
#' }
#'
#' \strong{Key points:}
#' \itemize{
#'   \item The declared value and the fitted estimate are \emph{variances}, not
#'         standard deviations. A value of 0.5 means \eqn{\sigma^2 = 0.5}, not
#'         \eqn{\sigma = 0.5}.
#'   \item Each diffusion parameter appears in \code{fit$theta} as
#'         \code{DIFF_<STATE>} (e.g. \code{DIFF_CENTRAL}). Standard errors and
#'         \code{fit$estimates} treat them as regular thetas.
#'   \item SDE models use finite differences for the EKF covariance
#'         propagation.
#'   \item SAEM is not supported with SDE models; a hard error is raised.
#'   \item \code{fit$uses_sde} is \code{TRUE} whenever a \code{[diffusion]}
#'         block was present.
#' }
#'
#' @section Unit conversion and scaling (\code{[scaling]} block):
#'
#' Add an optional \code{[scaling]} block to the \code{.ferx} model file to
#' convert model predictions before comparing them to the observations in DV.
#' This is useful when the model is parameterised in one unit (e.g. amounts in
#' mg) but DV is recorded in another (e.g. concentrations in ng/mL).
#'
#' The engine divides each prediction by the scale factor before computing
#' residuals, so the fitted sigma is on the DV scale.
#'
#' \strong{Form A -- scalar divisor:}
#' \preformatted{
#' [scaling]
#'   obs_scale = 1000   # divide every predicted value by 1000
#' }
#'
#' \strong{Form B -- expression divisor:}
#' \preformatted{
#' [scaling]
#'   obs_scale = V      # divide by the individual volume (theta/eta expression)
#' }
#' Expressions may reference thetas, etas, and \code{[individual_parameters]}
#' variables. Differentiated \strong{exactly} on the default
#' \code{gradient = auto} -- both the outer and inner loops serve an expression
#' scale analytically via the eta-quotient rule (ferx-core #486). Do not set
#' \code{gradient = fd} for it.
#'
#' \strong{Per-CMT scaling (Form A or B per compartment):}
#' \preformatted{
#' [scaling]
#'   obs_scale[CMT=1] = 1000   # CMT 1 in ng/mL, dose in ug
#'   obs_scale[CMT=2] = 1      # CMT 2 already in correct units
#' }
#'
#' \strong{Form C -- ODE readout expression:}
#'
#' For ODE models, \code{y = <expr>} defines the observation as a function of
#' ODE state variables, thetas, etas, and individual parameters. This replaces
#' the \code{obs_cmt=} argument in \code{[structural_model]} and is required
#' when the observation is not a raw compartment amount:
#' \preformatted{
#' [structural_model]
#'   ode(states = [depot, central], ...)   # no obs_cmt= here
#'
#' [scaling]
#'   y = central / V   # observation = central compartment / volume
#' }
#' Per-CMT Form C: \code{y[CMT=1] = central / V1}, \code{y[CMT=2] = central2 / V2}.
#'
#' \strong{Constraints:}
#' \itemize{
#'   \item \code{[scaling]} is not yet supported on SDE (\code{[diffusion]}) models.
#'   \item A \emph{uniform} Form B or Form C scale is differentiated exactly
#'         under the default \code{gradient = auto}. The per-CMT variants
#'         (\code{obs_scale[CMT=N]}, \code{y[CMT=N]}) fall back to finite
#'         differences \strong{silently} -- there is no error and no warning, so
#'         check \code{fit$gradient_used} if it matters.
#' }
#'
#' @section Parameter declaration syntax (\code{[parameters]} block):
#'
#' \strong{Theta (fixed effects):}
#' \preformatted{
#'   theta CL(0.134, 0.001, 10.0)          # (initial, lower, upper)
#'   theta CL(0.1, 0.001)                  # lower-bound only (no upper bound)
#'   theta CL(0.134, 0.001, 10.0) (FIX)    # FIX at end (traditional)
#'   theta CL (FIX) (0.134, 0.001, 10.0)   # FIX anywhere (flexible placement)
#' }
#'
#' \strong{Omega (inter-individual variability):}
#' \preformatted{
#'   omega ETA_CL ~ 0.07           # variance parameterisation (default)
#'   omega ETA_CL ~ 0.07 (FIX)    # fixed omega, FIX at end
#'   omega ETA_CL (FIX) ~ 0.07    # fixed omega, FIX before the tilde
#' }
#'
#' The same flexible \code{(FIX)} placement applies to \code{sigma} and
#' \code{kappa} (IOV) declarations.
#'
#' \strong{Unused-parameter warning:} a \code{warning} severity message with
#' category \code{"unused_parameter"} is emitted by the parser when a
#' parameter is declared in \code{[parameters]} but never referenced in
#' \code{[individual_parameters]} or \code{[error_model]}. This usually
#' indicates a commented-out expression or a typo in the parameter name.
#' Inspect \code{ferx_get_warnings(fit)} or check \code{ferx_model_validate()}
#' before fitting.
#'
#' @section Steady-state dosing (\code{SS} and \code{II} columns):
#'
#' Mark a dosing row as a steady-state dose by setting \code{SS = 1} (or
#' \code{SS = 2}) in the NONMEM CSV and providing the dosing interval
#' \code{II} (in the same time units as TIME).
#'
#' \itemize{
#'   \item \strong{\code{SS = 1}}: the subject is assumed to be at steady
#'         state before this dose. The engine computes the steady-state
#'         initial conditions analytically (for 1/2/3-cpt models) or via
#'         pulse expansion (for ODE models) and uses them as the starting
#'         compartment state.
#'   \item \strong{\code{SS = 2}}: the steady-state concentration is added
#'         to the current compartment state (superposition).
#'   \item \strong{\code{II}}: dosing interval (required when \code{SS > 0}).
#'         Must match the units of TIME. Rows with \code{SS = 0} or missing
#'         \code{SS} ignore \code{II}.
#' }
#'
#' No changes to the \code{.ferx} model file are needed. Steady-state
#' handling is purely data-driven and works for all PK model families
#' (1-cpt, 2-cpt, 3-cpt oral/IV/infusion, ODE-based).
#'
#' @section Specifying model and data:
#'
#' There are two equivalent ways to supply the model and data:
#'
#' \strong{1. Path strings (classic style):}
#' \preformatted{
#' ferx_fit("pk.ferx", "data.csv")
#' ferx_fit(model = "pk.ferx", data = "data.csv")
#' }
#'
#' \strong{2. \code{ferx_model} object (pipe style):}
#' \preformatted{
#' "data.csv" |> ferx_model("pk.ferx") |> ferx_fit()
#' }
#'
#' Both dispatch to the same Rust backend. The \code{ferx_model} form is
#' convenient when combined with \code{\link{ferx_model_set_section}} to modify
#' model options in the same chain (see Examples). The data path stored in the
#' \code{ferx_model} object can always be overridden by supplying \code{data}
#' explicitly to \code{ferx_fit()}.
#'
#' \strong{3. Data declared in the model file:} a \code{.ferx} file may carry a
#' \code{[data]} block with \code{path = <csv>} (resolved relative to the model
#' file's directory). When \code{data} is omitted, \code{ferx_fit()} falls back
#' to that path, so \code{ferx_fit("pk.ferx")} fits against the declared dataset.
#' An explicit \code{data} argument always wins.
#'
#' @section Controlling estimation:
#'
#' \strong{Estimation method} - pass a single method or a vector to chain
#' methods in sequence (each stage seeds the next with its converged
#' parameters):
#' \preformatted{
#' ferx_fit(m, d, method = "focei")
#' ferx_fit(m, d, method = c("saem", "focei"))  # SAEM warm-start, FOCEI polish
#' }
#'
#' \strong{Standard errors} - the covariance step is on by default:
#' \preformatted{
#' ferx_fit(m, d, covariance = TRUE)   # default: produces SE / %RSE
#' ferx_fit(m, d, covariance = FALSE)  # skip for speed during development
#'
#' # Tune the FD step for the Hessian (default 0.01):
#' ferx_fit(m, d, fd_hessian_step = 0.1)    # increase when Hessian is ill-conditioned
#' ferx_fit(m, d, fd_hessian_step = 1e-3)   # decrease for very smooth surfaces
#' }
#'
#' \strong{Parallelism} - cap the Rust thread pool:
#' \preformatted{
#' ferx_fit(m, d, threads = parallel::detectCores(logical = FALSE))
#' }
#'
#' \strong{LOQ censoring}:
#' \preformatted{
#' ferx_fit(m, d, bloq_method = "m3")    # M3 likelihood (CENS: 1=LLOQ, -1=ULOQ)
#' ferx_fit(m, d, bloq_method = "drop")  # disable M3; use rows as observed
#' }
#'
#' \strong{Gradient method} for the inner EBE loop:
#' \preformatted{
#' ferx_fit(m, d, gradient = "auto")  # default: engine-selected
#' ferx_fit(m, d, gradient = "fd")    # force finite differences
#' }
#'
#' \strong{Optimizer trace} - write a per-iteration CSV for convergence
#' diagnostics:
#' \preformatted{
#' fit <- ferx_fit(m, d, optimizer_trace = TRUE)
#' plot(fit)
#' }
#'
#' @section Fine-tuning with \code{settings}:
#'
#' The \code{settings} argument forwards a named list of low-level options to
#' the Rust backend without requiring a new wrapper release. Arguments with
#' dedicated parameters (e.g. \code{method}, \code{covariance}) cannot be
#' duplicated in \code{settings} - pass them via the named argument.
#'
#' For the complete settings reference - including which options apply to which
#' method, valid combinations, and a convergence troubleshooting guide - see
#' \href{https://ferx-nlme.github.io/model-dsl/fit-options.html}{ferx-nlme.github.io/model-dsl/fit-options}.
#'
#' \strong{Outer optimizer selection:}
#' \preformatted{
#' ferx_fit(m, d, settings = list(optimizer = "auto"))        # default
#' ferx_fit(m, d, settings = list(optimizer = "bobyqa"))      # derivative-free
#' ferx_fit(m, d, settings = list(optimizer = "slsqp"))       # gradient-based
#' ferx_fit(m, d, settings = list(optimizer = "lbfgs"))
#' ferx_fit(m, d, settings = list(optimizer = "bfgs"))
#' ferx_fit(m, d, settings = list(optimizer = "trust_region"))
#' ferx_fit(m, d, settings = list(optimizer = "mma"))
#' }
#'
#' \strong{Iteration cap and inner loop:}
#' \preformatted{
#' ferx_fit(m, d, settings = list(
#'   maxiter       = 500L,
#'   inner_maxiter = 100L,
#'   inner_tol     = 1e-6
#' ))
#' }
#'
#' \strong{SAEM tuning:}
#' \preformatted{
#' ferx_fit(m, d, method = "saem", settings = list(
#'   n_exploration = 200,
#'   n_convergence = 400,
#'   n_mh_steps    = 3,
#'   omega_burnin  = 20L,
#'   seed          = 42L
#' ))
#'
#' # HMC proposals:
#' ferx_fit(m, d, method = "saem", settings = list(n_leapfrog = 3L))
#' }
#'
#' \strong{SIR uncertainty (requires \code{sir = TRUE, covariance = TRUE}):}
#' \preformatted{
#' ferx_fit(m, d, sir = TRUE, settings = list(
#'   sir_samples   = 2000L,
#'   sir_resamples = 500L,
#'   sir_seed      = 1L
#' ))
#' }
#'
#' \strong{Global search before local refinement:}
#' \preformatted{
#' ferx_fit(m, d, settings = list(
#'   global_search  = TRUE,
#'   global_maxeval = 1000L
#' ))
#' }
#'
#' \strong{Trust-region CG budget:}
#' \preformatted{
#' ferx_fit(m, d, settings = list(
#'   optimizer          = "trust_region",
#'   steihaug_max_iters = 100L
#' ))
#' }
#'
#' \strong{Convergence tolerance for difficult subjects:}
#'
#' \code{max_unconverged_frac} (numeric in \[0, 1\], default \code{0.1}) is the
#' fraction of subjects allowed to have unconverged EBEs before the outer
#' optimizer \emph{rejects the trial step} (scoring it as \code{Inf}). It is an
#' in-loop guard on the search, not a post-hoc relaxation of the converged flag,
#' and it never raises an error; set it to \code{1.0} to disable the guard
#' entirely. \code{min_obs_for_convergence_check} (non-negative integer) is the
#' minimum number of observations a subject must have before its inner-loop
#' convergence counts toward that fraction; sparser subjects are excluded.
#' \preformatted{
#' ferx_fit(m, d, settings = list(
#'   max_unconverged_frac          = 0.1,
#'   min_obs_for_convergence_check = 2L
#' ))
#' }
#'
#' \strong{Stagnation guard (NLopt-based optimizers):}
#'
#' NLopt-based outer optimizers (BOBYQA, SLSQP, L-BFGS, MMA) short-circuit by
#' default once recent evaluations show no OFV improvement above 1e-3 over a
#' rolling window, letting them terminate quickly instead of burning through
#' \code{maxiter} at full inner-loop cost. Set \code{stagnation_guard = FALSE}
#' to disable this and let the optimizer run to its natural termination
#' criterion - useful when debugging or for problems with genuinely
#' slow-but-real OFV improvements.
#' \preformatted{
#' ferx_fit(m, d, settings = list(stagnation_guard = FALSE))
#' }
#'
#' \strong{Multi-start optimization:}
#'
#' Runs \code{n_starts} independent fits from perturbed initial values in
#' parallel (via rayon) and returns the converged result with the lowest OFV.
#' Default \code{n_starts = 1} disables multi-start (single run, no overhead).
#' Useful for models prone to local minima: Michaelis-Menten elimination,
#' full-block omega, or many covariate parameters. On an 8-core machine
#' \code{n_starts = 8} costs approximately the same wall-clock time as a single run.
#'
#' Start 0 always uses the exact initial values from the model file.
#' Starts 1..n apply a log-space perturbation of size \code{start_sigma}
#' (default 0.3) to each theta.
#'
#' \preformatted{
#' # 8 parallel starts (approx. same wall time as one run on 8 cores)
#' ferx_fit(m, d, settings = list(n_starts = 8L))
#'
#' # Wider perturbation for ridge-shaped surfaces (e.g. Michaelis-Menten):
#' ferx_fit(m, d, settings = list(n_starts = 8L, start_sigma = 0.5))
#'
#' # Pin the perturbation RNG seed for reproducibility (independent of the SAEM seed):
#' ferx_fit(m, d, settings = list(n_starts = 8L, multi_start_seed = 123L))
#' }
#'
#' @section Post-fit outputs and pipe extensions:
#'
#' \code{ferx_fit()} returns a \code{ferx_fit} object. All of the following
#' work both as standalone calls and as the tail of a \code{|>} pipe:
#'
#' \preformatted{
#' fit |> print()              # full parameter table
#' fit |> summary()            # compact diagnostic summary
#' fit |> ferx_model_inspect() # model structure auto-derived by the engine
#' fit |> plot()               # convergence trace (needs optimizer_trace = TRUE)
#' }
#'
#' Diagnostics data frame (PRED, IPRED, CWRES, ETAs, ...) lives in
#' \code{fit$sdtab}; tidy estimates, correlations, and ETA-covariate
#' correlations are plain fields (no accessor call needed):
#' \preformatted{
#' fit$sdtab
#' fit$ebe_etas
#' fit$individual_estimates
#' fit$estimates    # tidy data frame: theta / omega / sigma + SE / %RSE
#' fit$cor_matrix   # parameter correlation matrix (needs covariance = TRUE)
#' fit$eta_cov      # ETA-covariate correlation table
#' }
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' fit$theta
#' fit$ofv
#'
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # -- Classic style (path strings) ------------------------------------------
#'
#' # Minimal call: FOCEI with default optimizer, covariance step on
#' fit <- ferx_fit(ex$model, ex$data)
#'
#' # Named arguments (fully equivalent)
#' fit <- ferx_fit(model = ex$model, data = ex$data, method = "focei")
#'
#' # All common options at once
#' fit <- ferx_fit(ex$model, ex$data,
#'   method      = "focei",
#'   covariance  = TRUE,
#'   threads     = 4L,
#'   gradient    = "auto",
#'   verbose     = TRUE,
#'   scale_params = FALSE
#' )
#'
#' # SAEM warm-start, then FOCEI polish
#' fit <- ferx_fit(ex$model, ex$data, method = c("saem", "focei"))
#'
#' # Fine-tune via settings
#' fit <- ferx_fit(ex$model, ex$data,
#'   method   = "focei",
#'   settings = list(
#'     optimizer     = "slsqp",
#'     maxiter       = 500L,
#'     inner_maxiter = 100L,
#'     inner_tol     = 1e-6
#'   )
#' )
#'
#' # LOQ-censored observations (M3 method)
#' bloq <- ferx_example("warfarin_bloq")
#' fit  <- ferx_fit(bloq$model, bloq$data, method = "focei", bloq_method = "m3")
#'
#' # SIR parameter uncertainty
#' fit <- ferx_fit(ex$model, ex$data,
#'   sir      = TRUE,
#'   settings = list(sir_samples = 2000L, sir_resamples = 500L)
#' )
#'
#' # -- Pipe style (ferx_model object) -----------------------------------------
#'
#' # Basic pipe: data flows into ferx_model() which bundles the model path
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_fit(method = "focei", covariance = TRUE)
#'
#' # Modify a model section, then fit
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_model_set_section("fit_options", c(
#'     "  method     = focei",
#'     "  maxiter    = 500",
#'     "  covariance = true"
#'   )) |>
#'   ferx_fit()
#'
#' # Full pipeline including post-fit outputs
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_model_set_section("fit_options", c("  method = focei")) |>
#'   ferx_fit(covariance = TRUE, threads = 4L,
#'            settings = list(optimizer = "slsqp")) |>
#'   summary()
#'
#' # Inspect a section, then fit
#' ferx_model_get_section(ex$model, "parameters")
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_fit() |>
#'   (\(fit) fit$estimates)()
#'
#' # Override data stored in ferx_model at fit time
#' # (substitute the path to your own dataset for "other_cohort.csv")
#' m <- ferx_model(ex$data, ex$model)
#' ferx_fit(m, data = "other_cohort.csv", method = "focei")
#'
#' # -- Post-fit outputs -------------------------------------------------------
#'
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE, optimizer_trace = TRUE)
#'
#' summary(fit)              # compact diagnostic table
#' ferx_model_inspect(fit)   # model structure auto-derived by the engine
#' plot(fit)                 # convergence plot (optimizer_trace = TRUE required)
#'
#' fit$sdtab                 # per-observation diagnostics (PRED, IPRED, CWRES, etc.)
#' fit$ebe_etas              # per-subject empirical Bayes ETAs
#' fit$individual_estimates  # per-subject individual PK parameters
#' fit$estimates             # tidy data frame with SE and %RSE
#' fit$cor_matrix            # parameter correlation matrix
#' fit$eta_cov               # ETA-covariate correlation table
#' fit$eigenvalues           # sorted eigenvalues of parameter correlation matrix
#' fit$condition_number      # > 1000 flags potential ill-conditioning
#'
#' # -- SDE model (Extended Kalman Filter) -------------------------------------
#'
#' # The [diffusion] block in the .ferx file enables SDE mode.
#' # DIFF_CENTRAL is a within-subject process-noise variance on the central
#' # compartment.  SDE models use finite differences for this step.
#' sde <- ferx_example("warfarin_sde")
#' fit_sde <- ferx_fit(sde$model, sde$data)
#' fit_sde$uses_sde                  # TRUE
#' fit_sde$theta["DIFF_CENTRAL"]     # fitted diffusion variance
#' fit_sde$estimates                 # DIFF_CENTRAL appears in theta block
#'
#' # -- Multi-start (avoid local minima) ----------------------------------------
#' ex <- ferx_example("warfarin")
#' fit_ms <- ferx_fit(ex$model, ex$data, settings = list(n_starts = 4L))
#' fit_ms$ofv
#'
#' # -- NCA-based starting values (inits_from_nca) ------------------------------
#' # `inits_from_nca` is designed to mitigate stalling and local-minimum traps
#' # in **gradient-based estimation methods**, which navigate the likelihood
#' # surface using a local Jacobian / Hessian and so are sensitive to where the
#' # starting thetas sit. In ferx those are:
#' #   - method = "gn"         (FOCE Gauss-Newton)
#' #   - method = "gn_hybrid"  (FOCE-GN with SLSQP polish)
#' #   - method = "foce" / "focei" run with settings = list(optimizer =
#' #     "trust_region") or settings = list(optimizer = "slsqp") (the
#' #     default "auto" resolves to derivative-free "bobyqa" on these
#' #     poorly-started ODE/PD fits, which tolerates poor starts the best;
#' #     switch to slsqp / trust_region for gradient-based behaviour, then
#' #     NCA-derived starts help most)
#' # For these, NCA-derived starts often turn a non-converging or stagnating
#' # fit into a clean convergence.
#' #
#' # Stochastic / sampling-based methods ("saem", "imp") explore the space
#' # globally and are far less sensitive to starting values, so `inits_from_nca`
#' # gives a smaller benefit. The better tool for them is multi-start (run
#' # several fits from perturbed starts and keep the lowest OFV) - see the
#' # Multi-start section above (`settings = list(n_starts = 4L)`), or chain
#' # SAEM into FOCEI: `method = c("saem", "focei")`.
#'
#' # TRUE = default strategy ("nca_sweep"); rescues a trust_region fit that
#' # would otherwise stall on poor defaults.
#' fit <- ferx_fit(ex$model, ex$data,
#'   settings       = list(optimizer = "trust_region"),
#'   inits_from_nca = TRUE
#' )
#'
#' # Pick a strategy explicitly. "nca_ebe" handles large between-subject
#' # variability better than "nca_sweep" but falls back to it for ODE models.
#' fit <- ferx_fit(ex$model, ex$data, method = "gn",
#'                 inits_from_nca = "nca_ebe")
#'
#' # Fastest variant: NCA arithmetic only, no grid sweep.
#' fit <- ferx_fit(ex$model, ex$data, inits_from_nca = "nca")
#'
#' # Same flag from the model file's [fit_options] (the call-time arg wins
#' # and warns on conflict). Inspect what NCA produces without fitting via
#' # `ferx_inits_from_nca()` -- see ?ferx_inits_from_nca.
#' #
#' # Combine with multi-start when even good inits aren't enough - good starts
#' # PLUS a small perturbed ensemble is often the most robust setup.
#' fit <- ferx_fit(ex$model, ex$data,
#'   method         = "gn",
#'   inits_from_nca = "nca_sweep",
#'   settings       = list(n_starts = 4L)
#' )
#'
#' # -- Deep Compartment Model (covariate neural network) -----------------------
#' # Requires ferx-r built with the `nn` cargo feature.
#' # The [covariate_nn TYPICAL_PK] block replaces analytical covariate functions
#' # with a small MLP: WT + CRCL -> multiplicative modulator on TVCL/TVV1/etc.
#' # ferx_fit() call is identical to any other model.
#' dcm <- ferx_example("warfarin_dcm")
#' fit_dcm <- ferx_fit(dcm$model, dcm$data, method = "focei")
#'
#' # NN metadata is in fit$neural_networks (one sub-list per [covariate_nn] block)
#' fit_dcm$neural_networks[[1]]$name          # "TYPICAL_PK"
#' fit_dcm$neural_networks[[1]]$shape         # e.g. c(2, 8, 8, 5)
#' fit_dcm$neural_networks[[1]]$n_weights     # total weight + bias count
#' fit_dcm$neural_networks[[1]]$input_names   # c("WT", "CRCL")
#' fit_dcm$neural_networks[[1]]$output_names  # c("CL", "V1", "Q", "V2", "KA")
#'
#' # NN weight thetas are suppressed from the THETA table in print();
#' # a NEURAL NETWORKS summary block is shown instead.
#' print(fit_dcm)
#'
#' # See inst/examples/ex_warfarin_dcm.R for a full worked example including
#' # interpretability heuristics and comparison against a no-covariate baseline.
#' }
#'
#' @param ... Reserved for future use. Unrecognised arguments raise an error.
#' @seealso
#'   \itemize{
#'     \item \href{https://ferx-nlme.github.io/model-dsl/fit-options.html}{Fit
#'       options reference} - full documentation of every \code{[fit_options]}
#'       key and \code{settings} knob, organised by method with a
#'       convergence-troubleshooting guide and a method x option compatibility
#'       table.
#'     \item \code{\link{ferx_fit_async}} - non-blocking version for long runs.
#'     \item \code{\link{ferx_check_init}} - 5-iteration pilot to validate
#'       starting values before a full run.
#'     \item \code{\link{ferx_inits_from_nca}} - NCA-derived starting values.
#'     \item \code{\link{ferx_get_warnings}} - structured warnings from the fit.
#'     \item \code{fit$estimates} - tidy parameter table with SE / \%RSE.
#'     \item \code{\link{plot.ferx_fit}} - convergence trace plot.
#'   }
#' @family fitting
#' @export
ferx_fit <- function(model, data = NULL,
                     method = NULL,
                     covariance = NULL,
                     verbose = NULL,
                     bloq_method = NULL,
                     threads = NULL,
                     mu_referencing = NULL,
                     sir = NULL,
                     gradient = NULL,
                     optimizer_trace = FALSE,
                     scale_params = FALSE,
                     inits_from_nca = FALSE,
                     fd_hessian_step = NULL,
                     settings = NULL,
                     output = NULL,
                     include_data = FALSE,
                     ignore = NULL,
                     accept = NULL,
                     ignore_ids = NULL,
                     ...) {
  # `ferx_fit()` is a plain function (not an S3 generic) so that IDE argument
  # completion and inline help surface every argument - an S3 generic only
  # exposes `model`, `data`, and `...`, hiding the real options from the user
  # (see issue #52). Dispatch on the `model` type is handled inline below:
  # `model` may be a path to a `.ferx` file or a `ferx_model` object produced
  # by `ferx_model()` (the pipe style, where `data` defaults to the path stored
  # in the object).
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    if (is.null(data)) data <- .ferx_model_data_path(model$model)
    if (is.null(data)) {
      stop(
        "No data supplied. Pass `data`, or add a `[data]` block ",
        "(`path = ...`) to the model file."
      )
    }
    model <- model$model
  }
  # Guard against silently swallowed typos: any leftover `...` argument is an
  # unrecognised name, not a value forwarded anywhere. Use `...length()` /
  # `...names()` rather than `list(...)` so the promises are never forced -
  # forcing them could trigger side effects or surface an evaluation error
  # before this clearer "unused argument(s)" message (both base R >= 4.1.0,
  # which the package already assumes via its use of the `|>` pipe).
  if (...length() > 0L) {
    nms <- ...names()
    nms <- if (is.null(nms)) rep("", ...length()) else nms
    labels <- ifelse(nzchar(nms), nms, "<unnamed>")
    stop(
      "unused argument(s) passed to ferx_fit(): ",
      paste(labels, collapse = ", ")
    )
  }
  # ferx_data support: when data is a ferx_data object produced by
  # ferx_apply_selection(), extract the source path and merge any stored conditions
  # with the explicit ignore/accept/ignore_ids args.
  if (inherits(data, "ferx_data")) {
    fd <- data
    data <- attr(fd, "source_path")
    if (is.null(data)) {
      stop(
        "`data` is a ferx_data object built from a data.frame (no source file path). ",
        "The Rust engine requires a CSV file path. ",
        "Pass the original CSV file path to ferx_apply_selection() instead:\n",
        "  sel <- ferx_apply_selection(\"path/to/data.csv\", ignore = ...)\n",
        "  fit <- ferx_fit(model, sel)"
      )
    }
    ignore      <- unique(c(attr(fd, "ignore"),    ignore))
    accept      <- unique(c(attr(fd, "accept"),    accept))
    ignore_ids  <- unique(c(attr(fd, "ignore_ids"), ignore_ids))
  }
  # `gradient = NULL` (the default) defers to the model file's
  # `[fit_options] gradient`; an explicit value overrides it. The engine reads
  # "" as "keep the model-file value" (#558).
  if (is.null(gradient)) {
    gradient_arg <- ""
  } else {
    gradient_arg <- match.arg(gradient, c("auto", "ad", "fd"))
  }
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  stopifnot(file.exists(model), file.exists(data))
  # `covariance` / `verbose` / `mu_referencing` / `sir` default to NULL, meaning
  # "use the model file's [fit_options] value" (falling back to the engine
  # default when the model file is silent). They are forwarded to the engine as
  # sentinel strings: "" keeps the model-file value, "true"/"false" override it.
  # This stops an accepted R default from silently overruling the model file -
  # the same fix as `method` (#558). Pass the argument explicitly to override.
  .flag_arg <- function(x, nm) {
    if (is.null(x)) {
      return("")
    }
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be TRUE or FALSE", nm))
    }
    if (x) "true" else "false"
  }
  covariance_arg     <- .flag_arg(covariance, "covariance")
  verbose_arg        <- .flag_arg(verbose, "verbose")
  mu_referencing_arg <- .flag_arg(mu_referencing, "mu_referencing")
  sir_arg            <- .flag_arg(sir, "sir")
  # Read the model file's [fit_options] once, up front: it informs both the SIR
  # precondition below and the conflict warnings further down.
  model_file_opts <- .ferx_parse_model_fit_options(model)
  # `sir = TRUE` needs covariance estimation. Resolve the effective covariance
  # setting with the usual precedence (explicit R arg > model file > engine
  # default, which is on) and reject the request when covariance is actually
  # disabled - whether by the `covariance` argument or by the model file - rather
  # than letting the engine silently skip SIR with only a buried warning.
  cov_effective <- if (!is.null(covariance)) {
    isTRUE(covariance)
  } else if ("covariance" %in% names(model_file_opts)) {
    !identical(tolower(as.character(model_file_opts[["covariance"]])), "false")
  } else {
    TRUE
  }
  if (isTRUE(sir) && !cov_effective) {
    stop(
      "`sir = TRUE` requires covariance estimation, but covariance is disabled ",
      if (!is.null(covariance)) "via the `covariance` argument." else
        "in the model file [fit_options].",
      " Set covariance = TRUE."
    )
  }
  # `method = NULL` (the default) means "use whatever the model file's
  # [fit_options] selected, falling back to FOCEI" rather than forcing FOCEI from
  # the R side. We send an empty character vector to the engine; the Rust glue
  # then keeps the model-file method instead of overriding it (#558). A
  # user-supplied `method` always wins. The downstream `imp` / `bayes` chain
  # guards below no-op cleanly on `character(0)`.
  if (is.null(method)) {
    method <- character(0L)
  } else {
    if (!is.character(method) || length(method) == 0L) {
      stop("`method` must be a non-empty character vector")
    }
    method <- vapply(method, .normalize_method_token, character(1L), USE.NAMES = FALSE)
  }
  # `imp` may appear at most once. By default it is an estimator (NONMEM
  # METHOD=IMP) and may sit anywhere in the chain; with `imp_eval_only = TRUE`
  # (NONMEM EONLY=1) it only evaluates the marginal -2 log L and must be the
  # terminal stage. The engine rejects malformed chains too, but surfacing the
  # error R-side avoids a round-trip and gives a clearer message. (`imp_eval_only`
  # set via the model file's [fit_options] rather than `settings` is enforced
  # engine-side.)
  imp_positions <- which(method == "imp")
  if (length(imp_positions) > 0L) {
    if (length(imp_positions) > 1L) {
      stop(
        "`method` may include \"imp\" at most once; got ",
        length(imp_positions), " occurrences"
      )
    }
    imp_eval_only <- isTRUE(settings[["imp_eval_only"]])
    if (imp_eval_only && imp_positions != length(method)) {
      stop(
        "`\"imp\"` with `imp_eval_only = TRUE` must be the final stage of the ",
        "method chain; got method = c(\"",
        paste(method, collapse = "\", \""), "\"). ",
        "Drop `imp_eval_only` to run it as an estimator anywhere in the chain."
      )
    }
  }
  # `bayes` is a standalone MCMC estimator - it does not warm-start from or feed
  # another stage, and chaining it would run the (expensive) sampler and then
  # silently discard its posterior (the final stage's result wins). Reject
  # chains R-side rather than waste the run.
  if (any(method == "bayes") && length(method) > 1L) {
    stop(
      "`\"bayes\"` must be the only method (it runs standalone); got ",
      "method = c(\"", paste(method, collapse = "\", \""), "\")"
    )
  }
  if (is.null(bloq_method)) {
    bloq_arg <- ""
  } else {
    bloq_arg <- match.arg(tolower(bloq_method), c("drop", "m3"))
  }
  if (is.null(threads)) {
    threads_arg <- 0L
  } else {
    if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads != as.integer(threads) || threads < 0L) {
      stop("`threads` must be NULL or a non-negative integer scalar")
    }
    threads_arg <- as.integer(threads)
  }

  if (!is.logical(optimizer_trace) || length(optimizer_trace) != 1L || is.na(optimizer_trace)) {
    stop("`optimizer_trace` must be TRUE or FALSE")
  }
  if (!is.logical(scale_params) || length(scale_params) != 1L || is.na(scale_params)) {
    stop("`scale_params` must be TRUE or FALSE")
  }
  # Merge optimizer_trace into settings so apply_fit_option handles it on the
  # Rust side (it's already in framework_keys()).  User-supplied settings take
  # precedence if somehow duplicated, which Rust will reject as a duplicate key
  # - but we forbid that below via the RESERVED list.
  if (isTRUE(optimizer_trace)) {
    settings <- c(list(optimizer_trace = TRUE), settings)
  }
  # Only inject scale_params when TRUE - the Rust default is FALSE (off, since
  # ferx-core issue #99), so omitting it preserves the behaviour callers expect
  # when they don't pass the argument.
  if (isTRUE(scale_params)) {
    settings <- c(list(scale_params = TRUE), settings)
  }
  # inits_from_nca: TRUE/FALSE or one of "nca", "nca_sweep", "nca_ebe". TRUE is
  # an alias for "nca_sweep"; FALSE disables (the default). Merged into settings
  # so apply_fit_option handles it on the Rust side (it's in framework_keys()).
  # The dedicated arg always wins over a `settings` duplicate - including
  # FALSE, which strips any settings value rather than silently letting it
  # through.
  inits_from_nca_explicit <- "inits_from_nca" %in% names(match.call())[-1L]
  inits_value <- NULL
  if (is.logical(inits_from_nca)) {
    if (length(inits_from_nca) != 1L || is.na(inits_from_nca)) {
      stop("`inits_from_nca` must be TRUE/FALSE or one of \"nca\", \"nca_sweep\", \"nca_ebe\"")
    }
    if (isTRUE(inits_from_nca)) inits_value <- "nca_sweep"
  } else if (is.character(inits_from_nca) && length(inits_from_nca) == 1L && !is.na(inits_from_nca)) {
    inits_value <- match.arg(tolower(inits_from_nca), c("nca", "nca_sweep", "nca_ebe"))
  } else {
    stop("`inits_from_nca` must be TRUE/FALSE or one of \"nca\", \"nca_sweep\", \"nca_ebe\"")
  }
  if (inits_from_nca_explicit) {
    settings[["inits_from_nca"]] <- NULL
  }
  if (!is.null(inits_value)) {
    settings <- c(list(inits_from_nca = inits_value), settings)
  }
  # fd_hessian_step: dedicated arg wins over any `settings` duplicate.
  if (!is.null(fd_hessian_step)) {
    if (!is.numeric(fd_hessian_step) || length(fd_hessian_step) != 1L ||
        !is.finite(fd_hessian_step) || fd_hessian_step <= 0) {
      stop("`fd_hessian_step` must be a positive finite numeric scalar (e.g. 0.01)")
    }
    settings[["fd_hessian_step"]] <- NULL
    settings <- c(list(fd_hessian_step = fd_hessian_step), settings)
  }
  # Validate and normalise data-selection args.
  if (!is.null(ignore)) {
    if (!is.character(ignore)) stop("`ignore` must be a character vector of filter expressions")
    ignore <- as.character(ignore)
  }
  if (!is.null(accept)) {
    if (!is.character(accept)) stop("`accept` must be a character vector of filter expressions")
    accept <- as.character(accept)
  }
  if (!is.null(ignore_ids)) {
    if (!is.numeric(ignore_ids) && !is.character(ignore_ids))
      stop("`ignore_ids` must be a numeric or character vector of subject IDs")
    ignore_ids <- as.character(ignore_ids)
  }

  settings_parts <- .ferx_settings_to_strings(settings)
  # Append data-selection conditions after the helper (which enforces unique
  # names - selection keys are repeatable so we add them directly as parallel
  # vectors). Each element of `ignore` / `accept` becomes its own key entry;
  # Rust's apply_fit_option deduplicates exact-string repeats against any
  # conditions already declared in the [data_selection] model block.
  if (length(ignore) > 0L) {
    settings_parts$keys   <- c(settings_parts$keys,   rep("ignore", length(ignore)))
    settings_parts$values <- c(settings_parts$values, ignore)
  }
  if (length(accept) > 0L) {
    settings_parts$keys   <- c(settings_parts$keys,   rep("accept", length(accept)))
    settings_parts$values <- c(settings_parts$values, accept)
  }
  if (length(ignore_ids) > 0L) {
    settings_parts$keys   <- c(settings_parts$keys,   rep("ignore_subjects", length(ignore_ids)))
    settings_parts$values <- c(settings_parts$values, ignore_ids)
  }
  # Effective settings as sent to Rust, with merged defaults, for summary display.
  settings_used <- if (length(settings_parts$keys) > 0L) {
    setNames(lapply(settings_parts$values, .ferx_parse_setting_value), settings_parts$keys)
  } else {
    list()
  }

  # Read [fit_options] from the model file and detect conflicts with R call-time
  # arguments. Precedence: dedicated R args > settings list > model file. We
  # warn (not error) so the user is informed which value actually runs. Only
  # args the caller *explicitly* passed are flagged - accepting a default is
  # not an override of the model file. (`model_file_opts` was parsed above.)
  explicit_args <- names(match.call())[-1L]
  dedicated_explicit <- list(
    method         = if ("method"         %in% explicit_args) paste(method, collapse = ", "),
    covariance     = if ("covariance"     %in% explicit_args && nzchar(covariance_arg)) covariance_arg,
    bloq_method    = if ("bloq_method"    %in% explicit_args && !is.null(bloq_method)) bloq_arg,
    threads        = if ("threads"        %in% explicit_args && !is.null(threads)) as.character(threads),
    mu_referencing = if ("mu_referencing" %in% explicit_args && nzchar(mu_referencing_arg)) mu_referencing_arg,
    sir            = if ("sir"            %in% explicit_args && nzchar(sir_arg)) sir_arg,
    gradient       = if ("gradient"       %in% explicit_args && nzchar(gradient_arg)) gradient_arg,
    inits_from_nca  = if ("inits_from_nca"  %in% explicit_args) (if (is.null(inits_value)) "off" else inits_value),
    fd_hessian_step = if ("fd_hessian_step" %in% explicit_args && !is.null(fd_hessian_step))
      format(fd_hessian_step, scientific = FALSE, trim = TRUE, digits = 17)
  )
  .ferx_warn_fit_option_conflicts(model_file_opts, dedicated_explicit, settings_parts)

  fit_started_at <- Sys.time()
  raw <- ferx_rust_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    method = method,
    covariance = covariance_arg,
    verbose = verbose_arg,
    bloq_method = bloq_arg,
    threads = threads_arg,
    mu_referencing = mu_referencing_arg,
    sir = sir_arg,
    gradient = gradient_arg,
    settings_keys = settings_parts$keys,
    settings_values = settings_parts$values
  )

  if (length(raw) == 0) {
    stop("Model fitting failed. Check model and data files.")
  }

  # Structure the result
  result <- raw

  # Name the theta vector
  if (length(result$theta) > 0 && length(result$theta_names) > 0) {
    names(result$theta) <- result$theta_names
  }

  # Reshape omega into a matrix. A model with no random effects (n_eta = 0 - e.g.
  # a fixed-effects `[binary_model]`) returns an empty omega; keep it a 0x0 matrix
  # rather than a bare numeric(0), so downstream `nrow(result$omega)` is 0, not
  # NULL (the latter poisons the eta-metadata guards below with NA - #271).
  if (length(result$omega) > 0 && !is.null(result$omega_dim)) {
    d <- result$omega_dim
    result$omega <- matrix(result$omega, nrow = d, ncol = d)
    eta_nms_om <- result$eta_names
    omega_dim_nms <- if (!is.null(eta_nms_om) && length(eta_nms_om) == d) eta_nms_om else paste0("OMEGA(", seq_len(d), ",", seq_len(d), ")")
    rownames(result$omega) <- colnames(result$omega) <- omega_dim_nms
  } else {
    result$omega <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  }

  # Name SE vectors
  if (length(result$se_theta) > 0 && length(result$theta_names) > 0) {
    names(result$se_theta) <- result$theta_names
  }
  if (length(result$se_theta) == 0) result$se_theta <- NULL
  if (length(result$se_omega) == 0) result$se_omega <- NULL
  if (length(result$se_sigma) == 0) result$se_sigma <- NULL

  # SIR: NaN ess => not computed; flat [lo, hi, lo, hi, ...] => (n, 2) matrix
  if (is.null(result$sir_ess) || !is.finite(result$sir_ess)) {
    result$sir_ess <- NULL
  }
  reshape_ci <- function(v, row_names = NULL) {
    if (length(v) == 0) {
      return(NULL)
    }
    m <- matrix(v,
      ncol = 2, byrow = TRUE,
      dimnames = list(row_names, c("lower", "upper"))
    )
    m
  }
  result$sir_ci_theta <- reshape_ci(result$sir_ci_theta, result$theta_names)
  n_eta <- if (is.null(dim(result$omega))) NULL else nrow(result$omega)
  eta_names <- if (!is.null(n_eta)) {
    en <- result$eta_names
    if (!is.null(en) && length(en) == n_eta) en else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
  } else NULL
  result$sir_ci_omega <- reshape_ci(result$sir_ci_omega, eta_names)
  sn <- result$sigma_names
  sig_names <- if (!is.null(sn) && length(sn) == length(result$sigma)) sn else paste0("SIGMA(", seq_along(result$sigma), ")")
  result$sir_ci_sigma <- reshape_ci(result$sir_ci_sigma, sig_names)

  # Normalize trace_path: NULL/empty means no trace was written
  tp <- result$trace_path
  if (is.null(tp) || length(tp) == 0L || !nzchar(tp[[1L]])) {
    result$trace_path <- NULL
    result$trace <- NULL
  } else {
    .ferx_state$last_trace_path <- result$trace_path
    .ferx_state$last_trace_time <- fit_started_at
    .ferx_state$last_trace_model <- model
    # Read the trace CSV into the fit object itself, not just its path, so
    # it survives a save/load round-trip and callers don't need to re-read
    # a temp file that may since have been deleted (see ferx_trace()).
    result$trace <- tryCatch(.ferx_read_trace_csv(result$trace_path),
                              error = function(e) NULL)
  }

  # Normalize shrinkage: NaN ? NA (consistent with other optional numerics)
  if (!is.null(result$shrinkage_eta)) {
    result$shrinkage_eta[!is.finite(result$shrinkage_eta)] <- NA_real_
  }
  if (!is.null(result$shrinkage_eps) && !is.finite(result$shrinkage_eps)) {
    result$shrinkage_eps <- NA_real_
  }
  if (!is.null(result$dw_statistic) && !is.finite(result$dw_statistic)) {
    result$dw_statistic <- NA_real_
  }
  if (!is.null(result$iwres_lag1_r) && !is.finite(result$iwres_lag1_r)) {
    result$iwres_lag1_r <- NA_real_
  }

  # Normalize covariance_status: missing from older binaries ? "not_requested"
  if (is.null(result$covariance_status)) {
    result$covariance_status <- "not_requested"
  }

  # Model-selection surface (ferx-core #1177). `max_abs_correlation` arrives as
  # NaN when the fit has no covariance matrix to read one off, and the two
  # tri-state verdicts arrive absent when the engine recorded none. Normalize
  # both to NA so check_strictness() meets a single "no input" shape whether
  # the fit came from here or from ferx_load_fit().
  if (!is.null(result$max_abs_correlation) &&
        !is.finite(result$max_abs_correlation)) {
    result$max_abs_correlation <- NA_real_
  }
  result$left_init <- .fitrx_unwrap_opt_lgl(result$left_init)
  result$stalled_at_init <- .fitrx_unwrap_opt_lgl(result$stalled_at_init)

  # Reshape cov_matrix into a named square matrix (param ? param)
  d <- result$cov_matrix_dim %||% 0L
  if (!is.null(result$cov_matrix) && length(result$cov_matrix) > 0L && d > 0L) {
    m <- matrix(result$cov_matrix, nrow = d, ncol = d, byrow = TRUE)
    n_theta <- length(result$theta_names)
    n_eta <- result$omega_dim %||% 0L
    n_sigma <- length(result$sigma)
    n_omega_packed <- d - n_theta - n_sigma
    # Determine parameterisation: diagonal (n_omega_packed == n_eta) or block
    eta_nms <- if (!is.null(result$eta_names) && length(result$eta_names) == n_eta) result$eta_names else NULL
    omega_names <- if (n_omega_packed == n_eta) {
      if (!is.null(eta_nms)) eta_nms
      else paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
    } else {
      # Block lower-triangle: L(i,j) for i >= j, column-major
      nm <- character(n_omega_packed)
      k <- 0L
      for (i in seq_len(n_eta)) {
        for (j in seq_len(i)) {
          k <- k + 1L
          nm[k] <- if (!is.null(eta_nms)) sprintf("%s,%s", eta_nms[i], eta_nms[j]) else sprintf("OMEGA(%d,%d)", i, j)
        }
      }
      nm
    }
    sig_nms <- if (!is.null(result$sigma_names) && length(result$sigma_names) == n_sigma) result$sigma_names else NULL
    pnames <- c(
      result$theta_names,
      if (n_omega_packed > 0L) omega_names else character(0L),
      if (n_sigma > 0L) (if (!is.null(sig_nms)) sig_nms else paste0("SIGMA(", seq_len(n_sigma), ")")) else character(0L)
    )
    if (length(pnames) == d) rownames(m) <- colnames(m) <- pnames
    result$cov_matrix <- m
  } else {
    result$cov_matrix <- NULL
  }
  result$cov_matrix_dim <- NULL

  result <- .ferx_apply_cov_sentinels(result)

  # ETA normality (Shapiro-Wilk) - computed in R from per-subject EBEs.
  # Fold every flagged ETA into a single warning rather than one per ETA
  # (ferx-core#163).
  result$eta_normality <- .ferx_compute_eta_normality(result$ebe_etas)
  nn_msg <- .ferx_eta_normality_warning(result$eta_normality)
  if (!is.null(nn_msg)) {
    result$warnings <- c(result$warnings, nn_msg)
  }

  # Reshape omega_iov into a named matrix (NULL when no IOV)
  d_iov <- result$omega_iov_dim %||% 0L
  if (!is.null(result$omega_iov) && length(result$omega_iov) > 0L && d_iov > 0L) {
    m_iov <- matrix(result$omega_iov, nrow = d_iov, ncol = d_iov)
    if (length(result$kappa_names) == d_iov) {
      rownames(m_iov) <- colnames(m_iov) <- result$kappa_names
    }
    result$omega_iov <- m_iov
    if (length(result$kappa_names) > 0L && length(result$shrinkage_kappa) == length(result$kappa_names)) {
      names(result$shrinkage_kappa) <- result$kappa_names
    }
    if (length(result$se_kappa) == 0L) {
      result$se_kappa <- NULL
    } else {
      n_tri <- d_iov * (d_iov + 1L) / 2L
      if (length(result$se_kappa) == d_iov) {
        names(result$se_kappa) <- result$kappa_names
      } else if (length(result$se_kappa) == n_tri) {
        # block kappa: label lower-triangle elements as NAME (diagonal) or COV_i_j
        tri_names <- character(n_tri)
        idx <- 1L
        for (j in seq_len(d_iov)) {
          for (i in j:d_iov) {
            tri_names[idx] <- if (i == j) {
              result$kappa_names[i]
            } else {
              paste0("COV_", result$kappa_names[j], "_", result$kappa_names[i])
            }
            idx <- idx + 1L
          }
        }
        names(result$se_kappa) <- tri_names
      }
    }
    # Per-occasion shrinkage: set column names from kappa_names; NULL when
    # the Rust glue returned an empty/NULL frame (no IOV or unbalanced design).
    if (is.data.frame(result$shrinkage_kappa_by_occ) &&
        nrow(result$shrinkage_kappa_by_occ) > 0L) {
      kn <- result$kappa_names
      # Columns are: occ, <kappa1>, <kappa2>, ... already named by Rust glue
      # but guard against older binaries that may omit kappa column names.
      if (!is.null(kn) && length(kn) == ncol(result$shrinkage_kappa_by_occ) - 1L) {
        colnames(result$shrinkage_kappa_by_occ) <- c("occ", kn)
      }
    } else {
      result$shrinkage_kappa_by_occ <- NULL
    }
    # Sample-size-weighted IOV (ferx-core #1031). The engine leaves both vectors
    # empty unless some kappa carries `weight = <expr>`, so NULL them out on the
    # ordinary IOV path and keep them parallel to kappa_names otherwise.
    result[c("kappa_weights", "kappa_weight_typical")] <-
      .ferx_name_kappa_weights(result$kappa_weights, result$kappa_weight_typical,
                               result$kappa_names, d_iov)
  } else {
    result$omega_iov <- NULL
    result$se_kappa <- NULL
    result$shrinkage_kappa <- NULL
    result$shrinkage_kappa_by_occ <- NULL
    result$kappa_names <- NULL
    result$ebe_kappas <- NULL
    result$kappa_weights <- NULL
    result$kappa_weight_typical <- NULL
  }

  # Reshape omega_param_corr into a square matrix using omega_dim to guard
  # against unexpected vector lengths (NULL when diagonal or absent).
  if (!is.null(result$omega_param_corr) && length(result$omega_param_corr) > 0L) {
    d_pc <- result$omega_dim %||% 0L
    if (d_pc > 0L && length(result$omega_param_corr) == d_pc * d_pc) {
      result$omega_param_corr <- matrix(
        result$omega_param_corr, nrow = d_pc, ncol = d_pc, byrow = TRUE
      )
    } else {
      result$omega_param_corr <- NULL
    }
  } else {
    result$omega_param_corr <- NULL
  }

  # Reshape omega_iov_param_corr into a square matrix (NULL when diagonal or absent)
  if (!is.null(result$omega_iov_param_corr) && length(result$omega_iov_param_corr) > 0L) {
    d_ipc <- result$omega_iov_dim %||% 0L
    if (d_ipc > 0L && length(result$omega_iov_param_corr) == d_ipc * d_ipc) {
      result$omega_iov_param_corr <- matrix(
        result$omega_iov_param_corr, nrow = d_ipc, ncol = d_ipc, byrow = TRUE
      )
    } else {
      result$omega_iov_param_corr <- NULL
    }
  } else {
    result$omega_iov_param_corr <- NULL
  }

  # eta_log_transformed: empty vector -> NULL
  if (is.null(result$eta_log_transformed) || length(result$eta_log_transformed) == 0L) {
    result$eta_log_transformed <- NULL
  }

  # Parameter transform metadata - fall back gracefully for older ferx-core binaries
  # that don't populate these vectors (empty character() from FFI ? treat as absent).
  n_eta_fit   <- if (!is.null(result$omega)) nrow(result$omega) else 0L
  n_theta_fit <- length(result$theta)
  n_sigma_fit <- length(result$sigma)

  if (is.null(result$eta_param_types) || length(result$eta_param_types) != n_eta_fit) {
    result$eta_param_types <- rep("log_normal", n_eta_fit)
  }
  if (is.null(result$eta_linked_theta) || length(result$eta_linked_theta) != n_eta_fit) {
    result$eta_linked_theta <- rep("", n_eta_fit)
  }
  if (is.null(result$theta_transforms) || length(result$theta_transforms) != n_theta_fit) {
    result$theta_transforms <- rep("identity", n_theta_fit)
  }
  if (is.null(result$sigma_types) || length(result$sigma_types) != n_sigma_fit) {
    result$sigma_types <- rep("proportional", n_sigma_fit)
  }
  names(result$theta_transforms) <- names(result$theta)
  names(result$sigma_types)      <- names(result$sigma)

  # Clean up internal fields
  result$theta_names <- NULL
  result$omega_dim <- NULL
  result$omega_iov_dim <- NULL

  # Print mu-referencing detections as informational messages.
  for (eta_names in .ferx_mu_ref_detection_names(result$warnings)) {
    message("Mu-referencing detected for: ", eta_names)
  }

  # Store the dataset name (basename without extension) from the data path
  result$data_name <- tools::file_path_sans_ext(basename(data))

  # Fall back to the model file's basename when the .ferx file declares no name.
  # ferx-core returns "Unnamed" (not NULL/"") when no [model_name] is declared -
  # treat it the same as absent. See fit_result_to_list() in src/rust/src/lib.rs.
  # Long-term fix: engine should return "" so this sentinel check can be removed.
  if (is.null(result$model_name) || !nzchar(result$model_name) || identical(result$model_name, "Unnamed")) {
    result$model_name <- tools::file_path_sans_ext(basename(model))
  }

  # `result$model_structure` is now populated by the Rust engine in
  # `fit_result_to_list()` from the parsed CompiledModel (ferx-core#49), so it
  # reflects exactly what ferx-core ran. The R-side `.ferx_parse_structure()`
  # remains for the pre-fit `ferx_model_inspect(path)` workflow only.

  # Store the requested gradient method, and the method the engine actually
  # resolved to. The engine reports `gradient_method_inner` as one of
  # Map verbose engine labels to the short tokens used
  # by the `gradient` argument so the two values can be compared directly
  # (e.g. requested "auto" -> used "ad"). When the engine omitted the field
  # (older binaries return ""), `gradient_used` is NA.
  # `gradient` defaults to NULL ("defer to the model file / engine default").
  # Assigning that raw NULL would drop the field entirely (and with it the
  # RUN INFO line), so resolve it to the requested token: the explicit arg when
  # given, otherwise the model-file value, otherwise the engine default "auto".
  result$gradient <- if (nzchar(gradient_arg)) {
    gradient_arg
  } else {
    mf_grad <- if ("gradient" %in% names(model_file_opts)) {
      model_file_opts[["gradient"]]
    } else if ("gradient_method" %in% names(model_file_opts)) {
      model_file_opts[["gradient_method"]]
    } else {
      NULL
    }
    if (!is.null(mf_grad)) tolower(as.character(mf_grad)) else "auto"
  }
  result$gradient_used <- .ferx_short_gradient_label(result$gradient_method_inner)

  # Store effective settings (already serialised from settings_parts above)
  result$call_settings <- settings_used

  # Raw [fit_options] from the model file, for inspection and conflict auditing.
  result$model_file_settings <- if (length(model_file_opts) > 0L)
    as.list(model_file_opts) else list()

  # Stash inputs so ferx_save_fit() can embed `model.ferx` and (optionally)
  # `data.csv` into a portable .fitrx bundle. Read failures here must not
  # break the fit, so wrap in tryCatch.
  result$model_source <- tryCatch(
    paste(readLines(model, warn = FALSE), collapse = "\n"),
    error = function(e) ""
  )

  # Source-file provenance. The Rust binding emits the four fields
  # directly (paths as the caller-supplied strings, hashes via
  # `ferx_core::io::hash::sha256_file`). Always run `normalizePath` on
  # the path fields here, regardless of source: the old R-side behavior
  # was to store absolute paths (`normalizePath(data)`), and downstream
  # code (notably ferx_save_fit, where `fit$data_path` is opened from a
  # potentially-different working directory at save time) depends on
  # that. Falling back to NULL ? NA when the binding returns nothing
  # keeps older binaries working.
  empty_to_null <- function(x) {
    if (is.null(x) || length(x) == 0L || !nzchar(x[[1L]])) NULL else as.character(x)
  }
  normalize_or_na <- function(path_str) {
    if (is.null(path_str)) return(NA_character_)
    tryCatch(normalizePath(path_str, mustWork = FALSE),
             error = function(e) NA_character_)
  }
  result$model_path <- normalize_or_na(
    empty_to_null(result$model_path) %||% model
  )
  result$data_path <- normalize_or_na(
    empty_to_null(result$data_path) %||% data
  )
  result$model_hash <- empty_to_null(result$model_hash) %||% NA_character_
  result$data_hash <- empty_to_null(result$data_hash) %||% NA_character_

  # Assemble the structured-warning table. Core supplies severity/category for
  # every warning it emitted (including Durbin-Watson autocorrelation); the R
  # side appends only the diagnostics it computes itself (condition number,
  # ETA normality). No string re-parsing of core messages happens here.
  result$warnings_structured <- .ferx_assemble_structured_warnings(raw, result)

  # Reconstruct the per-iteration IMPMAP parameter trace as a data.frame.
  # The Rust side passes flat vectors + metadata; NULL when not collected.
  # Belt-and-suspenders: only keep it when the caller actually asked for it
  # (settings= arg or [fit_options] in the model file) - guards against the
  # trace leaking onto the fit object from an intermediate stage of a method
  # chain (e.g. c("impmap", "focei")) even when impmap_trace was never
  # requested for that chain.
  impmap_requested <- .ferx_impmap_trace_requested(settings_used, model_file_opts)
  if (!is.null(result$impmap_trace)) {
    result$impmap_trace <- if (impmap_requested) {
      .reconstruct_impmap_trace(result$impmap_trace)
    } else {
      NULL
    }
  }

  # Derived fields (cor_matrix, estimates, eta_cov) - shared with
  # ferx_load_fit() via .ferx_populate_derived_fields() so the two
  # construction paths can't diverge. Uses result$data_path (normalised
  # above), not the raw `data` argument, so a fresh fit and a loaded fit of
  # the same model resolve the dataset identically.
  result <- .ferx_populate_derived_fields(result)

  class(result) <- "ferx_fit"

  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L || is.na(output)) {
      stop("`output` must be a single .fitrx file path or NULL")
    }
    ferx_save_fit(result, output, include_data = isTRUE(include_data))
  }

  result
}

# Extract legacy mu-reference detection entries from the warnings vector.
# This intentionally ignores warnings such as "not mu-referenced".
.ferx_mu_ref_detection_names <- function(warnings) {
  mu_ref_warnings <- grep(
    "^(\\[[^]]+\\]\\s*)?mu-ref:",
    warnings %||% character(0),
    value = TRUE
  )
  sub("^(\\[[^]]+\\]\\s*)?mu-ref:\\s*", "", mu_ref_warnings)
}

# Map the engine's gradient-kind label to the short tokens used by the
# `gradient` argument so requested/used can be compared at a glance.
# Returns NA_character_ when the engine field is missing or unrecognized so
# downstream display code can fall back gracefully on older binaries.
.ferx_short_gradient_label <- function(kind) {
  if (is.null(kind) || !nzchar(kind)) return(NA_character_)
  switch(kind,
    "analytic (Dual2)"   = "analytic",
    "finite differences" = "fd",
    "N/A"                = "N/A",
    kind
  )
}

# Compute Shapiro-Wilk normality test for each ETA, one row per ID.
# Returns a data.frame with columns: eta, W, p_val, flag.
# Converts the raw FFI sentinels (empty vec, NaN) to idiomatic R (NULL) and
# appends a warning when condition_number > 1000 (NONMEM convention).
# Inf is kept as-is: it signals a non-positive eigenvalue, which the covariance
# step already warned about (covariance_status != "computed").
.ferx_apply_cov_sentinels <- function(result) {
  result$eigenvalues <- if (length(result$cov_eigenvalues) == 0L) NULL else result$cov_eigenvalues
  result$condition_number <- if (is.nan(result$cov_condition_number)) NULL else result$cov_condition_number
  result$cov_eigenvalues <- NULL
  result$cov_condition_number <- NULL
  if (!is.null(result$condition_number) && is.finite(result$condition_number) &&
        result$condition_number > 1000) {
    result$warnings <- c(
      result$warnings,
      sprintf(
        "High condition number (%.1f) -- parameter space may be ill-conditioned",
        result$condition_number
      )
    )
  }
  result
}

# `shapiro.test()` is hard-capped at 5000 observations; for larger N each
# ETA is systematically subsampled to 5000 values to keep the diagnostic
# from failing, and a single warning is emitted to flag that the result
# is approximate.
.ferx_compute_eta_normality <- function(ebe_etas) {
  if (is.null(ebe_etas) || !is.data.frame(ebe_etas)) {
    return(NULL)
  }
  # ebe_etas is purpose-built: ID + one column per BSV eta. The eta
  # declaration name is whatever the model file uses (ETA_CL, eta_CL,
  # RE_CL, ...), so a "^ETA" prefix filter would silently drop valid
  # columns from non-conventional models. Treat every non-ID column as
  # an eta.
  eta_cols <- setdiff(names(ebe_etas), "ID")
  if (length(eta_cols) == 0L) {
    return(NULL)
  }
  one_per <- ebe_etas[, eta_cols, drop = FALSE]
  sw_cap <- 5000L
  if (nrow(one_per) > sw_cap) {
    warning(
      "ETA normality: ", nrow(one_per), " subjects exceeds the ",
      sw_cap, "-sample limit of shapiro.test(); each ETA was subsampled ",
      "to ", sw_cap, " values. Shapiro-Wilk is also known to over-reject ",
      "normality at large N -- treat the result as a rough indicator and ",
      "prefer QQ-plots for diagnosis on large datasets.",
      call. = FALSE
    )
  }
  do.call(rbind, lapply(eta_cols, function(col) {
    vals <- one_per[[col]]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 3L) {
      return(data.frame(
        eta = col, W = NA_real_, p_val = NA_real_, flag = "",
        stringsAsFactors = FALSE
      ))
    }
    if (length(unique(vals)) == 1L) {
      return(data.frame(
        eta = col, W = NA_real_, p_val = NA_real_, flag = "",
        stringsAsFactors = FALSE
      ))
    }
    if (length(vals) > sw_cap) {
      vals <- vals[round(seq.int(1L, length(vals), length.out = sw_cap))]
    }
    sw <- shapiro.test(vals)
    flg <- if (sw$p.value < 0.05) "[!]" else ""
    data.frame(
      eta = col, W = round(as.numeric(sw$statistic), 3),
      p_val = round(sw$p.value, 4), flag = flg, stringsAsFactors = FALSE
    )
  }))
}

# Internal: canonicalise a single `method` token to its engine name, folding
# the documented case-insensitive aliases before `match.arg`.
.normalize_method_token <- function(m) {
  # `tolower` *before* the `[^a-z0-9]` strip - otherwise uppercase letters
  # get smashed to underscores ("IMP" -> "___") before tolower runs on
  # them, and the canonicalisation silently drops the case-insensitive
  # entry points the docs advertise.
  normalised <- gsub("[^a-z0-9]", "_", tolower(m))
  # Fold IMP aliases before `match.arg` - partial matching won't see
  # `importance_sampling` as a prefix of `imp`, so the documented
  # `c("focei", "importance_sampling")` form has to be mapped explicitly.
  if (normalised %in% c("importance_sampling", "importancesampling")) {
    normalised <- "imp"
  }
  # IMPMAP aliases. Fold before the IMP set would never match these (they
  # carry the `_map` suffix), so order is irrelevant; kept explicit for
  # the documented `"importance_sampling_map"` long form.
  if (normalised %in% c("importance_sampling_map", "importancesamplingmap")) {
    normalised <- "impmap"
  }
  if (normalised %in% c("bayesian", "mcmc")) {
    normalised <- "bayes"
  }
  # `method = "agq"` was removed (ferx-core #251): adaptive quadrature is now
  # `method = "laplace"` with `settings = list(n_agq = N)` (exact Hessian anchor),
  # or `method = "focei"` with `n_agq > 1` (Gauss-Newton anchor). Reject the old
  # tokens with a pointer rather than routing them to a method that no longer
  # exists. These `contains("gauss")` near-misses must be caught here before the
  # Gauss-Newton fold below. Mirrors ferx-core's `parse_method_token`.
  if (normalised %in% c(
    "agq", "aghq", "gauss_hermite", "adaptive_gaussian_quadrature"
  )) {
    stop(
      "method = \"agq\" has been removed: use method = \"laplace\" with ",
      "settings = list(n_agq = N) (n_agq = 1 is the Laplace approximation, ",
      "n_agq > 1 is adaptive Gauss-Hermite quadrature). For the ",
      "Gauss-Newton-anchored quadrature use method = \"focei\" with n_agq > 1.",
      call. = FALSE
    )
  }
  # `laplacian` is NONMEM's spelling; fold it onto the canonical `laplace`.
  if (normalised == "laplacian") {
    normalised <- "laplace"
  }
  # `match.arg` uses exact-match-first semantics, so listing both "imp" and
  # "impmap" is unambiguous (each exact token wins over the other's prefix).
  match.arg(
    normalised,
    c(
      "foce", "focei", "laplace", "saem", "gn", "gn_hybrid",
      "imp", "impmap", "bayes"
    )
  )
}

# Internal: rebuild the IMPMAP per-iteration parameter trace from the flat
# representation the Rust side returns (flat values + n_rows/n_cols/col_names)
# into a data.frame with one row per Monte-Carlo EM iteration.
.reconstruct_impmap_trace <- function(tr) {
  mat <- matrix(tr$flat, nrow = tr$n_rows, ncol = tr$n_cols, byrow = FALSE)
  df <- as.data.frame(mat)
  colnames(df) <- tr$col_names
  df
}

# Internal: was `impmap_trace = TRUE` actually requested by the caller, via
# either the `settings=` argument (already parsed into `settings_used`) or a
# `[fit_options] impmap_trace = true` line in the model file
# (`model_file_opts`, raw strings from the parser)?
.ferx_impmap_trace_requested <- function(settings_used, model_file_opts) {
  mf_val <- if ("impmap_trace" %in% names(model_file_opts)) {
    model_file_opts[["impmap_trace"]]
  } else {
    ""
  }
  isTRUE(settings_used[["impmap_trace"]]) ||
    identical(tolower(as.character(mf_val)), "true")
}

.FERX_BAYES_RHAT_THRESHOLD <- 1.01

# Styled "[!]" marker (yellow where colour is available, matching the shrinkage
# flag) when `max_rhat` exceeds the threshold; "" otherwise. Leading spaces
# separate it from the preceding text.
.ferx_bayes_rhat_flag <- function(max_rhat) {
  if (is.finite(max_rhat) && max_rhat > .FERX_BAYES_RHAT_THRESHOLD) {
    paste0("  ", .ferx_style("[!]", "yellow"))
  } else {
    ""
  }
}

#' @export
print.ferx_fit <- function(x, ...) {
  bar <- strrep("=", 60)
  cat(bar, "\n", sep = "")
  cat(" NONLINEAR MIXED EFFECTS MODEL ESTIMATION\n")
  cat(bar, "\n", sep = "")

  # Compact identification line: model and dataset on one line when both
  # present (a common case), else fall back to single-field lines.
  model_lbl <- if (!is.null(x$model_name) && nzchar(x$model_name)) x$model_name else NULL
  data_lbl  <- if (!is.null(x$data_name)  && nzchar(x$data_name))  x$data_name  else NULL
  if (!is.null(model_lbl) || !is.null(data_lbl)) {
    cat(sprintf(" Model: %-24s Dataset: %s\n",
                model_lbl %||% "(unnamed)",
                data_lbl %||% "(unnamed)"))
  }

  # Method | Gradient | Subjects | Obs - one pipe-separated line so the run's
  # shape fits in a single glance.
  sde_tag <- if (isTRUE(x$uses_sde)) " + SDE (EKF)" else ""
  method_str <- if (!is.null(x$method_chain) && length(x$method_chain) > 1) {
    paste(toupper(x$method_chain), collapse = " -> ")
  } else {
    toupper(x$method %||% "?")
  }
  method_str <- paste0(method_str, sde_tag)
  grad_str <- if (!is.null(x$gradient_used) && !is.na(x$gradient_used) &&
                    nzchar(x$gradient_used)) {
    toupper(x$gradient_used)
  } else {
    "?"
  }
  cat(sprintf(" Method: %s | Gradient: %s | Subjects: %s | Obs: %s\n",
              method_str, grad_str,
              format(x$n_subjects %||% NA),
              format(x$n_obs %||% NA)))

  # STATUS line - the most important diagnostic; colour green/red so the
  # eye locks on it immediately.
  status_lbl <- if (isTRUE(x$converged)) "CONVERGED" else "NOT CONVERGED"
  status_style <- if (isTRUE(x$converged)) "green" else "red"
  status_tail <- character(0)
  # Categories come from the engine; no R-side string parsing needed.
  ws <- x$warnings_structured
  if (!is.null(ws) && is.data.frame(ws) && all(c("severity", "category") %in% names(ws))) {
    crit_rows <- ws[ws$severity == "critical", , drop = FALSE]
    if (nrow(crit_rows) > 0L) {
      cats <- unique(crit_rows$category)
      status_tail <- c(status_tail, sprintf("(%s)", paste(cats, collapse = ", ")))
    }
  }
  if (!is.null(x$n_iterations)) status_tail <- c(status_tail, sprintf("%d iterations", x$n_iterations))
  if (!is.null(x$wall_time_secs)) status_tail <- c(status_tail, sprintf("%.1fs", x$wall_time_secs))
  status_tail_str <- if (length(status_tail) > 0L) paste0("   ", paste(status_tail, collapse = "   ")) else ""
  cat("\n STATUS: ", .ferx_style(status_lbl, status_style), status_tail_str, "\n", sep = "")

  # OFV / AIC / BIC on a single line.
  cat(sprintf(" OFV: %.4f    AIC: %.4f    BIC: %.4f\n", x$ofv, x$aic, x$bic))

  # DATA SELECTION block (only when [data_selection] rules were active).
  ex <- x$exclusions
  if (!is.null(ex)) {
    cat("\n", .ferx_style("DATA SELECTION", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    cat(sprintf(
      "  Records read: %d    Obs excl.: %d    Doses excl.: %d    Other excl.: %d\n",
      ex$n_records_total %||% 0L,
      ex$n_obs_excluded  %||% 0L,
      ex$n_dose_excluded %||% 0L,
      ex$n_other_excluded %||% 0L
    ))
    if (length(ex$excluded_subject_ids) > 0L)
      cat(sprintf("  Subjects excluded: %s\n",
                  paste(ex$excluded_subject_ids, collapse = ", ")))
    strip_pfx <- function(x) sub("^(ignore|accept): ", "", x)
    for (cond in ex$fired_ignore)
      cat(sprintf("  Fired ignore: %s\n", strip_pfx(cond)))
    for (cond in ex$fired_accept)
      cat(sprintf("  Fired accept: %s\n", strip_pfx(cond)))
  }

  if (!is.null(x$model_structure)) {
    cat("\n", .ferx_style("MODEL STRUCTURE (auto-derived)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    .ferx_print_structure(x$model_structure)
  }

  # THETA
  cat("\n", .ferx_style("THETA", "bold"), "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  cat(sprintf("%-16s %12s %12s %10s\n", "Parameter", "Estimate", "SE", "%RSE"))
  cat(strrep("-", 52), "\n", sep = "")
  theta_names <- names(x$theta)
  if (is.null(theta_names)) theta_names <- paste0("THETA", seq_along(x$theta))
  # Logical mask of NN-weight thetas to skip - they appear in the compact
  # NEURAL NETWORKS block below (Option E). 1-based indices to match `seq_along`.
  nn_skip <- logical(length(x$theta))
  if (!is.null(x$neural_networks)) {
    for (nn in x$neural_networks) {
      if (!is.null(nn$weights_offset) && !is.null(nn$n_weights) && nn$n_weights > 0) {
        nn_skip[seq(nn$weights_offset + 1L, nn$weights_offset + nn$n_weights)] <- TRUE
      }
    }
  }
  for (i in seq_along(x$theta)) {
    if (nn_skip[i]) next
    est       <- x$theta[i]
    transform <- if (!is.null(x$theta_transforms) && length(x$theta_transforms) >= i) x$theta_transforms[i] else "identity"
    if (!is.null(x$se_theta) && length(x$se_theta) >= i) {
      se_val  <- x$se_theta[i]
      rse     <- if (abs(est) > 1e-12) abs(se_val / est) * 100 else NaN
      se_str  <- sprintf("%.6f", se_val)
      # Highlight high-uncertainty estimates (%RSE > 50) - same convention
      # used by the [!] shrinkage flag and by ferx_get_warnings() severity colours.
      rse_raw <- sprintf("%.1f", rse)
      rse_str <- if (!is.nan(rse) && rse > 50) .ferx_style(rse_raw, "yellow") else rse_raw
    } else {
      se_val  <- NA_real_
      se_str  <- "N/A"
      rse_str <- "N/A"
    }
    scale_tag <- switch(transform,
      log              = "  [log scale]",
      logit            = "  [logit scale]",
      logit_probability = "  [logit scale]",
      ""
    )
    cat(sprintf("%-16s %12.6f %12s %10s%s\n", theta_names[i], est, se_str, rse_str, scale_tag))
    # Natural-scale typical value for log/logit transforms
    if (transform == "log") {
      tv <- exp(est)
      ci_str <- if (!is.na(se_val)) {
        sprintf("95%% CI: [%.4f, %.4f]", exp(est - 1.96 * se_val), exp(est + 1.96 * se_val))
      } else ""
      cat(sprintf("  %-14s %12.4f                    %s\n", "(typical)", tv, ci_str))
    } else if (transform %in% c("logit", "logit_probability")) {
      tv <- .ferx_inv_logit(est)
      ci_str <- if (!is.na(se_val)) {
        sprintf("95%% CI: [%.4f, %.4f]", .ferx_inv_logit(est - 1.96 * se_val), .ferx_inv_logit(est + 1.96 * se_val))
      } else ""
      cat(sprintf("  %-14s %12.4f                    %s\n", "(typical)", tv, ci_str))
    }
  }

  # NEURAL NETWORKS
  # Compact summary block when `[covariate_nn]` blocks are present (only
  # populated when ferx-r was built with `--features nn` and the model
  # declares NN blocks). The per-weight thetas (W_NN_*, B_NN_*) are still
  # in `x$theta`; this block exists so users get a readable overview
  # without scrolling past 100+ rows.
  if (!is.null(x$neural_networks) && length(x$neural_networks) > 0) {
    cat("\n", .ferx_style("NEURAL NETWORKS", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    for (nn in x$neural_networks) {
      shape_str <- paste(nn$shape, collapse = ", ")
      cat(sprintf(
        "%s  shape=[%s]  activation=%s/%s  n_weights=%d\n",
        nn$name, shape_str, nn$hidden_activation, nn$output_activation, nn$n_weights
      ))
      cat(sprintf("  inputs:  [%s]\n", paste(nn$input_names, collapse = ", ")))
      cat(sprintf("  outputs: [%s]\n", paste(nn$output_names, collapse = ", ")))
      if (length(nn$weights) > 0) {
        cat(sprintf(
          "  weights: min %.4f  max %.4f  mean %.4f  std %.4f\n",
          min(nn$weights), max(nn$weights), mean(nn$weights), stats::sd(nn$weights)
        ))
      }
    }
  }

  # OMEGA
  om <- x$omega
  if (is.null(dim(om))) om <- matrix(om, 1, 1)
  n_eta <- nrow(om)
  # A fixed-effects-only (naive-pooled) fit - ferx-core #989 - has no Omega to
  # report. Printing the header above an empty body reads as an estimation that
  # failed rather than one that was never requested; ferx-core's own console
  # output suppresses the same section. Everything below is already a no-op at
  # n_eta = 0 (the loops are over seq_len(n_eta), and has_offdiag stays FALSE).
  if (n_eta > 0L) {
    cat("\n", .ferx_style("OMEGA  (between-subject variability)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
  }
  has_offdiag <- FALSE
  for (i in seq_len(n_eta)) {
    var_ii   <- om[i, i]
    eta_type <- if (!is.null(x$eta_param_types) && length(x$eta_param_types) >= i) x$eta_param_types[i] else "log_normal"
    linked_theta_name <- if (!is.null(x$eta_linked_theta) && length(x$eta_linked_theta) >= i) x$eta_linked_theta[i] else ""
    se_val <- .omega_se_at(x$se_omega, n_eta, i, i)
    se_str <- if (!is.na(se_val)) sprintf("%.6f", se_val) else "N/A"

    extra <- if (eta_type == "log_normal") {
      cv <- if (var_ii > 0) sqrt(exp(var_ii) - 1) * 100 else 0
      sprintf("CV%% = %.1f", cv)
    } else if (eta_type == "additive") {
      sd_val  <- sqrt(max(var_ii, 0))
      tv_name <- if (nzchar(linked_theta_name)) linked_theta_name else NULL
      tv_val  <- if (!is.null(tv_name)) x$theta[tv_name] else NA_real_
      cv_part <- if (!is.null(tv_val) && !is.na(tv_val) && abs(tv_val) > 1e-12) {
        sprintf("  CV%% = %.1f (at %s=%.3g)", sd_val / abs(tv_val) * 100, tv_name, tv_val)
      } else ""
      sprintf("SD = %.4f%s", sd_val, cv_part)
    } else if (eta_type %in% c("logit", "logit_probability")) {
      # Natural-scale +/-1 SD range using linked theta (if known)
      tv_name <- if (nzchar(linked_theta_name)) linked_theta_name else NULL
      tv_val  <- if (!is.null(tv_name)) x$theta[tv_name] else NA_real_
      if (!is.null(tv_val) && !is.na(tv_val)) {
        sd_logit <- sqrt(max(var_ii, 0))
        lo <- .ferx_inv_logit(tv_val - sd_logit)
        hi <- .ferx_inv_logit(tv_val + sd_logit)
        sprintf("SD_logit = %.4f  +/-1SD = [%.3f, %.3f]", sd_logit, lo, hi)
      } else {
        sprintf("SD_logit = %.4f", sqrt(max(var_ii, 0)))
      }
    } else {
      "CV% = N/A"
    }
    type_tag <- switch(eta_type,
      log_normal        = "[log-normal]",
      additive          = "[additive]",
      logit             = "[logit]",
      logit_probability = "[logit]",
      custom            = "[custom]",
      ""
    )
    eta_label <- if (!is.null(x$eta_names) && length(x$eta_names) >= i && nzchar(x$eta_names[i])) x$eta_names[i] else sprintf("OMEGA(%d,%d)", i, i)
    sd_note <- if (!is.null(x$omega_init_as_sd) && length(x$omega_init_as_sd) >= i && isTRUE(x$omega_init_as_sd[i])) "  [initial specified as SD]" else ""
    cat(sprintf(
      "  %-24s %-13s = %.6f  %s  SE = %s%s\n",
      eta_label, type_tag, var_ii, extra, se_str, sd_note
    ))
    for (j in seq_len(i - 1L)) {
      if (abs(om[i, j]) > 1e-15) has_offdiag <- TRUE
    }
  }
  if (has_offdiag) {
    cat("  --- Correlations ---\n")
    for (i in seq_len(n_eta)) {
      for (j in seq_len(i - 1L)) {
        cov_ij <- om[i, j]
        param_corr <- if (!is.null(x$omega_param_corr)) {
          x$omega_param_corr[i, j]
        } else {
          var_i <- om[i, i]
          var_j <- om[j, j]
          if (var_i > 0 && var_j > 0) cov_ij / (sqrt(var_i) * sqrt(var_j)) else 0
        }
        corr_label <- if (!is.null(x$omega_param_corr)) "param corr" else "corr"
        se_ij <- .omega_se_at(x$se_omega, n_eta, i, j)
        se_part <- if (!is.na(se_ij)) sprintf("  SE = %.6f", se_ij) else ""
        has_nms <- !is.null(x$eta_names) && length(x$eta_names) >= i
        if (has_nms) {
          lbl_i <- if (nzchar(x$eta_names[i])) x$eta_names[i] else sprintf("OMEGA(%d,%d)", i, i)
          lbl_j <- if (nzchar(x$eta_names[j])) x$eta_names[j] else sprintf("OMEGA(%d,%d)", j, j)
          cat(sprintf(
            "  %s ~ %s : cov = %.6f  (%s = %.4f)%s\n",
            lbl_i, lbl_j, cov_ij, corr_label, param_corr, se_part
          ))
        } else {
          cat(sprintf(
            "  OMEGA(%d,%d) : cov = %.6f  (%s = %.4f)%s\n",
            i, j, cov_ij, corr_label, param_corr, se_part
          ))
        }
      }
    }
  }

  # OMEGA_IOV (Inter-Occasion Variability)
  if (!is.null(x$omega_iov)) {
    cat("\n", .ferx_style("OMEGA_IOV  (inter-occasion variability)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    m_iov <- x$omega_iov
    if (is.null(dim(m_iov))) m_iov <- matrix(m_iov, 1, 1)
    n_kap <- nrow(m_iov)
    kap_names <- x$kappa_names
    if (is.null(kap_names)) kap_names <- paste0("KAPPA", seq_len(n_kap))
    n_se <- length(x$se_kappa)
    n_tri <- n_kap * (n_kap + 1L) / 2L
    is_block_se <- (n_se == n_tri && n_kap > 1L)
    # For block kappa the SEs are packed as lower triangle (column-major).
    # Diagonal element for column j (1-indexed) sits at flat index j*n_kap - j*(j-1)/2.
    diag_se_idx <- function(j) j * n_kap - j * (j - 1L) / 2L - (n_kap - j)
    iov_has_offdiag <- FALSE
    for (i in seq_len(n_kap)) {
      var_ii   <- m_iov[i, i]
      kap_type <- if (!is.null(x$kappa_param_types) && length(x$kappa_param_types) >= i) x$kappa_param_types[i] else "log_normal"
      se_idx   <- if (is_block_se) diag_se_idx(i) else i
      se_str   <- if (!is.null(x$se_kappa) && n_se >= se_idx) sprintf("%.6f", x$se_kappa[se_idx]) else "N/A"
      shr_str  <- if (!is.null(x$shrinkage_kappa) && length(x$shrinkage_kappa) >= i) {
        sprintf("%.1f%%", x$shrinkage_kappa[i] * 100)
      } else "N/A"
      cv_str <- if (kap_type == "log_normal") {
        if (var_ii > 0) sprintf("%.1f", sqrt(exp(var_ii) - 1) * 100) else "0.0"
      } else if (kap_type == "additive") {
        sprintf("SD = %.4f", sqrt(max(var_ii, 0)))
      } else "N/A"
      cat(sprintf(
        "  %s = %.6f  (CV%% = %s)  SE = %s  Shrinkage = %s\n",
        kap_names[i], var_ii, cv_str, se_str, shr_str
      ))
      # Sample-size-weighted IOV (ferx-core #1031): the estimate above is the
      # *unweighted* gamma^2 - the quantity a published MBMA reports - so print
      # the effective SD at a typical weight next to it.
      wline <- .ferx_format_kappa_weight(x, i, var_ii, kap_names[i])
      if (!is.null(wline)) cat(wline, "\n", sep = "")
      for (j in seq_len(i - 1L)) {
        if (abs(m_iov[i, j]) > 1e-15) iov_has_offdiag <- TRUE
      }
    }
    if (iov_has_offdiag) {
      cat("  --- Correlations ---\n")
      for (i in seq_len(n_kap)) {
        for (j in seq_len(i - 1L)) {
          cov_ij <- m_iov[i, j]
          param_corr <- if (!is.null(x$omega_iov_param_corr)) {
            x$omega_iov_param_corr[i, j]
          } else {
            var_i <- m_iov[i, i]
            var_j <- m_iov[j, j]
            if (var_i > 0 && var_j > 0) cov_ij / (sqrt(var_i) * sqrt(var_j)) else 0
          }
          corr_label <- if (!is.null(x$omega_iov_param_corr)) "param corr" else "corr"
          cat(sprintf(
            "  %s ~ %s : cov = %.6f  (%s = %.4f)\n",
            kap_names[i], kap_names[j], cov_ij, corr_label, param_corr
          ))
        }
      }
    }
    if (!is.null(x$shrinkage_kappa_by_occ) &&
        is.data.frame(x$shrinkage_kappa_by_occ) &&
        nrow(x$shrinkage_kappa_by_occ) >= 1L) {
      cat("  --- Shrinkage by occasion ---\n")
      df_occ <- x$shrinkage_kappa_by_occ
      kap_cols <- setdiff(colnames(df_occ), "occ")
      for (row_i in seq_len(nrow(df_occ))) {
        parts <- vapply(kap_cols, function(cn) {
          sprintf("%s %.1f%%", cn, df_occ[[cn]][row_i] * 100)
        }, character(1L))
        cat(sprintf("  Occasion %d: %s\n", df_occ$occ[row_i], paste(parts, collapse = "  ")))
      }
    }
  }

  # SIGMA - sigma is on the SD scale (see ferx-core `docs/src/model-file/error-model.md`).
  # For proportional components CV% = sigma * 100; for additive components the
  # value is in observation units and no CV% applies. Variance = sigma^2 in
  # both cases, mirroring the YAML output (ferx-core#57).
  cat("\n", .ferx_style("SIGMA  (residual error)", "bold"), "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  for (i in seq_along(x$sigma)) {
    s   <- x$sigma[i]
    nm  <- if (length(x$sigma_names) >= i && nzchar(x$sigma_names[i])) {
      x$sigma_names[i]
    } else {
      sprintf("SIGMA(%d)", i)
    }
    typ <- if (length(x$sigma_types) >= i) x$sigma_types[i] else NA_character_
    se_str <- if (!is.null(x$se_sigma) && length(x$se_sigma) >= i) {
      sprintf("%.6f", x$se_sigma[i])
    } else {
      "N/A"
    }
    var_str  <- sprintf("var = %.6f", s * s)
    type_tag <- switch(as.character(typ),
      proportional = "[proportional]",
      additive     = "[additive]",
      ""
    )
    sd_note <- if (!is.null(x$sigma_init_as_sd) && length(x$sigma_init_as_sd) >= i && isTRUE(x$sigma_init_as_sd[i])) "  [initial specified as SD]" else ""
    if (identical(typ, "proportional")) {
      cat(sprintf("  %-16s %-14s = %.6f  (%s, CV%% = %.1f)  SE = %s%s\n",
                  nm, type_tag, s, var_str, s * 100, se_str, sd_note))
    } else if (identical(typ, "additive")) {
      cat(sprintf("  %-16s %-14s = %.6f  (SD = %.4f, %s)  SE = %s%s\n",
                  nm, type_tag, s, s, var_str, se_str, sd_note))
    } else {
      cat(sprintf("  %-16s = %.6f  (%s)  SE = %s%s\n",
                  nm, s, var_str, se_str, sd_note))
    }
  }

  # SIR uncertainty
  if (!is.null(x$sir_ess)) {
    cat("\n", .ferx_style("SIR  (95% CI from sampling importance resampling)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    cat(sprintf("Effective sample size: %.1f\n", x$sir_ess))
    print_ci <- function(m) {
      if (is.null(m)) {
        return(invisible())
      }
      for (i in seq_len(nrow(m))) {
        cat(sprintf("  %s : [%.6f, %.6f]\n", rownames(m)[i], m[i, 1], m[i, 2]))
      }
    }
    print_ci(x$sir_ci_theta)
    print_ci(x$sir_ci_omega)
    print_ci(x$sir_ci_sigma)
  }

  # Importance Sampling marginal log-likelihood (IMP terminal stage)
  if (!is.null(x$importance_sampling)) {
    is_r <- x$importance_sampling
    cat("\n", .ferx_style("IMPORTANCE SAMPLING  (marginal log-likelihood)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    cat(sprintf(
      "  -2 log L (IS): %.4f  (MC SE = %.4f, K = %d, df = %g)\n",
      is_r$minus2_log_likelihood,
      is_r$mc_standard_error,
      is_r$n_samples,
      is_r$proposal_df
    ))
    cat(sprintf(
      "  ESS / K: min = %.3f, median = %.3f\n",
      is_r$ess_min,
      is_r$ess_median
    ))
    if (identical(is_r$kappa_treatment, "fixed_at_mode")) {
      cat("  Note: kappa fixed at EBE (partial marginal - not fully comparable to NONMEM IMP)\n")
    }
    low_ids <- is_r$low_ess_subject_ids
    if (!is.null(low_ids) && length(low_ids) > 0L) {
      cat(sprintf(
        "  Low-ESS subjects (ESS/K below threshold): %d\n",
        length(low_ids)
      ))
    }
  }

  # BAYES posterior summary (method = "bayes", terminal stage)
  if (!is.null(x$bayes)) {
    b <- x$bayes
    cat("\n", .ferx_style("BAYES  (posterior summary)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    flag <- .ferx_bayes_rhat_flag(b$max_rhat)
    cat(sprintf(
      "  Chains: %d   Warmup: %d   Draws/chain: %d   Divergent: %d   Max R-hat: %.4f%s\n",
      as.integer(b$n_chains), as.integer(b$n_warmup),
      as.integer(b$n_draws_per_chain), as.integer(b$n_divergent),
      b$max_rhat, flag
    ))
    if (!is.null(b$param_names) && length(b$param_names) > 0L) {
      # Assemble the parallel posterior vectors into one table so the column
      # set is declared in a single place.
      tbl <- data.frame(
        param = b$param_names,
        mean = b$mean, sd = b$sd, q025 = b$q025, q975 = b$q975,
        rhat = b$rhat, ess = b$ess_bulk,
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "  %-14s %11s %10s %11s %11s %7s %8s\n",
        "Parameter", "Mean", "SD", "2.5%", "97.5%", "R-hat", "ESS"
      ))
      for (i in seq_len(nrow(tbl))) {
        cat(sprintf(
          "  %-14s %11.4g %10.4g %11.4g %11.4g %7.3f %8.0f\n",
          tbl$param[i], tbl$mean[i], tbl$sd[i],
          tbl$q025[i], tbl$q975[i], tbl$rhat[i], tbl$ess[i]
        ))
      }
    }
  }

  # SHRINKAGE - compact single-line layout: ETA_CL: 8.2%  ETA_V: 12.1%  ETA_KA: 34.7% [!]
  # The [!] flag (yellow when colour available) highlights values above 30%,
  # consistent with the %RSE highlighter in the THETA table.
  has_shrinkage <- (!is.null(x$shrinkage_eta) && any(!is.na(x$shrinkage_eta))) ||
    (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps))
  if (has_shrinkage) {
    cat("\n", .ferx_style("SHRINKAGE", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    parts <- character(0)
    if (!is.null(x$shrinkage_eta)) {
      for (k in seq_along(x$shrinkage_eta)) {
        sh  <- x$shrinkage_eta[k]
        if (is.na(sh)) next
        lbl <- if (!is.null(x$eta_names) && length(x$eta_names) >= k && nzchar(x$eta_names[k])) x$eta_names[k] else sprintf("ETA%d", k)
        val <- sprintf("%.1f%%", sh * 100)
        flag <- if (sh * 100 > 30) paste0(" ", .ferx_style("[!]", "yellow")) else ""
        parts <- c(parts, sprintf("%s: %s%s", lbl, val, flag))
      }
    }
    if (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps)) {
      val <- sprintf("%.1f%%", x$shrinkage_eps * 100)
      flag <- if (x$shrinkage_eps * 100 > 30) paste0(" ", .ferx_style("[!]", "yellow")) else ""
      parts <- c(parts, sprintf("EPS: %s%s", val, flag))
    }
    if (length(parts) > 0L) cat(" ", paste(parts, collapse = "   "), "\n", sep = "")
  }

  # DIAGNOSTICS - compact: Covariance + condition number + DW on one line where possible.
  cov_status <- if (!is.null(x$covariance_status)) x$covariance_status else "unknown"
  cov_str <- switch(cov_status,
    computed      = "computed",
    failed        = .ferx_style("FAILED", "red"),
    not_requested = "not requested",
    cov_status
  )
  cond_part <- if (!is.null(x$condition_number)) {
    if (is.infinite(x$condition_number)) "Cond: Inf" else sprintf("Cond: %.1f", x$condition_number)
  } else NULL
  dw <- x$dw_statistic
  dw_part <- if (!is.null(dw) && !is.na(dw)) {
    sprintf("DW: %.2f [%s]", dw, .dw_label(dw))
  } else NULL
  iwres_part <- if (!is.null(x$iwres_lag1_r) && !is.na(x$iwres_lag1_r)) {
    sprintf("IWRES lag-1 r: %.3f", x$iwres_lag1_r)
  } else NULL
  diag_parts <- c(sprintf("Covariance: %s", cov_str), cond_part, dw_part, iwres_part)
  cat("\n", .ferx_style("DIAGNOSTICS", "bold"), "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  cat(" ", paste(diag_parts, collapse = "   "), "\n", sep = "")

  # Run info - tighter; only show fields with useful content.
  cat("\n", .ferx_style("RUN INFO", "bold"), "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  if (!is.null(x$gradient)) {
    cat(sprintf(" Gradient (requested): %s", x$gradient))
    if (!is.null(x$gradient_used) && !is.na(x$gradient_used)) {
      cat(sprintf("   (used: %s)", x$gradient_used))
    }
    cat("\n")
  }
  if (!is.null(x$ferx_version)) {
    cat(" ferx v", as.character(utils::packageVersion("ferx")),
        " (core v", x$ferx_version, ")\n", sep = "")
  }

  mfs <- x$model_file_settings %||% list()
  cs  <- x$call_settings %||% list()
  if (length(mfs) > 0L || length(cs) > 0L) {
    cat("\n", .ferx_style("SETTINGS  (model file / call-time override)", "bold"), "\n", sep = "")
    cat(strrep("-", 60), "\n", sep = "")
    all_keys <- union(names(mfs), names(cs))
    for (nm in all_keys) {
      mval <- mfs[[nm]]
      cval <- cs[[nm]]
      if (!is.null(mval) && !is.null(cval) &&
          !identical(tolower(as.character(mval)), tolower(as.character(cval)))) {
        cat(sprintf("  %-28s %s  [model: %s] *\n", nm, cval, mval))
      } else if (!is.null(cval)) {
        cat(sprintf("  %-28s %s\n", nm, cval))
      } else {
        cat(sprintf("  %-28s %s  [model only]\n", nm, mval))
      }
    }
    if (any(vapply(all_keys, function(nm) {
      mval <- mfs[[nm]]
      cval <- cs[[nm]]
      !is.null(mval) && !is.null(cval) &&
        !identical(tolower(as.character(mval)), tolower(as.character(cval)))
    }, logical(1L)))) {
      cat("  (* model file value overridden by call-time argument)\n")
    }
  }

  # Warning summary - a compact tally and a call-to-action rather than a wall
  # of message strings. Severity comes from fit$warnings_structured (PR 2);
  # the legacy flat fit$warnings vector is used only as a count fallback.
  ws <- x$warnings_structured
  if (!is.null(ws) && is.data.frame(ws) && nrow(ws) > 0L) {
    n_crit <- sum(ws$severity == "critical")
    n_warn <- sum(ws$severity == "warning")
    n_info <- sum(ws$severity == "info")
    cat("\n", strrep("-", 60), "\n", sep = "")
    parts <- character(0)
    if (n_crit > 0) parts <- c(parts, .ferx_style(sprintf("%d critical", n_crit), "red"))
    if (n_warn > 0) parts <- c(parts, .ferx_style(sprintf("%d warning", n_warn), "yellow"))
    if (n_info > 0) parts <- c(parts, .ferx_style(sprintf("%d info", n_info), "dim"))
    cat(" ", paste(parts, collapse = "   "),
        "  --  call ", .ferx_style("ferx_get_warnings(fit)", "bold"),
        " for details\n", sep = "")
  } else if (length(x$warnings) > 0L) {
    cat("\n", strrep("-", 60), "\n", sep = "")
    cat(" ", length(x$warnings), " warning(s)  --  inspect fit$warnings\n", sep = "")
  }

  cat(bar, "\n", sep = "")
  invisible(x)
}

#' Summarize a ferx fit result
#'
#' Returns a compact list with the most-used diagnostic fields: model name,
#' dataset, estimation method (including the full chain when multiple methods
#' were chained), OFV/AIC/BIC, parameter estimates, standard errors, shrinkage,
#' effective settings, SIR ESS (when available), warnings, and run metadata.
#' Printing the returned object via \code{print()} or auto-printing formats
#' all fields in a human-readable layout.
#'
#' @param object A `ferx_fit` object returned by \code{\link{ferx_fit}}.
#' @param ... Ignored.
#' @return A `ferx_summary` list (invisibly). Fields include: \code{model_name},
#'   \code{data_name}, \code{method}, \code{method_chain}, \code{gradient},
#'   \code{gradient_used},
#'   \code{converged}, \code{ofv}, \code{aic}, \code{bic}, \code{n_subjects},
#'   \code{n_obs}, \code{n_parameters}, \code{n_iterations}, \code{theta},
#'   \code{se_theta}, \code{omega}, \code{se_omega}, \code{sigma},
#'   \code{se_sigma}, \code{shrinkage_eta}, \code{shrinkage_eps},
#'   \code{covariance_status}, \code{eigenvalues}, \code{condition_number}
#'   (see \code{\link{ferx_fit}} for definitions), \code{wall_time_secs},
#'   \code{ferx_version}, \code{call_settings}, \code{model_file_settings},
#'   \code{sir_ess}, \code{importance_sampling}, \code{warnings},
#'   \code{ebe_convergence_warnings}, \code{max_unconverged_subjects},
#'   \code{total_ebe_fallbacks}, \code{dw_statistic}, \code{iwres_lag1_r}.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' summary(fit)
#'
#' # Chained estimation: SAEM warm-start followed by FOCEI refinement
#' fit_chain <- ferx_fit(ex$model, ex$data, method = c("saem", "focei"))
#' summary(fit_chain)
#'
#' # With custom settings (shown in the Settings block)
#' fit_custom <- ferx_fit(ex$model, ex$data,
#'   settings = list(optimizer = "slsqp", maxiter = 200L))
#' summary(fit_custom)
#' }
#' @family diagnostics
#' @export
summary.ferx_fit <- function(object, ...) {
  x <- object
  s <- list(
    model_name = x$model_name %||% NA_character_,
    data_name = x$data_name %||% NA_character_,
    gradient = x$gradient %||% NA_character_,
    gradient_used = x$gradient_used %||% NA_character_,
    method = x$method,
    method_chain = x$method_chain,
    converged = x$converged,
    ofv = x$ofv,
    aic = x$aic,
    bic = x$bic,
    n_subjects = x$n_subjects,
    n_obs = x$n_obs,
    n_parameters = x$n_parameters,
    n_iterations = x$n_iterations,
    theta = x$theta,
    se_theta = x$se_theta,
    omega = x$omega,
    se_omega = x$se_omega,
    sigma = x$sigma,
    se_sigma = x$se_sigma,
    shrinkage_eta = x$shrinkage_eta,
    shrinkage_eps = x$shrinkage_eps,
    covariance_status = x$covariance_status,
    eigenvalues = x$eigenvalues,
    condition_number = x$condition_number,
    wall_time_secs = x$wall_time_secs,
    ferx_version = x$ferx_version,
    ebe_convergence_warnings = x$ebe_convergence_warnings,
    max_unconverged_subjects = x$max_unconverged_subjects,
    total_ebe_fallbacks = x$total_ebe_fallbacks,
    model_structure = x$model_structure,
    call_settings = x$call_settings %||% list(),
    model_file_settings = x$model_file_settings %||% list(),
    sir_ess = x$sir_ess,
    importance_sampling = x$importance_sampling,
    bayes = x$bayes,
    warnings = x$warnings,
    uses_sde = isTRUE(x$uses_sde),
    dw_statistic = x$dw_statistic,
    iwres_lag1_r = x$iwres_lag1_r,
    saem_n_subjects_hmc = x$saem_n_subjects_hmc
  )
  class(s) <- "ferx_summary"
  s
}

#' @export
print.ferx_summary <- function(x, ...) {
  bar <- strrep("-", 50)
  cat(bar, "\n", sep = "")
  cat("ferx NLME Fit Summary\n")
  cat(bar, "\n", sep = "")

  cat(sprintf("Model:     %s\n", x$model_name %||% "?"))
  cat(sprintf("Dataset:   %s\n", x$data_name %||% "?"))
  cat(sprintf("Converged: %s\n", if (isTRUE(x$converged)) "YES" else "NO"))

  sde_tag <- if (isTRUE(x$uses_sde)) " + SDE (EKF)" else ""
  if (!is.null(x$method_chain) && length(x$method_chain) > 1L) {
    cat("Method:    ", paste(toupper(x$method_chain), collapse = " -> "), sde_tag, "\n", sep = "")
  } else {
    cat(sprintf("Method:    %s%s\n", toupper(x$method %||% "?"), sde_tag))
  }
  used <- x$gradient_used
  if (is.null(used) || (length(used) == 1L && is.na(used))) {
    cat(sprintf("Gradient:  %s (requested)\n", x$gradient %||% "?"))
  } else {
    cat(sprintf("Gradient:  %s (requested) -> %s (used)\n",
                x$gradient %||% "?", used))
  }
  cat(sprintf("ferx v%s (core v%s)\n",
              as.character(utils::packageVersion("ferx")),
              x$ferx_version %||% "?"))

  if (!is.null(x$bayes)) {
    b <- x$bayes
    cat(sprintf(
      "Bayes:     %d chains, %d draws/chain, max R-hat = %.4f%s\n",
      as.integer(b$n_chains), as.integer(b$n_draws_per_chain), b$max_rhat,
      .ferx_bayes_rhat_flag(b$max_rhat)
    ))
  }

  if (length(x$model_file_settings) > 0L || length(x$call_settings) > 0L) {
    cat("\nSettings (model file / call-time override):\n")
    all_keys <- union(names(x$model_file_settings), names(x$call_settings))
    for (nm in all_keys) {
      mval <- x$model_file_settings[[nm]]
      cval <- x$call_settings[[nm]]
      if (!is.null(mval) && !is.null(cval) &&
          !identical(tolower(as.character(mval)), tolower(as.character(cval)))) {
        cat(sprintf("  %-28s %s  [model: %s] *\n", nm, cval, mval))
      } else if (!is.null(cval)) {
        cat(sprintf("  %-28s %s\n", nm, cval))
      } else {
        cat(sprintf("  %-28s %s  [model only]\n", nm, mval))
      }
    }
    if (any(vapply(all_keys, function(nm) {
      mval <- x$model_file_settings[[nm]]
      cval <- x$call_settings[[nm]]
      !is.null(mval) && !is.null(cval) &&
        !identical(tolower(as.character(mval)), tolower(as.character(cval)))
    }, logical(1L)))) {
      cat("  (* model file value overridden by call-time argument)\n")
    }
  }

  if (!is.null(x$model_structure)) {
    ms <- x$model_structure
    iiv_str <- if (length(ms$iiv) > 0L) paste(ms$iiv, collapse = ", ") else "none"
    # Via the shared helper so a sample-size-weighted kappa (ferx-core #1031)
    # is annotated here exactly as print.ferx_fit() annotates it.
    iov_lbl <- .ferx_iov_labels(ms)
    iov_str <- if (length(iov_lbl) > 0L) paste(iov_lbl, collapse = ", ") else "none"
    cat(sprintf("\nStructure: %s  |  IIV: %s  |  IOV: %s  |  Residual: %s\n",
                .ferx_format_structural(ms), iiv_str, iov_str, ms$residual))
  }

  cat(sprintf("\nOFV: %.4f  AIC: %.4f  BIC: %.4f\n", x$ofv, x$aic, x$bic))
  cat(sprintf(
    "Subjects: %d  Obs: %d  Params: %s  Iter: %d\n",
    x$n_subjects, x$n_obs,
    if (is.null(x$n_parameters) || is.na(x$n_parameters)) "?" else format(x$n_parameters),
    x$n_iterations
  ))

  if (!is.null(x$shrinkage_eta) && any(!is.na(x$shrinkage_eta))) {
    sh_lbls <- if (!is.null(x$eta_names) && length(x$eta_names) == length(x$shrinkage_eta)) x$eta_names else sprintf("ETA%d", seq_along(x$shrinkage_eta))
    sh_str <- paste(sprintf("%s=%.1f%%", sh_lbls, x$shrinkage_eta * 100), collapse = "  ")
    cat("Shrinkage:", sh_str, "\n")
  }
  if (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps)) {
    cat(sprintf("EPS shrinkage: %.1f%%\n", x$shrinkage_eps * 100))
  }

  if (!is.null(x$dw_statistic) && !is.na(x$dw_statistic)) {
    cat(sprintf("DW: %.2f [%s]\n", x$dw_statistic, .dw_label(x$dw_statistic)))
    if (!is.null(x$iwres_lag1_r) && !is.na(x$iwres_lag1_r)) {
      cat(sprintf("lag-1 r: %.3f\n", x$iwres_lag1_r))
    }
  }

  if (!is.null(x$sir_ess)) {
    cat(sprintf("SIR ESS:   %.1f\n", x$sir_ess))
  }
  if (!is.null(x$importance_sampling)) {
    cat(sprintf(
      "IS -2logL: %.4f  (MC SE %.4f, ESS/K min %.2f)\n",
      x$importance_sampling$minus2_log_likelihood,
      x$importance_sampling$mc_standard_error,
      x$importance_sampling$ess_min
    ))
  }

  cov_str <- switch(x$covariance_status %||% "unknown",
    computed = "computed",
    failed = "FAILED",
    not_requested = "not requested",
    x$covariance_status
  )
  cond_str <- if (!is.null(x$condition_number)) {
    if (is.infinite(x$condition_number)) "  Cond: Inf" else sprintf("  Cond: %.1f", x$condition_number)
  } else ""
  wall <- if (!is.null(x$wall_time_secs)) sprintf("%.1fs", x$wall_time_secs) else "?"
  cat(sprintf("Covariance: %s%s  |  Wall time: %s\n", cov_str, cond_str, wall))

  warn_msgs <- x$warnings
  if (!is.null(x$ebe_convergence_warnings) && x$ebe_convergence_warnings > 0L) {
    warn_msgs <- c(warn_msgs, sprintf(
      "%d EBE convergence warning(s) (max unconverged subjects: %s)",
      x$ebe_convergence_warnings,
      x$max_unconverged_subjects %||% "?"
    ))
  }
  if (!is.null(x$total_ebe_fallbacks) && x$total_ebe_fallbacks > 0L) {
    warn_msgs <- c(warn_msgs, sprintf("%d EBE fallback(s)", x$total_ebe_fallbacks))
  }
  if (length(warn_msgs) > 0L) {
    cat("\nWarnings:\n")
    for (w in warn_msgs) cat("  *", w, "\n")
  }

  cat(bar, "\n", sep = "")
  invisible(x)
}

# Emit a warning() for each model file [fit_options] key that disagrees with
# an explicit R call-time value. `dedicated_explicit` is a named list keyed by
# R argument name; entries are NULL for args the user did not pass explicitly
# (those silently defer to the model file). `settings_parts` is the output of
# .ferx_settings_to_strings(): every key in it was supplied by the caller.
.ferx_warn_fit_option_conflicts <- function(model_file_opts,
                                            dedicated_explicit,
                                            settings_parts) {
  if (length(model_file_opts) == 0L) return(invisible(NULL))

  warn <- function(key, model_val, call_val) {
    call_str  <- as.character(call_val)
    model_str <- as.character(model_val)
    if (!identical(tolower(model_str), tolower(call_str))) {
      warning(
        "Model file [fit_options] sets `", key, " = ", model_str,
        "` but ferx_fit() argument overrides it with `", call_str, "`.",
        " The call-time value will be used.",
        call. = FALSE
      )
    }
  }

  # Aliases the model file may use for each dedicated R argument.
  key_aliases <- list(
    method         = "method",
    covariance     = "covariance",
    bloq_method    = c("bloq_method", "bloq"),
    threads        = "threads",
    mu_referencing = "mu_referencing",
    sir            = "sir",
    gradient       = c("gradient", "gradient_method"),
    inits_from_nca = "inits_from_nca"
  )
  mf_names <- names(model_file_opts)
  for (arg in names(dedicated_explicit)) {
    call_val <- dedicated_explicit[[arg]]
    if (is.null(call_val)) next
    for (alias in key_aliases[[arg]]) {
      if (alias %in% mf_names) warn(alias, model_file_opts[[alias]], call_val)
    }
  }

  # Settings list keys are always explicit by construction.
  for (i in seq_along(settings_parts$keys)) {
    mkey <- tolower(settings_parts$keys[[i]])
    if (mkey %in% mf_names) {
      warn(mkey, model_file_opts[[mkey]], settings_parts$values[[i]])
    }
  }
  invisible(NULL)
}

# Reverse of the stringify step in .ferx_settings_to_strings: parse a single
# string value back to its natural R type for human-readable summary display.
.ferx_parse_setting_value <- function(v) {
  if (v == "true") return(TRUE)
  if (v == "false") return(FALSE)
  if (v == "null") return(NA)
  num <- suppressWarnings(as.numeric(v))
  if (!is.na(num)) return(num)
  v
}

# Convert a named-list of settings into two parallel character vectors for the
# Rust FFI. Each value is stringified in a Rust-parser-friendly form (logicals
# ? "true"/"false", numerics with full precision, NULL/NA ? "null"). Strict
# validation happens in Rust (apply_fit_option); here we only enforce shape.
.ferx_settings_to_strings <- function(settings) {
  if (is.null(settings)) {
    return(list(keys = character(0), values = character(0)))
  }
  if (!is.list(settings)) {
    stop("`settings` must be NULL or a named list")
  }
  if (length(settings) == 0L) {
    return(list(keys = character(0), values = character(0)))
  }
  nms <- names(settings)
  if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms) != 0L) {
    stop("`settings` must be a uniquely-named list (all entries must have a non-empty name)")
  }
  stringify <- function(key, v) {
    if (is.null(v) || (length(v) == 1L && is.na(v))) {
      return("null")
    }
    if (length(v) != 1L) {
      stop(sprintf("`settings$%s` must be a length-1 scalar", key))
    }
    if (is.logical(v)) {
      return(if (isTRUE(v)) "true" else "false")
    }
    if (is.numeric(v)) {
      if (!is.finite(v)) {
        stop(sprintf("`settings$%s` must be a finite number", key))
      }
      return(format(v, scientific = FALSE, trim = TRUE, digits = 17))
    }
    if (is.character(v)) {
      return(v)
    }
    stop(sprintf("`settings$%s` has unsupported type `%s`", key, class(v)[1L]))
  }
  values <- vapply(
    seq_along(settings),
    function(i) stringify(nms[i], settings[[i]]),
    character(1L)
  )
  list(keys = nms, values = values)
}

# Parse the [fit_options] block of a .ferx file into a named character vector.
# Keys are lower-cased; values are the raw trimmed strings from the file
# (before Rust type coercion). Returns an empty character(0) when the section
# is absent or empty. Uses .ferx_extract_blocks() so comments are stripped.
.ferx_parse_model_fit_options <- function(path) {
  blocks <- .ferx_extract_blocks(path)
  lines  <- blocks[["fit_options"]] %||% character(0)
  if (length(lines) == 0L) return(setNames(character(0), character(0)))
  parts <- lapply(strsplit(lines, "=", fixed = TRUE), function(x) {
    if (length(x) < 2L) return(NULL)
    list(key = tolower(trimws(x[[1L]])),
         val = trimws(paste(x[-1L], collapse = "=")))
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) return(setNames(character(0), character(0)))
  keys <- vapply(parts, `[[`, character(1L), "key")
  vals <- vapply(parts, `[[`, character(1L), "val")
  setNames(vals, keys)
}
