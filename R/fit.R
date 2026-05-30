#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend
#' powered by automatic differentiation (Enzyme).
#'
#' @param model Path to a \code{.ferx} model file, or a \code{\link{ferx_model}} object.
#' @param data Path to a NONMEM-format CSV file. Required columns: ID, TIME,
#'   DV, EVID, AMT, CMT. Optional columns recognised by the engine: RATE
#'   (infusion rate), MDV (missing-DV flag), II (dosing interval, required when
#'   SS > 0), SS (steady-state flag: 1 = pre-dose at steady state, 2 = add SS
#'   concentration to current state), CENS (BLOQ flag for M3 method), OCC
#'   (occasion index for IOV), and any covariate columns referenced in the
#'   model. See the steady-state section below for SS/II details.
#' @param method Estimation method(s). Either a single string or a character
#'   vector of methods to run in sequence. Each stage is seeded with the
#'   previous stage's converged parameters, and only the final stage produces
#'   the reported covariance/diagnostics. Supported methods: \code{"foce"},
#'   \code{"focei"}, \code{"saem"}, \code{"gn"} (Gauss-Newton / BHHH),
#'   \code{"gn_hybrid"} (Gauss-Newton followed by a FOCEI polish step), or
#'   \code{"imp"} (also accepted as \code{"importance_sampling"} or
#'   \code{"importance-sampling"}; importance-sampling marginal
#'   log-likelihood). The \code{"imp"} stage is a diagnostic terminal
#'   stage that does not update parameters and must be the *last* entry
#'   of the chain, e.g. \code{c("focei", "imp")} or
#'   \code{c("saem", "imp")}. Example chain: \code{c("saem", "focei")}.
#'   SAEM fully supports inter-occasion variability (IOV / kappa) models.
#' @param covariance Logical; compute the covariance step for standard errors
#' @param verbose Logical; print progress during estimation
#' @param bloq_method Handling of observations below the lower limit of
#'   quantification. \code{NULL} (default) keeps whatever the model file
#'   specified; \code{"m3"} enables Beal's M3 likelihood (requires a
#'   \code{CENS} column in the data, with \code{DV} carrying the LLOQ value
#'   on \code{CENS=1} rows); \code{"drop"} disables M3.
#' @param threads Number of worker threads for the per-subject parallel loops
#'   in the Rust backend (inner EBE search, SAEM, SIR). \code{NULL} (default)
#'   leaves the backend's thread pool at its default of one worker per logical
#'   CPU. Pass an integer to cap the pool - useful on shared machines, in CI,
#'   or to avoid SMT-induced contention (try
#'   \code{parallel::detectCores(logical = FALSE)}). The setting is per-call,
#'   so successive fits in the same R session can use different values.
#' @param mu_referencing Logical. If \code{TRUE} (default), automatically
#'   detects mu-referencing from the model structure for faster and more
#'   accurate convergence. Applies to all estimation methods. Set to
#'   \code{FALSE} to disable for comparison purposes. Detection works
#'   automatically for standard parameterizations such as
#'   \code{PARAM = THETA * exp(ETA)}; unusual parameterizations fall back
#'   silently to zero-centred ETA initialisation with no error. No changes
#'   to the \code{.ferx} model file are needed. Check \code{fit$warnings}
#'   to see which ETAs were detected.
#' @param gradient Inner-loop (per-subject EBE) gradient method.
#'   One of \code{"auto"} (default), \code{"ad"}, or \code{"fd"}.
#'
#'   The inner optimizer is BFGS; what differs is how the gradient of the
#'   individual NLL w.r.t. ETA is computed:
#'   \itemize{
#'     \item \code{"ad"}: reverse-mode automatic differentiation via
#'       Enzyme. One forward + one reverse pass per gradient call,
#'       regardless of the number of etas. Requires the package to have
#'       been compiled with autodiff support (see
#'       \code{ferx_rust_autodiff_enabled()}) and the model to have an
#'       analytical PK path. ODE-based models currently have no AD
#'       implementation and will silently fall back to FD even if
#'       \code{"ad"} is requested.
#'     \item \code{"fd"}: central finite differences, \code{2 * n_eta}
#'       forward evaluations per gradient call. Always available.
#'     \item \code{"auto"} (default): use AD whenever the package was
#'       compiled with it and the model has an analytical PK path,
#'       otherwise FD. This is the right choice for almost all uses.
#'   }
#'
#'   \strong{When AD wins.} AD's per-gradient cost is roughly independent
#'   of the number of etas, while FD scales linearly. For typical
#'   analytical PK fits we have measured AD ~5x faster per gradient call
#'   on 1-cpt oral (\code{n_eta = 3}) and ~1.5-1.8x faster per call on 2-
#'   and 3-cpt infusion. The advantage grows with more random effects,
#'   more observations per subject, and (when implemented) ODE forward
#'   models. On small problems the \emph{wall-clock} gap is often small
#'   because non-gradient work (NLopt steps, population likelihood
#'   reductions, parallel scheduling) dominates the total fit time.
#'
#'   \strong{When to set \code{"fd"} explicitly.} Mainly for validation
#'   (cross-check that AD and FD converge to the same OFV on a new model),
#'   for reproducibility against a known FD baseline, or to sidestep an
#'   Enzyme codegen pathology on an unusual model structure. Both methods
#'   converge to the same optimum within line-search tolerance on
#'   well-conditioned problems.
#'
#'   Set \code{FERX_TIME_GRADIENTS=1} in the environment to print
#'   per-gradient-call timings at the end of a fit, which is the easiest
#'   way to check which method is faster on your specific model/data.
#' @param sir Logical; run Sampling Importance Resampling after the fit to
#'   produce non-parametric parameter uncertainty intervals. Requires
#'   \code{covariance = TRUE}. Tuning knobs (\code{sir_samples},
#'   \code{sir_resamples}, \code{sir_seed}) still flow through
#'   \code{settings}. Default \code{FALSE}.
#' @param optimizer_trace Logical. If \code{TRUE}, write a per-iteration CSV
#'   trace to a temporary file and store its path in \code{fit$trace_path}.
#'   Pass the result to \code{\link{ferx_trace}} or
#'   \code{\link{ferx_plot_trace}} to inspect optimizer progress. Default
#'   \code{FALSE}.
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
#' @param settings Optional named list of estimation-method-specific options
#'   forwarded to the Rust \code{FitOptions}. Use this to tune knobs that do
#'   not have a dedicated \code{ferx_fit()} argument, without needing a new
#'   wrapper release for each option. Recognized keys include SAEM
#'   (\code{n_exploration}, \code{n_convergence}, \code{n_mh_steps},
#'   \code{adapt_interval}, \code{omega_burnin}: initial iterations during
#'   which Omega is held fixed while the sampler warms up (default 20,
#'   guards against Omega collapse on sparse data); \code{seed}; and
#'   \code{n_leapfrog} /
#'   \code{saem_n_leapfrog}: integer leapfrog steps for HMC proposals,
#'   default 0 = Metropolis-Hastings; positive value uses one HMC proposal
#'   per subject per iteration, requires AD build, see examples below),
#'   SIR tuning (\code{sir_samples},
#'   \code{sir_resamples}, \code{sir_seed}), Importance Sampling
#'   (\code{is_samples}, \code{is_proposal_df}, \code{is_seed},
#'   \code{is_low_ess_threshold}; used only when an \code{"imp"} stage is
#'   in the method chain), Gauss-Newton (\code{gn_lambda}),
#'   optimizer selection (\code{optimizer} - one of \code{"bobyqa"} (default),
#'   \code{"slsqp"}, \code{"lbfgs"}, \code{"nlopt_lbfgs"}, \code{"mma"},
#'   \code{"bfgs"}, \code{"trust_region"}; also \code{global_search},
#'   \code{global_maxeval}), the outer optimizer iteration cap
#'   (\code{maxiter}), the inner (per-subject EBE) loop
#'   (\code{inner_maxiter}, \code{inner_tol}), the Steihaug CG budget
#'   for \code{optimizer = "trust_region"} (\code{steihaug_max_iters}),
#'   the NLopt outer-loop stagnation guard
#'   (\code{stagnation_guard} - logical, default \code{TRUE}; set
#'   \code{FALSE} to let SLSQP / L-BFGS run to their own xtol/ftol or to
#'   \code{maxiter} instead of short-circuiting on a numerically-flat
#'   OFV plateau; only consumed by FOCE / FOCEI / GN-hybrid),
#'   the reconverged-gradient schedule
#'   (\code{reconverge_gradient_interval} - integer, default \code{0}; how
#'   often to re-solve each subject's inner EBE loop during the population
#'   gradient instead of holding it fixed. \code{0} never reconverges (the
#'   cheap fixed-EBE gradient); \code{1} reconverges every gradient
#'   evaluation; \code{N} reconverges every \code{N}-th. The fixed-EBE
#'   gradient can stall \code{slsqp} above the derivative-free
#'   (\code{bobyqa}) optimum on ill-conditioned non-IOV fits; reconverging
#'   recovers the full surface at ~5-6x the per-gradient cost. IOV models
#'   always reconverge and ignore this; only consumed by FOCE / FOCEI /
#'   GN-hybrid),
#'   and multi-start optimization (\code{n_starts}, \code{start_sigma},
#'   \code{multi_start_seed}).
#'   Values that duplicate a dedicated argument
#'   (\code{method}, \code{covariance}, \code{verbose},
#'   \code{bloq_method}, \code{threads}, \code{sir}) are rejected - pass
#'   them via the dedicated argument. Unknown keys and malformed values
#'   also raise an error.
#'   \strong{Precedence}: dedicated \code{ferx_fit()} arguments win over
#'   \code{settings}, which in turn win over the model file's
#'   \code{[fit_options]} block. A \code{warning()} is issued whenever a
#'   call-time value overrides a different value from \code{[fit_options]}.
#'   Inspect \code{fit$model_file_settings} and \code{fit$call_settings} to
#'   audit the full picture.
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
#' @return A list with components:
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
#'   \item{se_omega}{Standard errors for omega diagonal}
#'   \item{se_sigma}{Standard errors for sigma (on the SD scale, like
#'     \code{sigma} itself)}
#'   \item{sdtab}{Data frame with ID, TIME, DV, PRED, IPRED, CWRES, IWRES,
#'     EBE_OFV, N_OBS; OCC if any subject carries an occasion column; CENS
#'     if any rows are below LLOQ. Per-subject ETAs live in
#'     \code{ebe_etas}; per-subject parameter values in
#'     \code{individual_estimates}.}
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
#'     and ETA-normality checks. Inspect with \code{\link{ferx_warnings}}.}
#'   \item{sir_ess}{SIR effective sample size (NULL if SIR not run)}
#'   \item{sir_ci_theta, sir_ci_omega, sir_ci_sigma}{SIR 95\% CI matrices
#'     with columns \code{lower} and \code{upper} (NULL if SIR not run)}
#'   \item{importance_sampling}{Importance-sampling marginal log-likelihood
#'     diagnostics, populated only when the method chain ends with an
#'     \code{"imp"} stage (\code{NULL} otherwise). A list with
#'     \code{minus2_log_likelihood} (the IS estimate of \eqn{-2 \log L},
#'     comparable across models in the same nesting family),
#'     \code{mc_standard_error} (Monte-Carlo SE; scales as
#'     \eqn{1/\sqrt{is\_samples}}), \code{n_samples},
#'     \code{proposal_df}, \code{ess_min} and \code{ess_median}
#'     (per-subject effective sample size as a fraction of
#'     \code{is_samples}), \code{kappa_treatment}
#'     (\code{"not_applicable"}, \code{"fixed_at_mode"} for partial-marginal
#'     IOV fits, or \code{"marginalized"}), and parallel vectors
#'     \code{low_ess_subject_ids} / \code{low_ess_subject_frac} flagging
#'     subjects whose ESS fraction fell below \code{is_low_ess_threshold}.}
#'   \item{trace_path}{Path to the optimizer trace CSV, or \code{NULL} when
#'     \code{optimizer_trace = FALSE}. Pass to \code{\link{ferx_trace}}
#'     or \code{\link{ferx_plot_trace}}.}
#'   \item{ebe_convergence_warnings}{Number of outer iterations in which too
#'     many EBEs were unconverged (step was rejected by the guard).}
#'   \item{max_unconverged_subjects}{Worst-case number of unconverged subjects
#'     in a single outer iteration.}
#'   \item{total_ebe_fallbacks}{Total Nelder-Mead fallback invocations across
#'     all subjects and iterations.}
#'   \item{covariance_status}{String: \code{"computed"}, \code{"failed"}, or
#'     \code{"not_requested"}.}
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
#'     caller (one of \code{"auto"}, \code{"ad"}, \code{"fd"}). The
#'     \emph{resolved} method may differ - see \code{gradient_used} and the
#'     \code{gradient} argument.}
#'   \item{gradient_used}{The inner-loop gradient method the engine actually
#'     used: \code{"ad"}, \code{"fd"}, or \code{"N/A"} (derivative-free /
#'     sampling step). When \code{gradient = "auto"} this records which
#'     branch resolved at fit time; for ODE models with \code{gradient = "ad"}
#'     this shows the silent fallback to \code{"fd"}. \code{NA} on older
#'     engine binaries that did not populate this field.}
#'   \item{gradient_method_inner}{Engine label for the inner-loop gradient
#'     ("Enzyme AD", "finite differences", or "N/A"). Prefer
#'     \code{gradient_used} for display.}
#'   \item{gradient_method_outer}{Engine label for the population-level
#'     optimizer's gradient ("Enzyme AD", "finite differences", or "N/A").}
#'   \item{ferx_version}{ferx-core library version string.}
#'   \item{cov_matrix}{Full parameter covariance matrix as a named numeric
#'     matrix (params ? params). Row/column names use declared variable names
#'     (\code{"TVCL"}, \code{"ETA_CL"}, \code{"EPS_PROP"}); fallback is
#'     \code{"OMEGA(1,1)"} / \code{"SIGMA(1)"}. \code{NULL} when covariance
#'     step was not run or failed. Use \code{\link{ferx_cor_matrix}} to
#'     inspect correlations.}
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
#'     \code{"1-cpt oral"}, \code{"ODE"}, or \code{NULL} when the structural
#'     form is not one of the known analytical PK families), \code{iiv}
#'     (omega names), \code{iov} (kappa names), \code{residual} (error type
#'     string). Use \code{\link{ferx_model_inspect}} to view this before or
#'     after fitting.}
#'   \item{shrinkage_kappa}{Shrinkage values for kappa EBEs (\code{NULL} if no IOV).}
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
#'     rows in \code{ferx_estimates()} and \code{print.ferx_fit()}.}
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
#'   \item{saem_n_subjects_hmc}{Integer count of subjects whose SAEM E-step
#'     proposals were generated by HMC (i.e. \code{n_leapfrog > 0} and the AD
#'     build is active). Subjects that fell back to MH (e.g. ODE models) are
#'     not counted. \code{NA_integer_} when the method is not SAEM, when
#'     \code{n_leapfrog = 0} (MH only), or when AD is unavailable.}
#'   \item{neural_networks}{Named list of sub-lists, one per \code{[covariate_nn]}
#'     block declared in the model file. Empty list when the \code{nn} cargo
#'     feature is off or no NN blocks are present. Each sub-list contains:
#'     \code{name} (block name), \code{shape} (integer vector of layer widths),
#'     \code{hidden_activation}, \code{output_activation}, \code{n_weights}
#'     (total weight + bias count), \code{weights_offset} (0-based start index
#'     into \code{fit$theta}), \code{input_names}, \code{output_names}, and
#'     \code{weights} (numeric vector of trained values).}
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
#'         \code{ferx_estimates()} treat them as regular thetas.
#'   \item Autodifferentiation (Enzyme AD) is incompatible with the EKF
#'         covariance propagation - the engine automatically falls back to
#'         finite differences and emits a warning.
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
#' \strong{Form A -- scalar divisor (safe with AD):}
#' \preformatted{
#' [scaling]
#'   obs_scale = 1000   # divide every predicted value by 1000
#' }
#'
#' \strong{Form B -- expression divisor (forces \code{gradient = fd}):}
#' \preformatted{
#' [scaling]
#'   obs_scale = V      # divide by the individual volume (theta/eta expression)
#' }
#' Expressions may reference thetas, etas, and \code{[individual_parameters]}
#' variables. AD gradients are not supported for Form B; add
#' \code{gradient = fd} to \code{[fit_options]} or the engine will error.
#'
#' \strong{Per-CMT scaling (Form A or B per compartment):}
#' \preformatted{
#' [scaling]
#'   obs_scale[CMT=1] = 1000   # CMT 1 in ng/mL, dose in ug
#'   obs_scale[CMT=2] = 1      # CMT 2 already in correct units
#' }
#'
#' \strong{Form C -- ODE readout expression (forces \code{gradient = fd}):}
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
#'   \item Forms B, per-CMT, and Form C force finite-difference gradients; add
#'         \code{gradient = fd} to \code{[fit_options]} explicitly, or the
#'         engine raises an informative error.
#' }
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
#' convenient when combined with \code{\link{ferx_set_section}} to modify model
#' options in the same chain (see Examples). The data path stored in the
#' \code{ferx_model} object can always be overridden by supplying \code{data}
#' explicitly to \code{ferx_fit()}.
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
#' }
#'
#' \strong{Parallelism} - cap the Rust thread pool:
#' \preformatted{
#' ferx_fit(m, d, threads = parallel::detectCores(logical = FALSE))
#' }
#'
#' \strong{BLOQ handling}:
#' \preformatted{
#' ferx_fit(m, d, bloq_method = "m3")    # Beal M3 likelihood (needs CENS column)
#' ferx_fit(m, d, bloq_method = "drop")  # discard BLOQ rows
#' }
#'
#' \strong{Gradient method} for the inner EBE loop:
#' \preformatted{
#' ferx_fit(m, d, gradient = "auto")  # default: AD when available, else FD
#' ferx_fit(m, d, gradient = "ad")    # force automatic differentiation
#' ferx_fit(m, d, gradient = "fd")    # force finite differences
#' }
#'
#' \strong{Optimizer trace} - write a per-iteration CSV for convergence
#' diagnostics:
#' \preformatted{
#' fit <- ferx_fit(m, d, optimizer_trace = TRUE)
#' ferx_plot_trace(fit)
#' }
#'
#' @section Fine-tuning with \code{settings}:
#'
#' The \code{settings} argument forwards a named list of low-level options to
#' the Rust backend without requiring a new wrapper release. Arguments with
#' dedicated parameters (e.g. \code{method}, \code{covariance}) cannot be
#' duplicated in \code{settings} - pass them via the named argument.
#'
#' \strong{Outer optimizer selection:}
#' \preformatted{
#' ferx_fit(m, d, settings = list(optimizer = "bobyqa"))      # default
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
#' # HMC proposals (requires AD build; falls back to MH otherwise):
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
#' \code{max_unconverged_frac} (numeric in \[0, 1\]) flags the fit as converged
#' even if up to this fraction of subjects failed their inner EBE loop, instead
#' of raising an error - useful for large datasets with a handful of difficult
#' subjects. \code{min_obs_for_convergence_check} (non-negative integer) is the
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
#' fit |> ferx_estimates()     # tidy data frame: theta / omega / sigma + SE / %RSE
#' fit |> ferx_cor_matrix()    # parameter correlation matrix (needs covariance = TRUE)
#' fit |> ferx_model_inspect() # model structure auto-derived by the engine
#' fit |> ferx_plot_trace()    # convergence trace (needs optimizer_trace = TRUE)
#' }
#'
#' Diagnostics data frame (PRED, IPRED, CWRES, ETAs, ...) lives in
#' \code{fit$sdtab} and can be used directly:
#' \preformatted{
#' fit$sdtab
#' fit$ebe_etas
#' fit$individual_estimates
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
#' # BLOQ (M3 method)
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
#'   ferx_set_section("fit_options", c(
#'     "  method     = focei",
#'     "  maxiter    = 500",
#'     "  covariance = true"
#'   )) |>
#'   ferx_fit()
#'
#' # Full pipeline including post-fit outputs
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_set_section("fit_options", c("  method = focei")) |>
#'   ferx_fit(covariance = TRUE, threads = 4L,
#'            settings = list(optimizer = "slsqp")) |>
#'   summary()
#'
#' # Inspect then continue (ferx_get_section returns the ferx_model invisibly)
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_get_section("parameters") |>
#'   ferx_fit() |>
#'   ferx_estimates()
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
#' ferx_estimates(fit)       # tidy data frame with SE and %RSE
#' ferx_cor_matrix(fit)      # parameter correlation matrix
#' ferx_model_inspect(fit)   # model structure auto-derived by the engine
#' ferx_plot_trace(fit)      # convergence plot (optimizer_trace = TRUE required)
#'
#' fit$sdtab                 # per-observation diagnostics (PRED, IPRED, CWRES, etc.)
#' fit$ebe_etas              # per-subject empirical Bayes ETAs
#' fit$individual_estimates  # per-subject individual PK parameters
#' fit$eigenvalues           # sorted eigenvalues of parameter correlation matrix
#' fit$condition_number      # > 1000 flags potential ill-conditioning
#'
#' # Covariance diagnostics in a pipe
#' ferx_fit(ex$model, ex$data, covariance = TRUE) |> ferx_cor_matrix()
#'
#' # -- SDE model (Extended Kalman Filter) -------------------------------------
#'
#' # The [diffusion] block in the .ferx file enables SDE mode.
#' # DIFF_CENTRAL is a within-subject process-noise variance on the central
#' # compartment.  Autodiff is forced to finite differences automatically.
#' sde <- ferx_example("warfarin_sde")
#' fit_sde <- ferx_fit(sde$model, sde$data)
#' fit_sde$uses_sde                  # TRUE
#' fit_sde$theta["DIFF_CENTRAL"]     # fitted diffusion variance
#' ferx_estimates(fit_sde)           # DIFF_CENTRAL appears in theta block
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
#' #     default optimizer "bobyqa" is derivative-free and tolerates poor
#' #     starts the best of the three; switch to slsqp / trust_region for
#' #     gradient-based behaviour, then NCA-derived starts help most)
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
#' @family fitting
#' @export
ferx_fit <- function(model, data = NULL,
                     method = "focei",
                     covariance = TRUE,
                     verbose = TRUE,
                     bloq_method = NULL,
                     threads = NULL,
                     mu_referencing = TRUE,
                     sir = FALSE,
                     gradient = c("auto", "ad", "fd"),
                     optimizer_trace = FALSE,
                     scale_params = FALSE,
                     inits_from_nca = FALSE,
                     settings = NULL,
                     output = NULL,
                     include_data = FALSE,
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
    if (is.null(data)) {
      stop(
        "No data path supplied. Either pass `data` to ferx_fit() ",
        "or include it in ferx_model()."
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
  gradient <- match.arg(gradient)
  if (is.null(data)) {
    stop("`data` is required. Pass a path to a NONMEM CSV file.")
  }
  stopifnot(file.exists(model), file.exists(data))
  if (!is.logical(covariance) || length(covariance) != 1L || is.na(covariance)) {
    stop("`covariance` must be TRUE or FALSE")
  }
  if (!is.logical(mu_referencing) || length(mu_referencing) != 1L || is.na(mu_referencing)) {
    stop("`mu_referencing` must be TRUE or FALSE")
  }
  if (!is.logical(sir) || length(sir) != 1L || is.na(sir)) {
    stop("`sir` must be TRUE or FALSE")
  }
  if (sir && !covariance) {
    stop("`sir = TRUE` requires `covariance = TRUE`")
  }
  if (!is.character(method) || length(method) == 0L) {
    stop("`method` must be a non-empty character vector")
  }
  method <- vapply(
    method,
    function(m) {
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
      match.arg(
        normalised,
        c("foce", "focei", "saem", "gn", "gn_hybrid", "imp")
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
  # `imp` is a diagnostic terminal stage - engine rejects malformed chains too,
  # but surfacing the error R-side avoids a round-trip and gives a clearer
  # message anchored to the R argument.
  imp_positions <- which(method == "imp")
  if (length(imp_positions) > 0L) {
    if (length(imp_positions) > 1L) {
      stop(
        "`method` may include \"imp\" at most once; got ",
        length(imp_positions), " occurrences"
      )
    }
    if (length(method) == 1L) {
      stop(
        "`method = \"imp\"` is invalid: importance sampling is a diagnostic ",
        "terminal stage and must follow another method, e.g. ",
        "`c(\"focei\", \"imp\")` or `c(\"saem\", \"imp\")`"
      )
    }
    if (imp_positions != length(method)) {
      stop(
        "`\"imp\"` must be the final stage of the method chain; got ",
        "method = c(\"", paste(method, collapse = "\", \""), "\")"
      )
    }
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
  settings_parts <- .ferx_settings_to_strings(settings)
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
  # not an override of the model file.
  model_file_opts <- .ferx_parse_model_fit_options(model)
  explicit_args <- names(match.call())[-1L]
  dedicated_explicit <- list(
    method         = if ("method"         %in% explicit_args) paste(method, collapse = ", "),
    covariance     = if ("covariance"     %in% explicit_args) tolower(as.character(covariance)),
    bloq_method    = if ("bloq_method"    %in% explicit_args && !is.null(bloq_method)) bloq_arg,
    threads        = if ("threads"        %in% explicit_args && !is.null(threads)) as.character(threads),
    mu_referencing = if ("mu_referencing" %in% explicit_args) tolower(as.character(mu_referencing)),
    sir            = if ("sir"            %in% explicit_args) tolower(as.character(sir)),
    gradient       = if ("gradient"       %in% explicit_args) gradient,
    inits_from_nca = if ("inits_from_nca" %in% explicit_args) (if (is.null(inits_value)) "off" else inits_value)
  )
  .ferx_warn_fit_option_conflicts(model_file_opts, dedicated_explicit, settings_parts)

  fit_started_at <- Sys.time()
  raw <- ferx_rust_fit(
    model_path = normalizePath(model),
    data_path = normalizePath(data),
    method = method,
    covariance = covariance,
    verbose = verbose,
    bloq_method = bloq_arg,
    threads = threads_arg,
    mu_referencing = mu_referencing,
    sir = sir,
    gradient = gradient,
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

  # Reshape omega into a matrix
  if (length(result$omega) > 0 && !is.null(result$omega_dim)) {
    d <- result$omega_dim
    result$omega <- matrix(result$omega, nrow = d, ncol = d)
    eta_nms_om <- result$eta_names
    omega_dim_nms <- if (!is.null(eta_nms_om) && length(eta_nms_om) == d) eta_nms_om else paste0("OMEGA(", seq_len(d), ",", seq_len(d), ")")
    rownames(result$omega) <- colnames(result$omega) <- omega_dim_nms
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
  } else {
    .ferx_state$last_trace_path <- result$trace_path
    .ferx_state$last_trace_time <- fit_started_at
    .ferx_state$last_trace_model <- model
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

  # ETA normality (Shapiro-Wilk) - computed in R from per-subject EBEs
  result$eta_normality <- .ferx_compute_eta_normality(result$ebe_etas)
  # Push normality warnings into the warnings vector
  if (!is.null(result$eta_normality)) {
    for (i in seq_len(nrow(result$eta_normality))) {
      row <- result$eta_normality[i, ]
      if (nzchar(row$flag) && !is.na(row$p_val)) {
        result$warnings <- c(
          result$warnings,
          sprintf("%s Shapiro-Wilk p=%.4f - distribution may be non-normal", row$eta, row$p_val)
        )
      }
    }
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
  } else {
    result$omega_iov <- NULL
    result$se_kappa <- NULL
    result$shrinkage_kappa <- NULL
    result$kappa_names <- NULL
    result$ebe_kappas <- NULL
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

  # Print mu-referencing detections as informational messages
  mu_ref_warnings <- grep("mu-ref", result$warnings, value = TRUE)
  for (w in mu_ref_warnings) {
    eta_names <- sub("^mu-ref:\\s*", "", w)
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
  # "Enzyme AD", "finite differences", or "N/A"; map to the short tokens used
  # by the `gradient` argument so the two values can be compared directly
  # (e.g. requested "auto" -> used "ad"). When the engine omitted the field
  # (older binaries return ""), `gradient_used` is NA.
  result$gradient <- gradient
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

  class(result) <- "ferx_fit"

  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L || is.na(output)) {
      stop("`output` must be a single .fitrx file path or NULL")
    }
    ferx_save_fit(result, output, include_data = isTRUE(include_data))
  }

  result
}

# Assemble the structured-warning data frame for a fit result. Core supplies
# the severity/category/message triples it classified (parallel vectors on the
# raw FFI list); the R side appends only diagnostics it computes itself and
# that core never sees: the condition-number flag and the per-ETA Shapiro-Wilk
# normality flags. Returns a data frame with columns severity, category,
# message - always a data frame (zero rows when there are no warnings).
.ferx_assemble_structured_warnings <- function(raw, result) {
  df <- data.frame(
    severity      = as.character(raw$warnings_severity      %||% character(0)),
    category      = as.character(raw$warnings_category      %||% character(0)),
    message       = as.character(raw$warnings_message       %||% character(0)),
    source_method = as.character(raw$warnings_source_method %||% character(0)),
    stringsAsFactors = FALSE
  )
  extra <- list()
  # Condition number (NONMEM convention: > 1000 flags ill-conditioning).
  if (!is.null(result$condition_number) && is.finite(result$condition_number) &&
        result$condition_number > 1000) {
    extra[[length(extra) + 1L]] <- data.frame(
      severity      = "critical",
      category      = "condition_number",
      message       = sprintf(
        "High condition number (%.1f) -- parameter space may be ill-conditioned",
        result$condition_number
      ),
      source_method = "",
      stringsAsFactors = FALSE
    )
  }
  # Per-ETA Shapiro-Wilk normality flags (computed R-side from the EBEs).
  if (!is.null(result$eta_normality)) {
    for (i in seq_len(nrow(result$eta_normality))) {
      row <- result$eta_normality[i, ]
      if (!is.na(row$flag) && nzchar(row$flag) && !is.na(row$p_val)) {
        extra[[length(extra) + 1L]] <- data.frame(
          severity      = "warning",
          category      = "eta_normality",
          message       = sprintf(
            "%s Shapiro-Wilk p=%.4f - distribution may be non-normal",
            row$eta, row$p_val
          ),
          source_method = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(extra) > 0L) df <- rbind(df, do.call(rbind, extra))
  df
}

# Map the engine's gradient-kind label to the short tokens used by the
# `gradient` argument so requested/used can be compared at a glance.
# Returns NA_character_ when the engine field is missing or unrecognized so
# downstream display code can fall back gracefully on older binaries.
.ferx_short_gradient_label <- function(kind) {
  if (is.null(kind) || !nzchar(kind)) return(NA_character_)
  switch(kind,
    "Enzyme AD"          = "ad",
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
      # used by the [!] shrinkage flag and by ferx_warnings() severity colours.
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
  cat("\n", .ferx_style("OMEGA  (between-subject variability)", "bold"), "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")
  om <- x$omega
  if (is.null(dim(om))) om <- matrix(om, 1, 1)
  n_eta <- nrow(om)
  has_offdiag <- FALSE
  for (i in seq_len(n_eta)) {
    var_ii   <- om[i, i]
    eta_type <- if (!is.null(x$eta_param_types) && length(x$eta_param_types) >= i) x$eta_param_types[i] else "log_normal"
    linked_theta_name <- if (!is.null(x$eta_linked_theta) && length(x$eta_linked_theta) >= i) x$eta_linked_theta[i] else ""
    se_str <- if (!is.null(x$se_omega) && length(x$se_omega) >= i) {
      sprintf("%.6f", x$se_omega[i])
    } else {
      "N/A"
    }

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
        has_nms <- !is.null(x$eta_names) && length(x$eta_names) >= i
        if (has_nms) {
          lbl_i <- if (nzchar(x$eta_names[i])) x$eta_names[i] else sprintf("OMEGA(%d,%d)", i, i)
          lbl_j <- if (nzchar(x$eta_names[j])) x$eta_names[j] else sprintf("OMEGA(%d,%d)", j, j)
          cat(sprintf(
            "  %s ~ %s : cov = %.6f  (%s = %.4f)\n",
            lbl_i, lbl_j, cov_ij, corr_label, param_corr
          ))
        } else {
          cat(sprintf(
            "  OMEGA(%d,%d) : cov = %.6f  (%s = %.4f)\n",
            i, j, cov_ij, corr_label, param_corr
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
        "  --  call ", .ferx_style("ferx_warnings(fit)", "bold"),
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
#'   settings = list(optimizer = "slsqp", max_iter = 200L))
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
    warnings = x$warnings,
    uses_sde = isTRUE(x$uses_sde),
    dw_statistic = x$dw_statistic,
    iwres_lag1_r = x$iwres_lag1_r,
    saem_n_subjects_hmc = x$saem_n_subjects_hmc
  )
  class(s) <- "ferx_summary"
  s
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
    iov_str <- if (length(ms$iov) > 0L) paste(ms$iov, collapse = ", ") else "none"
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

#' Fit a model asynchronously
#'
#' Submits \code{\link{ferx_fit}} to a background process and returns
#' immediately, leaving the R session free. Call \code{\link{ferx_collect}}
#' on the returned handle to block-wait for the result and see live trace
#' progress.
#'
#' In RStudio the fit is submitted via \code{rstudioapi::jobRunScript()} and
#' appears in the Jobs pane. In Positron, plain R, or any other interactive
#' session it uses \code{callr::r_bg()}. Both backends return the same
#' \code{ferx_job} handle and the same \code{ferx_collect()} call retrieves
#' the result.
#'
#' The function forces \code{optimizer_trace = TRUE} internally so it has a
#' parseable progress channel. When you did not explicitly request the trace,
#' the temp CSV is deleted after \code{ferx_collect()} returns; pass
#' \code{optimizer_trace = TRUE} explicitly to keep it.
#'
#' @param model,data,... Forwarded to \code{\link{ferx_fit}}. See its
#'   documentation for the full argument list.
#' @param tail_n Integer. Number of recent trace rows to keep visible per
#'   poll in \code{\link{ferx_collect}}. Default 6.
#' @param poll_interval Numeric. Seconds between polls in
#'   \code{\link{ferx_collect}}. Default 0.5.
#'
#' @return A \code{ferx_job} handle. Pass it to \code{\link{ferx_collect}}
#'   to retrieve the \code{ferx_fit} result. In non-interactive sessions
#'   (\code{Rscript}, knitr, batch mode) the function falls back to a
#'   synchronous \code{\link{ferx_fit}} call and returns the result directly.
#'
#' @section Non-interactive use:
#' When called in a non-interactive session, \code{ferx_fit_async()} silently
#' falls back to \code{\link{ferx_fit}} since background jobs only add value
#' in an interactive console.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#' handle <- ferx_fit_async(ex$model, ex$data, method = "focei")
#' # R prompt is free immediately
#' print(handle)
#' fit <- ferx_collect(handle)
#' summary(fit)
#' }
#'
#' @family fitting
#' @export
ferx_fit_async <- function(model, data = NULL, ...,
                           tail_n = 6L,
                           poll_interval = 0.5) {
  if (!interactive()) {
    return(ferx_fit(model, data, ...))
  }
  if (!is.numeric(tail_n) || length(tail_n) != 1L || tail_n < 1L) {
    stop("`tail_n` must be a positive integer scalar")
  }
  if (!is.numeric(poll_interval) || length(poll_interval) != 1L ||
        poll_interval <= 0) {
    stop("`poll_interval` must be a positive numeric scalar")
  }
  tail_n <- as.integer(tail_n)
  dots <- list(...)
  user_wanted_trace <- isTRUE(dots$optimizer_trace)
  model_label <- if (is.character(model)) basename(model) else "model"

  use_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    rstudioapi::hasFun("jobRunScript")

  if (use_rstudio) {
    handle <- .ferx_async_rstudio(model, data, dots, model_label,
                                   tail_n, poll_interval, user_wanted_trace)
  } else {
    if (!requireNamespace("callr", quietly = TRUE)) {
      stop(
        "ferx_fit_async() requires the `callr` package. ",
        "Install with: install.packages(\"callr\")"
      )
    }
    handle <- .ferx_async_callr(model, data, dots, model_label,
                                 tail_n, poll_interval, user_wanted_trace)
  }

  cat(sprintf(
    "Fitting %s in background [%s]. Call ferx_collect(handle) when ready.\n",
    model_label, handle$backend
  ))
  invisible(handle)
}

.ferx_async_callr <- function(model, data, dots, model_label,
                               tail_n, poll_interval, user_wanted_trace) {
  bg <- callr::r_bg(
    func = function(model, data, dots) {
      dots$optimizer_trace <- TRUE
      do.call(ferx::ferx_fit, c(list(model = model, data = data), dots))
    },
    args = list(model = model, data = data, dots = dots),
    package = TRUE,
    supervise = TRUE
  )
  structure(
    list(
      backend          = "callr",
      model_label      = model_label,
      bg               = bg,
      tail_n           = tail_n,
      poll_interval    = poll_interval,
      user_wanted_trace = user_wanted_trace
    ),
    class = "ferx_job"
  )
}

.ferx_async_rstudio <- function(model, data, dots, model_label,
                                 tail_n, poll_interval, user_wanted_trace) {
  rds_path     <- tempfile("ferx_result_", fileext = ".rds")
  sidecar_path <- tempfile("ferx_trace_loc_", fileext = ".txt")

  script_path <- tempfile("ferx_job_", fileext = ".R")
  args_path   <- tempfile("ferx_args_", fileext = ".rds")
  saveRDS(list(model = model, data = data, dots = dots), args_path)

  # The job script deletes args_path itself immediately after readRDS so the
  # file is gone before callr::r_bg() is called - no race between cleanup and
  # the inner bg process starting up.
  writeLines(
    c(
      "local({",
      sprintf("  saved <- readRDS(%s)", deparse(args_path)),
      sprintf("  unlink(%s)", deparse(args_path)),
      "  model        <- saved$model",
      "  data         <- saved$data",
      "  dots         <- saved$dots",
      sprintf("  rds_path     <- %s", deparse(rds_path)),
      sprintf("  sidecar_path <- %s", deparse(sidecar_path)),
      "  dots$optimizer_trace <- TRUE",
      "  bg <- callr::r_bg(",
      "    func = function(model, data, dots) {",
      "      do.call(ferx::ferx_fit, c(list(model = model, data = data), dots))",
      "    },",
      "    args = list(model = model, data = data, dots = dots),",
      "    package = TRUE, supervise = TRUE",
      "  )",
      "  trace_written <- FALSE",
      "  while (bg$is_alive()) {",
      "    if (!trace_written) {",
      "      lines <- c(",
      "        tryCatch(bg$read_output_lines(), error = function(e) character(0)),",
      "        tryCatch(bg$read_error_lines(),  error = function(e) character(0))",
      "      )",
      "      path <- ferx:::.ferx_find_trace_from_lines(lines)",
      "      if (!is.null(path)) {",
      "        writeLines(path, sidecar_path)",
      "        trace_written <- TRUE",
      "      }",
      "    }",
      "    Sys.sleep(0.2)",
      "  }",
      "  result <- bg$get_result()",
      "  saveRDS(result, rds_path)",
      "})"
    ),
    script_path
  )

  # If jobRunScript() throws, clean up both temp files before propagating.
  job_id <- tryCatch(
    rstudioapi::jobRunScript(
      path       = script_path,
      name       = model_label,
      workingDir = getwd()
    ),
    error = function(e) {
      unlink(c(script_path, args_path))
      stop(e)
    }
  )

  structure(
    list(
      backend           = "rstudio",
      model_label       = model_label,
      job_id            = job_id,
      rds_path          = rds_path,
      sidecar_path      = sidecar_path,
      script_path       = script_path,
      tail_n            = tail_n,
      poll_interval     = poll_interval,
      user_wanted_trace = user_wanted_trace
    ),
    class = "ferx_job"
  )
}

#' Retrieve the result of a background fit
#'
#' Blocks until the background fit started by \code{\link{ferx_fit_async}}
#' completes, optionally printing live optimizer trace progress to the console.
#' Pressing \code{Ctrl+C} cancels the wait and kills the background process;
#' \code{NULL} is returned.
#'
#' @param handle A \code{ferx_job} handle returned by
#'   \code{\link{ferx_fit_async}}.
#' @param verbose Logical. Whether to print optimizer trace progress while
#'   waiting. Default \code{TRUE}.
#'
#' @return The \code{ferx_fit} object produced by the background fit.
#'   Returns \code{NULL} if the user interrupts.
#'
#' @section Errors:
#' If the background fit itself errors (e.g. bad model file, engine failure),
#' \code{ferx_collect()} re-throws the error in the calling session.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#' handle <- ferx_fit_async(ex$model, ex$data, method = "focei")
#' fit <- ferx_collect(handle)
#' summary(fit)
#' }
#'
#' @family fitting
#' @export
ferx_collect <- function(handle, verbose = TRUE) {
  if (!inherits(handle, "ferx_job")) {
    stop("`handle` must be a ferx_job object returned by ferx_fit_async()")
  }
  if (handle$backend == "callr") {
    result <- .ferx_collect_callr(handle, verbose)
  } else {
    result <- .ferx_collect_rstudio(handle, verbose)
  }
  if (is.null(result)) return(invisible(NULL))
  if (!handle$user_wanted_trace && !is.null(result$trace_path) &&
        file.exists(result$trace_path)) {
    unlink(result$trace_path)
    result$trace_path <- NULL
  }
  result
}

.ferx_collect_callr <- function(handle, verbose) {
  bg            <- handle$bg
  tail_n        <- handle$tail_n
  poll_interval <- handle$poll_interval
  trace_path    <- NULL
  state         <- list(last_iter = -1L, lines_printed = 0L)

  tryCatch({
    while (bg$is_alive()) {
      # Always drain both streams to prevent the bg process from blocking
      # on a full OS pipe buffer, regardless of whether we display output.
      lines <- c(
        tryCatch(bg$read_output_lines(), error = function(e) character(0)),
        tryCatch(bg$read_error_lines(),  error = function(e) character(0))
      )
      if (verbose && is.null(trace_path)) {
        trace_path <- .ferx_find_trace_from_lines(lines)
      }
      if (verbose && !is.null(trace_path) && file.exists(trace_path)) {
        state <- .ferx_print_trace_tail(trace_path, tail_n, state)
      }
      Sys.sleep(poll_interval)
    }
    bg$get_result()
  }, interrupt = function(e) {
    if (bg$is_alive()) bg$kill()
    message("\nFit interrupted by user.")
    NULL
  })
}

.ferx_collect_rstudio <- function(handle, verbose) {
  job_id        <- handle$job_id
  rds_path      <- handle$rds_path
  sidecar_path  <- handle$sidecar_path
  tail_n        <- handle$tail_n
  poll_interval <- handle$poll_interval
  trace_path    <- NULL
  state         <- list(last_iter = -1L, lines_printed = 0L)

  on_done <- function() {
    if (!is.null(handle$script_path) && file.exists(handle$script_path))
      unlink(handle$script_path)
  }

  # Declare outside the repeat so it's always in scope after the loop,
  # including when the loop body runs zero times (job already done).
  job_state <- "unknown"

  tryCatch({
    repeat {
      job_state <- tryCatch(
        rstudioapi::jobGetState(job_id),
        error = function(e) "unknown"
      )
      if (job_state %in% c("succeeded", "failed")) break

      if (verbose) {
        if (is.null(trace_path) && file.exists(sidecar_path)) {
          trace_path <- trimws(readLines(sidecar_path, warn = FALSE)[1])
          if (!nzchar(trace_path)) trace_path <- NULL
        }
        if (!is.null(trace_path) && file.exists(trace_path)) {
          state <- .ferx_print_trace_tail(trace_path, tail_n, state)
        }
      }
      Sys.sleep(poll_interval)
    }
    on_done()
    if (job_state == "failed" || !file.exists(rds_path)) {
      stop("Background fit failed. Check the Jobs pane for details.")
    }
    result <- readRDS(rds_path)
    unlink(rds_path)
    if (file.exists(sidecar_path)) unlink(sidecar_path)
    result
  }, interrupt = function(e) {
    tryCatch(rstudioapi::jobRemove(job_id), error = function(e) NULL)
    on_done()
    if (file.exists(sidecar_path)) unlink(sidecar_path)
    message("\nFit interrupted by user.")
    NULL
  })
}

#' @export
print.ferx_job <- function(x, ...) {
  status <- if (x$backend == "callr") {
    if (tryCatch(x$bg$is_alive(), error = function(e) FALSE)) "running" else "done"
  } else {
    tryCatch(rstudioapi::jobGetState(x$job_id), error = function(e) "unknown")
  }
  backend_note <- if (x$backend == "rstudio") {
    "rstudio (visible in Jobs pane)"
  } else {
    "callr"
  }
  cat(sprintf("<ferx_job> [%s] %s\n", status, x$model_label))
  cat(sprintf("Backend : %s\n", backend_note))
  cat("Call ferx_collect(handle) to retrieve the result.\n")
  invisible(x)
}

# Parse a trace CSV path from a character vector of stdout/stderr lines.
# ferx-core emits "[ferx] optimizer trace -> /path/to/file.csv" on stderr.
#
# The path may contain spaces (e.g. when TMPDIR resolves under
# "/Users/Some User/...") so we anchor on ".csv" at end-of-line rather
# than refusing any whitespace inside the match. Greedy capture from the
# first "/" through the trailing ".csv" on the announcement line, with a
# trailing "\r" stripped to tolerate Windows line endings.
# NOTE: "/.+\.csv$" assumes the path is the last token on the line. If
# ferx-core ever appends text after the path (e.g. "(appended)"), this
# match will fail silently and fall through to the fallback.
.ferx_find_trace_from_lines <- function(lines) {
  if (length(lines) == 0L) return(NULL)
  parse_one <- function(line) {
    line <- sub("\r$", "", line)
    m <- regexpr("/.+\\.csv$", line, perl = TRUE)
    if (m == -1L) return(NA_character_)
    regmatches(line, m)
  }
  # Primary: lines that look like the trace announcement.
  cand <- grep("optimizer trace", lines, value = TRUE, fixed = TRUE)
  if (length(cand) > 0L) {
    parsed <- vapply(cand, parse_one, character(1L), USE.NAMES = FALSE)
    parsed <- parsed[!is.na(parsed)]
    if (length(parsed) > 0L) return(parsed[[length(parsed)]])
  }
  # Fallback: any line that mentions "trace" and ends in ".csv".
  cand2 <- grep("trace", lines, value = TRUE, ignore.case = TRUE)
  if (length(cand2) == 0L) return(NULL)
  parsed <- vapply(cand2, parse_one, character(1L), USE.NAMES = FALSE)
  parsed <- parsed[!is.na(parsed)]
  if (length(parsed) == 0L) NULL else parsed[[length(parsed)]]
}

# Keep the old name as an internal alias so existing tests that mock the
# bg-object interface still work; .ferx_find_trace_from_bg() now delegates
# to the lines-based helper.
.ferx_find_trace_from_bg <- function(bg) {
  stdout_lines <- tryCatch(bg$read_output_lines(), error = function(e) character(0))
  stderr_lines <- tryCatch(bg$read_error_lines(), error = function(e) character(0))
  .ferx_find_trace_from_lines(c(stdout_lines, stderr_lines))
}

# Tail the trace CSV in place: each poll reads the file, prints any new rows
# beyond `state$last_iter`, and updates the state. Append-only display so it
# works across terminals without relying on ANSI cursor movement.
.ferx_print_trace_tail <- function(path, tail_n, state) {
  tr <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(tr) || nrow(tr) == 0L) return(state)
  last <- tr$iter[nrow(tr)]
  if (is.null(last) || !is.finite(last) || last <= state$last_iter) return(state)

  if (state$last_iter < 0L) {
    cat(sprintf("  %5s  %12s  %12s  %12s\n",
                "iter", "OFV", "grad_norm", "optimizer"))
  }
  new_rows <- tr[tr$iter > state$last_iter, , drop = FALSE]
  if (nrow(new_rows) > tail_n) {
    new_rows <- utils::tail(new_rows, tail_n)
  }
  for (i in seq_len(nrow(new_rows))) {
    row <- new_rows[i, ]
    gn <- row$grad_norm
    gn_str <- if (is.null(gn) || is.na(gn)) "         -  " else sprintf("%12.4g", gn)
    opt_str <- row$optimizer %||% ""
    cat(sprintf("  %5d  %12.4f  %s  %12s\n",
                as.integer(row$iter), row$ofv, gn_str, opt_str))
  }
  state$last_iter <- last
  state$lines_printed <- state$lines_printed + nrow(new_rows)
  state
}
