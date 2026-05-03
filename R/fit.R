#' Fit a nonlinear mixed effects model
#'
#' Fits a NLME model using FOCE or FOCEI estimation with a Rust backend
#' powered by automatic differentiation (Enzyme).
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV file (ID, TIME, DV, EVID, AMT, CMT, ...)
#' @param method Estimation method(s). Either a single string or a character
#'   vector of methods to run in sequence. Each stage is seeded with the
#'   previous stage's converged parameters, and only the final stage produces
#'   the reported covariance/diagnostics. Supported methods: \code{"foce"},
#'   \code{"focei"}, \code{"saem"}, \code{"gn"} (Gauss-Newton / BHHH), or
#'   \code{"gn_hybrid"} (Gauss-Newton followed by a FOCEI polish step).
#'   Example chain: \code{c("saem", "focei")}.
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
#'   CPU. Pass an integer to cap the pool — useful on shared machines, in CI,
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
#'   Pass the result to \code{\link{ferx_read_trace}} or
#'   \code{\link{ferx_plot_trace}} to inspect optimizer progress. Default
#'   \code{FALSE}.
#' @param scale_params Logical. If \code{TRUE} (default), apply a per-coordinate
#'   scaling layer on top of the existing log/Cholesky parameterization so
#'   that every transformed parameter the outer optimizer sees is O(1).
#'
#'   \strong{Why it helps.} Population PK parameters often differ by orders
#'   of magnitude on the natural scale (e.g. \code{CL = 0.0015}, \code{V =
#'   200}, \code{Ka = 0.8}). Even after the log transform, the packed
#'   optimizer vector mixes large negatives (\code{log(0.0015) ≈ -6.5}) with
#'   large positives (\code{log(200) ≈ 5.3}), which leaves the Hessian
#'   poorly conditioned and slows down NLopt / BFGS / Gauss-Newton. With
#'   \code{scale_params = TRUE} each coordinate is divided by
#'   \code{|x0[i]|} (when \code{|x0[i]| > 0.1}, otherwise by \code{1.0}), so
#'   the optimizer works in a near-unit-magnitude space. The scale is
#'   computed once from the initial point and held fixed for the entire run.
#'
#'   \strong{When to set \code{FALSE}.} Mathematically the scaling is
#'   transparent — the OFV, estimates, standard errors, and diagnostics
#'   should match the unscaled fit to numerical tolerance. Set
#'   \code{scale_params = FALSE} to:
#'   \itemize{
#'     \item reproduce the bit-exact iteration trajectory of a pre-scaling
#'       fit (e.g. when comparing against a baseline produced by an older
#'       \code{ferx} build),
#'     \item debug suspected scaling-related issues by toggling the layer
#'       off without changing anything else,
#'     \item validate that scaled and unscaled fits land on the same
#'       optimum on a new model.
#'   }
#'
#'   Applies to all outer optimizers: NLopt FOCE/FOCEI (BOBYQA, SLSQP,
#'   L-BFGS, MMA), the hand-rolled BFGS, Gauss-Newton / BHHH, and the SAEM
#'   M-step.
#' @param settings Optional named list of estimation-method-specific options
#'   forwarded to the Rust \code{FitOptions}. Use this to tune knobs that do
#'   not have a dedicated \code{ferx_fit()} argument, without needing a new
#'   wrapper release for each option. Recognized keys include SAEM
#'   (\code{n_exploration}, \code{n_convergence}, \code{n_mh_steps},
#'   \code{adapt_interval}, \code{seed}), SIR tuning (\code{sir_samples},
#'   \code{sir_resamples}, \code{sir_seed}), Gauss-Newton (\code{gn_lambda}),
#'   optimizer selection (\code{optimizer} — one of \code{"bobyqa"} (default),
#'   \code{"slsqp"}, \code{"lbfgs"}, \code{"nlopt_lbfgs"}, \code{"mma"},
#'   \code{"bfgs"}, \code{"trust_region"}; also \code{global_search},
#'   \code{global_maxeval}), the outer optimizer iteration cap
#'   (\code{maxiter}), the inner (per-subject EBE) loop
#'   (\code{inner_maxiter}, \code{inner_tol}), and the Steihaug CG budget
#'   for \code{optimizer = "trust_region"} (\code{steihaug_max_iters}).
#'   Values that duplicate a dedicated argument
#'   (\code{method}, \code{covariance}, \code{verbose},
#'   \code{bloq_method}, \code{threads}, \code{sir}) are rejected — pass
#'   them via the dedicated argument. Unknown keys and malformed values
#'   also raise an error. Settings apply on top of the model file's
#'   \code{[fit_options]} block.
#'
#' @return A list with components:
#'   \item{converged}{Logical; did the optimizer converge}
#'   \item{method}{Estimation method used}
#'   \item{n_iterations}{Number of outer iterations run}
#'   \item{ofv}{Objective function value (-2 log-likelihood)}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{theta}{Named numeric vector of fixed effect estimates}
#'   \item{omega}{Between-subject variability covariance matrix}
#'   \item{sigma}{Residual error parameter estimates}
#'   \item{se_theta}{Standard errors for theta (NULL if covariance step failed)}
#'   \item{se_omega}{Standard errors for omega diagonal}
#'   \item{se_sigma}{Standard errors for sigma}
#'   \item{sdtab}{Data frame with ID, TIME, DV, PRED, IPRED, CWRES, IWRES, ETA1..n, OCC (if IOV)}
#'   \item{warnings}{Character vector of warnings}
#'   \item{sir_ess}{SIR effective sample size (NULL if SIR not run)}
#'   \item{sir_ci_theta, sir_ci_omega, sir_ci_sigma}{SIR 95\% CI matrices
#'     with columns \code{lower} and \code{upper} (NULL if SIR not run)}
#'   \item{trace_path}{Path to the optimizer trace CSV, or \code{NULL} when
#'     \code{optimizer_trace = FALSE}. Pass to \code{\link{ferx_read_trace}}
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
#'   \item{wall_time_secs}{Total wall-clock time for the fit in seconds.}
#'   \item{model_name}{Model name from the \code{.ferx} file. Falls back to
#'     the model file's basename (without extension) when the file declares
#'     no name.}
#'   \item{data_name}{Dataset name, derived as the basename of the \code{data}
#'     path with its extension stripped.}
#'   \item{gradient}{The inner-loop gradient method as requested by the
#'     caller (one of \code{"auto"}, \code{"ad"}, \code{"fd"}). The
#'     \emph{resolved} method may differ — see the \code{gradient} argument.}
#'   \item{ferx_version}{ferx-nlme library version string.}
#'   \item{cov_matrix}{Full parameter covariance matrix as a named numeric
#'     matrix (params × params). \code{NULL} when covariance step was not run
#'     or failed. Use \code{\link{ferx_cor_matrix}} to inspect correlations.}
#'   \item{eta_normality}{Data frame with Shapiro-Wilk normality test for each
#'     ETA: columns \code{eta}, \code{W}, \code{p_val}, \code{flag}. A
#'     \code{[!]} flag appears when \code{p < 0.05}. \code{NULL} when fewer
#'     than 3 subjects or no ETA columns in \code{sdtab}.}
#'   \item{omega_iov}{IOV variance matrix for kappa parameters (\code{NULL} if no IOV).}
#'   \item{kappa_names}{Names of kappa (IOV) parameters (\code{NULL} if no IOV).}
#'   \item{se_kappa}{Standard errors for kappa parameters: length \code{d}
#'     (diagonal kappa) or \code{d*(d+1)/2} (block_kappa, lower-triangle order).
#'     \code{NULL} if covariance step not run or no IOV.}
#'   \item{shrinkage_kappa}{Shrinkage values for kappa EBEs (\code{NULL} if no IOV).}
#'   \item{ebe_kappas}{Data frame with columns \code{ID}, \code{OCC}, and one
#'     column per kappa parameter containing per-subject per-occasion kappa
#'     EBEs. \code{ID} carries the original subject identifier (matching
#'     \code{sdtab}); \code{OCC} carries the labeled occasion in first-seen
#'     order. \code{NULL} when no IOV.}
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Default: derivative-free BOBYQA — robust on the FOCE surface
#' fit <- ferx_fit(ex$model, ex$data)
#'
#' # Gradient-based SLSQP — faster on smooth, well-behaved problems
#' fit_slsqp <- ferx_fit(ex$model, ex$data,
#'                       settings = list(optimizer = "slsqp"))
#'
#' # Second-order trust region with a tuned CG budget
#' fit_tr <- ferx_fit(ex$model, ex$data,
#'                    settings = list(optimizer          = "trust_region",
#'                                    steihaug_max_iters = 100L))
#'
#' # Fine-tune inner (per-subject EBE) loop via `settings`
#' fit_fast <- ferx_fit(ex$model, ex$data,
#'                      settings = list(inner_maxiter = 100L,
#'                                      inner_tol     = 1e-6))
#'
#' # Chain SAEM to FOCEI (SAEM explores, FOCEI polishes):
#' result <- ferx_fit("warfarin.ferx", "warfarin.csv",
#'                    method = c("saem", "focei"))
#'
#' # Compare with mu-referencing off
#' fit_no_mu <- ferx_fit("warfarin.ferx", "warfarin.csv", mu_referencing = FALSE)
#'
#' # Check which ETAs were detected
#' result$warnings
#'
#' # Likelihood-based BLOQ handling (M3):
#' bloq <- ferx_example("warfarin_bloq")
#' result <- ferx_fit(bloq$model, bloq$data, method = "focei", bloq_method = "m3")
#'
#' # Tune SAEM phase lengths via `settings`:
#' result <- ferx_fit("warfarin.ferx", "warfarin.csv",
#'                    method  = "saem",
#'                    settings = list(n_exploration = 200,
#'                                    n_convergence = 400,
#'                                    seed = 42L))
#' }
#'
#' @export
ferx_fit <- function(model, data,
                     method = "foce",
                     covariance = TRUE,
                     verbose = TRUE,
                     bloq_method = NULL,
                     threads = NULL,
                     mu_referencing = TRUE,
                     sir = FALSE,
                     gradient = c("auto", "ad", "fd"),
                     optimizer_trace = FALSE,
                     scale_params = TRUE,
                     max_unconverged_frac = NULL,
                     min_obs_for_convergence_check = NULL,
                     settings = NULL) {
  gradient <- match.arg(gradient)
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
      match.arg(
        tolower(gsub("[^a-z0-9]", "_", m)),
        c("foce", "focei", "saem", "gn", "gn_hybrid")
      )
    },
    character(1L),
    USE.NAMES = FALSE
  )
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
  # — but we forbid that below via the RESERVED list.
  if (isTRUE(optimizer_trace)) {
    settings <- c(list(optimizer_trace = TRUE), settings)
  }
  # Only inject scale_params when FALSE — the Rust default is true, so omitting
  # it preserves the behaviour callers expect when they don't pass the argument.
  if (isFALSE(scale_params)) {
    settings <- c(list(scale_params = FALSE), settings)
  }
  if (!is.null(max_unconverged_frac)) {
    if (!is.numeric(max_unconverged_frac) || length(max_unconverged_frac) != 1L ||
        !is.finite(max_unconverged_frac) || max_unconverged_frac < 0 || max_unconverged_frac > 1) {
      stop("`max_unconverged_frac` must be a numeric scalar between 0 and 1")
    }
    settings <- c(list(max_unconverged_frac = max_unconverged_frac), settings)
  }
  if (!is.null(min_obs_for_convergence_check)) {
    if (!is.numeric(min_obs_for_convergence_check) ||
        length(min_obs_for_convergence_check) != 1L ||
        !is.finite(min_obs_for_convergence_check) ||
        min_obs_for_convergence_check != as.integer(min_obs_for_convergence_check) ||
        min_obs_for_convergence_check < 0L) {
      stop("`min_obs_for_convergence_check` must be a non-negative integer scalar")
    }
    settings <- c(list(min_obs_for_convergence_check = as.integer(min_obs_for_convergence_check)),
                  settings)
  }
  settings_parts <- .ferx_settings_to_strings(settings)

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
    if (length(v) == 0) return(NULL)
    m <- matrix(v, ncol = 2, byrow = TRUE,
                dimnames = list(row_names, c("lower", "upper")))
    m
  }
  result$sir_ci_theta <- reshape_ci(result$sir_ci_theta, result$theta_names)
  n_eta <- if (is.null(dim(result$omega))) NULL else nrow(result$omega)
  eta_names <- if (!is.null(n_eta)) paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")") else NULL
  result$sir_ci_omega <- reshape_ci(result$sir_ci_omega, eta_names)
  sig_names <- paste0("SIGMA(", seq_along(result$sigma), ")")
  result$sir_ci_sigma <- reshape_ci(result$sir_ci_sigma, sig_names)

  # Normalize trace_path: NULL/empty means no trace was written
  tp <- result$trace_path
  if (is.null(tp) || length(tp) == 0L || !nzchar(tp[[1L]])) {
    result$trace_path <- NULL
  } else {
    .ferx_state$last_trace_path  <- result$trace_path
    .ferx_state$last_trace_time  <- fit_started_at
    .ferx_state$last_trace_model <- model
  }

  # Normalize shrinkage: NaN → NA (consistent with other optional numerics)
  if (!is.null(result$shrinkage_eta)) {
    result$shrinkage_eta[!is.finite(result$shrinkage_eta)] <- NA_real_
  }
  if (!is.null(result$shrinkage_eps) && !is.finite(result$shrinkage_eps)) {
    result$shrinkage_eps <- NA_real_
  }

  # Normalize covariance_status: missing from older binaries → "not_requested"
  if (is.null(result$covariance_status)) {
    result$covariance_status <- "not_requested"
  }

  # Reshape cov_matrix into a named square matrix (param × param)
  d <- result$cov_matrix_dim %||% 0L
  if (!is.null(result$cov_matrix) && length(result$cov_matrix) > 0L && d > 0L) {
    m <- matrix(result$cov_matrix, nrow = d, ncol = d, byrow = TRUE)
    n_theta  <- length(result$theta_names)
    n_eta    <- result$omega_dim %||% 0L
    n_sigma  <- length(result$sigma)
    n_omega_packed <- d - n_theta - n_sigma
    # Determine parameterisation: diagonal (n_omega_packed == n_eta) or block
    omega_names <- if (n_omega_packed == n_eta) {
      paste0("OMEGA(", seq_len(n_eta), ",", seq_len(n_eta), ")")
    } else {
      # Block lower-triangle: L(i,j) for i >= j, column-major
      nm <- character(n_omega_packed)
      k  <- 0L
      for (i in seq_len(n_eta)) {
        for (j in seq_len(i)) {
          k <- k + 1L
          nm[k] <- sprintf("OMEGA_L(%d,%d)", i, j)
        }
      }
      nm
    }
    pnames <- c(
      result$theta_names,
      if (n_omega_packed > 0L) omega_names else character(0L),
      if (n_sigma > 0L) paste0("SIGMA(", seq_len(n_sigma), ")") else character(0L)
    )
    if (length(pnames) == d) rownames(m) <- colnames(m) <- pnames
    result$cov_matrix <- m
  } else {
    result$cov_matrix <- NULL
  }
  result$cov_matrix_dim <- NULL

  # ETA normality (Shapiro-Wilk) — computed in R from sdtab EBEs
  result$eta_normality <- .ferx_compute_eta_normality(result$sdtab)
  # Push normality warnings into the warnings vector
  if (!is.null(result$eta_normality)) {
    for (i in seq_len(nrow(result$eta_normality))) {
      row <- result$eta_normality[i, ]
      if (nzchar(row$flag) && !is.na(row$p_val)) {
        result$warnings <- c(
          result$warnings,
          sprintf("%s Shapiro-Wilk p=%.4f — distribution may be non-normal", row$eta, row$p_val)
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
    if (length(result$kappa_names) > 0L && length(result$shrinkage_kappa) == length(result$kappa_names))
      names(result$shrinkage_kappa) <- result$kappa_names
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

  # Fall back to the model file's basename when the .ferx file declares no name
  if (is.null(result$model_name) || !nzchar(result$model_name)) {
    result$model_name <- tools::file_path_sans_ext(basename(model))
  }

  # Store the requested gradient method
  result$gradient <- gradient

  class(result) <- "ferx_fit"
  result
}

# Compute Shapiro-Wilk normality test for each ETA, one row per ID.
# Returns a data.frame with columns: eta, W, p_val, flag.
.ferx_compute_eta_normality <- function(sdtab) {
  if (is.null(sdtab) || !is.data.frame(sdtab)) return(NULL)
  eta_cols <- grep("^ETA", names(sdtab), value = TRUE)
  if (length(eta_cols) == 0L) return(NULL)
  id_col   <- if ("ID" %in% names(sdtab)) "ID" else names(sdtab)[1L]
  one_per  <- sdtab[!duplicated(sdtab[[id_col]]), eta_cols, drop = FALSE]
  do.call(rbind, lapply(eta_cols, function(col) {
    vals <- one_per[[col]]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 3L) {
      return(data.frame(eta = col, W = NA_real_, p_val = NA_real_, flag = "",
                        stringsAsFactors = FALSE))
    }
    sw  <- shapiro.test(vals)
    flg <- if (sw$p.value < 0.05) "[!]" else ""
    data.frame(eta = col, W = round(as.numeric(sw$statistic), 3),
               p_val = round(sw$p.value, 4), flag = flg, stringsAsFactors = FALSE)
  }))
}

#' @export
print.ferx_fit <- function(x, ...) {
  bar <- strrep("=", 60)
  cat(bar, "\n", sep = "")
  cat("NONLINEAR MIXED EFFECTS MODEL ESTIMATION\n")
  cat(bar, "\n\n", sep = "")

  if (!is.null(x$model_name) && nzchar(x$model_name)) {
    cat("Model:     ", x$model_name, "\n", sep = "")
  }
  if (!is.null(x$data_name) && nzchar(x$data_name)) {
    cat("Dataset:   ", x$data_name, "\n", sep = "")
  }
  cat("Converged: ", if (isTRUE(x$converged)) "YES" else "NO", "\n", sep = "")
  if (!is.null(x$method_chain) && length(x$method_chain) > 1) {
    cat("Estimation chain:  ", paste(x$method_chain, collapse = " -> "), "\n", sep = "")
  } else {
    cat("Estimation method: ", x$method, "\n", sep = "")
  }
  if (!is.null(x$n_iterations)) {
    cat("Iterations: ", x$n_iterations, "\n", sep = "")
  }

  cat("\n--- Objective Function ---\n")
  cat(sprintf("OFV:  %.4f\n", x$ofv))
  cat(sprintf("AIC:  %.4f\n", x$aic))
  cat(sprintf("BIC:  %.4f\n", x$bic))

  n_par <- if (is.null(x$n_parameters)) NA_integer_ else x$n_parameters
  cat(sprintf(
    "\nSubjects: %d  Observations: %d  Parameters: %s\n",
    x$n_subjects, x$n_obs,
    if (is.na(n_par)) "?" else format(n_par)
  ))

  # THETA
  cat("\n--- THETA Estimates ---\n")
  cat(sprintf("%-16s %12s %12s %10s\n", "Parameter", "Estimate", "SE", "%RSE"))
  cat(strrep("-", 52), "\n", sep = "")
  theta_names <- names(x$theta)
  if (is.null(theta_names)) theta_names <- paste0("THETA", seq_along(x$theta))
  for (i in seq_along(x$theta)) {
    est <- x$theta[i]
    if (!is.null(x$se_theta) && length(x$se_theta) >= i) {
      se_val <- x$se_theta[i]
      rse <- if (abs(est) > 1e-12) abs(se_val / est) * 100 else NaN
      se_str <- sprintf("%.6f", se_val)
      rse_str <- sprintf("%.1f", rse)
    } else {
      se_str <- "N/A"; rse_str <- "N/A"
    }
    cat(sprintf("%-16s %12.6f %12s %10s\n", theta_names[i], est, se_str, rse_str))
  }

  # OMEGA
  cat("\n--- OMEGA Estimates ---\n")
  om <- x$omega
  if (is.null(dim(om))) om <- matrix(om, 1, 1)
  n_eta <- nrow(om)
  has_offdiag <- FALSE
  for (i in seq_len(n_eta)) {
    var_ii <- om[i, i]
    cv_pct <- if (var_ii > 0) sqrt(var_ii) * 100 else 0
    se_str <- if (!is.null(x$se_omega) && length(x$se_omega) >= i) {
      sprintf("%.6f", x$se_omega[i])
    } else {
      "N/A"
    }
    cat(sprintf(
      "  OMEGA(%d,%d) = %.6f  (CV%% = %.1f)  SE = %s\n",
      i, i, var_ii, cv_pct, se_str
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
        var_i <- om[i, i]; var_j <- om[j, j]
        corr <- if (var_i > 0 && var_j > 0) cov_ij / (sqrt(var_i) * sqrt(var_j)) else 0
        cat(sprintf(
          "  OMEGA(%d,%d) = %.6f  (corr = %.4f)\n",
          i, j, cov_ij, corr
        ))
      }
    }
  }

  # OMEGA_IOV (Inter-Occasion Variability)
  if (!is.null(x$omega_iov)) {
    cat("\n--- OMEGA_IOV Estimates (Inter-Occasion Variability) ---\n")
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
      var_ii <- m_iov[i, i]
      cv_pct <- if (var_ii > 0) sqrt(var_ii) * 100 else 0
      se_idx <- if (is_block_se) diag_se_idx(i) else i
      se_str <- if (!is.null(x$se_kappa) && n_se >= se_idx) {
        sprintf("%.6f", x$se_kappa[se_idx])
      } else {
        "N/A"
      }
      shr_str <- if (!is.null(x$shrinkage_kappa) && length(x$shrinkage_kappa) >= i) {
        sprintf("%.1f%%", x$shrinkage_kappa[i] * 100)
      } else {
        "N/A"
      }
      cat(sprintf(
        "  %s = %.6f  (CV%% = %.1f)  SE = %s  Shrinkage = %s\n",
        kap_names[i], var_ii, cv_pct, se_str, shr_str
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
          var_i <- m_iov[i, i]; var_j <- m_iov[j, j]
          corr <- if (var_i > 0 && var_j > 0) cov_ij / (sqrt(var_i) * sqrt(var_j)) else 0
          cat(sprintf(
            "  %s ~ %s : cov = %.6f  (corr = %.4f)\n",
            kap_names[i], kap_names[j], cov_ij, corr
          ))
        }
      }
    }
  }

  # SIGMA
  cat("\n--- SIGMA Estimates ---\n")
  for (i in seq_along(x$sigma)) {
    se_str <- if (!is.null(x$se_sigma) && length(x$se_sigma) >= i) {
      sprintf("%.6f", x$se_sigma[i])
    } else {
      "N/A"
    }
    cat(sprintf("  SIGMA(%d) = %.6f  SE = %s\n", i, x$sigma[i], se_str))
  }

  # SIR uncertainty
  if (!is.null(x$sir_ess)) {
    cat("\n--- SIR Uncertainty (95% CI) ---\n")
    cat(sprintf("Effective sample size: %.1f\n", x$sir_ess))
    print_ci <- function(m) {
      if (is.null(m)) return(invisible())
      for (i in seq_len(nrow(m))) {
        cat(sprintf("  %s : [%.6f, %.6f]\n", rownames(m)[i], m[i, 1], m[i, 2]))
      }
    }
    print_ci(x$sir_ci_theta)
    print_ci(x$sir_ci_omega)
    print_ci(x$sir_ci_sigma)
  }

  # Shrinkage
  has_shrinkage <- (!is.null(x$shrinkage_eta) && any(!is.na(x$shrinkage_eta))) ||
                   (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps))
  if (has_shrinkage) {
    cat("\n--- Shrinkage ---\n")
    if (!is.null(x$shrinkage_eta)) {
      for (k in seq_along(x$shrinkage_eta)) {
        sh <- x$shrinkage_eta[k]
        if (!is.na(sh)) cat(sprintf("  ETA%d shrinkage: %.1f%%\n", k, sh * 100))
      }
    }
    if (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps)) {
      cat(sprintf("  EPS shrinkage:  %.1f%%\n", x$shrinkage_eps * 100))
    }
  }

  # Run info
  cov_status <- if (!is.null(x$covariance_status)) x$covariance_status else "unknown"
  cov_str <- switch(cov_status,
    computed      = "computed",
    failed        = "FAILED",
    not_requested = "not requested",
    cov_status
  )
  cat("\n--- Run Info ---\n")
  cat("  Covariance:", cov_str, "\n")
  if (!is.null(x$gradient)) {
    cat("  Gradient (requested): ", x$gradient, "\n", sep = "")
  }
  if (!is.null(x$wall_time_secs)) {
    cat(sprintf("  Wall time:  %.1fs\n", x$wall_time_secs))
  }
  if (!is.null(x$ferx_version)) {
    cat("  ferx v", x$ferx_version, "\n", sep = "")
  }

  if (length(x$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in x$warnings) cat("  *", w, "\n")
  }

  cat(bar, "\n", sep = "")
  invisible(x)
}

#' Summarize a ferx fit result
#'
#' Returns a compact list with the most-used diagnostic fields: OFV/AIC/BIC,
#' parameter estimates, standard errors, shrinkage, and run metadata.
#'
#' @param object A `ferx_fit` object returned by \code{\link{ferx_fit}}.
#' @param ... Ignored.
#' @return A `ferx_summary` list (invisibly). Print method formats the output.
#' @export
summary.ferx_fit <- function(object, ...) {
  x <- object
  s <- list(
    model_name      = x$model_name %||% NA_character_,
    data_name       = x$data_name %||% NA_character_,
    gradient        = x$gradient %||% NA_character_,
    method          = x$method,
    converged       = x$converged,
    ofv             = x$ofv,
    aic             = x$aic,
    bic             = x$bic,
    n_subjects      = x$n_subjects,
    n_obs           = x$n_obs,
    n_parameters    = x$n_parameters,
    n_iterations    = x$n_iterations,
    theta           = x$theta,
    se_theta        = x$se_theta,
    omega           = x$omega,
    se_omega        = x$se_omega,
    sigma           = x$sigma,
    se_sigma        = x$se_sigma,
    shrinkage_eta   = x$shrinkage_eta,
    shrinkage_eps   = x$shrinkage_eps,
    covariance_status = x$covariance_status,
    wall_time_secs  = x$wall_time_secs,
    ferx_version    = x$ferx_version,
    ebe_convergence_warnings = x$ebe_convergence_warnings,
    max_unconverged_subjects = x$max_unconverged_subjects,
    total_ebe_fallbacks      = x$total_ebe_fallbacks
  )
  class(s) <- "ferx_summary"
  s
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @export
print.ferx_summary <- function(x, ...) {
  cat(sprintf("ferx %s — %s  [gradient: %s (requested)]\n",
              x$ferx_version %||% "?", toupper(x$method), x$gradient %||% "?"))
  cat(sprintf("Model: %s  |  Dataset: %s  |  Converged: %s\n",
              x$model_name %||% "?",
              x$data_name %||% "?",
              if (isTRUE(x$converged)) "YES" else "NO"))
  cat(sprintf("OFV: %.4f  AIC: %.4f  BIC: %.4f\n", x$ofv, x$aic, x$bic))
  cat(sprintf("Subjects: %d  Obs: %d  Params: %d  Iter: %d\n",
              x$n_subjects, x$n_obs, x$n_parameters %||% NA, x$n_iterations))

  if (!is.null(x$shrinkage_eta) && any(!is.na(x$shrinkage_eta))) {
    sh_str <- paste(sprintf("ETA%d=%.1f%%", seq_along(x$shrinkage_eta),
                            x$shrinkage_eta * 100), collapse = "  ")
    cat("Shrinkage:", sh_str, "\n")
  }
  if (!is.null(x$shrinkage_eps) && !is.na(x$shrinkage_eps)) {
    cat(sprintf("EPS shrinkage: %.1f%%\n", x$shrinkage_eps * 100))
  }

  cov_str <- switch(x$covariance_status %||% "unknown",
    computed = "computed", failed = "FAILED", not_requested = "not requested",
    x$covariance_status)
  wall <- if (!is.null(x$wall_time_secs)) sprintf("%.1fs", x$wall_time_secs) else "?"
  cat(sprintf("Covariance: %s  |  Wall time: %s\n", cov_str, wall))

  invisible(x)
}

# Convert a named-list of settings into two parallel character vectors for the
# Rust FFI. Each value is stringified in a Rust-parser-friendly form (logicals
# → "true"/"false", numerics with full precision, NULL/NA → "null"). Strict
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
    if (is.null(v) || (length(v) == 1L && is.na(v))) return("null")
    if (length(v) != 1L) {
      stop(sprintf("`settings$%s` must be a length-1 scalar", key))
    }
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.numeric(v)) {
      if (!is.finite(v)) {
        stop(sprintf("`settings$%s` must be a finite number", key))
      }
      return(format(v, scientific = FALSE, trim = TRUE, digits = 17))
    }
    if (is.character(v)) return(v)
    stop(sprintf("`settings$%s` has unsupported type `%s`", key, class(v)[1L]))
  }
  values <- vapply(
    seq_along(settings),
    function(i) stringify(nms[i], settings[[i]]),
    character(1L)
  )
  list(keys = nms, values = values)
}
