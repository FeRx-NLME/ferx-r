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
#'   \item{omega_iov}{IOV variance matrix for kappa parameters (NULL if no IOV)}
#'   \item{kappa_names}{Names of kappa (IOV) parameters (NULL if no IOV)}
#'   \item{se_kappa}{Standard errors for kappa parameters: length \code{d} (diagonal kappa) or \code{d*(d+1)/2} (block_kappa, lower-triangle order). NULL if covariance step not run or no IOV.}
#'   \item{shrinkage_kappa}{Shrinkage values for kappa EBEs (NULL if no IOV)}
#'   \item{ebe_kappas}{Data frame with columns ID, OCC, and one column per kappa parameter
#'     containing per-subject per-occasion kappa EBEs (NULL if no IOV)}
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

  settings_parts <- .ferx_settings_to_strings(settings)

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

  class(result) <- "ferx_fit"
  result
}

#' @export
print.ferx_fit <- function(x, ...) {
  bar <- strrep("=", 60)
  cat(bar, "\n", sep = "")
  cat("NONLINEAR MIXED EFFECTS MODEL ESTIMATION\n")
  cat(bar, "\n\n", sep = "")

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

  if (length(x$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in x$warnings) cat("  *", w, "\n")
  }

  cat(bar, "\n", sep = "")
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
