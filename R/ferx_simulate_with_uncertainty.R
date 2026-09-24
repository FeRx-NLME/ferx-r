#' Simulate with parameter-uncertainty propagation
#'
#' For each parameter set drawn from the uncertainty distribution, the
#' per-subject random-effect / residual-error simulator runs
#' \code{n_sim_per_draw} times. The result includes both individual
#' variability (etas, epsilons) and parameter uncertainty - useful for
#' uncertainty-aware VPCs, dose-recommendation intervals, and any analysis
#' where treating the ML estimates as fixed would understate variability.
#'
#' Two uncertainty sources are supported:
#' \describe{
#'   \item{\code{method = "asymptotic"}}{Multivariate normal around the ML
#'     estimate in the engine's packed (log-theta, Cholesky-omega, log-sigma)
#'     parameter space, using \code{fit$cov_matrix}. Requires \code{fit}
#'     to come from a \code{ferx_fit()} call with \code{covariance = TRUE}.}
#'   \item{\code{method = "sir"}}{Sample with replacement from
#'     \code{fit$sir_resamples}. Requires the fit to have been run with
#'     \code{sir = TRUE} and \code{sir_keep_samples = TRUE} (passed via
#'     \code{settings}).}
#' }
#'
#' @section Bounded thetas under \code{method = "asymptotic"}:
#' A theta with a non-negative lower bound is packed as \code{log(theta)}, so
#' its asymptotic draws are log-normal: always above 0, with no ceiling. A
#' draw past a declared upper bound is rejected and the whole parameter set
#' redrawn, which keeps the draws inside the bound but truncates their
#' distribution.
#'
#' For a \code{"logit_probability"} theta - one declared on \eqn{(0, 1)} and
#' used as \code{inv_logit(logit(THETA) + ETA)} - the draws are therefore not
#' the logit-normal ones the fit implies. With an upper bound below 1 they are
#' truncated, and so pulled low. With an upper bound at or above 1 (an
#' undeclared bound defaults to \code{1e9}) a draw at or above 1 is simulated,
#' and the model's \code{logit()} clamps it so that every subject in that draw
#' gets a probability of 1. When more than 0.1% of draws are expected to reach
#' 1, the function warns (class \code{"ferx_logit_probability_draws"}), naming
#' the theta and the expected share.
#'
#' \code{method = "sir"} does not share the bias: its draws come from the same
#' proposal but are resampled by likelihood weight, so a value reaches the
#' pool only in proportion to how well it fits. Neither does a probability
#' declared on the logit scale - \code{inv_logit(LOGIT_F + ETA_F)} with a
#' negative lower bound on \code{LOGIT_F} - since that theta is packed, and so
#' drawn, on the logit scale.
#'
#' @param model Path to a .ferx model file
#' @param data Path to a NONMEM-format CSV (provides population structure).
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
#' @param fit A \code{ferx_fit} result. Must carry either \code{cov_matrix}
#'   (asymptotic) or \code{sir_resamples} (SIR) depending on \code{method}.
#' @param n_uncertainty_draws Number of parameter sets to draw from the
#'   uncertainty distribution
#' @param n_sim_per_draw Number of eta/eps replicates per parameter draw
#' @param method Either \code{"asymptotic"} (default) or \code{"sir"}
#' @param seed Random seed for reproducibility
#'
#' @return A data.frame with columns: DRAW, SIM, ID, TIME, CMT, IPRED, DV_SIM,
#'   OBSERVED. \code{CMT} is the observation compartment; \code{OBSERVED} is
#'   \code{NA} throughout (this path is Gaussian-only - a drug-driven ODE-TTE
#'   endpoint is not supported here, and \code{OBSERVED} never echoes the input
#'   \code{DV}; see \code{\link{ferx_simulate}} for what each column holds).
#'   Row count: \code{n_uncertainty_draws * n_sim_per_draw * n_obs}.
#'
#' @inheritSection ferx_simulate Errors raised by the engine
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, covariance = TRUE)
#'
#' # Asymptotic (default): fast, MVN draws around the ML estimate
#' sims <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10
#' )
#' head(sims)         # SIM, ID, TIME, IPRED, DV_SIM
#' range(sims$DV_SIM) # sanity check on simulated range
#'
#' # SIR method: requires sir = TRUE + sir_keep_samples = TRUE at fit time
#' fit_sir <- ferx_fit(ex$model, ex$data,
#'                     covariance = TRUE, sir = TRUE,
#'                     settings = list(sir_keep_samples = TRUE))
#' sims_sir <- ferx_simulate_with_uncertainty(
#'   ex$model, ex$data, fit_sir,
#'   n_uncertainty_draws = 200, n_sim_per_draw = 10,
#'   method = "sir"
#' )
#'
#' # Summarise: 90% prediction interval per time point
#' pi90 <- aggregate(DV_SIM ~ TIME, data = sims,
#'                   FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
#' }
#'
#' @family simulation
#' @export
ferx_simulate_with_uncertainty <- function(model, data, fit,
                                           n_uncertainty_draws = 100L,
                                           n_sim_per_draw = 1L,
                                           method = c("asymptotic", "sir"),
                                           seed = 42L) {
  stopifnot(file.exists(model), file.exists(data))
  method <- match.arg(method)

  # Check the cheap, user-facing arg constraints before the fit-object
  # validation so callers passing `n_uncertainty_draws = 0` see the obvious
  # error rather than a follow-on complaint about `cov_matrix` / `sir_resamples`.
  if (n_uncertainty_draws < 1L) {
    stop("`n_uncertainty_draws` must be >= 1.")
  }
  if (n_sim_per_draw < 1L) {
    stop("`n_sim_per_draw` must be >= 1.")
  }

  fit_pieces <- validate_fit_for_params(fit)
  unc_pieces <- validate_fit_for_uncertainty(fit, method)

  # A refusal is an R error, classed like a refused `ferx_fit()` (#385).
  res <- .ferx_engine_call(
    ferx_rust_simulate_with_uncertainty(
      model_path           = normalizePath(model),
      data_path            = normalizePath(data),
      theta                = fit_pieces$theta,
      omega_flat           = fit_pieces$omega_flat,
      omega_dim            = fit_pieces$omega_dim,
      sigma                = fit_pieces$sigma,
      omega_iov_flat       = fit_pieces$omega_iov_flat,
      omega_iov_dim        = fit_pieces$omega_iov_dim,
      method               = method,
      cov_matrix_flat      = unc_pieces$cov_matrix_flat,
      cov_matrix_dim       = unc_pieces$cov_matrix_dim,
      sir_resamples_flat   = unc_pieces$sir_resamples_flat,
      sir_resamples_n      = unc_pieces$sir_resamples_n,
      sir_resamples_dim    = unc_pieces$sir_resamples_dim,
      residual_rho         = fit_pieces$residual_rho,
      n_uncertainty_draws  = as.integer(n_uncertainty_draws),
      n_sim_per_draw       = as.integer(n_sim_per_draw),
      seed                 = as.integer(seed)
    ),
    model, data
  )

  # Same `simulation_warnings` channel `ferx_simulate()` uses - here it carries
  # the design-point count (a kept empty-DV record; see the `data` note above),
  # which otherwise diverges silently from what `ferx_fit()` scored.
  res <- .ferx_surface_sim_warnings(res, "ferx_simulate_with_uncertainty")

  if (method == "asymptotic") {
    .ferx_warn_logit_probability_draws(fit, "ferx_simulate_with_uncertainty")
  }
  res
}

# Internal: share of asymptotic draws that put each `logit_probability` theta
# at or above 1 (#373).
#
# The engine draws in its packed space, where a theta with a non-negative lower
# bound is packed as `log(theta)`, and `fit$cov_matrix` is the covariance in
# that space. So a `logit_probability` theta - declared on (0, 1) - is drawn
# log-normally: `exp(log(theta) + sd * Z)` stays above 0 but has no ceiling at
# 1. The share that reaches 1 is `P(log(theta) + sd * Z >= 0)`, i.e.
# `pnorm(log(theta) / sd)`.
#
# Returns a named numeric vector (theta name -> share), with one entry per
# `logit_probability` theta that carries a positive, finite variance; empty
# when there is none (no such theta, a FIXed one, or no covariance matrix).
.ferx_logit_probability_draw_share <- function(fit) {
  tf    <- fit$theta_transforms
  theta <- fit$theta
  cov   <- fit$cov_matrix
  if (is.null(tf) || is.null(theta) || !is.matrix(cov)) return(numeric(0))
  idx <- which(as.character(tf) == "logit_probability")
  # `theta` occupies the leading segment of the packed vector, so theta i is
  # row/column i of `cov_matrix`.
  idx <- idx[idx <= length(theta) & idx <= nrow(cov)]
  if (length(idx) == 0L) return(numeric(0))
  est <- as.numeric(theta[idx])
  sd  <- sqrt(pmax(diag(cov)[idx], 0))
  ok  <- is.finite(est) & est > 0 & est < 1 & is.finite(sd) & sd > 0
  if (!any(ok)) return(numeric(0))
  nms <- names(theta)[idx]
  if (is.null(nms)) nms <- paste0("THETA", idx)
  stats::setNames(stats::pnorm(log(est[ok]) / sd[ok]), nms[ok])
}

# Internal: warn when asymptotic draws of a `logit_probability` theta leave
# (0, 1) often enough to matter (#373). What becomes of such a draw is decided
# by the engine, not here, and neither outcome is a draw from the intended
# distribution:
#   - declared upper bound below 1 (the bundled `bioavailability` example):
#     the whole parameter vector is rejected and redrawn, so the theta's
#     distribution is truncated - the draws stay in (0, 1) but are pulled low;
#   - declared upper bound at or above 1 (an undeclared bound defaults to 1e9):
#     the draw is simulated, and the model's `logit()` clamps its argument to
#     `1 - 1e-15`, so every subject in that draw gets a probability of ~1.
# Sampling such a theta on the logit scale needs a change in the engine's
# sampler; until then the warning names the theta, the share, and the ways
# round it.
.ferx_warn_logit_probability_draws <- function(fit, fn, threshold = 1e-3) {
  share <- .ferx_logit_probability_draw_share(fit)
  share <- share[share >= threshold]
  if (length(share) == 0L) return(invisible(share))
  what <- paste0(names(share), " ", sprintf("%.1f%%", 100 * share),
                 collapse = ", ")
  msg <- paste0(
    fn, "(method = \"asymptotic\"): expected share of draws at or above 1 ",
    "for logit_probability theta ", what, ". The engine draws thetas on ",
    "the log scale from `fit$cov_matrix`, so nothing holds such a theta below ",
    "1: a draw ",
    "past a declared upper bound is rejected (truncating the distribution), ",
    "and one at or above 1 is simulated with the probability clamped to 1. ",
    "Use method = \"sir\", which resamples the proposal by likelihood ",
    "weight, or declare the parameter on the logit scale - F = inv_logit(LOGIT_F + ETA_F) ",
    "with a negative lower bound on LOGIT_F - so it is drawn there."
  )
  warning(warningCondition(msg, class = "ferx_logit_probability_draws",
                           call = NULL))
  invisible(share)
}

# Internal: pull the uncertainty payload out of a ferx_fit result, validate it
# matches the requested method, and return flat representations for FFI. The
# empty branches give the Rust side something to ignore - it switches on
# `method` and only inspects the relevant arrays.
validate_fit_for_uncertainty <- function(fit, method) {
  if (method == "asymptotic") {
    cov <- fit$cov_matrix
    if (is.null(cov) || length(cov) == 0L) {
      stop("`fit$cov_matrix` is empty - re-fit with `covariance = TRUE` for ",
           "asymptotic uncertainty.")
    }
    # `cov` is a square R matrix after `process_fit_result()`. Flatten
    # row-major to match the engine's `DMatrix::from_row_slice` reader on
    # the Rust side - same convention used for `omega` above. The
    # transpose is mathematically a no-op for any well-formed covariance
    # (which must be symmetric), but the explicit `t()` keeps the FFI
    # contract uniform with `omega_flat` and prevents a future footgun
    # if the Rust reader ever changes.
    if (!is.matrix(cov) || nrow(cov) != ncol(cov)) {
      stop("`fit$cov_matrix` must be a square matrix.")
    }
    list(
      cov_matrix_flat    = as.numeric(t(cov)),
      cov_matrix_dim     = as.integer(nrow(cov)),
      sir_resamples_flat = numeric(0),
      sir_resamples_n    = 0L,
      sir_resamples_dim  = 0L
    )
  } else {
    resamples <- fit$sir_resamples
    n <- fit$sir_resamples_n
    d <- fit$sir_resamples_dim
    if (is.null(resamples) || length(resamples) == 0L ||
        is.null(n) || n == 0L || is.null(d) || d == 0L) {
      stop("`fit$sir_resamples` is empty - re-fit with `sir = TRUE` and ",
           "`sir_keep_samples = TRUE` in `settings` for SIR uncertainty.")
    }
    list(
      cov_matrix_flat    = numeric(0),
      cov_matrix_dim     = 0L,
      sir_resamples_flat = as.numeric(resamples),
      sir_resamples_n    = as.integer(n),
      sir_resamples_dim  = as.integer(d)
    )
  }
}
