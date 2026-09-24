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
#' @section Bounded thetas and \code{logit_probability} draws:
#' The engine draws in its packed parameter space. A theta with a non-negative
#' lower bound is packed as \code{log(theta)}, so its asymptotic draws are
#' log-normal: always above 0, with no ceiling. A draw past a declared upper
#' bound is rejected and the whole parameter set redrawn, which keeps the
#' draws inside the bound but truncates their distribution.
#'
#' For a \code{"logit_probability"} theta - one declared on \eqn{(0, 1)} and
#' used as \code{inv_logit(logit(THETA) + ETA)} - the draws are therefore not
#' the logit-normal ones the fit implies:
#' \itemize{
#'   \item With a declared upper bound of 1 or below (the bundled
#'     \code{bioavailability} example uses 0.999), draws past it are rejected,
#'     so the distribution is truncated and pulled low.
#'   \item With a declared upper bound above 1 (an undeclared bound defaults to
#'     \code{1e9}), a draw above 1 is simulated, and the model's
#'     \code{logit()} clamps it so that every subject in that draw gets a
#'     probability of 1.
#' }
#' \code{method = "asymptotic"} warns (class
#' \code{"ferx_logit_probability_draws"}) when more than 0.1% of draws are
#' expected to be rejected or clamped this way, naming the theta, its declared
#' upper bound and the share. \code{method = "sir"} starts from the same
#' log-packed proposal: its likelihood weighting corrects the truncation, but
#' with an upper bound above 1 a clamped draw still gets a finite weight and
#' can reach the resample pool. So under \code{method = "sir"} the function
#' checks the pool itself and warns when any pooled draw is above 1.
#'
#' Declaring the probability on the logit scale avoids all of this:
#' \code{F = inv_logit(LOGIT_F + ETA_F)} with a negative lower bound on
#' \code{LOGIT_F} is packed, and so drawn, on the logit scale. Drawing a
#' \code{logit_probability} theta on the logit scale is an engine change
#' (\href{https://github.com/FeRx-NLME/ferx-core/issues/1548}{ferx-core #1548}).
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

  .ferx_warn_logit_probability_draws(
    fit, .ferx_theta_packing(model), method, "ferx_simulate_with_uncertainty"
  )
  res
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

# Internal: the model's declared theta bounds and transforms, as the engine
# parsed them (#373). NULL when the model cannot be read, so a diagnostic built
# on it degrades to nothing rather than failing a simulation that ran.
.ferx_theta_packing <- function(model) {
  if (!is.character(model) || length(model) != 1L || !file.exists(model)) {
    return(NULL)
  }
  info <- tryCatch(ferx_rust_theta_packing(normalizePath(model)),
                   error = function(e) NULL)
  if (is.null(info) || length(info$names) == 0L) return(NULL)
  info
}

# Internal: how often the uncertainty draws of each `logit_probability` theta
# leave (0, 1), and what the engine does with them (#373).
#
# The engine draws in its packed space: `log(theta)` when the declared lower
# bound is non-negative, the natural scale otherwise. `fit$cov_matrix` is the
# covariance there, so an asymptotic draw of a probability has no ceiling at 1.
# Then:
#   - declared upper <= 1: a draw past the bound is rejected and redrawn
#     (`outcome = "rejected"`), truncating the distribution;
#   - declared upper > 1: a draw in (1, upper] is simulated with the model's
#     `logit()` clamping it to a probability of 1 (`outcome = "clamped"`).
# Under `method = "asymptotic"` the share is the Gaussian tail of the proposal;
# under `method = "sir"` it is counted in the retained resample pool, which is
# where a clamped draw that won a likelihood weight ends up.
#
# `packing` is `.ferx_theta_packing()`'s list. Returns a data frame with one
# row per exposed theta (theta, upper, share, n, outcome); zero rows when there
# is nothing to report.
.ferx_logit_probability_draw_share <- function(fit, packing, method = "asymptotic") {
  empty <- data.frame(theta = character(0), upper = numeric(0),
                      share = numeric(0), n = integer(0),
                      outcome = character(0), stringsAsFactors = FALSE)
  theta <- fit$theta
  if (is.null(packing) || is.null(theta)) return(empty)
  n_theta <- min(length(theta), length(packing$transform))
  idx <- which(packing$transform[seq_len(n_theta)] == "logit_probability")
  if (length(idx) == 0L) return(empty)

  rows <- lapply(idx, function(i) {
    lower    <- packing$lower[i]
    upper    <- packing$upper[i]
    log_pack <- isTRUE(lower >= 0)
    clamped  <- isTRUE(upper > 1)
    if (method == "sir") {
      pool <- .ferx_sir_pool(fit)
      if (is.null(pool) || ncol(pool) < i) return(NULL)
      # Resampled draws already sit inside the box, so only clamping can show.
      if (!clamped) return(NULL)
      draws <- if (log_pack) exp(pool[, i]) else pool[, i]
      k <- sum(draws > 1)
      share <- k / length(draws)
      n <- length(draws)
    } else {
      cov <- fit$cov_matrix
      if (!is.matrix(cov) || nrow(cov) < i) return(NULL)
      est <- as.numeric(theta[i])
      sd  <- sqrt(max(cov[i, i], 0))
      if (!is.finite(est) || est <= 0 || !is.finite(sd) || sd <= 0) return(NULL)
      # P(draw > c) for the proposal centred on the estimate.
      above <- function(c) {
        if (!is.finite(c)) return(0)
        if (log_pack) stats::pnorm((log(est) - log(c)) / sd)
        else stats::pnorm((est - c) / sd)
      }
      # Clamped draws are the ones in (1, upper]; beyond that they are rejected.
      share <- if (clamped) above(1) - above(upper) else above(upper)
      n <- NA_integer_
    }
    nm <- names(theta)[i]
    if (is.null(nm) || !nzchar(nm)) nm <- packing$names[i]
    data.frame(theta = nm, upper = upper, share = share, n = n,
               outcome = if (clamped) "clamped" else "rejected",
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

# Internal: the retained SIR resample pool as an n x d matrix of packed
# vectors, or NULL.
.ferx_sir_pool <- function(fit) {
  v <- fit$sir_resamples
  n <- fit$sir_resamples_n
  d <- fit$sir_resamples_dim
  if (is.null(v) || is.null(n) || is.null(d) || n < 1L || d < 1L ||
      length(v) != n * d) {
    return(NULL)
  }
  matrix(as.numeric(v), nrow = n, ncol = d, byrow = TRUE)
}

# Internal: warn when uncertainty draws of a `logit_probability` theta leave
# (0, 1) often enough to matter (#373). Under `method = "asymptotic"` the
# threshold is an expected share of 0.1%; under `method = "sir"` any clamped
# draw in the pool is reported, since it is a draw that will be simulated.
.ferx_warn_logit_probability_draws <- function(fit, packing, method, fn,
                                               threshold = 1e-3) {
  ex <- .ferx_logit_probability_draw_share(fit, packing, method)
  keep <- if (method == "sir") ex$share > 0 else ex$share >= threshold
  ex <- ex[keep, , drop = FALSE]
  if (nrow(ex) == 0L) return(invisible(ex))
  what <- vapply(seq_len(nrow(ex)), function(r) {
    e <- ex[r, ]
    pct <- sprintf("%.1f%%", 100 * e$share)
    bound <- format(e$upper, digits = 4)
    if (method == "sir") {
      sprintf(paste0("%s: %d of %d pooled SIR draws (%s) are above 1 ",
                     "(declared upper bound %s), so logit() clamps each to a ",
                     "probability of 1 for every subject in that draw"),
              e$theta, as.integer(round(e$share * e$n)), e$n, pct, bound)
    } else if (e$outcome == "clamped") {
      sprintf(paste0("%s: about %s of draws land above 1 (declared upper ",
                     "bound %s), so logit() clamps each to a probability of 1 ",
                     "for every subject in that draw"), e$theta, pct, bound)
    } else {
      sprintf(paste0("%s: about %s of draws exceed the declared upper bound ",
                     "%s and are rejected, truncating the distribution below ",
                     "it"), e$theta, pct, bound)
    }
  }, character(1))
  msg <- paste0(
    fn, "(method = \"", method, "\"): draws of a logit_probability theta ",
    "leave (0, 1). The engine draws it on the log scale, not the logit scale ",
    "(ferx-core #1548).\n  ", paste(what, collapse = "\n  "), "\n",
    "Declare the probability on the logit scale - ",
    "F = inv_logit(LOGIT_F + ETA_F) with a negative lower bound on LOGIT_F - ",
    "so it is drawn there."
  )
  warning(warningCondition(msg, class = "ferx_logit_probability_draws",
                           call = NULL))
  invisible(ex)
}
