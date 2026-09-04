# Strictness gate for candidate ranking (ferx-core #1177).
#
# The gates themselves are evaluated here rather than in the engine: the R fit
# object is a list, not a `FitResult`, so `ferx_core::check_strictness` cannot
# be pointed at one. What the engine does supply - at fit time, and through
# `.fitrx` - is every ingredient the gates need that R could not re-derive:
# `max_abs_correlation` (over the natural-scale covariance), `near_boundary`
# and its warning text, and `stalled_at_init`. Keep the gate list and the
# defaults in step with `ferx_core::Strictness`.

# The engine reports the status lowercase on a fresh fit ("computed") and
# CamelCase on one rebuilt from a bundle ("Computed"). Fold both.
.ferx_covariance_status_token <- function(status) {
  s <- as.character(status %||% "not_requested")
  s <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", s)
  tolower(s)
}

# Kish effective sample size below which a SIR fallback counts as having
# delivered no uncertainty. Mirrors ferx_core::MIN_SIR_FALLBACK_ESS.
.FERX_MIN_SIR_FALLBACK_ESS <- 2

# Relative displacement below which a free parameter counts as "still at its
# initial value". Mirrors ferx_core::INIT_STALL_REL_TOL; used for the message
# only, since the verdict itself comes from the engine.
.FERX_INIT_STALL_REL_TOL <- 0.01

# A threshold gate is off when its argument is NULL or NA.
.ferx_threshold <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (length(x) != 1L) stop("`", what, "` must be a single number, NA, or NULL.")
  if (is.na(x)) return(NULL)
  if (!is.numeric(x)) stop("`", what, "` must be a single number, NA, or NULL.")
  as.numeric(x)
}

#' Check a fit against the gates a candidate must pass to be ranked
#'
#' Before a criterion such as \code{\link{ferx_bic}} is trusted for model
#' selection, the fit behind it has to be sound. This applies the pyDarwin gate
#' set - convergence, optionally an uncertainty step, a condition-number and a
#' parameter-correlation ceiling, no estimate pinned to a bound, and no run
#' that never left its initial estimates - and reports which gates failed.
#'
#' The two threshold gates are evaluated only when their input exists: a fit
#' without a covariance matrix has no condition number to test. Such a gate is
#' reported in \code{skipped} rather than passed silently or failed; combine it
#' with \code{require_covariance = TRUE} to make it mandatory.
#'
#' Note that \code{reject_init_stall} assumes a \strong{cold start} from the
#' model file's initial estimates, where "did not move" means the OFV is the
#' OFV of the initial estimates. A candidate warm-started from its parent's
#' final estimates legitimately converges close to where it began, and nothing
#' on the fit distinguishes the two - turn the gate off for such a search and
#' rely on \code{require_converged}.
#'
#' @param fit A \code{ferx_fit} object.
#' @param require_converged Fail a fit with \code{converged = FALSE}. This
#'   already covers an internal runaway-guard hit, which demotes convergence.
#' @param require_covariance Fail unless the covariance step produced
#'   uncertainty: \code{"computed"} with a stored matrix, or a SIR fallback
#'   whose importance weights did not collapse onto a single draw (effective
#'   sample size at least 2).
#' @param max_condition_number Fail when \code{fit$cov_condition_number}
#'   exceeds this. \code{NULL} or \code{NA} turns the gate off.
#' @param max_correlation Fail when any parameter correlation exceeds this in
#'   absolute value. \code{NULL} or \code{NA} turns the gate off. The
#'   correlation is read from the natural-scale covariance the engine records
#'   (\code{fit$max_abs_correlation}), which for a \code{block_omega} model is
#'   not the largest off-diagonal of \code{fit$cor_matrix}.
#' @param reject_on_boundary Fail a fit with a theta pinned to a declared
#'   optimizer bound.
#' @param reject_init_stall Fail a fit that never left its initial estimates.
#' @return A list with \code{passed} (single logical, \code{TRUE} when
#'   \code{failures} is empty), \code{failures} (character vector, one entry
#'   per failed gate) and \code{skipped} (character vector, one entry per gate
#'   that was enabled but had no input to test on this fit).
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn")
#' v <- check_strictness(fit)
#' v$passed
#' v$failures
#'
#' # Rank everything, report everything:
#' check_strictness(fit,
#'   require_converged = FALSE, max_condition_number = NULL,
#'   max_correlation = NULL, reject_on_boundary = FALSE,
#'   reject_init_stall = FALSE
#' )
#' @family utilities
#' @seealso \code{\link{ferx_bic}} for the criterion these gates protect, and
#'   \code{\link{check_diagnostics}} for the wider post-fit diagnostic sweep.
#' @export
check_strictness <- function(fit,
                             require_converged = TRUE,
                             require_covariance = FALSE,
                             max_condition_number = 1000,
                             max_correlation = 0.95,
                             reject_on_boundary = TRUE,
                             reject_init_stall = TRUE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }
  max_cn <- .ferx_threshold(max_condition_number, "max_condition_number")
  max_r <- .ferx_threshold(max_correlation, "max_correlation")

  failures <- character()
  skipped <- character()

  if (isTRUE(require_converged) && !isTRUE(fit$converged)) {
    failures <- c(failures, "did not converge (`converged = FALSE`)")
  }

  if (isTRUE(require_covariance)) {
    status <- .ferx_covariance_status_token(fit$covariance_status)
    reason <- switch(status,
      "computed" = if (is.null(fit$cov_matrix)) {
        "covariance step reported computed but stored no matrix"
      },
      "not_requested" = "covariance step not requested (`covariance = FALSE`)",
      "failed" = "covariance step failed",
      # Accepted when the user opted into SIR as the uncertainty source for a
      # non-PD Hessian and it delivered - which a run whose weights collapsed
      # onto one draw did not, whatever it returned.
      "sir_fallback" = {
        ess <- suppressWarnings(as.numeric(fit$sir_ess %||% NA_real_))
        if (length(ess) != 1L || !is.finite(ess)) {
          "SIR fallback recorded no effective sample size to judge it by"
        } else if (ess < .FERX_MIN_SIR_FALLBACK_ESS) {
          sprintf(
            paste0(
              "SIR fallback collapsed: effective sample size %.2f is below %g ",
              "(every credible interval has zero width)"
            ),
            ess, .FERX_MIN_SIR_FALLBACK_ESS
          )
        }
      },
      paste0("covariance status not recognised: ", status)
    )
    if (!is.null(reason)) failures <- c(failures, reason)
  }

  if (!is.null(max_cn)) {
    cn <- suppressWarnings(as.numeric(fit$cov_condition_number %||% NA_real_))
    # The engine reports NaN both for "no covariance matrix" and for a matrix
    # whose condition number could not be formed; `cov_matrix` separates them.
    if (is.null(fit$cov_matrix) || length(cn) != 1L || is.na(cn)) {
      if (is.null(fit$cov_matrix)) {
        skipped <- c(
          skipped,
          "condition number: no covariance matrix (covariance step not run or failed)"
        )
      } else {
        failures <- c(failures, "condition number is not a number")
      }
    } else if (cn > max_cn) {
      failures <- c(failures, sprintf(
        "condition number %s exceeds %s", format(cn, digits = 5), format(max_cn, digits = 5)
      ))
    }
  }

  if (!is.null(max_r)) {
    r <- suppressWarnings(as.numeric(fit$max_abs_correlation %||% NA_real_))
    if (length(r) != 1L || !is.finite(r)) {
      skipped <- c(
        skipped,
        "parameter correlation: no covariance matrix with two or more free parameters"
      )
    } else if (r > max_r) {
      failures <- c(failures, sprintf(
        "parameter correlation |r| = %.4f exceeds %.4f", r, max_r
      ))
    }
  }

  if (isTRUE(reject_on_boundary) && isTRUE(fit$near_boundary)) {
    # The engine names the pinned parameters from the warning's structured
    # payload, which does not cross into R; its message text stands in.
    what <- fit$boundary_estimate_message %||% ""
    if (!nzchar(what)) what <- "see fit$warnings"
    failures <- c(failures, paste0("estimate pinned to a declared bound: ", what))
  }

  if (isTRUE(reject_init_stall)) {
    stalled <- fit$stalled_at_init
    if (is.null(stalled) || length(stalled) != 1L || is.na(stalled)) {
      skipped <- c(
        skipped,
        paste0(
          "init stall: fit carries no initial estimates (or no optimizer ",
          "verdict) to compare against"
        )
      )
    } else if (isTRUE(stalled)) {
      failures <- c(failures, sprintf(
        paste0(
          "stalled at the initial estimates: no free parameter moved more ",
          "than %g%% of its initial value"
        ),
        .FERX_INIT_STALL_REL_TOL * 100
      ))
    }
  }

  list(
    passed = length(failures) == 0L,
    failures = failures,
    skipped = skipped
  )
}
