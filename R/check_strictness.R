#' Eligibility gates for a candidate fit
#'
#' \code{fit$converged} is a single flag. It does not distinguish a genuine
#' optimum from a run that never left its initial estimates, a parameter pinned
#' to a declared bound, an ill-conditioned covariance step or a near-singular
#' correlation matrix. Under automation all of those are model-selection
#' errors: a candidate that never moved is ranked on an OFV that says nothing
#' about the model. \code{check_strictness()} evaluates the pyDarwin-style
#' gates and returns the reasons, so a search report can say why a candidate
#' was excluded.
#'
#' The defaults are pyDarwin's posture. Set every argument to \code{FALSE} /
#' \code{NULL} for a "rank everything, report everything" run.
#'
#' The two threshold gates are evaluated only when their input exists: a fit
#' without a covariance matrix has no condition number to test. Such a gate is
#' reported under \code{skipped} rather than being failed or silently passed;
#' combine it with \code{require_covariance = TRUE} to make it mandatory.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} or
#'   \code{\link{ferx_load_fit}}.
#' @param require_converged Fail a fit with \code{converged = FALSE}. This
#'   already covers an internal runaway-guard hit, which demotes
#'   \code{converged}.
#' @param require_covariance Fail unless the covariance step ran and produced
#'   uncertainty. A SIR fallback counts as uncertainty (it is what the user
#'   configured that fallback to produce) unless its importance weights
#'   collapsed onto a single draw, i.e. an effective sample size below 2, in
#'   which case every interval has zero width and it fails like a failed step.
#' @param max_condition_number Fail when \code{fit$condition_number} (largest
#'   over smallest eigenvalue of the free-parameter correlation matrix) exceeds
#'   this. \code{NULL} disables the gate.
#' @param max_correlation Fail when any off-diagonal of the covariance matrix's
#'   correlation form exceeds this in absolute value, on the natural theta /
#'   OMEGA / SIGMA scale. \code{NULL} disables the gate.
#' @param reject_on_boundary Fail a fit with a theta pinned to a declared
#'   bound - the same predicate \code{\link{ferx_bootstrap}} applies as
#'   \code{skip_estimate_near_boundary}.
#' @param reject_init_stall Fail a fit that never left its initial estimates.
#'   This gate assumes a cold start from the model file's initial estimates.
#'   A candidate warm-started from its parent's final estimates legitimately
#'   converges within a tolerance of where it began and the fit carries nothing
#'   that separates the two, so turn this off for such a search and rely on
#'   \code{require_converged}.
#' @return A list with elements:
#'   \describe{
#'     \item{passed}{\code{TRUE} when no gate failed.}
#'     \item{failures}{Character vector, one entry per failed gate. Empty when
#'       the fit passed.}
#'     \item{skipped}{Character vector, one entry per gate that had no input to
#'       judge by.}
#'   }
#' @section Differences from the engine:
#' \code{ferx_core::model_selection::check_strictness()} names the two
#' parameters behind a failed correlation gate; the R fit carries the maximum
#' as a scalar, so the message here reports the value only. A condition number
#' that is present but \code{NaN} is reported as skipped rather than failed,
#' since an absent one reaches R as \code{NA} too.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' v <- check_strictness(fit)
#' if (!v$passed) print(v$failures)
#'
#' # Rank everything, report nothing as ineligible.
#' check_strictness(fit, require_converged = FALSE, max_condition_number = NULL,
#'                  max_correlation = NULL, reject_on_boundary = FALSE,
#'                  reject_init_stall = FALSE)
#' }
#' @seealso \code{\link{ferx_bic}} for the criteria a search ranks the eligible
#'   candidates on.
#' @family diagnostics
#' @export
check_strictness <- function(fit,
                             require_converged = TRUE,
                             require_covariance = FALSE,
                             max_condition_number = 1000,
                             max_correlation = 0.95,
                             reject_on_boundary = TRUE,
                             reject_init_stall = TRUE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit())",
         call. = FALSE)
  }
  # Validate the switches before reading the fit. `isTRUE()` silently treats
  # anything that is not a length-1 TRUE as "gate off", so a malformed search
  # configuration ("TRUE", 1, NA) would disable a gate and rank an ineligible
  # candidate rather than erroring.
  require_converged  <- .ferx_gate_flag(require_converged, "require_converged")
  require_covariance <- .ferx_gate_flag(require_covariance, "require_covariance")
  reject_on_boundary <- .ferx_gate_flag(reject_on_boundary, "reject_on_boundary")
  reject_init_stall  <- .ferx_gate_flag(reject_init_stall, "reject_init_stall")

  failures <- character(0)
  skipped  <- character(0)

  if (require_converged && !isTRUE(fit$converged)) {
    failures <- c(failures, "did not converge (converged = FALSE)")
  }

  if (require_covariance) {
    reason <- .ferx_covariance_gate(fit)
    if (!is.null(reason)) failures <- c(failures, reason)
  }

  max_cn <- .ferx_gate_threshold(max_condition_number, "max_condition_number")
  if (!is.null(max_cn)) {
    # ferx_fit() and ferx_covariance() rename the wire field to
    # `condition_number` (and clear the wire one) via
    # .ferx_apply_cov_sentinels(); ferx_load_fit() does not, and leaves the
    # wire name in place. Read whichever the fit carries.
    cn <- .ferx_gate_value(fit$condition_number %||% fit$cov_condition_number)
    if (is.null(cn)) {
      skipped <- c(skipped, paste(
        "condition number: no covariance matrix",
        "(covariance step not run or failed)"))
    } else if (cn > max_cn) {
      failures <- c(failures, sprintf(
        "condition number %.4e exceeds %.4e", cn, max_cn))
    }
  }

  max_r <- .ferx_gate_threshold(max_correlation, "max_correlation")
  if (!is.null(max_r)) {
    r <- .ferx_gate_value(fit$max_abs_correlation)
    if (is.null(r)) {
      skipped <- c(skipped, paste(
        "parameter correlation: no covariance matrix with two or more",
        "free parameters"))
    } else if (r > max_r) {
      failures <- c(failures, sprintf(
        "parameter correlation |r| = %.4f exceeds %.4f", r, max_r))
    }
  }

  if (reject_on_boundary) {
    # Same "no verdict" shape as the init-stall gate below: the predicate reads
    # structured warning details the .fitrx wire drops, so a bundle written
    # before the verdict existed carries NA rather than FALSE. Reporting that
    # as a pass would silently clear the gate for every older bundle.
    boundary <- fit$estimate_near_boundary
    if (is.null(boundary) || length(boundary) != 1L || is.na(boundary)) {
      skipped <- c(skipped, paste(
        "boundary: fit carries no boundary verdict (saved before the",
        "predicate existed)"))
    } else if (isTRUE(as.logical(boundary))) {
      failures <- c(failures, sprintf(
        "estimate pinned to a declared bound: %s", .ferx_boundary_detail(fit)))
    }
  }

  if (reject_init_stall) {
    stalled <- fit$stalled_at_init
    if (is.null(stalled) || length(stalled) != 1L || is.na(stalled)) {
      skipped <- c(skipped, paste(
        "init stall: fit carries no initial estimates (or no optimizer",
        "verdict) to compare against"))
    } else if (isTRUE(as.logical(stalled))) {
      failures <- c(failures, sprintf(paste(
        "stalled at the initial estimates: no free parameter moved more",
        "than %.0f%% of its initial value"), .ferx_init_stall_rel_tol * 100))
    }
  }

  list(passed = length(failures) == 0L,
       failures = failures,
       skipped = skipped)
}

# ferx-core's INIT_STALL_REL_TOL. Only used to word the failure; the verdict
# itself is the engine's (fit$stalled_at_init), never recomputed here.
.ferx_init_stall_rel_tol <- 1e-2

# ferx-core's MIN_SIR_FALLBACK_ESS: below this the importance weights have
# collapsed onto a single draw and every credible interval has zero width.
.ferx_min_sir_fallback_ess <- 2

# A single logical gate switch. Anything else is an error rather than a silent
# "gate off" - a search configuration that disables a gate by typo ranks
# candidates the gate exists to exclude.
.ferx_gate_flag <- function(x, arg) {
  if (length(x) != 1L || !is.logical(x) || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE", arg), call. = FALSE)
  }
  x
}

# A gate threshold, or NULL when the caller disabled the gate. An explicit NA
# disables it too - a threshold nothing can exceed is not a gate - but a value
# that merely *coerces* to NA ("typo") is a malformed argument, not a request
# to switch the gate off.
.ferx_gate_threshold <- function(x, arg) {
  if (is.null(x)) return(NULL)
  if (length(x) != 1L) {
    stop(sprintf("`%s` must be a single number or NULL", arg), call. = FALSE)
  }
  if (is.na(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) {
    stop(sprintf("`%s` must be a single number or NULL, not %s",
                 arg, dQuote(as.character(x), FALSE)), call. = FALSE)
  }
  v
}

# The value a threshold gate tests, or NULL when the fit has none to judge by.
# The FFI ships an absent value as NaN and a .fitrx round-trip can turn it into
# NA, so both read as "no input" here.
.ferx_gate_value <- function(x) {
  v <- suppressWarnings(as.numeric(x %||% NA_real_))
  if (length(v) != 1L || is.na(v)) return(NULL)
  v
}

# Normalise a covariance status to a comparison key: ferx_fit() carries the
# engine's snake_case token ("not_requested"), while ferx_load_fit() maps it to
# a CamelCase label ("NotRequested").
.ferx_cov_status_key <- function(status) {
  s <- as.character(status %||% "")
  if (length(s) != 1L || !nzchar(s)) return("unknown")
  tolower(gsub("_", "", s, fixed = TRUE))
}

# The reason require_covariance fails a fit, or NULL when it is satisfied.
.ferx_covariance_gate <- function(fit) {
  switch(.ferx_cov_status_key(fit$covariance_status),
    computed = if (length(fit$cov_matrix %||% numeric(0)) == 0L) {
      "covariance step reported computed but stored no matrix"
    } else {
      NULL
    },
    notrequested = "covariance step not requested (covariance = FALSE)",
    failed = "covariance step failed",
    sirfallback = {
      ess <- .ferx_gate_value(fit$sir_ess)
      if (is.null(ess)) {
        "SIR fallback recorded no effective sample size to judge it by"
      } else if (ess < .ferx_min_sir_fallback_ess) {
        sprintf(paste("SIR fallback collapsed: effective sample size %.2f is",
                      "below %g (every credible interval has zero width)"),
                ess, .ferx_min_sir_fallback_ess)
      } else {
        NULL
      }
    },
    "covariance step produced no uncertainty"
  )
}

# What was pinned to a bound, for the failure message. The engine records the
# offending parameters in a structured warning; a fit read back from a .fitrx
# keeps only the message text, and an older bundle may keep neither.
.ferx_boundary_detail <- function(fit) {
  ws <- fit$warnings_structured
  if (is.data.frame(ws) && nrow(ws) > 0L && !is.null(ws$category)) {
    hit <- ws$message[ws$category == "boundary_estimate"]
    hit <- hit[nzchar(as.character(hit))]
    if (length(hit) > 0L) return(as.character(hit)[[1L]])
  }
  "see fit$warnings"
}
