#' Kaplan-Meier mean covariate VPC for a time-to-event model
#'
#' Checks whether a single-event time-to-event (TTE) model captures how a
#' baseline covariate drives the hazard (Hooker and Karlsson, PAGE 21 (2012)
#' Abstr 2564). At every step of the Kaplan-Meier curve it takes the mean of
#' the covariate over the subjects still at risk, in the data and in each of
#' \code{n_sim} replicates of the observed design, and returns the observed
#' curve beside the pointwise band of the replicate curves, per stratum.
#' Everything is computed by ferx-core (\code{vpc_kmmc()}, ferx-core #1886);
#' this function only carries the result into R. It never refits.
#'
#' @section Reading it:
#' A covariate that raises the hazard makes its at-risk mean drift as the
#' high-hazard subjects leave the risk set. A model that lacks the effect
#' drifts too little, and the observed curve leaves the band. Every replicate
#' keeps the data's covariates, so the curves differ only in who stays at
#' risk: at time 0 the observed value and the whole band are the population
#' mean, and the signal is the gap that opens later.
#'
#' @section Censoring in the replicates:
#' Under \code{censoring = "observed"} every replicate subject draws its
#' censoring time from its own stratum's reverse Kaplan-Meier censoring
#' distribution, unconditionally. Unlike \code{\link{ferx_vpc_tte}}, a subject
#' observed censored does not keep its own time: censored subjects are
#' selected for a low hazard, so keeping their times would remove the same
#' covariate values at the same times in every replicate as in the data, and
#' the band would hug the observed curve. The draw assumes censoring is
#' independent of the covariate within a stratum, so stratify on what drives
#' dropout. The observed records are read on \code{[0, H]} only.
#'
#' @section What it refuses:
#' Every input \code{\link{ferx_vpc_tte}} refuses, as an R error of class
#' \code{ferx_engine_error} with \code{e$code == "E_VPC_TTE_INPUT"}, and also:
#' a repeated-events endpoint (use \code{\link{ferx_vpc_tte}}, which draws
#' its RTTE VPC); a subject on the endpoint that lacks the covariate or whose
#' baseline value is missing; a \code{level} no such subject has; and
#' stratifying on the covariate itself. A covariate that changes within a
#' subject uses each subject's first value, with a
#' \code{W_VPC_KMMC_TIME_VARYING_COVARIATE} warning.
#'
#' @inheritParams ferx_vpc_tte
#' @param covariate The baseline covariate whose at-risk mean is taken: each
#'   subject's first non-missing value.
#' @param level \code{NULL} (default) for the arithmetic mean (for a 0/1
#'   covariate, the share coded 1), or one value: the curve is then the share
#'   of at-risk subjects at exactly that value, for a multi-level code whose
#'   mean means nothing.
#' @param min_at_risk A curve is undefined where fewer subjects than this are
#'   at risk; applied to the observed curve and to each replicate on its own.
#'
#' @return An object of class \code{ferx_vpc_kmmc}, a list with:
#'   \describe{
#'     \item{\code{observed}}{The observed curve, one row per row of each
#'       stratum's Kaplan-Meier curve: \code{STRATUM}, \code{time},
#'       \code{n_risk}, \code{n_event}, \code{n_censor} and \code{mean}, the
#'       covariate's mean (or share at \code{level}) over the \code{n_risk}
#'       subjects at risk at \code{time} (event or censoring time
#'       \eqn{\ge} \code{time}). The curve is left-continuous: a row's value
#'       holds from the previous row's time up to and including its own.}
#'     \item{\code{band}}{The replicate band, one row per grid time per
#'       stratum: \code{STRATUM}, \code{time}, \code{obs} (the observed curve
#'       at \code{time}; \code{NA} where fewer than \code{min_at_risk}
#'       subjects are at risk), \code{lo}, \code{med}, \code{hi} and
#'       \code{n_used}.}
#'     \item{\code{strata}}{One row per stratum, as in
#'       \code{\link{ferx_vpc_tte}}.}
#'     \item{\code{covariate}, \code{level}, \code{min_at_risk}, \code{cmt},
#'       \code{horizon}, \code{n_sim}, \code{probs}, \code{censoring}}{The
#'       settings used and the horizon \code{H}.}
#'     \item{\code{warnings}}{The engine's warnings, also raised as an R
#'       warning.}
#'   }
#'
#' @examples
#' ex <- ferx_example("tte_weibull")
#' d <- utils::read.csv(ex$data)
#' d$WT <- 50 + 3 * (d$ID %% 10)
#' kv <- ferx_vpc_kmmc(ex$model, d, covariate = "WT", n_sim = 50L, seed = 1L)
#' kv
#' plot(kv)
#'
#' @seealso \code{\link{plot.ferx_vpc_kmmc}}, \code{\link{ferx_vpc_tte}}
#' @family simulation
#' @export
ferx_vpc_kmmc <- function(model, data = NULL, covariate, level = NULL,
                          min_at_risk = 10L, n_sim = 200L, seed = 42L,
                          fit = NULL, strata = NULL, by_dose = FALSE,
                          censoring = c("observed", "horizon"), cmt = NULL,
                          horizon = NULL, probs = c(0.05, 0.5, 0.95),
                          grid = "auto") {
  censoring <- match.arg(censoring)
  if (missing(covariate) || !is.character(covariate) || length(covariate) != 1L ||
      is.na(covariate) || !nzchar(covariate)) {
    stop("`covariate` must be a single covariate name.", call. = FALSE)
  }
  if (!is.null(level) &&
      (!is.numeric(level) || length(level) != 1L || !is.finite(level))) {
    stop("`level` must be NULL or a single finite number.", call. = FALSE)
  }
  .ferx_vpc_tte_count(min_at_risk, "min_at_risk", min = 1)
  args <- .ferx_vpc_args(model, data, fit, n_sim, seed, strata, by_dose,
                         censoring, cmt, horizon, probs, grid)
  extra <- list(
    covariate = covariate,
    # NaN is the glue's "the mean" (extendr refuses an NA scalar).
    level = if (is.null(level)) NaN else as.numeric(level),
    min_at_risk = as.integer(min_at_risk)
  )
  # The engine points a repeated-events endpoint at its own `vpc_rtte`, which
  # ferx_vpc_tte() calls for such an endpoint.
  res <- .ferx_engine_spelling(
    .ferx_vpc_call(args, fit, ferx_rust_vpc_kmmc, ferx_rust_vpc_kmmc_from_fit,
                   extra, "ferx_vpc_kmmc"),
    "E_VPC_TTE_INPUT", "`vpc_rtte`", "`ferx_vpc_tte()`"
  )
  structure(
    list(
      observed = .ferx_vpc_stack(res$strata, "observed"),
      band = .ferx_vpc_stack(res$strata, "band"),
      strata = .ferx_vpc_strata(res$strata),
      covariate = res$covariate,
      level = as.numeric(res$level),
      min_at_risk = res$min_at_risk,
      cmt = res$cmt,
      horizon = res$horizon,
      n_sim = res$n_sim,
      probs = res$probs,
      censoring = res$censoring,
      warnings = as.character(attr(res, "simulation_warnings", exact = TRUE))
    ),
    class = "ferx_vpc_kmmc"
  )
}

# The y-axis label of a KMMC curve: the mean, or the share at `level`.
.ferx_vpc_kmmc_label <- function(x) {
  if (is.na(x$level)) {
    paste("Mean", x$covariate, "at risk")
  } else {
    sprintf("Share at risk with %s = %s", x$covariate, format(x$level))
  }
}

#' @export
print.ferx_vpc_kmmc <- function(x, ...) {
  pct <- paste0(format(100 * x$probs, trim = TRUE), "%")
  cat("Kaplan-Meier mean covariate VPC (ferx_vpc_kmmc)\n")
  cat("  ", .ferx_vpc_kmmc_label(x), "; at least ", x$min_at_risk,
      " at risk\n", sep = "")
  cat("  endpoint cmt ", x$cmt, ", horizon ", format(x$horizon), ", ",
      x$n_sim, " replicates, ", x$censoring, " censoring\n", sep = "")
  cat("  band ", pct[1], "-", pct[3], " with the ", pct[2], " line\n",
      sep = "")
  for (i in seq_len(nrow(x$strata))) {
    s <- as.character(x$strata$STRATUM[i])
    obs <- x$observed[x$observed$STRATUM == s, , drop = FALSE]
    cat(sprintf("  %s: %d subjects, %d events\n", s, x$strata$n_subjects[i],
                sum(obs$n_event)))
  }
  .ferx_vpc_print_warnings(x$warnings)
  invisible(x)
}

#' Plot a Kaplan-Meier mean covariate VPC
#'
#' Draws each stratum's observed Kaplan-Meier mean covariate curve over its
#' replicate band from \code{\link{ferx_vpc_kmmc}}, one panel per stratum, with
#' base graphics. The band is a step ribbon: the band at grid time \code{t[k]}
#' holds until \code{t[k + 1]}. The observed curve stops where fewer than
#' \code{min_at_risk} subjects are at risk.
#'
#' @param x A \code{ferx_vpc_kmmc} object.
#' @param y Unused.
#' @param ... Passed to \code{graphics::plot()}.
#' @return \code{x}, invisibly.
#' @examples
#' ex <- ferx_example("tte_weibull")
#' d <- utils::read.csv(ex$data)
#' d$WT <- 50 + 3 * (d$ID %% 10)
#' plot(ferx_vpc_kmmc(ex$model, d, covariate = "WT", n_sim = 50L, seed = 1L))
#' @export
plot.ferx_vpc_kmmc <- function(x, y = NULL, ...) {
  panels <- levels(x$strata$STRATUM)
  n <- length(panels)
  nc <- ceiling(sqrt(n))
  old <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (s in panels) {
    o <- x$observed[x$observed$STRATUM == s & x$observed$n_risk >= x$min_at_risk, ,
                    drop = FALSE]
    # Left-continuous (core's `KmmcCurve::at`): a row's mean holds up to and
    # including its time, and the first row's from time 0, so the curve is
    # drawn as an upper step.
    start <- if (nrow(o)) o$mean[1L] else NA_real_
    .ferx_vpc_step_panel(x$band[x$band$STRATUM == s, , drop = FALSE],
                         o$time, o$mean, x$horizon, start = start,
                         ylab = .ferx_vpc_kmmc_label(x),
                         main = if (n > 1L) s else "", step = "S", ...)
  }
  invisible(x)
}
