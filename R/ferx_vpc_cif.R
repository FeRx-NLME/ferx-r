#' Cumulative-incidence visual predictive check for competing risks
#'
#' Simulates the observed design \code{n_sim} times and returns, per stratum,
#' the observed Aalen-Johansen cumulative incidence of every cause and the
#' observed all-cause Kaplan-Meier curve, each beside the pointwise band of the
#' same curve over the replicates. Every time-to-event (TTE) endpoint of the
#' model is a cause. Everything is computed by ferx-core (\code{vpc_cif()},
#' ferx-core #1887); this function only carries the result into R. It never
#' refits.
#'
#' @section Why not one Kaplan-Meier VPC per cause:
#' One minus the Kaplan-Meier curve of one cause, with the other causes as
#' censorings, is not that cause's cumulative incidence: it estimates the
#' incidence in a world without the other causes, and overstates it. And the
#' all-cause curve cannot see a model that swaps the causes' rates while
#' keeping their sum; the cumulative incidences can. So
#' \code{\link{ferx_vpc_tte}} refuses a model with several TTE endpoints unless
#' \code{cmt} names one, and this function draws all of them.
#'
#' @section The data:
#' Every TTE endpoint must be single-event. A subject with any cause record
#' needs exactly one record on every cause, all at the same time (its first
#' event or its censoring), with at most one of them an event (\code{DV = 1});
#' the other causes are \code{DV = 0} there. That is the layout of
#' \code{ferx_example("tte_competing_risks")}. The simulator competes only the
#' causes a subject has records on, so a subject with fewer would be simulated
#' under fewer causes than the data imply. A model with one TTE endpoint is
#' accepted: its cumulative incidence is one minus the Kaplan-Meier curve of
#' \code{\link{ferx_vpc_tte}}.
#'
#' @section The curves and the bands:
#' The observed curves match \code{survival::survfit}'s multi-state fit
#' (\code{Surv(time, factor(status))} with the censoring as the first level)
#' row for row: \code{surv} is the all-cause Kaplan-Meier curve, and each
#' \code{cif_<cmt>} column the probability of having had that cause by
#' \code{time}. At every time the observed cumulative incidences and
#' \code{surv} add up to 1, but the band medians (or limits) of the causes and
#' of the survival need not: each band is pointwise, over the replicates still
#' under follow-up at that time, as in \code{\link{ferx_vpc_tte}}.
#'
#' @section Censoring in the replicates:
#' As in \code{\link{ferx_vpc_tte}}, on the all-cause time: under
#' \code{censoring = "observed"} a subject observed censored on every cause
#' keeps its own time, and a subject with an event draws a censoring time from
#' its stratum's reverse Kaplan-Meier curve given that it exceeds the event
#' time; a replicate's first event past its censoring time becomes a censoring
#' there. \code{"horizon"} censors every replicate at \code{H} only. The
#' observed records are read on \code{[0, H]} only.
#'
#' @section What it refuses:
#' An input the curves would be wrong on is an R error of class
#' \code{ferx_engine_error} with \code{e$code == "E_VPC_TTE_INPUT"}: a model
#' without a TTE endpoint; a repeated-events endpoint (use
#' \code{\link{ferx_vpc_tte}}, which draws its RTTE VPC); a subject with
#' records on some causes but not on all, more than one record on a cause,
#' cause records at different times, or events on two causes; left truncation;
#' an interval-censored record; a duplicated subject id; no record on any
#' cause; and the band, grid and stratification inputs
#' \code{\link{ferx_vpc_tte}} refuses.
#'
#' @inheritParams ferx_vpc_tte
#' @param model Path to a .ferx model file whose \code{[event_model]}
#'   endpoints are the competing causes.
#' @param data Path to a NONMEM-format CSV, or a data.frame, holding the
#'   \strong{observed} records: per subject one row on every cause's
#'   \code{CMT}, at the same time, \code{DV = 1} on the cause that occurred (if
#'   any) and \code{DV = 0} on the others. When omitted, the data \code{fit}
#'   was fitted to, otherwise the model file's \code{[data]} block.
#' @param horizon The time the replicates are simulated to. \code{NULL}
#'   (default) is the latest time recorded on a cause.
#' @param grid The times the bands are evaluated at: \code{"auto"} (0, the
#'   stratum's observed all-cause times and 101 evenly spaced points up to
#'   \code{H}), \code{"all_steps"} (also every step of every replicate curve),
#'   or a numeric vector of times.
#'
#' @return An object of class \code{ferx_vpc_cif}, a list with:
#'   \describe{
#'     \item{\code{observed}}{The observed curves, one row per row of each
#'       stratum's all-cause Kaplan-Meier curve: \code{STRATUM}, \code{time},
#'       \code{n_risk}, \code{n_event} (of any cause), \code{n_censor},
#'       \code{surv} (all-cause survival), then per cause \code{cif_<cmt>}
#'       (its cumulative incidence just after \code{time}) and
#'       \code{n_event_<cmt>} (its events at \code{time}), in \code{cmts}
#'       order.}
#'     \item{\code{band}}{The causes' bands, one row per cause per grid time
#'       per stratum: \code{STRATUM}, \code{cmt}, \code{time}, \code{obs} (the
#'       observed cumulative incidence at \code{time}; \code{NA} past the
#'       observed curve's last time), \code{lo}, \code{med}, \code{hi} (the
#'       band at \code{probs}; \code{NA} where too few replicates are still
#'       under follow-up) and \code{n_used}, the number of replicate curves
#'       defined there.}
#'     \item{\code{survival_band}}{The all-cause survival band, laid out as
#'       \code{band} without \code{cmt}.}
#'     \item{\code{strata}}{One row per stratum, as in
#'       \code{\link{ferx_vpc_tte}}.}
#'     \item{\code{cmts}}{The causes' compartments, ascending.}
#'     \item{\code{horizon}, \code{n_sim}, \code{probs}, \code{censoring}}{The
#'       horizon \code{H} the replicates were simulated to, and the settings
#'       used.}
#'     \item{\code{warnings}}{The simulator's diagnostics and the data
#'       reader's, also raised as an R warning.}
#'   }
#'
#' @examples
#' ex <- ferx_example("tte_competing_risks")
#' cv <- ferx_vpc_cif(ex$model, ex$data, n_sim = 50L, seed = 1L)
#' cv
#' head(cv$observed)
#' plot(cv)
#' plot(cv, survival = TRUE)
#'
#' @seealso \code{\link{plot.ferx_vpc_cif}}, \code{\link{ferx_vpc_tte}},
#'   \code{\link{ferx_predict_survival}}
#' @family simulation
#' @export
ferx_vpc_cif <- function(model, data = NULL, n_sim = 200L, seed = 42L,
                         fit = NULL, strata = NULL, by_dose = FALSE,
                         censoring = c("observed", "horizon"), horizon = NULL,
                         probs = c(0.05, 0.5, 0.95), grid = "auto") {
  censoring <- match.arg(censoring)
  args <- .ferx_vpc_args(model, data, fit, n_sim, seed, strata, by_dose,
                         censoring, NULL, horizon, probs, grid)
  # Every TTE endpoint is a cause: the glue takes no `cmt`.
  args$shared$cmt <- NULL
  # The engine points a repeated-events endpoint at its own `vpc_rtte`, which
  # ferx_vpc_tte() calls for such an endpoint.
  res <- .ferx_engine_spelling(
    .ferx_vpc_call(args, fit, ferx_rust_vpc_cif, ferx_rust_vpc_cif_from_fit,
                   list(), "ferx_vpc_cif"),
    "E_VPC_TTE_INPUT", "`vpc_rtte`", "`ferx_vpc_tte()`"
  )
  structure(
    list(
      observed = .ferx_vpc_stack(res$strata, "observed"),
      band = .ferx_vpc_stack(res$strata, "band"),
      survival_band = .ferx_vpc_stack(res$strata, "survival_band"),
      strata = .ferx_vpc_strata(res$strata),
      cmts = res$cmts,
      horizon = res$horizon,
      n_sim = res$n_sim,
      probs = res$probs,
      censoring = res$censoring,
      warnings = as.character(attr(res, "simulation_warnings", exact = TRUE))
    ),
    class = "ferx_vpc_cif"
  )
}

#' @export
print.ferx_vpc_cif <- function(x, ...) {
  pct <- paste0(format(100 * x$probs, trim = TRUE), "%")
  cat("Cumulative-incidence VPC for competing risks (ferx_vpc_cif)\n")
  cat("  causes cmt ", paste(x$cmts, collapse = ", "), ", horizon ",
      format(x$horizon), ", ", x$n_sim, " replicates, ", x$censoring,
      " censoring\n", sep = "")
  cat("  band ", pct[1], "-", pct[3], " with the ", pct[2], " line\n",
      sep = "")
  for (i in seq_len(nrow(x$strata))) {
    s <- as.character(x$strata$STRATUM[i])
    obs <- x$observed[x$observed$STRATUM == s, , drop = FALSE]
    per_cause <- vapply(x$cmts, function(k) {
      sprintf("%d on cmt %d", sum(obs[[paste0("n_event_", k)]]), k)
    }, character(1))
    cat(sprintf("  %s: %d subjects, events %s\n", s, x$strata$n_subjects[i],
                paste(per_cause, collapse = ", ")))
  }
  .ferx_vpc_print_warnings(x$warnings)
  invisible(x)
}

#' Plot a cumulative-incidence VPC
#'
#' Draws each cause's observed cumulative incidence over its replicate band
#' from \code{\link{ferx_vpc_cif}}, one panel per stratum and cause, with base
#' graphics; optionally also each stratum's all-cause Kaplan-Meier curve over
#' its band. The curves are right-continuous steps, and so are the bands: the
#' band at grid time \code{t[k]} holds until \code{t[k + 1]}. A panel with no
#' band at any time says so.
#'
#' @param x A \code{ferx_vpc_cif} object.
#' @param y Unused.
#' @param cmts The causes to draw; \code{NULL} (default) for every cause.
#' @param survival \code{TRUE} to add one panel per stratum with the all-cause
#'   survival and its band.
#' @param ... Passed to \code{graphics::plot()}.
#' @return \code{x}, invisibly.
#' @examples
#' ex <- ferx_example("tte_competing_risks")
#' cv <- ferx_vpc_cif(ex$model, ex$data, n_sim = 50L, seed = 1L)
#' plot(cv, cmts = 2)
#' @export
plot.ferx_vpc_cif <- function(x, y = NULL, cmts = NULL, survival = FALSE,
                              ...) {
  if (!is.null(cmts) &&
      (!is.numeric(cmts) || !length(cmts) || anyNA(cmts) ||
       !all(cmts %in% x$cmts))) {
    stop("`cmts` must be NULL or causes of this VPC (cmt ",
         paste(x$cmts, collapse = ", "), ").", call. = FALSE)
  }
  if (!is.logical(survival) || length(survival) != 1L || is.na(survival)) {
    stop("`survival` must be TRUE or FALSE.", call. = FALSE)
  }
  draw <- if (is.null(cmts)) x$cmts else x$cmts[x$cmts %in% cmts]
  panels <- list()
  for (s in levels(x$strata$STRATUM)) {
    tag <- if (nrow(x$strata) > 1L) paste0(s, ": ") else ""
    obs <- x$observed[x$observed$STRATUM == s, , drop = FALSE]
    for (k in draw) {
      panels[[length(panels) + 1L]] <- list(
        band = x$band[x$band$STRATUM == s & x$band$cmt == k, , drop = FALSE],
        time = obs$time, y = obs[[paste0("cif_", k)]], start = 0,
        ylab = "Cumulative incidence", main = sprintf("%scause cmt %d", tag, k)
      )
    }
    if (survival) {
      panels[[length(panels) + 1L]] <- list(
        band = x$survival_band[x$survival_band$STRATUM == s, , drop = FALSE],
        time = obs$time, y = obs$surv, start = 1,
        ylab = "All-cause survival", main = paste0(tag, "all causes")
      )
    }
  }
  n <- length(panels)
  nc <- ceiling(sqrt(n))
  old <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (p in panels) {
    .ferx_vpc_step_panel(p$band, p$time, p$y, x$horizon, start = p$start,
                         ylab = p$ylab, main = p$main,
                         note = .ferx_vpc_cif_no_band_label, ...)
  }
  invisible(x)
}

# The note a CIF panel with no band carries. Unlike its siblings, core's
# `vpc_cif` raises no W_VPC_EMPTY_STRATUM (ferx-core at 13e6c8bf), so the
# shared note's "see $warnings" would point at nothing (#778 review, finding
# 1); this one says the reason itself. Back to the shared note once core
# raises the warning for `vpc_cif` too.
.ferx_vpc_cif_no_band_label <- paste(
  "No simulated band here: too few replicates under follow-up at these grid",
  "times (raise n_sim, or use grid times inside the horizon)"
)
