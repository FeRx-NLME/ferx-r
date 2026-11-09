#' Kaplan-Meier visual predictive check for a time-to-event model
#'
#' Simulates the observed design \code{n_sim} times and returns, per stratum,
#' the observed Kaplan-Meier curve beside the pointwise band of the replicate
#' Kaplan-Meier curves, for one single-event time-to-event (TTE) endpoint.
#' Everything is computed by ferx-core (\code{vpc_tte()}, ferx-core #1886);
#' this function only carries the result into R. It never refits.
#'
#' @section The curve and the band:
#' The observed curve uses \code{survival::survfit}'s tie convention (at a
#' tied time, events count before censorings), so \code{x$observed} matches
#' \code{survfit(Surv(time, event) ~ 1)} row for row, with the Nelson-Aalen
#' cumulative hazard in \code{cumhaz}. The band at each grid time is R's
#' default (type 7) quantile over the replicate curves still under follow-up
#' there. A replicate past its last follow-up time is left out rather than
#' counted as 0 or 1, and a grid time with fewer than
#' \code{ceiling(1 / min(p, 1 - p)) - 1} replicates left (19 for a 5\%-95\%
#' band) has an \code{NA} band: its tail quantile would be an interpolation
#' between the two most extreme curves. So \code{n_sim} must be at least one
#' more than that, 20 for the default band.
#'
#' @section Censoring in the replicates:
#' The replicates are simulated to the horizon \code{H} (\code{horizon}, or
#' the latest time recorded on the endpoint). Under the default
#' \code{censoring = "observed"} each subject then gets a censoring time, and
#' a replicate event after it becomes a censoring there: a subject observed
#' censored keeps its own time, and a subject observed with an event draws
#' one from the reverse Kaplan-Meier curve of the observed censorings in its
#' own stratum, given that it exceeds the event time. The replicates are then
#' censored like the data, which is what makes the observed curve comparable
#' with the band. \code{censoring = "horizon"} censors every replicate at
#' \code{H} only, which gives too narrow a band when the data are censored
#' before \code{H}. The draws are seeded by \code{seed}.
#'
#' @section What it refuses:
#' An input the curve or the band would be wrong on is an R error of class
#' \code{ferx_engine_error} with \code{e$code == "E_VPC_TTE_INPUT"}: a model
#' without a TTE endpoint; several TTE endpoints and no \code{cmt} (under
#' competing risks one minus the Kaplan-Meier curve of one cause is not its
#' cumulative incidence; with \code{cmt} set, the curve is that cause's
#' survival with the other causes as censoring); a \code{cmt} that is not a TTE
#' endpoint; a repeated-events endpoint; left truncation (delayed entry); an
#' interval-censored record; more than one record per subject on the endpoint;
#' no record on it at all; band probabilities not strictly increasing inside
#' (0, 1), or too few replicates for them; a grid time that is negative; and a
#' stratification covariate a subject lacks, or with more than 10 values (bin
#' it into a categorical covariate first).
#'
#' @param model Path to a .ferx model file with an \code{[event_model]}
#'   endpoint.
#' @param data Path to a NONMEM-format CSV, or a data.frame, holding the
#'   \strong{observed} TTE records (\code{DV = 1} event, \code{DV = 0}
#'   right-censored), one per subject on the endpoint. When omitted, the data
#'   \code{fit} was fitted to, otherwise the model file's \code{[data]} block.
#' @param n_sim Number of simulated replicates; at least 20 for the default
#'   \code{probs}.
#' @param seed Random seed for the replicates and the censoring draws.
#' @param fit Optional \code{ferx_fit} result. When given, the replicates are
#'   simulated at the fitted parameters, as in \code{\link{ferx_simulate}};
#'   otherwise at the model file's initial values.
#' @param strata \code{NULL} (one stratum, \code{"all"}) or the name of a
#'   covariate: subjects are grouped by its exact value at their first record.
#' @param by_dose \code{TRUE} to stratify by the amount of each subject's
#'   first dose instead (0 for a subject with no dose). Cannot be combined with
#'   \code{strata}.
#' @param censoring \code{"observed"} (default) or \code{"horizon"}; see the
#'   section above.
#' @param cmt The TTE endpoint's compartment. Required when the model has more
#'   than one TTE endpoint; \code{NULL} takes the only one.
#' @param horizon The time the replicates are simulated to. \code{NULL}
#'   (default) is the latest time recorded on the endpoint.
#' @param probs Three band probabilities, strictly increasing inside (0, 1):
#'   the lower bound, the middle line and the upper bound.
#' @param grid The times the band is evaluated at: \code{"auto"} (0, the
#'   stratum's observed times and 101 evenly spaced points up to \code{H}),
#'   \code{"all_steps"} (also every step of every replicate curve: exact step
#'   ribbons, on a large grid), or a numeric vector of times.
#'
#' @return An object of class \code{ferx_vpc_tte}, a list with:
#'   \describe{
#'     \item{\code{observed}}{The observed Kaplan-Meier curve, one row per
#'       distinct time per stratum: \code{STRATUM}, \code{time},
#'       \code{n_risk}, \code{n_event}, \code{n_censor}, \code{surv} and the
#'       Nelson-Aalen \code{cumhaz}.}
#'     \item{\code{band}}{The replicate band, one row per grid time per
#'       stratum: \code{STRATUM}, \code{time}, \code{obs} (the observed curve
#'       read at \code{time}; \code{NA} past its last time unless it reached
#'       0), \code{lo}, \code{med}, \code{hi} (the band at \code{probs};
#'       \code{NA} where too few replicates are still under follow-up) and
#'       \code{n_used}, the number of replicate curves defined there.}
#'     \item{\code{strata}}{One row per stratum: \code{STRATUM} (its label:
#'       \code{"all"}, \code{"NAME=value"} or \code{"dose=value"}),
#'       \code{value} (the covariate or dose value; \code{NA} without strata)
#'       and \code{n_subjects}.}
#'     \item{\code{cmt}, \code{horizon}, \code{n_sim}, \code{probs},
#'       \code{censoring}}{The endpoint checked, the horizon \code{H} the
#'       replicates were simulated to, and the settings used.}
#'     \item{\code{warnings}}{The simulator's diagnostics and the data
#'       reader's, also raised as an R warning.}
#'   }
#'   \code{STRATUM} is a factor whose levels are in the engine's stratum order
#'   (ascending value).
#'
#' @examples
#' ex <- ferx_example("tte_weibull")
#' vpc <- ferx_vpc_tte(ex$model, ex$data, n_sim = 50L, seed = 1L)
#' vpc
#' plot(vpc)
#'
#' @seealso \code{\link{plot.ferx_vpc_tte}}, \code{\link{ferx_predict_survival}}
#' @family simulation
#' @export
ferx_vpc_tte <- function(model, data = NULL, n_sim = 200L, seed = 42L,
                         fit = NULL, strata = NULL, by_dose = FALSE,
                         censoring = c("observed", "horizon"), cmt = NULL,
                         horizon = NULL, probs = c(0.05, 0.5, 0.95),
                         grid = "auto") {
  censoring <- match.arg(censoring)
  if (is.null(data)) data <- .ferx_default_data(model, fit)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  data <- .ferx_data_arg(data)
  stopifnot(file.exists(model), file.exists(data))

  .ferx_vpc_tte_count(n_sim, "n_sim", min = 1)
  .ferx_vpc_tte_count(seed, "seed", min = 0)
  if (!is.logical(by_dose) || length(by_dose) != 1L || is.na(by_dose)) {
    stop("`by_dose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(strata) &&
      (!is.character(strata) || length(strata) != 1L || is.na(strata) ||
       !nzchar(strata))) {
    stop("`strata` must be NULL or a single covariate name.", call. = FALSE)
  }
  if (!is.null(strata) && by_dose) {
    stop("Pass either `strata` (a covariate) or `by_dose = TRUE`, not both.",
         call. = FALSE)
  }
  stratify <- if (by_dose) "dose" else if (is.null(strata)) "none" else "covariate"
  if (!is.null(cmt)) .ferx_vpc_tte_count(cmt, "cmt", min = 1)
  if (!is.null(horizon) &&
      (length(horizon) != 1L || !is.numeric(horizon) || !is.finite(horizon) ||
       horizon <= 0)) {
    stop("`horizon` must be a single finite positive number (or NULL).",
         call. = FALSE)
  }
  # The engine checks the probabilities' values (E_VPC_TTE_INPUT); only the
  # shape is the glue's.
  if (!is.numeric(probs) || length(probs) != 3L) {
    stop("`probs` must be a numeric vector of three probabilities.",
         call. = FALSE)
  }
  if (is.numeric(grid)) {
    if (length(grid) == 0L || anyNA(grid)) {
      stop("A numeric `grid` must hold at least one time and no NA.",
           call. = FALSE)
    }
    grid_mode <- "times"
    grid_times <- as.numeric(grid)
  } else if (is.character(grid) && length(grid) == 1L &&
             grid %in% c("auto", "all_steps")) {
    grid_mode <- grid
    grid_times <- numeric(0)
  } else {
    stop("`grid` must be \"auto\", \"all_steps\" or a numeric vector of times.",
         call. = FALSE)
  }

  shared <- list(
    n_sim = as.integer(n_sim),
    seed = as.integer(seed),
    probs = as.numeric(probs),
    grid = grid_mode,
    grid_times = grid_times,
    censoring = censoring,
    stratify = stratify,
    strata_covariate = if (is.null(strata)) "" else strata,
    cmt = if (is.null(cmt)) 0L else as.integer(cmt),
    # -1 is the "unset" sentinel the glue reads as the latest recorded time.
    horizon = if (is.null(horizon)) -1 else as.numeric(horizon)
  )
  paths <- list(model_path = normalizePath(model),
                data_path = normalizePath(data))
  # A refusal is an R error, classed like a refused `ferx_fit()` (#385).
  res <- if (is.null(fit)) {
    .ferx_engine_call(do.call(ferx_rust_vpc_tte, c(paths, shared)), model, data)
  } else {
    fit_pieces <- validate_fit_for_params(fit)
    .ferx_engine_call(
      do.call(ferx_rust_vpc_tte_from_fit, c(paths, fit_pieces[c(
        "theta", "omega_flat", "omega_dim", "sigma", "omega_iov_flat",
        "omega_iov_dim", "residual_rho", "residual_rho_pairs",
        "mixture_overrides", "fit_bindings", "iov_occasion"
      )], shared)),
      model, data
    )
  }
  res <- .ferx_surface_sim_warnings(res, "ferx_vpc_tte")
  .ferx_vpc_tte_object(res)
}

# A single whole number >= `min`, for the integer arguments the glue takes.
.ferx_vpc_tte_count <- function(x, arg, min) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) ||
      x < min || x > .Machine$integer.max) {
    stop("`", arg, "` must be a single whole number >= ", min, ".",
         call. = FALSE)
  }
}

# The glue's per-stratum lists as long tables with a STRATUM factor. Only a
# reshape: every number is core's.
.ferx_vpc_tte_object <- function(res) {
  strata <- res$strata
  labels <- vapply(strata, function(s) s$label, character(1))
  tag <- function(field) {
    out <- do.call(rbind, lapply(strata, function(s) {
      d <- s[[field]]
      cbind(data.frame(STRATUM = rep(s$label, nrow(d)),
                       stringsAsFactors = FALSE), d)
    }))
    out$STRATUM <- factor(out$STRATUM, levels = labels)
    rownames(out) <- NULL
    out
  }
  structure(
    list(
      observed = tag("observed"),
      band = tag("band"),
      strata = data.frame(
        STRATUM = factor(labels, levels = labels),
        value = vapply(strata, function(s) as.numeric(s$value), numeric(1)),
        n_subjects = vapply(strata, function(s) s$n_subjects, integer(1))
      ),
      cmt = res$cmt,
      horizon = res$horizon,
      n_sim = res$n_sim,
      probs = res$probs,
      censoring = res$censoring,
      warnings = as.character(attr(res, "simulation_warnings", exact = TRUE))
    ),
    class = "ferx_vpc_tte"
  )
}

#' @export
print.ferx_vpc_tte <- function(x, ...) {
  pct <- paste0(format(100 * x$probs, trim = TRUE), "%")
  cat("Kaplan-Meier VPC (ferx_vpc_tte)\n")
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
  if (length(x$warnings)) {
    cat("  ", length(x$warnings), " warning",
        if (length(x$warnings) > 1L) "s" else "",
        " (see `$warnings`)\n", sep = "")
  }
  invisible(x)
}

#' Plot a Kaplan-Meier VPC
#'
#' Draws each stratum's observed Kaplan-Meier curve over the replicate band
#' from \code{\link{ferx_vpc_tte}}. With the \pkg{vpc} package installed the
#' plot is \code{vpc::plot_vpc()}'s time-to-event plot (a ggplot); otherwise
#' it is drawn with base graphics, one panel per stratum. Either way the band
#' is a step ribbon: the band at grid time \code{t[k]} holds until
#' \code{t[k + 1]}.
#'
#' @param x A \code{ferx_vpc_tte} object.
#' @param y Unused.
#' @param engine \code{"auto"} (\pkg{vpc} when installed, else base
#'   graphics), \code{"vpc"} or \code{"base"}.
#' @param ... Passed to \code{vpc::plot_vpc()} (for example \code{xlab},
#'   \code{title}, \code{smooth}) or to \code{graphics::plot()}.
#' @return With \pkg{vpc}, the ggplot object; otherwise \code{x}, invisibly.
#' @export
plot.ferx_vpc_tte <- function(x, y = NULL, engine = c("auto", "vpc", "base"),
                              ...) {
  engine <- match.arg(engine)
  has_vpc <- requireNamespace("vpc", quietly = TRUE)
  if (engine == "vpc" && !has_vpc) {
    stop("engine = \"vpc\" needs the vpc package.", call. = FALSE)
  }
  if (engine == "vpc" || (engine == "auto" && has_vpc)) {
    return(.ferx_vpc_tte_plot_vpc(x, ...))
  }
  .ferx_vpc_tte_plot_base(x, ...)
}

# The band and the observed curve in the vpcdb shape vpc::vpc_tte() builds
# (vpc 1.2.x): `sim_km` has one row per grid interval [bin_min, bin_max) with
# qmin / qmed / qmax, `obs_km` the observed steps.
.ferx_vpc_tte_vpcdb <- function(x) {
  sim_km <- do.call(rbind, lapply(levels(x$band$STRATUM), function(s) {
    b <- x$band[x$band$STRATUM == s, , drop = FALSE]
    b <- b[order(b$time), , drop = FALSE]
    nxt <- c(b$time[-1L], x$horizon)
    keep <- !is.na(b$lo) & !is.na(b$med) & !is.na(b$hi)
    data.frame(strat = rep(s, sum(keep)), bin = which(keep),
               bin_min = b$time[keep], bin_max = nxt[keep],
               bin_mid = b$time[keep],
               qmin = b$lo[keep], qmed = b$med[keep], qmax = b$hi[keep],
               stringsAsFactors = FALSE)
  }))
  # The observed curve starts at S(0) = 1, as survfit's plot does.
  obs_km <- do.call(rbind, lapply(levels(x$observed$STRATUM), function(s) {
    o <- x$observed[x$observed$STRATUM == s, , drop = FALSE]
    data.frame(strat = s, time = c(0, o$time), surv = c(1, o$surv),
               stringsAsFactors = FALSE)
  }))
  cens <- x$observed[x$observed$n_censor > 0L, , drop = FALSE]
  cens_dat <- data.frame(strat = as.character(cens$STRATUM), time = cens$time,
                         y = cens$surv, stringsAsFactors = FALSE)
  stratified <- nrow(x$strata) > 1L
  levs <- levels(x$strata$STRATUM)
  sim_km$strat <- factor(sim_km$strat, levels = levs)
  obs_km$strat <- factor(obs_km$strat, levels = levs)
  cens_dat$strat <- factor(cens_dat$strat, levels = levs)
  if (stratified) cens_dat$STRATUM <- cens_dat$strat
  vpc::as_vpcdb(
    sim_km = sim_km, obs = obs_km, obs_km = obs_km, cens_dat = cens_dat,
    stratify = if (stratified) "STRATUM", stratify_pars = if (stratified) "STRATUM",
    rtte = FALSE, type = "time-to-event", facet = "wrap", scales = "fixed",
    # vpc's own TTE labeller reads a numeric stratum from the label and shows
    # NA for ours ("ARM=1"); the label is already what the panel should say.
    labeller = ggplot2::labeller(STRATUM = ggplot2::label_value),
    xlab = "Time", ylab = "Survival",
    show = list(sim_median = TRUE, bin_sep = FALSE, obs_ci = FALSE,
                sim_km = FALSE)
  )
}

.ferx_vpc_tte_plot_vpc <- function(x, smooth = FALSE, ...) {
  vpc::plot_vpc(.ferx_vpc_tte_vpcdb(x), smooth = smooth, ...)
}

.ferx_vpc_tte_plot_base <- function(x, ...) {
  panels <- levels(x$strata$STRATUM)
  n <- length(panels)
  nc <- ceiling(sqrt(n))
  old <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (s in panels) {
    b <- x$band[x$band$STRATUM == s, , drop = FALSE]
    b <- b[order(b$time), , drop = FALSE]
    o <- x$observed[x$observed$STRATUM == s, , drop = FALSE]
    graphics::plot(c(0, x$horizon), c(0, 1), type = "n", xlab = "Time",
                   ylab = "Survival", main = if (n > 1L) s else "", ...)
    keep <- !is.na(b$lo) & !is.na(b$hi)
    nxt <- c(b$time[-1L], x$horizon)
    # Opaque, so neighbouring rectangles leave no seam where they meet.
    graphics::rect(b$time[keep], b$lo[keep], nxt[keep], b$hi[keep],
                   col = "#C8D9EA", border = "#C8D9EA")
    graphics::lines(b$time[keep], b$med[keep], type = "s", lty = 2,
                    col = "steelblue", lwd = 2)
    graphics::lines(c(0, o$time), c(1, o$surv), type = "s", lwd = 2)
    cens <- o$n_censor > 0L
    graphics::points(o$time[cens], o$surv[cens], pch = "|")
  }
  invisible(x)
}
