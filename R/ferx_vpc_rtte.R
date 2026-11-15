# The repeated-events (RTTE) VPC: what ferx_vpc_tte() returns on an endpoint
# declared `type = rtte` (ferx-core #1886 PR2, `vpc_rtte`). The glue reads every
# number off core's strata; this file only stacks the strata, prints and plots.

# The glue's per-stratum lists as long tables with a STRATUM factor.
.ferx_vpc_rtte_object <- function(res) {
  strata <- res$strata
  orders <- .ferx_vpc_stack(strata, "orders")
  structure(
    list(
      mcf = .ferx_vpc_stack(strata, "mcf"),
      mcf_band = .ferx_vpc_stack(strata, "mcf_band"),
      absolute = .ferx_vpc_stack(strata, "absolute"),
      absolute_band = .ferx_vpc_stack(strata, "absolute_band"),
      gap = .ferx_vpc_stack(strata, "gap"),
      gap_band = .ferx_vpc_stack(strata, "gap_band"),
      orders = orders,
      strata = .ferx_vpc_strata(strata),
      cmt = res$cmt,
      horizon = res$horizon,
      n_sim = res$n_sim,
      probs = res$probs,
      censoring = res$censoring,
      max_order = res$max_order,
      warnings = as.character(attr(res, "simulation_warnings", exact = TRUE))
    ),
    class = "ferx_vpc_rtte"
  )
}

#' @export
print.ferx_vpc_rtte <- function(x, ...) {
  pct <- paste0(format(100 * x$probs, trim = TRUE), "%")
  cat("Repeated-events VPC (ferx_vpc_tte, RTTE endpoint)\n")
  cat("  endpoint cmt ", x$cmt, ", horizon ", format(x$horizon), ", ",
      x$n_sim, " replicates, ", x$censoring, " censoring\n", sep = "")
  cat("  band ", pct[1], "-", pct[3], " with the ", pct[2],
      " line; event orders 1-", x$max_order, "\n", sep = "")
  for (i in seq_len(nrow(x$strata))) {
    s <- as.character(x$strata$STRATUM[i])
    m <- x$mcf[x$mcf$STRATUM == s, , drop = FALSE]
    cat(sprintf("  %s: %d subjects, %d events\n", s, x$strata$n_subjects[i],
                sum(m$n_event)))
  }
  .ferx_vpc_print_warnings(x$warnings)
  invisible(x)
}

#' Plot a repeated-events VPC
#'
#' Draws the curves of a repeated-events (RTTE) VPC from
#' \code{\link{ferx_vpc_tte}}, each observed curve over its replicate band: the
#' mean cumulative function, and the Kaplan-Meier curves by event order on
#' absolute and on gap time. One panel per stratum and curve, with base
#' graphics. The band is a step ribbon: the band at grid time \code{t[k]} holds
#' until \code{t[k + 1]}.
#'
#' @param x A \code{ferx_vpc_rtte} object.
#' @param y Unused.
#' @param which Which families to draw: any of \code{"mcf"} (mean cumulative
#'   function), \code{"absolute"} (time to the \eqn{k}-th event) and
#'   \code{"gap"} (time from the \eqn{(k-1)}-th event to the \eqn{k}-th).
#' @param orders The event orders to draw for \code{"absolute"} and
#'   \code{"gap"}; \code{NULL} (default) for every order returned.
#' @param ... Passed to \code{graphics::plot()}.
#' @return \code{x}, invisibly.
#' @examples
#' \donttest{
#' # See ferx_vpc_tte() for building an RTTE VPC; then
#' # plot(rv, which = "mcf")
#' # plot(rv, which = "gap", orders = 2:3)
#' }
#' @export
plot.ferx_vpc_rtte <- function(x, y = NULL, which = c("mcf", "absolute", "gap"),
                               orders = NULL, ...) {
  which <- match.arg(which, several.ok = TRUE)
  if (!is.null(orders) &&
      (!is.numeric(orders) || anyNA(orders) || any(orders < 1))) {
    stop("`orders` must be NULL or positive event orders.", call. = FALSE)
  }
  keep <- function(k) is.null(orders) || k %in% orders
  panels <- list()
  for (s in levels(x$strata$STRATUM)) {
    tag <- if (nrow(x$strata) > 1L) paste0(s, ": ") else ""
    if ("mcf" %in% which) {
      panels[[length(panels) + 1L]] <- list(
        band = x$mcf_band[x$mcf_band$STRATUM == s, , drop = FALSE],
        obs = x$mcf[x$mcf$STRATUM == s, , drop = FALSE],
        y = "mcf", start = 0, ylab = "Mean cumulative function",
        main = paste0(tag, "MCF")
      )
    }
    for (clock in intersect(c("absolute", "gap"), which)) {
      band <- x[[paste0(clock, "_band")]]
      obs <- x[[clock]]
      ks <- x$orders$order[x$orders$STRATUM == s & x$orders$clock == clock]
      for (k in ks[vapply(ks, keep, logical(1))]) {
        panels[[length(panels) + 1L]] <- list(
          band = band[band$STRATUM == s & band$order == k, , drop = FALSE],
          obs = obs[obs$STRATUM == s & obs$order == k, , drop = FALSE],
          y = "surv", start = 1,
          ylab = if (clock == "absolute") "P(fewer than k events)" else "Gap survival",
          main = sprintf("%s%s time, event %d", tag,
                         if (clock == "absolute") "Absolute" else "Gap", k)
        )
      }
    }
  }
  if (!length(panels)) {
    stop("Nothing to draw for `which` and `orders`.", call. = FALSE)
  }
  n <- length(panels)
  nc <- ceiling(sqrt(n))
  old <- graphics::par(mfrow = c(ceiling(n / nc), nc), mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (p in panels) {
    .ferx_vpc_step_panel(p$band, p$obs$time, p$obs[[p$y]], x$horizon,
                         start = p$start, ylab = p$ylab, main = p$main, ...)
  }
  invisible(x)
}

# One VPC panel with base graphics: the band (`time`, `lo`, `med`, `hi`) as an
# opaque step ribbon over [0, horizon), its middle line dashed, and the
# observed step curve from (0, `start`): `step = "s"` for a right-continuous
# curve (Kaplan-Meier, MCF), `"S"` for a left-continuous one (KMMC). The y
# range covers both.
.ferx_vpc_step_panel <- function(band, obs_time, obs_y, horizon, start, ylab,
                                 main, step = "s", ...) {
  b <- .ferx_vpc_tte_bins(band, horizon)
  ys <- c(start, obs_y, b$lo, b$hi)
  ys <- ys[is.finite(ys)]
  ylim <- if (length(ys)) range(ys) else c(0, 1)
  graphics::plot(c(0, horizon), ylim, type = "n", xlab = "Time", ylab = ylab,
                 main = main, ...)
  if (nrow(b)) {
    graphics::rect(b$time, b$lo, b$bin_max, b$hi,
                   col = "#C8D9EA", border = "#C8D9EA")
    graphics::lines(c(b$time, utils::tail(b$bin_max, 1L)),
                    c(b$med, utils::tail(b$med, 1L)), type = "s", lty = 2,
                    col = "steelblue", lwd = 2)
  }
  if (length(obs_time)) {
    graphics::lines(c(0, obs_time), c(start, obs_y), type = step, lwd = 2)
  }
  invisible()
}
