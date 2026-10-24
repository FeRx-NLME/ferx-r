#' Plot the distribution of each random-effect parameter
#'
#' One panel per individual parameter that carries an ETA: a histogram of the
#' individual estimates (\code{fit$individual_estimates}) with an EBE rug, the
#' population density of the parameter that the fitted model implies, and,
#' where the parameter is not plainly log-normal, a dashed log-normal
#' reference with the same \eqn{\omega}, so the effect of the scale or shape
#' can be seen. With \code{qq = TRUE} each parameter also gets a Q-Q panel of
#' its EBE \eqn{\eta} against \eqn{N(0, \omega)}.
#'
#' The implied density is a change of variables from \eqn{\eta \sim N(0,
#' \omega)} through the parameter's scale as the engine classified it
#' (\code{fit$eta_param_info}): \eqn{P = TV e^{h(\eta)}} for log-normal,
#' \eqn{P = TV + h(\eta)} for additive and
#' \eqn{P = \mathrm{logit}^{-1}(\mathrm{logit}(TV) + h(\eta))} for the logit
#' scales, with \eqn{h} the identity unless the engine reports a shape
#' transform (\code{boxcox}, \code{tdist}, \code{johndraper}) and its
#' estimated shape. Each subject's typical value is read back from its
#' individual estimate and EBE, so with covariates on the typical value the
#' density is the mixture over subjects. A parameter on the \code{"custom"}
#' scale, or one whose scale is unknown, is drawn as a histogram only.
#'
#' EBEs shrink towards 0, so the histogram and the Q-Q panel understate the
#' spread when \eqn{\eta}-shrinkage is high; each panel quotes it. With
#' \code{conddist = TRUE} (a fit run with the \code{conddist} option) the
#' histogram and the Q-Q panel use draws from each subject's conditional
#' distribution, \eqn{N(\mathrm{COND\_MEAN}, \mathrm{COND\_SD}^2)} as
#' \code{\link{ferx_conddist}} reports it, instead of the point EBEs; those
#' do not shrink.
#'
#' @param fit A \code{ferx_fit} object, as \code{\link{ferx_fit}} returns it
#'   (a fit loaded with \code{\link{ferx_load_fit}} carries no
#'   \code{individual_estimates}).
#' @param param Optional character vector of individual parameter names
#'   (columns of \code{fit$individual_estimates}) to plot. \code{NULL} (the
#'   default) plots every parameter that carries an ETA.
#' @param qq Logical; add a Q-Q panel under each parameter. Default
#'   \code{FALSE}.
#' @param conddist Logical; use conditional-distribution draws in place of
#'   the EBEs. Default \code{FALSE}. Requires \code{fit$cond_dist}.
#' @param n_draws Draws per subject when \code{conddist = TRUE}. Default 20.
#' @param breaks Passed to \code{\link[graphics]{hist}}.
#' @param ... Ignored.
#'
#' @return Invisibly, a named list with one element per plotted parameter:
#'   \code{eta}, \code{param_type}, \code{shape} (\code{NULL} or a list of
#'   \code{kind}, \code{value}), \code{omega}, \code{shrinkage}, the plotted
#'   \code{values}, and \code{density}, a data frame of \code{x},
#'   \code{density} and \code{reference} (\code{NA} where no density or
#'   reference is drawn).
#'
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "foce", covariance = FALSE)
#' plot_param_dist(fit)
#' plot_param_dist(fit, param = "CL", qq = TRUE)
#'
#' @family diagnostics
#' @export
plot_param_dist <- function(fit, param = NULL, qq = FALSE, conddist = FALSE,
                            n_draws = 20L, breaks = "Sturges", ...) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object.", call. = FALSE)
  }
  if (!is.data.frame(fit$individual_estimates)) {
    stop(
      "The fit carries no `individual_estimates` to plot. A fit loaded with ",
      "ferx_load_fit() does not keep them; plot the fit ferx_fit() returned.",
      call. = FALSE
    )
  }
  specs <- .ferx_param_dist_specs(fit)
  if (!is.null(param)) {
    unknown <- setdiff(param, names(specs))
    if (length(unknown) > 0L) {
      stop(sprintf(
        "No ETA-bearing individual parameter named %s; available: %s.",
        paste(sQuote(unknown, FALSE), collapse = ", "),
        if (length(specs)) paste(names(specs), collapse = ", ") else "none"
      ), call. = FALSE)
    }
    specs <- specs[param]
  }
  if (length(specs) == 0L) {
    stop("The fit has no individual parameter with an ETA to plot.", call. = FALSE)
  }
  draws <- if (isTRUE(conddist)) .ferx_param_dist_draws(fit, n_draws) else NULL

  panels <- lapply(specs, .ferx_param_dist_panel, fit = fit, draws = draws)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  n <- length(panels)
  ncol <- min(n, 3L)
  nrow <- ceiling(n / ncol)
  if (isTRUE(qq)) {
    # Each parameter's Q-Q panel sits directly under its histogram.
    idx <- matrix(0L, nrow = 2L * nrow, ncol = ncol)
    for (k in seq_len(n)) {
      r <- (k - 1L) %/% ncol
      cc <- (k - 1L) %% ncol + 1L
      idx[2L * r + 1L, cc] <- 2L * k - 1L
      idx[2L * r + 2L, cc] <- 2L * k
    }
    graphics::layout(idx)
  } else {
    par(mfrow = c(nrow, ncol))
  }
  par(mar = c(4.5, 4, 3, 1))
  for (p in panels) {
    .ferx_param_dist_draw(p, breaks = breaks)
    if (isTRUE(qq)) .ferx_param_dist_qq(p)
  }
  invisible(lapply(panels, function(p) p[c(
    "eta", "param_type", "shape", "omega", "shrinkage", "values", "density"
  )]))
}

# -- Specs: which ETA enters which parameter, on what scale --

# One spec per individual parameter that carries exactly one classified ETA,
# from `fit$eta_param_info` (#478). A fit without it, or whose entries name
# no parameter (a `.fitrx` written before #478), falls back to `eta_param_types`, pairing ETA_X with a column X of
# `individual_estimates` when one exists.
.ferx_param_dist_specs <- function(fit) {
  ie <- fit$individual_estimates
  eta_names <- as.character(fit$eta_names %||% character())
  if (!is.data.frame(ie) || length(eta_names) == 0L) return(list())
  info <- fit$eta_param_info
  if (is.data.frame(info) && nrow(info) > 0L && any(nzchar(info$individual_param))) {
    rows <- data.frame(
      eta = as.character(info$eta_name),
      param = as.character(info$individual_param),
      type = as.character(info$param_type),
      stringsAsFactors = FALSE
    )
    rows$shape_kind <- if (!is.null(info$shape_kind)) as.character(info$shape_kind) else NA_character_
    rows$shape_value <- if (!is.null(info$shape_value)) as.numeric(info$shape_value) else NA_real_
  } else {
    types <- as.character(fit$eta_param_types %||% character())
    if (length(types) != length(eta_names)) types <- rep(NA_character_, length(eta_names))
    rows <- data.frame(
      eta = eta_names, param = sub("^ETA_?", "", eta_names), type = types,
      shape_kind = NA_character_, shape_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[nzchar(rows$param) & rows$param %in% setdiff(names(ie), "ID") &
                 rows$eta %in% eta_names, , drop = FALSE]
  # A parameter two ETAs enter has no one-ETA density.
  multi <- names(which(table(rows$param) > 1L))
  rows <- rows[!rows$param %in% multi, , drop = FALSE]
  specs <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    shape <- if (!is.na(r$shape_kind) && nzchar(r$shape_kind) && is.finite(r$shape_value)) {
      list(kind = r$shape_kind, value = r$shape_value)
    }
    list(param = r$param, eta = r$eta, type = r$type, shape = shape)
  })
  names(specs) <- rows$param
  specs
}

# -- Shape transforms h(eta) and their derivative --
#
# These mirror ferx-core's `shape_g` (src/parser/eta_shape.rs) for drawing
# only; every reported number (CV%, quantiles) comes from the engine.
.ferx_shape_h <- function(eta, kind = NULL, s = NA_real_) {
  if (is.null(kind)) return(eta)
  switch(kind,
    boxcox = if (abs(s) < 1e-8) eta else expm1(s * eta) / s,
    tdist = {
      e2 <- eta^2
      eta * (1 + (e2 + 1) / (4 * s) + (5 * e2^2 + 16 * e2 + 3) / (96 * s^2) +
               (3 * e2^3 + 19 * e2^2 + 17 * e2 - 15) / (384 * s^3))
    },
    johndraper = {
      u <- log1p(abs(eta))
      sign(eta) * (if (abs(s) < 1e-8) u else expm1(s * u) / s)
    },
    stop(sprintf("Unknown shape transform %s.", sQuote(kind, FALSE)), call. = FALSE)
  )
}

.ferx_shape_dh <- function(eta, kind = NULL, s = NA_real_) {
  if (is.null(kind)) return(rep(1, length(eta)))
  switch(kind,
    boxcox = exp(s * eta),
    tdist = {
      e2 <- eta^2
      1 + (3 * e2 + 1) / (4 * s) + (25 * e2^2 + 48 * e2 + 3) / (96 * s^2) +
        (21 * e2^3 + 95 * e2^2 + 51 * e2 - 15) / (384 * s^3)
    },
    johndraper = (abs(eta) + 1)^(s - 1),
    stop(sprintf("Unknown shape transform %s.", sQuote(kind, FALSE)), call. = FALSE)
  )
}

# The link that puts a parameter on the scale its ETA is added on.
.ferx_param_link <- function(type) {
  switch(type,
    log_normal = list(fun = log, inv = exp, deriv = function(p) 1 / p),
    additive = list(fun = identity, inv = identity, deriv = function(p) rep(1, length(p))),
    logit = ,
    logit_probability = list(
      fun = stats::qlogis, inv = stats::plogis,
      deriv = function(p) 1 / (p * (1 - p))
    ),
    NULL
  )
}

# Density of P where link(P) = base + h(eta), eta ~ N(0, omega), averaged
# over the `base` values (one per subject), evaluated at `x`.
.ferx_param_density <- function(x, base, omega, type, shape = NULL) {
  link <- .ferx_param_link(type)
  if (is.null(link) || !is.finite(omega) || omega <= 0) return(rep(NA_real_, length(x)))
  base <- base[is.finite(base)]
  if (length(base) == 0L) return(rep(NA_real_, length(x)))
  sd <- sqrt(omega)
  eta <- seq(-8 * sd, 8 * sd, length.out = 4001L)
  h <- .ferx_shape_h(eta, shape$kind, shape$value %||% NA_real_)
  dh <- .ferx_shape_dh(eta, shape$kind, shape$value %||% NA_real_)
  ok <- is.finite(h) & is.finite(dh) & dh > 0
  eta <- eta[ok]; h <- h[ok]; dh <- dh[ok]
  phi <- stats::dnorm(eta, 0, sd)
  out <- numeric(length(x))
  for (b in unique(base)) {
    w <- sum(base == b)
    p <- link$inv(b + h)
    f <- phi * abs(link$deriv(p)) / dh
    keep <- is.finite(p) & is.finite(f) & !duplicated(p)
    if (sum(keep) < 2L) next
    out <- out + w * stats::approx(p[keep], f[keep], xout = x, rule = 1)$y
  }
  out[is.na(out)] <- 0
  out / length(base)
}

# -- Panel data --

.ferx_param_dist_panel <- function(spec, fit, draws) {
  ie <- fit$individual_estimates
  ebe <- fit$ebe_etas
  i_eta <- match(spec$eta, fit$eta_names)
  omega <- if (is.matrix(fit$omega) && i_eta <= nrow(fit$omega)) fit$omega[i_eta, i_eta] else NA_real_
  shrink <- if (length(fit$shrinkage_eta) >= i_eta) fit$shrinkage_eta[[i_eta]] else NA_real_
  pv <- ie[[spec$param]]
  eta_i <- if (is.data.frame(ebe) && spec$eta %in% names(ebe)) {
    ebe[[spec$eta]][match(as.character(ie$ID), as.character(ebe$ID))]
  } else {
    rep(NA_real_, length(pv))
  }
  type <- if (is.na(spec$type)) "unknown" else spec$type
  link <- .ferx_param_link(type)
  shape_s <- spec$shape$value %||% NA_real_
  # Each subject's typical value on the link scale: link(P_i) - h(eta_i).
  base <- if (!is.null(link)) {
    suppressWarnings(link$fun(pv) - .ferx_shape_h(eta_i, spec$shape$kind, shape_s))
  } else {
    rep(NA_real_, length(pv))
  }

  values <- pv
  eta_vals <- eta_i
  if (!is.null(draws)) {
    d <- draws[draws$ETA == spec$eta, , drop = FALSE]
    b <- base[match(as.character(d$ID), as.character(ie$ID))]
    eta_vals <- d$eta
    values <- if (!is.null(link)) {
      link$inv(b + .ferx_shape_h(d$eta, spec$shape$kind, shape_s))
    } else {
      rep(NA_real_, nrow(d))
    }
  }

  fin <- c(values[is.finite(values)], pv[is.finite(pv)])
  dens <- NULL
  if (!is.null(link) && is.finite(omega) && omega > 0 && any(is.finite(base))) {
    lo <- min(fin); hi <- max(fin)
    tv <- link$inv(stats::median(base, na.rm = TRUE))
    q <- link$inv(stats::median(base, na.rm = TRUE) +
                    .ferx_shape_h(c(-3, 3) * sqrt(omega), spec$shape$kind, shape_s))
    rng <- range(c(lo, hi, q, tv), finite = TRUE)
    pad <- 0.05 * diff(rng)
    if (type %in% c("log_normal", "logit", "logit_probability")) rng[1] <- max(rng[1] - pad, 0)
    else rng[1] <- rng[1] - pad
    rng[2] <- if (type %in% c("logit", "logit_probability")) min(rng[2] + pad, 1) else rng[2] + pad
    x <- seq(rng[1], rng[2], length.out = 512L)
    x <- x[x > 0 | !type %in% c("log_normal", "logit", "logit_probability")]
    if (type %in% c("logit", "logit_probability")) x <- x[x < 1]
    f <- .ferx_param_density(x, base, omega, type, spec$shape)
    ref <- rep(NA_real_, length(x))
    plain_ln <- identical(type, "log_normal") && is.null(spec$shape)
    if (!plain_ln) {
      tv_nat <- link$inv(base)
      if (all(tv_nat[is.finite(tv_nat)] > 0)) {
        ref <- .ferx_param_density(x, log(tv_nat), omega, "log_normal")
      }
    }
    dens <- data.frame(x = x, density = f, reference = ref)
  } else if (length(fin) > 0L) {
    dens <- data.frame(x = range(fin), density = NA_real_, reference = NA_real_)
  }

  list(
    param = spec$param, eta = spec$eta, param_type = type, shape = spec$shape,
    omega = omega, shrinkage = shrink, values = values, ebe_values = pv,
    eta_values = eta_vals, density = dens, conddist = !is.null(draws)
  )
}

# Conditional-distribution draws, N(COND_MEAN, COND_SD^2) per subject x ETA.
.ferx_param_dist_draws <- function(fit, n_draws) {
  cd <- ferx_conddist(fit)
  n_draws <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 1L) {
    stop("`n_draws` must be a positive integer.", call. = FALSE)
  }
  k <- rep(seq_len(nrow(cd)), each = n_draws)
  data.frame(
    ID = cd$ID[k], ETA = as.character(cd$ETA[k]),
    eta = stats::rnorm(length(k), cd$COND_MEAN[k], cd$COND_SD[k]),
    stringsAsFactors = FALSE
  )
}

.ferx_param_dist_tag <- function(p) {
  scale <- switch(p$param_type,
    log_normal = "log-normal", additive = "additive",
    logit = , logit_probability = "logit", custom = "custom scale",
    "scale unknown"
  )
  if (!is.null(p$shape)) {
    sym <- switch(p$shape$kind, boxcox = "lambda", johndraper = "lambda", tdist = "nu", "shape")
    label <- switch(p$shape$kind, boxcox = "box-cox", johndraper = "john-draper", tdist = "t", p$shape$kind)
    scale <- sprintf("%s, %s %s=%s", scale, label, sym, format(signif(p$shape$value, 3)))
  }
  scale
}

.ferx_param_dist_shrink_note <- function(p) {
  if (isTRUE(p$conddist)) return("conditional-distribution draws")
  if (is.finite(p$shrinkage)) {
    sprintf("EBEs; eta-shrinkage %.0f%%%s", 100 * p$shrinkage,
            if (p$shrinkage >= 0.2) " (spread understated)" else "")
  } else {
    "EBEs (shrink towards 0)"
  }
}

.ferx_param_dist_draw <- function(p, breaks) {
  vals <- p$values[is.finite(p$values)]
  d <- p$density
  xlim <- range(c(vals, d$x), finite = TRUE)
  h <- if (length(vals) > 0L) {
    graphics::hist(vals, breaks = breaks, plot = FALSE)
  }
  ymax <- max(c(h$density, d$density, d$reference, 0), na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  plot(NA, xlim = xlim, ylim = c(0, 1.05 * ymax),
       xlab = p$param, ylab = "Density",
       main = sprintf("%s  [%s]", p$param, .ferx_param_dist_tag(p)))
  graphics::mtext(.ferx_param_dist_shrink_note(p), side = 3, line = 0.2, cex = 0.7)
  if (!is.null(h)) {
    graphics::rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1L], h$density,
                   col = "grey80", border = "grey40")
  }
  graphics::rug(p$ebe_values[is.finite(p$ebe_values)])
  leg <- character(); lty <- integer(); col <- character()
  if (any(is.finite(d$density))) {
    graphics::lines(d$x, d$density, lwd = 2, col = "steelblue4")
    leg <- c(leg, "implied population"); lty <- c(lty, 1L); col <- c(col, "steelblue4")
  }
  if (any(is.finite(d$reference))) {
    graphics::lines(d$x, d$reference, lwd = 1.5, lty = 2, col = "firebrick")
    leg <- c(leg, "log-normal, same omega"); lty <- c(lty, 2L); col <- c(col, "firebrick")
  }
  if (length(leg)) {
    graphics::legend("topright", legend = leg, lty = lty, col = col, lwd = 1.5,
                     bty = "n", cex = 0.75)
  }
}

.ferx_param_dist_qq <- function(p) {
  e <- p$eta_values[is.finite(p$eta_values)]
  if (length(e) < 2L || !is.finite(p$omega) || p$omega <= 0) {
    plot.new()
    graphics::mtext(sprintf("%s: no Q-Q (no EBEs or omega)", p$eta), side = 3)
    return(invisible())
  }
  e <- sort(e)
  q <- stats::qnorm(stats::ppoints(length(e)), 0, sqrt(p$omega))
  lim <- range(c(e, q))
  plot(q, e, xlim = lim, ylim = lim, pch = 16, cex = 0.6,
       col = "steelblue4",
       xlab = sprintf("N(0, omega) quantile"), ylab = sprintf("%s quantile", p$eta),
       main = sprintf("%s Q-Q", p$eta))
  graphics::abline(0, 1, lty = 2)
  graphics::mtext(.ferx_param_dist_shrink_note(p), side = 3, line = 0.2, cex = 0.7)
}
