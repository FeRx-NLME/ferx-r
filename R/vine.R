#' Evaluate the fitted vine-copula density of the random effects
#'
#' For a vine SAEM fit (\code{method = "saem"} with \code{omega_dist = "vine"}),
#' evaluates the fitted joint density of the between-subject random effects (the
#' copula-distribution of IIV) over a grid, and returns it as a tidy data frame
#' suitable for contour / heat-map plotting. The fitted vine is reconstructed
#' from \code{fit$vine_copula}, so this works on a freshly fitted object and on
#' one reloaded from a \code{.fitrx} bundle.
#'
#' With two random effects the returned grid is the full bivariate density. With
#' more than two, the chosen pair is varied while the remaining random effects
#' are held at their marginal means (a density slice); pick the pair with
#' \code{etas}.
#'
#' @param fit A \code{ferx_fit} from a vine SAEM run (\code{fit$vine_copula} not
#'   \code{NULL}).
#' @param etas Optional character vector naming the random effect(s) to vary:
#'   two names for a bivariate grid, one for a 1-D marginal slice. Defaults to
#'   the first two etas (or the single eta when the model has only one).
#' @param n Number of grid points per varied dimension (default 50).
#' @param range_sd Half-width of each grid axis in marginal standard deviations
#'   around the marginal mean (default 3.5).
#'
#' @return A data frame. For a bivariate grid: one column per varied eta (named
#'   after the etas) plus \code{density}, with \code{n * n} rows. For a 1-D
#'   slice: the eta column plus \code{density}, with \code{n} rows.
#'
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin_vine")
#' fit <- ferx_fit(ex$model, ex$data, method = "saem",
#'                 settings = list(omega_dist = "vine"))
#' d <- ferx_vine_density(fit)
#' # contour: with ggplot2 -> geom_contour(aes(ETA_CL, ETA_V, z = density))
#' }
#'
#' @family simulation
#' @export
ferx_vine_density <- function(fit, etas = NULL, n = 50L, range_sd = 3.5) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object.", call. = FALSE)
  }
  vc <- fit$vine_copula
  if (is.null(vc)) {
    stop(
      "This fit has no vine-copula structure. Fit with method = \"saem\" and ",
      "omega_dist = \"vine\" to obtain one.",
      call. = FALSE
    )
  }
  if (!is.numeric(n) || length(n) != 1L || n < 2L) {
    stop("`n` must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(range_sd) || length(range_sd) != 1L || range_sd <= 0) {
    stop("`range_sd` must be a single positive number.", call. = FALSE)
  }

  eta_names <- as.character(vc$marginals$eta)
  means     <- as.numeric(vc$marginals$mean)
  sds       <- as.numeric(vc$marginals$sd)
  d         <- length(eta_names)

  # Resolve which dimension(s) to vary.
  if (is.null(etas)) {
    etas <- if (d >= 2L) eta_names[1:2] else eta_names[1]
  }
  if (!is.character(etas) || length(etas) < 1L || length(etas) > 2L) {
    stop("`etas` must name one or two random effects.", call. = FALSE)
  }
  idx <- match(etas, eta_names)
  if (anyNA(idx)) {
    stop(
      "Unknown eta(s): ", paste(etas[is.na(idx)], collapse = ", "),
      ". Available: ", paste(eta_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  n <- as.integer(n)
  axis <- function(k) {
    seq(means[k] - range_sd * sds[k], means[k] + range_sd * sds[k], length.out = n)
  }
  grid_x <- axis(idx[1])
  two_d  <- length(idx) == 2L
  grid_y <- if (two_d) axis(idx[2]) else numeric(0)

  flat <- .ferx_vine_flat(vc)
  res <- ferx_rust_vine_density(
    marg_mean   = flat$marg_mean,
    marg_sd     = flat$marg_sd,
    pair_tree   = flat$pair_tree,
    pair_label  = flat$pair_label,
    pair_family = flat$pair_family,
    param_tree  = flat$param_tree,
    param_label = flat$param_label,
    param_name  = flat$param_name,
    param_value = flat$param_value,
    eta_base    = means,                 # non-varied dims held at marginal mean
    idx_x       = as.integer(idx[1] - 1L),
    idx_y       = if (two_d) as.integer(idx[2] - 1L) else -1L,
    grid_x      = as.numeric(grid_x),
    grid_y      = as.numeric(grid_y)
  )

  if (is.null(res)) {
    stop("Vine density evaluation returned nothing (malformed vine summary).",
         call. = FALSE)
  }

  if (two_d) {
    out <- data.frame(
      x = as.numeric(res$x),
      y = as.numeric(res$y),
      density = as.numeric(res$density),
      stringsAsFactors = FALSE
    )
    names(out) <- c(eta_names[idx[1]], eta_names[idx[2]], "density")
  } else {
    out <- data.frame(
      x = as.numeric(res$x),
      density = as.numeric(res$density),
      stringsAsFactors = FALSE
    )
    names(out) <- c(eta_names[idx[1]], "density")
  }
  out
}

# Flatten a fit$vine_copula list (marginals / pairs / params data frames) into
# the parallel vectors the Rust glue consumes for vine density + simulation.
# Returns empty vectors when `vc` is NULL (Gaussian fit), so callers can pass
# the result straight through to the FFI unconditionally.
.ferx_vine_flat <- function(vc) {
  if (is.null(vc)) {
    return(list(
      marg_mean = numeric(), marg_sd = numeric(),
      pair_tree = integer(), pair_label = character(), pair_family = character(),
      param_tree = integer(), param_label = character(), param_name = character(),
      param_value = numeric()
    ))
  }
  list(
    marg_mean   = as.numeric(vc$marginals$mean),
    marg_sd     = as.numeric(vc$marginals$sd),
    pair_tree   = as.integer(vc$pairs$tree),
    pair_label  = as.character(vc$pairs$pair),
    pair_family = as.character(vc$pairs$family),
    param_tree  = as.integer(vc$params$tree),
    param_label = as.character(vc$params$pair),
    param_name  = as.character(vc$params$parameter),
    param_value = as.numeric(vc$params$estimate)
  )
}
