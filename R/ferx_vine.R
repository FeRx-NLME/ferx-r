# Regular-vine random-effect distributions (`omega_dist = "vine"`, ferx-core #225).
#
# `.ferx_vine_from_raw()` normalizes the FFI `vine_copula` list into tidy R
# structures shared identically by `ferx_fit()` and `ferx_load_fit()`.
# `ferx_vine_density()` evaluates exact 1-D marginal, 2-D slice/joint grids, or
# user-supplied points via `ferx_rust_vine_density()`.

.ferx_vine_from_raw <- function(raw) {
  if (is.null(raw) || !length(raw)) return(NULL)

  as_named_num <- function(x) {
    if (is.null(x) || !length(x)) return(stats::setNames(numeric(0), character(0)))
    nms <- names(x)
    v <- as.numeric(unlist(x, use.names = FALSE))
    if (!is.null(nms) && length(nms) == length(v)) names(v) <- nms
    v
  }

  # 1. Tidy `marginals` data.frame (with list-columns `params` and `se` for
  #    multi-parameter zero-mean families: normal, student_t, skew_normal, skew_t).
  m <- raw$marginals %||% list()
  m_eta <- as.character(unlist(m$eta, use.names = FALSE) %||% character(0))
  n_m <- length(m_eta)
  m_family <- as.character(unlist(m$family, use.names = FALSE) %||% rep(NA_character_, n_m))
  m_mean <- as.numeric(unlist(m$mean, use.names = FALSE) %||% rep(0, n_m))
  m_sd <- as.numeric(unlist(m$sd, use.names = FALSE) %||% rep(NA_real_, n_m))
  m_skew <- if (is.null(m$skewness)) rep(NA_real_, n_m) else as.numeric(.fitrx_unwrap_nullable_vec(m$skewness, NA_real_))
  m_kurt <- if (is.null(m$excess_kurtosis)) rep(NA_real_, n_m) else as.numeric(.fitrx_unwrap_nullable_vec(m$excess_kurtosis, NA_real_))
  m_params <- lapply(m$params %||% vector("list", n_m), as_named_num)
  m_se <- lapply(m$se %||% vector("list", n_m), as_named_num)

  marginals_df <- data.frame(
    eta = m_eta,
    family = m_family,
    mean = m_mean,
    sd = m_sd,
    skewness = m_skew,
    excess_kurtosis = m_kurt,
    stringsAsFactors = FALSE
  )
  marginals_df$params <- I(m_params)
  marginals_df$se <- I(m_se)

  # 2. Nested `trees` list + flat `pairs` data.frame + flat `params` data.frame.
  raw_trees <- raw$trees %||% list()
  trees_out <- vector("list", length(raw_trees))

  pair_rows <- list()
  param_rows <- list()

  for (i in seq_len(n_m)) {
    pv <- m_params[[i]]
    sv <- m_se[[i]]
    if (length(pv) > 0L) {
      pnms <- names(pv) %||% paste0("param", seq_along(pv))
      for (k in seq_along(pv)) {
        param_rows[[length(param_rows) + 1L]] <- data.frame(
          component = "marginal",
          target = m_eta[[i]],
          tree = 0L,
          family = m_family[[i]],
          rotation = 0L,
          parameter = pnms[[k]],
          estimate = unname(pv[[k]]),
          se = if (length(sv) >= k) unname(sv[[k]]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  for (ti in seq_along(raw_trees)) {
    tr <- raw_trees[[ti]]
    tree_idx <- as.integer(unlist(tr$tree, use.names = FALSE) %||% ti)[[1L]]
    raw_pairs <- tr$pairs %||% list()
    norm_pairs <- vector("list", length(raw_pairs))

    for (pi in seq_along(raw_pairs)) {
      p <- raw_pairs[[pi]]
      lbl <- as.character(unlist(p$label, use.names = FALSE) %||% "")[[1L]]
      fam <- as.character(unlist(p$family, use.names = FALSE) %||% "independence")[[1L]]
      rot <- as.integer(unlist(p$rotation, use.names = FALSE) %||% 0L)[[1L]]
      pv <- as_named_num(p$params)
      sv <- as_named_num(p$se)
      tau <- as.numeric(unlist(p$kendall_tau, use.names = FALSE) %||% 0)[[1L]]
      td_l <- as.numeric(unlist(p$tail_dep_lower, use.names = FALSE) %||% 0)[[1L]]
      td_u <- as.numeric(unlist(p$tail_dep_upper, use.names = FALSE) %||% 0)[[1L]]

      norm_pairs[[pi]] <- list(
        label = lbl,
        family = fam,
        rotation = rot,
        params = pv,
        se = sv,
        kendall_tau = tau,
        tail_dep_lower = td_l,
        tail_dep_upper = td_u
      )

      p1 <- if (length(pv) >= 1L) unname(pv[[1L]]) else NA_real_
      s1 <- if (length(sv) >= 1L) unname(sv[[1L]]) else NA_real_
      p2 <- if (length(pv) >= 2L) unname(pv[[2L]]) else NA_real_
      s2 <- if (length(sv) >= 2L) unname(sv[[2L]]) else NA_real_
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        tree = tree_idx,
        label = lbl,
        family = fam,
        rotation = rot,
        param1 = p1,
        se1 = s1,
        param2 = p2,
        se2 = s2,
        kendall_tau = tau,
        tail_dep_lower = td_l,
        tail_dep_upper = td_u,
        stringsAsFactors = FALSE
      )

      if (length(pv) > 0L) {
        pnms <- names(pv) %||% paste0("param", seq_along(pv))
        for (k in seq_along(pv)) {
          param_rows[[length(param_rows) + 1L]] <- data.frame(
            component = "copula",
            target = lbl,
            tree = tree_idx,
            family = fam,
            rotation = rot,
            parameter = pnms[[k]],
            estimate = unname(pv[[k]]),
            se = if (length(sv) >= k) unname(sv[[k]]) else NA_real_,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    trees_out[[ti]] <- list(
      tree = tree_idx,
      pairs = norm_pairs
    )
  }

  pairs_df <- if (length(pair_rows) > 0L) {
    do.call(rbind, pair_rows)
  } else {
    data.frame(
      tree = integer(0), label = character(0), family = character(0),
      rotation = integer(0), param1 = numeric(0), se1 = numeric(0),
      param2 = numeric(0), se2 = numeric(0), kendall_tau = numeric(0),
      tail_dep_lower = numeric(0), tail_dep_upper = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(pairs_df) <- NULL

  params_df <- if (length(param_rows) > 0L) {
    do.call(rbind, param_rows)
  } else {
    data.frame(
      component = character(0), target = character(0), tree = integer(0),
      family = character(0), rotation = integer(0), parameter = character(0),
      estimate = numeric(0), se = numeric(0), stringsAsFactors = FALSE
    )
  }
  rownames(params_df) <- NULL

  n_free <- as.integer(unlist(raw$n_free_params, use.names = FALSE) %||% nrow(params_df))[[1L]]
  json_str <- as.character(unlist(raw$json, use.names = FALSE) %||% "")[[1L]]

  list(
    marginals = marginals_df,
    trees = trees_out,
    pairs = pairs_df,
    params = params_df,
    n_free_params = n_free,
    json = json_str
  )
}

# Shared formatted lines for `print.ferx_fit()` and `print.ferx_summary()`.
.ferx_vine_report_lines <- function(fit) {
  vc <- fit$vine_copula
  if (is.null(vc)) return(NULL)

  lines <- character(0)
  m <- vc$marginals
  if (is.data.frame(m) && nrow(m) > 0L) {
    lines <- c(lines, "  Marginals (zero-mean constrained):")
    for (i in seq_len(nrow(m))) {
      pv <- m$params[[i]]
      sv <- m$se[[i]]
      p_str <- if (length(pv) > 0L) {
        pnms <- names(pv) %||% paste0("p", seq_along(pv))
        parts <- vapply(seq_along(pv), function(k) {
          se_k <- if (length(sv) >= k && is.finite(sv[[k]])) sprintf(" (SE %.4f)", sv[[k]]) else ""
          sprintf("%s=%.4f%s", pnms[[k]], pv[[k]], se_k)
        }, character(1L))
        paste(parts, collapse = ", ")
      } else {
        "fixed"
      }
      extra <- character(0)
      if (!is.na(m$skewness[[i]]) && abs(m$skewness[[i]]) > 1e-4) {
        extra <- c(extra, sprintf("skew=%.3f", m$skewness[[i]]))
      }
      if (!is.na(m$excess_kurtosis[[i]]) && abs(m$excess_kurtosis[[i]]) > 1e-4) {
        extra <- c(extra, sprintf("ex.kurt=%.3f", m$excess_kurtosis[[i]]))
      }
      extra_str <- if (length(extra) > 0L) paste0("  [", paste(extra, collapse = ", "), "]") else ""
      lines <- c(lines, sprintf(
        "    %-12s %-12s SD = %.4f  (%s)%s",
        m$eta[[i]], m$family[[i]], m$sd[[i]], p_str, extra_str
      ))
    }
  }

  p <- vc$pairs
  if (is.data.frame(p) && nrow(p) > 0L) {
    lines <- c(lines, "  Regular-vine pair-copulas:")
    for (i in seq_len(nrow(p))) {
      tr_pairs <- vc$trees
      # Format named parameters directly from the structured tree entry when available.
      p_fmt <- ""
      for (tr in tr_pairs) {
        if (identical(as.integer(tr$tree), as.integer(p$tree[[i]]))) {
          for (pe in tr$pairs) {
            if (identical(pe$label, p$label[[i]])) {
              pv <- pe$params
              sv <- pe$se
              if (length(pv) > 0L) {
                pnms <- names(pv) %||% paste0("p", seq_along(pv))
                parts <- vapply(seq_along(pv), function(k) {
                  se_k <- if (length(sv) >= k && is.finite(sv[[k]])) sprintf(" (SE %.4f)", sv[[k]]) else ""
                  sprintf("%s=%.4f%s", pnms[[k]], pv[[k]], se_k)
                }, character(1L))
                p_fmt <- paste0(" [", paste(parts, collapse = ", "), "]")
              }
            }
          }
        }
      }
      lines <- c(lines, sprintf(
        "    T%d  %-24s %-13s tau = %+.4f  lambda_L = %.3f  lambda_U = %.3f%s",
        p$tree[[i]], p$label[[i]], p$family[[i]],
        p$kendall_tau[[i]], p$tail_dep_lower[[i]], p$tail_dep_upper[[i]], p_fmt
      ))
    }
  }
  lines
}

#' Evaluate the fitted regular-vine random-effect density
#'
#' Evaluates the joint random-effect density \eqn{p_{\text{vine}}(\eta)} (or a
#' 1-D zero-mean marginal density / 2-D joint slice) from a model fitted with
#' \code{settings = list(omega_dist = "vine")}.
#'
#' When \code{at} is omitted, \code{ferx_vine_density()} builds an evenly spaced
#' grid spanning \code{[-range_sd * sd, +range_sd * sd]} on the selected random
#' effects (\code{etas}):
#' \itemize{
#'   \item \strong{1 random effect selected} (\code{length(etas) == 1L}, or a
#'     single-ETA model): evaluates the exact 1-D marginal density
#'     \eqn{f_j(\eta_j)} when \code{fixed_etas} is \code{NULL}, or the unnormalized
#'     1-D slice through \code{fixed_etas} when \code{fixed_etas} is supplied.
#'   \item \strong{2 random effects selected} (the default when the model has at
#'     least 2 random effects): evaluates the \code{n x n} grid over the two
#'     selected ETAs, holding any remaining random effects at \code{0} (or at the
#'     named values in \code{fixed_etas}). On a 2-ETA model this is the exact
#'     joint density \eqn{p_{\text{vine}}(\eta_1, \eta_2)}, which integrates to 1
#'     over \eqn{\mathbb{R}^2}.
#' }
#'
#' When \code{at} is supplied (a numeric matrix or data frame with one column per
#' model random effect), \code{ferx_vine_density()} evaluates the exact joint
#' density at every row of \code{at}.
#'
#' @param fit A \code{ferx_fit} object fitted with \code{omega_dist = "vine"}.
#' @param etas Optional character vector of 1 or 2 ETA names (from
#'   \code{fit$eta_names}) or integer indices selecting the grid axes when
#'   \code{at} is \code{NULL}. Defaults to the first two random effects (or the
#'   sole random effect when \code{nrow(fit$omega) == 1L}).
#' @param n Integer grid resolution per axis (default \code{50L}, at least
#'   \code{2L}).
#' @param range_sd Positive numeric half-width in units of each selected ETA's
#'   fitted marginal SD (default \code{3.5}).
#' @param at Optional numeric matrix or data frame of evaluation points with one
#'   column per random effect in \code{fit$eta_names} (named columns are aligned
#'   to \code{fit$eta_names}). Overrides \code{etas}, \code{n}, and
#'   \code{range_sd}.
#' @param fixed_etas Optional named numeric vector giving fixed conditioning
#'   values for random effects not in \code{etas} when \code{at} is \code{NULL}
#'   (defaults to \code{0} for unspecified coordinates).
#' @param log Logical; if \code{TRUE}, include only log-density in the
#'   \code{density} column? No: \code{density} always holds \eqn{p(\eta)} and
#'   \code{log_density} always holds \eqn{\log p(\eta)}; when \code{log = TRUE}
#'   \code{density} is computed via \code{exp(log_density)} without extra FFI
#'   evaluations. Kept for convenience when scripting.
#' @return A \code{data.frame} containing the evaluated ETA coordinates alongside
#'   \code{density} and \code{log_density} columns.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("vine_tail_iv")
#' fit <- ferx_fit(ex$model, ex$data)
#' grid <- ferx_vine_density(fit, n = 60L)
#' head(grid)
#' }
#' @family diagnostics
#' @export
ferx_vine_density <- function(fit, etas = NULL, n = 50L, range_sd = 3.5,
                              at = NULL, fixed_etas = NULL, log = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a `ferx_fit` object.", call. = FALSE)
  }
  vc <- fit$vine_copula
  if (is.null(vc) || !nzchar(vc$json %||% "")) {
    stop("`fit` does not carry a fitted vine copula (`fit$vine_copula` is NULL). ",
         "Fit the model with `method = \"saem\"` and `settings = list(omega_dist = \"vine\")`.",
         call. = FALSE)
  }

  eta_names <- as.character(vc$marginals$eta)
  d <- length(eta_names)
  if (d == 0L) {
    stop("`fit` has zero random effects.", call. = FALSE)
  }

  if (!is.null(at)) {
    if (is.data.frame(at)) {
      if (!is.null(colnames(at)) && all(eta_names %in% colnames(at))) {
        at <- at[, eta_names, drop = FALSE]
      }
      mat <- as.matrix(at)
    } else if (is.matrix(at)) {
      if (!is.null(colnames(at)) && all(eta_names %in% colnames(at))) {
        mat <- at[, eta_names, drop = FALSE]
      } else {
        mat <- at
      }
    } else if (is.numeric(at) && length(at) == d) {
      mat <- matrix(at, nrow = 1L, ncol = d)
    } else {
      stop("`at` must be a numeric matrix or data.frame with ", d,
           " columns matching `fit$eta_names`.", call. = FALSE)
    }
    storage.mode(mat) <- "double"
    if (ncol(mat) != d) {
      stop("`at` has ", ncol(mat), " columns, expected ", d, " (`",
           paste(eta_names, collapse = "`, `"), "`).", call. = FALSE)
    }
    if (any(!is.finite(mat))) {
      stop("`at` must contain finite numeric values.", call. = FALSE)
    }
    nr <- nrow(mat)
    eta_flat <- as.numeric(t(mat))
    lp <- ferx_rust_vine_density(vc$json, eta_flat, nr, d, 0L, TRUE)
    out <- as.data.frame(mat, stringsAsFactors = FALSE)
    colnames(out) <- eta_names
    out$density <- exp(lp)
    out$log_density <- lp
    if (isTRUE(log)) out$value <- lp else out$value <- out$density
    out$value <- NULL
    return(out)
  }

  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be a single integer >= 2.", call. = FALSE)
  }
  if (length(range_sd) != 1L || !is.numeric(range_sd) || !is.finite(range_sd) || range_sd <= 0) {
    stop("`range_sd` must be a single positive finite number.", call. = FALSE)
  }

  # Resolve `etas` to 1-based indices into `eta_names`.
  idx <- if (is.null(etas)) {
    seq_len(min(2L, d))
  } else if (is.character(etas)) {
    m <- match(etas, eta_names)
    if (anyNA(m)) {
      bad <- etas[is.na(m)]
      stop("Unknown ETA name(s) in `etas`: ", paste(bad, collapse = ", "),
           ". Available: ", paste(eta_names, collapse = ", "), ".", call. = FALSE)
    }
    m
  } else if (is.numeric(etas)) {
    m <- as.integer(etas)
    if (anyNA(m) || any(m < 1L | m > d)) {
      stop("Integer `etas` indices must lie in 1..", d, ".", call. = FALSE)
    }
    m
  } else {
    stop("`etas` must be NULL, a character vector of ETA names, or integer indices.",
         call. = FALSE)
  }
  idx <- unique(idx)
  if (length(idx) < 1L || length(idx) > 2L) {
    stop("`etas` must select 1 or 2 random effects for grid evaluation (got ",
         length(idx), "); pass `at =` to evaluate arbitrary ", d, "-D points.",
         call. = FALSE)
  }

  sds <- as.numeric(vc$marginals$sd)
  sds[!is.finite(sds) | sds <= 0] <- 0.1

  # Parse optional `fixed_etas`.
  base_row <- numeric(d)
  if (!is.null(fixed_etas)) {
    if (!is.numeric(fixed_etas) || is.null(names(fixed_etas)) || any(!is.finite(fixed_etas))) {
      stop("`fixed_etas` must be a named finite numeric vector.", call. = FALSE)
    }
    fm <- match(names(fixed_etas), eta_names)
    if (anyNA(fm)) {
      stop("Unknown ETA name(s) in `fixed_etas`: ",
           paste(names(fixed_etas)[is.na(fm)], collapse = ", "), ".", call. = FALSE)
    }
    base_row[fm] <- as.numeric(fixed_etas)
  }

  if (length(idx) == 1L) {
    j <- idx[[1L]]
    grid_x <- seq(-range_sd * sds[[j]], range_sd * sds[[j]], length.out = n)
    if (is.null(fixed_etas)) {
      # Exact 1-D marginal density f_j(eta_j), integrating to 1 over R.
      lp <- ferx_rust_vine_density(vc$json, grid_x, n, 1L, j, TRUE)
    } else {
      mat <- matrix(rep(base_row, each = n), nrow = n, ncol = d)
      mat[, j] <- grid_x
      lp <- ferx_rust_vine_density(vc$json, as.numeric(t(mat)), n, d, 0L, TRUE)
    }
    out <- data.frame(grid_x, density = exp(lp), log_density = lp, stringsAsFactors = FALSE)
    names(out)[[1L]] <- eta_names[[j]]
    return(out)
  }

  j1 <- idx[[1L]]
  j2 <- idx[[2L]]
  g1 <- seq(-range_sd * sds[[j1]], range_sd * sds[[j1]], length.out = n)
  g2 <- seq(-range_sd * sds[[j2]], range_sd * sds[[j2]], length.out = n)
  eg <- expand.grid(g1, g2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  nr <- nrow(eg)
  mat <- matrix(rep(base_row, each = nr), nrow = nr, ncol = d)
  mat[, j1] <- eg[[1L]]
  mat[, j2] <- eg[[2L]]
  lp <- ferx_rust_vine_density(vc$json, as.numeric(t(mat)), nr, d, 0L, TRUE)

  out <- data.frame(
    eg[[1L]],
    eg[[2L]],
    density = exp(lp),
    log_density = lp,
    stringsAsFactors = FALSE
  )
  names(out)[1:2] <- eta_names[c(j1, j2)]
  out
}
