#' GAM-based covariate pre-screening
#'
#' Screens all declared covariates against each ETA using independent
#' generalised additive model (GAM) regressions. For each ETA x covariate pair
#' the function fits \code{eta ~ f(cov)} and reports the AIC improvement over
#' the null model \code{eta ~ 1}. Covariates are ranked by
#' \code{delta_aic = AIC_null - AIC_best}; a positive value means the covariate
#' improves the null model for that ETA.
#'
#' This is the R equivalent of Xpose4's \code{xpose.gam()} (Jonsson and
#' Karlsson, \emph{Pharm Res} 1999). Like Xpose4, it uses independent
#' regressions (not stepwise backfitting), which is appropriate for a
#' pre-screening role where speed and interpretability matter most.
#'
#' The functional form candidates for continuous covariates are:
#' \itemize{
#'   \item Linear: \code{eta ~ 1 + x}
#'   \item Natural cubic spline: \code{eta ~ 1 + ns(x, df)} for each
#'     \code{df} in \code{spline_df} (when \code{n > df + 1})
#' }
#' For categorical covariates: one-hot encoding with the lowest observed level
#' as the reference.
#'
#' The AIC formula used is \code{n * log(RSS/n) + 2*p}, which gives the same
#' delta_aic as R's \code{AIC(lm())} and Xpose4's \code{gam::gam()}. The
#' formula is identical to the one used by the Rust implementation in
#' \code{ferx_tools::gam::gam_screen()}, so results are directly comparable.
#'
#' @section Shrinkage caveat:
#' EBE-based covariate screening is only informative when ETA shrinkage is low
#' (< 30\%). At high shrinkage the EBEs regress toward zero and the
#' relationship between the estimated ETA and a covariate is attenuated.
#' A warning is emitted for each ETA whose shrinkage exceeds
#' \code{shrinkage_warn}.
#'
#' @param fit A \code{ferx_fit} object with \code{fit$ebe_etas} and
#'   \code{fit$covtab} populated (the model must declare a
#'   \code{[covariates]} block).
#' @param etas Character vector of ETA names to screen. \code{NULL} (default)
#'   screens all ETAs present in \code{fit$ebe_etas}.
#' @param covariates Character vector of covariate names to screen.
#'   \code{NULL} (default) uses all covariates declared in
#'   \code{fit$covariate_types} (or inferred from \code{fit$covtab}).
#' @param spline_df Integer vector of natural-spline degrees of freedom to
#'   try for continuous covariates. Default \code{c(2L, 3L)} matches
#'   Xpose4's defaults (\code{smoother3 = "ns", arg3 = "df=2"} and
#'   \code{smoother4 = "ns", arg4 = "df=3"}).
#' @param include_linear Logical. Include the linear form as a candidate.
#'   Default \code{TRUE}.
#' @param shrinkage_warn Numeric fraction in \code{[0, 1]}. Warn when an
#'   ETA's shrinkage (from \code{fit$shrinkage_eta}) exceeds this threshold.
#'   Default \code{0.30}.
#' @return A data frame (returned invisibly; also printed) with one row per
#'   ETA x covariate pair and columns:
#'   \describe{
#'     \item{eta_name}{ETA name (e.g. \code{"ETA_CL"}).}
#'     \item{covariate}{Covariate name.}
#'     \item{delta_aic}{\code{AIC_null - AIC_best}. Positive means the
#'       covariate improves the null model.}
#'     \item{best_form}{Winning functional form: \code{"Linear"},
#'       \code{"Spline(df=2)"}, \code{"Spline(df=3)"}, or
#'       \code{"Categorical"}.}
#'     \item{aic}{AIC of the best model.}
#'     \item{aic_null}{AIC of the null (intercept-only) model.}
#'     \item{r_squared}{R-squared of the best model.}
#'     \item{shrinkage}{ETA shrinkage from the fit (fraction).}
#'   }
#'   Returns \code{NULL} with a message when the fit has no ETA or covariate
#'   data. Returns an empty data frame when no pairs can be screened.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("two_cpt_oral_cov")
#' fit <- ferx_fit(ex$model, ex$data, method = "focei", covariance = FALSE)
#' ferx_gam_screen(fit)
#' }
#' @family diagnostics
#' @importFrom splines ns
#' @export
ferx_gam_screen <- function(fit,
                             etas        = NULL,
                             covariates  = NULL,
                             spline_df   = c(2L, 3L),
                             include_linear = TRUE,
                             shrinkage_warn = 0.30) {

  # -- Guard: need ebe_etas ---------------------------------------------------
  if (is.null(fit$ebe_etas) || !is.data.frame(fit$ebe_etas)) {
    message("No ETAs (`fit$ebe_etas`) - the model declares no ",
            "inter-individual variability.")
    return(invisible(NULL))
  }

  # -- Guard: need covtab -----------------------------------------------------
  if (is.null(fit$covtab) || !is.data.frame(fit$covtab)) {
    message("No covariate table (`fit$covtab`) - the model needs a ",
            "[covariates] block.")
    return(invisible(NULL))
  }

  # -- ETA columns ------------------------------------------------------------
  ebe       <- fit$ebe_etas
  eta_id    <- if ("ID" %in% names(ebe)) "ID" else names(ebe)[1L]
  all_etas  <- setdiff(names(ebe), eta_id)
  if (!is.null(etas)) {
    all_etas <- intersect(etas, all_etas)
  }
  if (!length(all_etas)) {
    message("No ETA columns found in `fit$ebe_etas`.")
    return(invisible(NULL))
  }
  ebe[[eta_id]] <- as.character(ebe[[eta_id]])

  # -- Covariate names + types ------------------------------------------------
  covtab  <- fit$covtab
  cov_id  <- if ("ID" %in% names(covtab)) "ID" else names(covtab)[1L]

  types <- fit$covariate_types
  if (!is.null(types) && length(types)) {
    cov_names <- names(types)
  } else {
    cov_names <- setdiff(names(covtab), c(cov_id, "TIME", "EVID"))
    types <- stats::setNames(
      ifelse(vapply(covtab[cov_names], is.numeric, logical(1L)),
             "continuous", "categorical"),
      cov_names
    )
  }
  if (!is.null(covariates)) {
    cov_names <- intersect(covariates, cov_names)
  }
  cov_names <- intersect(cov_names, names(covtab))
  if (!length(cov_names)) {
    message("No declared covariates found in `fit$covtab`.")
    return(invisible(NULL))
  }

  # -- Aggregate covariates to one value per subject --------------------------
  # Median for continuous, most frequent level for categorical.
  mode_fn <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NA)
    ux <- unique(v)
    ux[which.max(tabulate(match(v, ux)))]
  }
  ids <- unique(covtab[[cov_id]])
  percov <- stats::setNames(
    data.frame(as.character(ids), stringsAsFactors = FALSE), cov_id
  )
  for (cov in cov_names) {
    if (identical(types[[cov]], "categorical")) {
      agg <- tapply(covtab[[cov]], covtab[[cov_id]], mode_fn)
    } else {
      agg <- tapply(covtab[[cov]], covtab[[cov_id]],
                    function(v) stats::median(v, na.rm = TRUE))
    }
    percov[[cov]] <- agg[as.character(ids)]
  }

  # -- Shrinkage lookup -------------------------------------------------------
  shrink_vec  <- fit$shrinkage_eta
  eta_names_r <- if (!is.null(names(shrink_vec))) names(shrink_vec) else NULL

  get_shrinkage <- function(eta_name) {
    if (is.null(shrink_vec) || !length(shrink_vec)) return(NA_real_)
    if (!is.null(eta_names_r)) {
      idx <- match(eta_name, eta_names_r)
      if (!is.na(idx)) return(shrink_vec[[idx]])
    }
    # Fall back to positional lookup using all_etas order.
    idx <- match(eta_name, all_etas)
    if (!is.na(idx) && idx <= length(shrink_vec)) shrink_vec[[idx]] else NA_real_
  }

  # -- OLS AIC helper ---------------------------------------------------------
  # AIC = n * log(RSS/n) + 2*p. Same formula as the Rust implementation and
  # equivalent delta-AIC to R's AIC(lm()).
  ols_aic <- function(y, X) {
    n   <- length(y)
    p   <- ncol(X)
    fit_ols <- lm.fit(X, y)
    rss <- sum(fit_ols$residuals^2)
    aic <- n * log(rss / n) + 2 * p
    sst <- sum((y - mean(y))^2)
    r2  <- if (sst < 1e-20) 0 else 1 - rss / sst
    list(aic = aic, r2 = r2)
  }

  # -- Screen -----------------------------------------------------------------
  merged <- merge(percov, ebe, by.x = cov_id, by.y = eta_id)

  rows <- vector("list", length(all_etas) * length(cov_names))
  k    <- 0L

  for (eta in all_etas) {

    shrinkage <- get_shrinkage(eta)
    if (!is.na(shrinkage) && is.finite(shrinkage) &&
        shrinkage > shrinkage_warn) {
      warning(sprintf(
        "%s: shrinkage %.1f%% exceeds the %.0f%% threshold; ",
        eta, shrinkage * 100, shrinkage_warn * 100
      ), "EBE-based covariate screening may be unreliable.",
      call. = FALSE)
    }

    y_all   <- as.numeric(merged[[eta]])
    valid_y <- is.finite(y_all)
    y_sub   <- y_all[valid_y]
    n_sub   <- length(y_sub)
    if (n_sub < 3L) next

    # Null AIC on the valid-ETA subset.
    null_res  <- ols_aic(y_sub, matrix(1, nrow = n_sub, ncol = 1L))
    aic_null  <- null_res$aic

    for (cov in cov_names) {
      k <- k + 1L
      x_all    <- merged[[cov]]
      valid    <- valid_y & !is.na(x_all)
      y        <- y_all[valid]
      x        <- x_all[valid]
      n        <- sum(valid)
      if (n < 3L) {
        rows[[k]] <- NULL
        next
      }

      # Per-pair null AIC (same subset as alternative).
      pn_res   <- ols_aic(y, matrix(1, nrow = n, ncol = 1L))
      aic_loc  <- pn_res$aic

      best_aic  <- Inf
      best_r2   <- 0
      best_form <- NA_character_

      type <- types[[cov]]

      if (identical(type, "categorical")) {
        levs <- sort(unique(x))
        if (length(levs) >= 2L) {
          dummies <- outer(x, levs[-1L], `==`) * 1.0
          X_cat   <- cbind(1, dummies)
          res     <- ols_aic(y, X_cat)
          if (res$aic < best_aic) {
            best_aic  <- res$aic
            best_r2   <- res$r2
            best_form <- "Categorical"
          }
        }
      } else {
        x_num <- as.numeric(x)
        if (include_linear) {
          res <- ols_aic(y, cbind(1, x_num))
          if (res$aic < best_aic) {
            best_aic  <- res$aic
            best_r2   <- res$r2
            best_form <- "Linear"
          }
        }
        for (df in spline_df) {
          if (n <= df + 1L) next
          basis <- tryCatch(
            splines::ns(x_num, df = df),
            error = function(e) NULL
          )
          if (is.null(basis)) next
          X_spl <- cbind(1, basis)
          res   <- ols_aic(y, X_spl)
          if (res$aic < best_aic) {
            best_aic  <- res$aic
            best_r2   <- res$r2
            best_form <- sprintf("Spline(df=%d)", df)
          }
        }
      }

      if (is.infinite(best_aic)) {
        rows[[k]] <- NULL
        next
      }

      rows[[k]] <- data.frame(
        eta_name  = eta,
        covariate = cov,
        delta_aic = aic_loc - best_aic,
        best_form = best_form,
        aic       = best_aic,
        aic_null  = aic_loc,
        r_squared = best_r2,
        shrinkage = shrinkage,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, Filter(Negate(is.null), rows))

  if (is.null(result) || nrow(result) == 0L) {
    message("No covariate pairs could be screened.")
    return(invisible(data.frame(
      eta_name = character(), covariate = character(),
      delta_aic = numeric(), best_form = character(),
      aic = numeric(), aic_null = numeric(),
      r_squared = numeric(), shrinkage = numeric(),
      stringsAsFactors = FALSE
    )))
  }

  # Order by ETA (preserving input order), then by delta_aic descending.
  eta_order <- match(result$eta_name, all_etas)
  result    <- result[order(eta_order, -result$delta_aic), ]
  rownames(result) <- NULL

  print(result)
  invisible(result)
}
