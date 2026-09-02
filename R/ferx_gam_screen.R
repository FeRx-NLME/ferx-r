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
#' All numerical computation (OLS, AIC, spline basis construction) is
#' performed in Rust via \code{ferx_rust_gam_screen()}. This wrapper handles
#' input validation, per-subject covariate aggregation, and result formatting.
#'
#' The AIC formula used is \code{n * log(RSS/n) + 2*p}, which gives the same
#' delta_aic as R's \code{AIC(lm())} and Xpose4's \code{gam::gam()}.
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
#'   Xpose4's defaults.
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
  ebe      <- fit$ebe_etas
  eta_id   <- if ("ID" %in% names(ebe)) "ID" else names(ebe)[1L]
  # `eta_cols_all` stays the unfiltered column order: `fit$shrinkage_eta` is
  # positional against it, so subsetting before the shrinkage lookup would
  # silently shift every ETA's shrinkage (see below).
  eta_cols_all <- setdiff(names(ebe), eta_id)
  all_etas     <- eta_cols_all
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
  # Aggregate on the ordinal subject index, not the raw ID: a dataset that
  # reuses a subject ID in a non-contiguous block has two distinct subjects
  # sharing that ID, and keying on ID would collapse them into one row (and
  # later double-weight that row against both subjects' ETAs). Same treatment
  # as `ferx_cov_screen()`.
  mode_fn <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NA)
    ux <- unique(v)
    ux[which.max(tabulate(match(v, ux)))]
  }
  subj <- .ferx_subject_index(covtab[[cov_id]])
  ids  <- covtab[[cov_id]][!duplicated(subj)]  # one ID per subject, subject order
  keys <- as.character(seq_along(ids))
  percov <- stats::setNames(
    data.frame(as.character(ids), stringsAsFactors = FALSE), cov_id
  )
  for (cov in cov_names) {
    if (identical(types[[cov]], "categorical")) {
      agg <- tapply(covtab[[cov]], subj, mode_fn)
    } else {
      agg <- tapply(covtab[[cov]], subj,
                    function(v) stats::median(v, na.rm = TRUE))
    }
    percov[[cov]] <- agg[keys]
  }

  # -- Align EBEs and per-subject covariates ----------------------------------
  # `percov` and `ebe` are each one row per subject in subject order, so align
  # positionally - an ID merge would cross-join two subjects that share a
  # reused ID. Fall back to an ID match when the counts disagree (and index
  # with match() rather than merge() for the same reason).
  # The two frames are kept apart rather than cbind()-ed: a covariate that
  # shares its name with an ETA column would otherwise give the combined frame
  # two columns of that name, and `[[` would silently return the wrong one.
  if (nrow(percov) == nrow(ebe)) {
    cov_rows <- seq_len(nrow(percov))
    eta_rows <- seq_len(nrow(ebe))
  } else {
    pidx     <- match(ebe[[eta_id]], percov[[cov_id]])
    eta_rows <- which(!is.na(pidx))
    cov_rows <- pidx[eta_rows]
  }
  n_subj <- length(eta_rows)

  if (n_subj == 0L) {
    message("No subjects in common between `fit$covtab` and `fit$ebe_etas`.")
    return(invisible(NULL))
  }

  # -- Shrinkage vector -------------------------------------------------------
  # `fit$shrinkage_eta` is an unnamed numeric vector, one entry per random
  # effect in declaration order, with the labels carried separately in
  # `fit$eta_names` (see `print.ferx_fit`). So the name lookup has to consult
  # `fit$eta_names` too, and the positional fallback has to index the
  # *unfiltered* ETA column order - indexing `all_etas` would report the first
  # ETA's shrinkage for whatever ETA `etas` happened to select.
  shrink_vec  <- fit$shrinkage_eta
  eta_names_r <- names(shrink_vec)
  if (is.null(eta_names_r) &&
      length(fit$eta_names) == length(shrink_vec)) {
    eta_names_r <- fit$eta_names
  }
  shrinkage_v <- vapply(all_etas, function(eta_name) {
    if (is.null(shrink_vec) || !length(shrink_vec)) return(NA_real_)
    if (!is.null(eta_names_r)) {
      idx <- match(eta_name, eta_names_r)
      if (!is.na(idx)) return(shrink_vec[[idx]])
    }
    idx <- match(eta_name, eta_cols_all)
    if (!is.na(idx) && idx <= length(shrink_vec)) shrink_vec[[idx]] else NA_real_
  }, numeric(1L))

  # -- Build column-major flat arrays for Rust --------------------------------
  # ETAs: each column is all per-subject values for one ETA.
  eta_flat <- unlist(lapply(all_etas, function(e) as.numeric(ebe[[e]][eta_rows])),
                     use.names = FALSE)

  # Covariates: coerce to numeric. For categorical covariates stored as
  # character, convert to factor integer codes (1, 2, ...) so Rust receives
  # a numeric vector with distinct integer levels for one-hot encoding.
  cov_flat <- unlist(lapply(cov_names, function(cov) {
    vals <- percov[[cov]][cov_rows]
    if (identical(types[[cov]], "categorical")) {
      as.numeric(as.factor(vals))
    } else {
      as.numeric(vals)
    }
  }), use.names = FALSE)

  cov_kinds_str <- unname(
    ifelse(types[cov_names] == "categorical", "categorical", "continuous")
  )

  # -- Delegate all computation to Rust ---------------------------------------
  raw <- ferx_rust_gam_screen(
    eta_names      = all_etas,
    eta_flat       = eta_flat,
    n_subjects     = n_subj,
    shrinkage      = shrinkage_v,
    cov_names      = cov_names,
    cov_flat       = cov_flat,
    cov_kinds      = cov_kinds_str,
    spline_df      = as.integer(spline_df),
    include_linear = isTRUE(include_linear),
    shrinkage_warn = shrinkage_warn
  )

  # Emit any shrinkage warnings collected by Rust.
  for (w in raw$warnings) warning(w, call. = FALSE)

  # -- Assemble result data frame ---------------------------------------------
  result <- data.frame(
    eta_name  = raw$eta_name,
    covariate = raw$covariate,
    delta_aic = raw$delta_aic,
    best_form = raw$best_form,
    aic       = raw$aic,
    aic_null  = raw$aic_null,
    r_squared = raw$r_squared,
    shrinkage = raw$shrinkage,
    stringsAsFactors = FALSE
  )

  if (nrow(result) == 0L) {
    message("No covariate pairs could be screened.")
    return(invisible(result))
  }

  # Order by ETA (preserving input order), then by delta_aic descending.
  eta_order <- match(result$eta_name, all_etas)
  result    <- result[order(eta_order, -result$delta_aic), ]
  rownames(result) <- NULL

  print(result)
  invisible(result)
}
