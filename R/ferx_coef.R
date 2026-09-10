# Name-based accessors for a fit's parameter estimates and standard errors.
# `fit$estimates` carries both a `param` column and (since #299) row names, so
# `fit$estimates["TVCL", "estimate"]` works - but a typo there still returns NA
# rather than erroring, because that is what `[` does to a data frame. These
# accessors are the loud version: an unknown name is an error naming the
# closest matches, so a mistyped coefficient can never be read as a missing one.
# They resolve against the row names, which are block-qualified where a
# declared name collides, so a colliding name is never resolved by position.

#' Pull a parameter estimate or its standard error by name
#'
#' Name-based accessors for \code{fit$estimates}. \code{ferx_coef()} returns the
#' point estimate on the scale it was estimated on (thetas as declared - so on
#' the log scale for a \code{log}-transformed theta; omega and kappa as
#' variances); \code{ferx_se()} returns the matching standard error from the
#' covariance step.
#'
#' Unlike \code{fit$estimates[param, ]}, an unrecognised name is an
#' \strong{error} listing the closest available names, not a silent
#' \code{NA} - a mistyped parameter cannot be mistaken for an unestimated one.
#'
#' @param fit A \code{ferx_fit} object from \code{\link{ferx_fit}} (or
#'   \code{\link{ferx_load_fit}}).
#' @param param Character vector of parameter names, as they appear in the row
#'   names of \code{fit$estimates}: theta names (\code{"TVCL"}), eta names for
#'   the omega diagonal (\code{"ETA_CL"}), sigma names (\code{"EPS_PROP"}) and
#'   kappa names for an IOV diagonal (\code{"KAPPA_CL"}). Unnamed parameters
#'   use the fallback labels \code{"THETA1"}, \code{"OMEGA(1,1)"},
#'   \code{"SIGMA(1)"}, \code{"KAPPA1"}. A name declared in more than one
#'   block is qualified by that block - a \code{CL} declared as both a theta
#'   and an eta is addressed as \code{"CL.theta"} and \code{"CL.omega"}, and
#'   the bare \code{"CL"} is an error rather than one of the two chosen
#'   silently. \code{NULL} (default) returns every parameter, in table order.
#' @return A named numeric vector, one element per requested name, in the order
#'   requested. For \code{ferx_se()}, \code{NA} where the covariance step did
#'   not run or produced no standard error for that parameter; a parameter
#'   declared \code{FIX} carries the engine's exact \code{0}, not \code{NA}.
#' @seealso \code{fit$estimates} for the full table (SE, relative standard
#'   error, confidence intervals, natural-scale back-transforms) and
#'   \code{fit$cor_matrix} for parameter correlations.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn")
#' ferx_coef(fit)                     # every parameter
#' ferx_coef(fit, "TVCL")             # one theta
#' ferx_se(fit, c("TVCL", "TVV"))     # standard errors
#' @family fitting
#' @export
ferx_coef <- function(fit, param = NULL) {
  .ferx_estimates_column(fit, param, "estimate", "ferx_coef")
}

#' @rdname ferx_coef
#' @export
ferx_se <- function(fit, param = NULL) {
  .ferx_estimates_column(fit, param, "se", "ferx_se")
}

# Shared body: validate, resolve names against `fit$estimates$param`, and pull
# one column as a named numeric vector.
.ferx_estimates_column <- function(fit, param, column, fn) {
  if (!inherits(fit, "ferx_fit")) {
    stop(sprintf("`%s()` needs a `ferx_fit` object (from ferx_fit() or ferx_load_fit()).", fn),
         call. = FALSE)
  }
  est <- fit$estimates
  if (is.null(est) || !is.data.frame(est) || nrow(est) == 0L) {
    stop(sprintf("`%s()`: this fit carries no `estimates` table.", fn), call. = FALSE)
  }

  # Row names are the addressable keys: identical to `param` except where a
  # name is declared in more than one block, where each is qualified by its
  # block so neither owns the bare name (see .ferx_estimate_keys()).
  available <- rownames(est)
  if (is.null(available)) available <- est$param
  if (is.null(param)) {
    idx <- seq_len(nrow(est))
  } else {
    if (!is.character(param)) {
      stop(sprintf("`%s()`: `param` must be a character vector of parameter names (or NULL).", fn),
           call. = FALSE)
    }
    idx     <- match(param, available)
    unknown <- param[is.na(idx)]
    if (length(unknown) > 0L) {
      # A bare name that collides across blocks resolves to no single row. Say
      # so and name the qualified keys, rather than silently handing back the
      # first - the silent-wrong-coefficient failure these accessors exist to
      # prevent.
      ambiguous <- unknown[unknown %in% est$param]
      if (length(ambiguous) > 0L) {
        stop(sprintf("`%s()`: %s declared in more than one block, so it does not name a single parameter.\nUse one of: %s.",
                     fn,
                     paste(sprintf("\"%s\" is", ambiguous), collapse = ", "),
                     paste(sprintf("\"%s\"", available[est$param %in% ambiguous]), collapse = ", ")),
             call. = FALSE)
      }
      stop(sprintf("`%s()`: unknown parameter%s %s.\n%s",
                   fn,
                   if (length(unknown) > 1L) "s" else "",
                   paste(sprintf("\"%s\"", unknown), collapse = ", "),
                   .ferx_param_hint(unknown, available)),
           call. = FALSE)
    }
  }

  out        <- as.numeric(est[[column]][idx])
  names(out) <- available[idx]
  if (column == "se" && length(out) > 0L && all(is.na(out))) {
    warning("This fit carries no standard errors for the requested parameters: ",
            "the covariance step was not run (`covariance = FALSE`), or it failed. ",
            "See `fit$covariance_status`.", call. = FALSE)
  }
  out
}

# "Did you mean" line for an unknown parameter name: the closest available
# names by edit distance, else the full list when it is short enough to print.
# Case-insensitive, because ETA_CL / eta_cl is the common slip.
.ferx_param_hint <- function(unknown, available) {
  near <- unique(unlist(lapply(unknown, function(u) {
    d <- utils::adist(u, available, ignore.case = TRUE)[1L, ]
    # Allow one edit per three characters, at least one - enough for a typo or
    # a case slip, not enough to propose an unrelated parameter.
    keep <- which(d <= max(1L, floor(nchar(u) / 3L)))
    available[keep[order(d[keep])]][seq_len(min(3L, length(keep)))]
  }), use.names = FALSE))
  near <- near[!is.na(near)]
  if (length(near) > 0L) {
    return(paste0("Did you mean: ", paste(near, collapse = ", "), "?"))
  }
  if (length(available) <= 30L) {
    return(paste0("Available: ", paste(available, collapse = ", "), "."))
  }
  paste0("Available (first 30 of ", length(available), "): ",
         paste(utils::head(available, 30L), collapse = ", "), ", ...")
}
