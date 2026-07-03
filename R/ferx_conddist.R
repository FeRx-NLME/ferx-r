#' Per-subject SAEM conditional distribution of random effects
#'
#' Returns the per-subject / per-eta conditional distribution summary
#' produced by the opt-in SAEM conditional-distribution pass
#' (\code{[fit_options] conddist = true}). This is the shrinkage-unbiased
#' analogue of saemix's \code{conddist.saemix} / Monolix's "Conditional
#' Distribution" task: rather than only the conditional \strong{mode} (the
#' EBE on \code{fit$ebe_etas}), it reports the conditional \strong{mean} and
#' \strong{SD} of each subject's eta, i.e. \code{p(eta_i | y_i; theta_hat)}.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} (or
#'   \code{\link{ferx_load_fit}}), fit with \code{method = "saem"} and
#'   \code{settings = list(conddist = TRUE)}.
#' @return A data frame (classed \code{ferx_conddist}) with one row per
#'   subject x eta: \code{ID}, \code{ETA}, \code{COND_MEAN}, \code{COND_SD},
#'   \code{COND_MODE}. The distribution-based eta-shrinkage (parallel to
#'   \code{fit$eta_names}) and the retained-draw bookkeeping are attached as
#'   attributes \code{"shrinkage"}, \code{"nsamp"}, \code{"burnin"}.
#' @family diagnostics
#' @export
#' @seealso \code{\link{ferx_fit}}'s \code{conddist} fit option (see
#'   \code{[fit_options]} in the model-file docs); \code{\link{ferx_get_warnings}}
#'   for other post-fit diagnostics.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin_saem")
#' fit <- ferx_fit(ex$model, ex$data, method = "saem",
#'                  settings = list(conddist = TRUE))
#' cd  <- ferx_conddist(fit)
#' head(cd)
#' attr(cd, "shrinkage")
#' }
ferx_conddist <- function(fit) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object.")
  }
  cd <- fit$cond_dist
  if (is.null(cd) || is.null(cd$data)) {
    stop(
      "This fit has no conditional-distribution results. Re-fit with ",
      "method = \"saem\" and settings = list(conddist = TRUE) (see ?ferx_fit)."
    )
  }

  out <- cd$data
  attr(out, "shrinkage") <- if (length(cd$shrinkage) == length(fit$eta_names)) {
    stats::setNames(cd$shrinkage, fit$eta_names)
  } else {
    cd$shrinkage
  }
  attr(out, "nsamp")  <- cd$nsamp
  attr(out, "burnin") <- cd$burnin
  class(out) <- c("ferx_conddist", class(out))
  out
}

#' @export
print.ferx_conddist <- function(x, ...) {
  nsamp  <- attr(x, "nsamp")
  burnin <- attr(x, "burnin")
  provenance <- if (is.na(nsamp) || is.na(burnin)) {
    "draws retained/burn-in unknown (loaded from a .fitrx bundle)"
  } else {
    sprintf("%d draws retained, %d burn-in", nsamp, burnin)
  }
  cat(sprintf("SAEM conditional distribution  (%s)\n", provenance))
  shrink <- attr(x, "shrinkage")
  if (!is.null(shrink) && length(shrink) > 0L) {
    cat("Distribution-based eta-shrinkage:\n")
    for (nm in names(shrink)) {
      cat(sprintf("  %-12s %6.1f%%\n", nm, 100 * shrink[[nm]]))
    }
  }
  cat("\n")
  print.data.frame(as.data.frame(x))
  invisible(x)
}
