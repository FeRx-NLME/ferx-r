#' BIC variants for model selection
#'
#' The four BIC conventions an automated model search ranks candidates on, as
#' implemented by \code{ferx_core::model_selection::bic()} and by
#' \code{pharmpy.modeling.calculate_bic}. \code{fit$bic} is the classical
#' \code{OFV + p * log(n_obs)}, which is the \code{"fixed"} variant here;
#' Pharmpy's \code{modelsearch} and \code{iivsearch} rank on \code{"mixed"} and
#' \code{"iiv"} instead, because penalising a random-effects structure on the
#' observation count systematically favours the wrong model.
#'
#' The penalty added to the OFV is:
#'
#' \describe{
#'   \item{mixed}{Delattre et al. (2014):
#'     \code{n_random * log(n_subjects) + n_fixed * log(n_obs)}, where a free
#'     theta counts as random when it enters an individual parameter carrying
#'     an ETA or KAPPA, every free OMEGA / OMEGA_IOV element is random, and
#'     SIGMA is fixed-class unless the residual error itself carries an ETA.}
#'   \item{iiv}{\code{n_omega * log(n_subjects)} - free BSV OMEGA elements
#'     only, for comparing variability structures.}
#'   \item{random}{\code{n_parameters * log(n_subjects)}.}
#'   \item{fixed}{\code{n_parameters * log(n_obs)} - reproduces \code{fit$bic}.}
#' }
#'
#' The class tally comes from \code{fit$bic_inputs}, which the engine fills
#' from the same packed parameter mask it counts \code{n_parameters} from, so
#' the variants can be recomputed for a saved fit without the model or data in
#' hand.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} or
#'   \code{\link{ferx_load_fit}}.
#' @param type Which penalty convention to apply: \code{"mixed"} (default),
#'   \code{"fixed"}, \code{"iiv"} or \code{"random"}.
#' @return A single numeric: the BIC under \code{type}. \code{NA_real_} when
#'   the tally cannot support the penalty - a fit saved before the tally
#'   existed, a tally that disagrees with \code{n_parameters}, or a zero
#'   subject / record count under a convention that needs its logarithm.
#'   (ferx-core's \code{bic()} reports \code{NaN} in exactly those cases.) A
#'   fit with no free parameter returns its OFV under every convention.
#' @examples
#' \dontrun{
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data)
#' ferx_bic(fit)                  # mixed, the search default
#' ferx_bic(fit, "fixed")         # == fit$bic
#' }
#' @seealso \code{\link{check_strictness}} for the eligibility gates a search
#'   applies before ranking on one of these.
#' @family diagnostics
#' @export
ferx_bic <- function(fit, type = c("mixed", "fixed", "iiv", "random")) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit())",
         call. = FALSE)
  }
  type <- match.arg(type)

  inp <- .ferx_bic_inputs(fit)
  if (is.null(inp)) return(NA_real_)

  n_par <- suppressWarnings(as.integer(fit$n_parameters %||% NA_integer_))
  # A tally that does not add up to `n_parameters` is a tally from a different
  # parameterisation (or an all-zero one off a pre-#1177 bundle); ranking on it
  # would be a wrong penalty rather than a missing one.
  if (length(n_par) != 1L || is.na(n_par) || .ferx_bic_n_free(inp) != n_par) {
    return(NA_real_)
  }

  ln_subj <- .ferx_pos_log(fit$n_subjects)
  ln_obs  <- .ferx_pos_log(inp$n_obs)
  # A class with no members never needs its logarithm to exist, so a fit with
  # no free parameter is its OFV whatever n_subjects / n_obs say.
  term <- function(count, ln) if (count == 0L) 0 else count * ln

  penalty <- switch(type,
    fixed  = term(n_par, ln_obs),
    random = term(n_par, ln_subj),
    iiv    = term(inp$omega, ln_subj),
    mixed  = term(.ferx_bic_n_random_class(inp), ln_subj) +
             term(.ferx_bic_n_fixed_class(inp), ln_obs)
  )
  if (!is.finite(penalty)) return(NA_real_)

  ofv <- suppressWarnings(as.numeric(fit$ofv %||% NA_real_))
  if (length(ofv) != 1L) return(NA_real_)
  ofv + penalty
}

# log(n) for a count that must be positive to be usable; NA_real_ otherwise, so
# a term that needs it propagates NA into the penalty.
.ferx_pos_log <- function(n) {
  n <- suppressWarnings(as.numeric(n %||% NA_real_))
  if (length(n) != 1L || is.na(n) || n <= 0) return(NA_real_)
  log(n)
}

# Normalise fit$bic_inputs (an FFI list, or the same list read back from a
# .fitrx) to integer counts. NULL when the fit carries no tally at all.
.ferx_bic_inputs <- function(fit) {
  raw <- fit$bic_inputs
  if (is.null(raw) || !length(raw)) return(NULL)
  count <- function(nm) {
    v <- suppressWarnings(as.integer(raw[[nm]] %||% NA_integer_))
    if (length(v) != 1L || is.na(v) || v < 0L) NA_integer_ else v
  }
  out <- list(
    n_obs        = count("n_obs"),
    theta_random = count("theta_random"),
    theta_fixed  = count("theta_fixed"),
    omega        = count("omega"),
    kappa        = count("kappa"),
    sigma        = count("sigma"),
    sigma_random = isTRUE(as.logical(raw$sigma_random %||% FALSE))
  )
  if (anyNA(unlist(out[c("n_obs", "theta_random", "theta_fixed",
                         "omega", "kappa", "sigma")]))) {
    return(NULL)
  }
  out
}

.ferx_bic_n_free <- function(inp) {
  inp$theta_random + inp$theta_fixed + inp$omega + inp$kappa + inp$sigma
}

# Free parameters the mixed BIC penalises on log(n_subjects).
.ferx_bic_n_random_class <- function(inp) {
  inp$theta_random + inp$omega + inp$kappa + if (inp$sigma_random) inp$sigma else 0L
}

# Free parameters the mixed BIC penalises on log(n_obs).
.ferx_bic_n_fixed_class <- function(inp) {
  inp$theta_fixed + if (inp$sigma_random) 0L else inp$sigma
}
