# BIC variants for candidate ranking (ferx-core #1177).
#
# The engine records the free-parameter tally by Delattre class on every fit
# (`fit$bic_inputs`) and ships the three non-"fixed" conventions alongside it.
# The arithmetic is mirrored here so a fit that carries the tally but not the
# variants - a `.fitrx` bundle, which persists `bic_inputs` but not the derived
# values - is still rankable. Keep this in step with `ferx_core::bic`.

# `count * ln`, where a class with no members does not need its log to exist.
# NULL means "no log available", which propagates to NaN.
.ferx_bic_term <- function(count, ln) {
  if (count == 0L) return(0)
  if (is.null(ln)) return(NULL)
  count * ln
}

# One BIC variant from `fit$bic_inputs`. Mirrors `ferx_core::bic()`, NaN cases
# included: a tally that does not add up to `n_parameters` (an all-zero tally
# on a bundle saved before the field existed is the common case), or a zero
# subject / record count under a convention that needs its log.
.ferx_bic_variant <- function(fit, type) {
  inp <- fit$bic_inputs
  if (is.null(inp)) return(NaN)
  counts <- c("theta_random", "theta_fixed", "omega", "kappa", "sigma", "n_obs")
  if (!all(counts %in% names(inp))) return(NaN)
  n <- vapply(inp[counts], function(x) as.integer(x %||% 0L), integer(1L))
  n_par <- as.integer(fit$n_parameters %||% NA_integer_)
  n_free <- sum(n[c("theta_random", "theta_fixed", "omega", "kappa", "sigma")])
  if (is.na(n_par) || n_free != n_par) return(NaN)

  n_subjects <- as.integer(fit$n_subjects %||% 0L)
  ln_subj <- if (!is.na(n_subjects) && n_subjects > 0L) log(n_subjects) else NULL
  ln_obs <- if (n[["n_obs"]] > 0L) log(n[["n_obs"]]) else NULL

  sigma_random <- isTRUE(inp$sigma_random)
  n_random <- n[["theta_random"]] + n[["omega"]] + n[["kappa"]] +
    if (sigma_random) n[["sigma"]] else 0L
  n_fixed <- n[["theta_fixed"]] + if (sigma_random) 0L else n[["sigma"]]

  penalty <- switch(type,
    fixed  = .ferx_bic_term(n_par, ln_obs),
    random = .ferx_bic_term(n_par, ln_subj),
    iiv    = .ferx_bic_term(n[["omega"]], ln_subj),
    mixed  = {
      r <- .ferx_bic_term(n_random, ln_subj)
      f <- .ferx_bic_term(n_fixed, ln_obs)
      if (is.null(r) || is.null(f)) NULL else r + f
    },
    stop("Unknown BIC type: ", type)
  )
  if (is.null(penalty)) NaN else as.numeric(fit$ofv) + penalty
}

# Fill `bic_mixed` / `bic_iiv` / `bic_random` on a fit that does not carry
# them. `ferx_fit()` gets them from the engine; a fit rebuilt from a `.fitrx`
# bundle gets them here. Called from `.ferx_populate_derived_fields()`.
.ferx_populate_bic_variants <- function(result) {
  for (type in c("mixed", "iiv", "random")) {
    field <- paste0("bic_", type)
    if (is.null(result[[field]])) {
      result[[field]] <- .ferx_bic_variant(result, type)
    }
  }
  result
}

#' BIC of a fit under one of the four penalty conventions
#'
#' The four variants of \code{pharmpy.modeling.calculate_bic}, which differ in
#' which free parameters are penalised on \code{log(n_subjects)} and which on
#' \code{log(n_obs)}. \code{"mixed"} is the Delattre et al. (2014) BIC that
#' Pharmpy's \code{iivsearch} and \code{modelsearch} rank candidates on, and is
#' the default here for the same reason.
#'
#' \tabular{ll}{
#'   \code{"mixed"}  \tab \code{OFV + n_random * log(n_subjects) + n_fixed * log(n_obs)} \cr
#'   \code{"iiv"}    \tab \code{OFV + n_omega * log(n_subjects)} \cr
#'   \code{"random"} \tab \code{OFV + n_parameters * log(n_subjects)} \cr
#'   \code{"fixed"}  \tab \code{OFV + n_parameters * log(n_obs)} \cr
#' }
#'
#' Class membership follows Pharmpy: every free omega (BSV) and kappa (IOV)
#' element is random; a free theta is random when it enters an individual
#' parameter carrying an eta or kappa and fixed otherwise; sigma is fixed
#' unless the residual error itself carries an eta (\code{iiv_on_ruv}), which
#' moves it into the random class. The counts are recorded by the engine at fit
#' time and are readable as \code{fit$bic_inputs}.
#'
#' @param fit A \code{ferx_fit} object.
#' @param type One of \code{"mixed"} (default), \code{"iiv"}, \code{"random"},
#'   \code{"fixed"}.
#' @return A single numeric. \code{"fixed"} returns the fit's own \code{bic},
#'   which the engine computed at fit time and which is the same quantity.
#'   The other three are \code{NaN} for a fit whose parameter tally cannot
#'   support the penalty - in practice a \code{.fitrx} bundle saved before the
#'   tally was recorded, which loads with an all-zero \code{bic_inputs}.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' ferx_bic(fit)                  # Delattre (mixed) BIC
#' ferx_bic(fit, "iiv")           # variability-structure criterion
#' fit$bic_inputs                 # the counts behind the penalties
#' @references Delattre M, Lavielle M, Poursat MA (2014). A note on BIC in
#'   mixed-effects models. \emph{Electronic Journal of Statistics} 8(1):456-475.
#' @family utilities
#' @seealso \code{\link{check_strictness}} for the gates a candidate should
#'   pass before its criterion is trusted.
#' @export
ferx_bic <- function(fit, type = c("mixed", "iiv", "random", "fixed")) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).")
  }
  type <- match.arg(type)
  # "fixed" is the engine's own `bic`, and stays right even on a bundle whose
  # tally predates ferx-core #1177 - which the recomputation could not.
  if (identical(type, "fixed")) return(as.numeric(fit$bic))
  stored <- fit[[paste0("bic_", type)]]
  if (!is.null(stored) && length(stored) == 1L) return(as.numeric(stored))
  .ferx_bic_variant(fit, type)
}
