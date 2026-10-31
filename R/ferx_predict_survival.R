#' Survival-function predictions for time-to-event endpoints
#'
#' For a model with one or more \code{[event_model]} (TTE) endpoints, evaluates
#' the survival function \eqn{S(t)}, cumulative hazard \eqn{H(t)} and hazard
#' \eqn{h(t)} on a user-supplied time grid, for every subject and every TTE CMT.
#' Population-level (\eqn{\eta = 0}) unless a \code{fit} is supplied, in which
#' case its fitted \code{theta} is used. With \code{type = "individual"} each
#' subject's curves are evaluated at its own empirical Bayes estimates from
#' \code{fit} instead. Each subject also carries its median and mean survival
#' time (analytic where available).
#'
#' For a model with multiple TTE CMTs (competing risks), each row also carries the
#' cause-specific cumulative incidence \code{cif} (\eqn{F_k(t)}) and the all-cause
#' survival \code{survival_all} (\eqn{S_\mathrm{all}(t)}), with
#' \eqn{\sum_k F_k(t) + S_\mathrm{all}(t) = 1}. For a single endpoint
#' \code{cif} \eqn{= 1 -} \code{survival} and \code{survival_all} \eqn{=}
#' \code{survival}.
#'
#' @param model Path to a .ferx model file containing at least one
#'   \code{[event_model]} block.
#' @param data Path to a NONMEM-format CSV, or a data.frame (see \code{\link{ferx_fit}}). When omitted, the data
#'   \code{fit} was fitted to (when \code{fit} is given and its data can still
#'   be found), otherwise the model file's \code{[data]} block
#'   (\code{path = ...}).
#' @param times Numeric vector of times at which to evaluate \eqn{S(t)},
#'   \eqn{H(t)}, \eqn{h(t)}.
#' @param fit Optional \code{ferx_fit} result. When provided, predictions use
#'   \code{fit$theta} instead of the model file's initial theta.
#' @param type \code{"population"} (the default) evaluates every subject at
#'   \eqn{\eta = 0}. \code{"individual"} evaluates each subject at its own
#'   EBEs, \code{fit$ebe_etas} (requires \code{fit}): the individual curves
#'   \eqn{\hat S_i(t)}, \eqn{\hat H_i(t)}, \eqn{\hat h_i(t)} and CIF. Every
#'   subject in \code{data} needs a row in \code{fit$ebe_etas}, matched by
#'   \code{ID}, so \code{data} is the fitted data or a subset of it. IOV
#'   kappas are taken as 0, which is exact since no hazard can read one. The
#'   EBEs are point estimates, so a subject with few observations (high
#'   shrinkage) gets a curve close to the population one. A \code{[mixture]}
#'   model is refused (\code{E_SUPPLIED_ETAS}). For an analytic hazard and a
#'   single-event record without delayed entry, \code{cum_hazard} at the
#'   record's time is the \eqn{H(T_i)} the fit's likelihood evaluated
#'   (ferx-core #1885).
#'
#' @return A data.frame with columns: ID, CMT, TIME, survival, cum_hazard,
#'   hazard, cif, survival_all, median_survival, mean_survival (one row per
#'   subject \eqn{\times} TTE CMT \eqn{\times} time).
#'
#'   The frame carries a \code{simulation_warnings} attribute (a character
#'   vector, empty for a clean run), re-emitted as a single R warning when
#'   non-empty, as \code{\link{ferx_predict}()} does: the engine's parse
#'   warnings (an unreferenced theta, say), its data-reader findings (a
#'   covariate missing for some subjects, a dose that never lands; those a
#'   compartment-free TTE model makes moot are withheld), the model/data
#'   checks, and the ODE solver diagnostics of this pass for a joint PK-TTE
#'   hazard. Ahead of them comes any occasion note: a model-file
#'   \code{iov_occasion} rule that replaces the \code{iov_column} labels, or,
#'   with \code{fit}, a fit that records no occasion rule and so takes the
#'   model file's.
#'
#' @inheritSection ferx_simulate Errors raised by the engine
#' @inheritSection ferx_simulate Theta level blocks
#'
#' @examples
#' \dontrun{
#' # `model` must contain at least one [event_model] (TTE) block.
#' preds <- ferx_predict_survival(
#'   "tte_model.ferx", "tte_data.csv",
#'   times = seq(0, 24, by = 2)
#' )
#' head(preds)
#'
#' # Each subject's curve at its fitted EBEs.
#' fit <- ferx_fit("tte_model.ferx", "tte_data.csv")
#' indiv <- ferx_predict_survival(
#'   "tte_model.ferx", times = seq(0, 24, by = 2),
#'   fit = fit, type = "individual"
#' )
#' }
#'
#' @family simulation
#' @export
ferx_predict_survival <- function(model, data = NULL, times, fit = NULL,
                                  type = c("population", "individual")) {
  type <- match.arg(type)
  if (type == "individual" && is.null(fit)) {
    stop("ferx_predict_survival: type = \"individual\" needs `fit`, whose ",
         "EBEs (`fit$ebe_etas`) give each subject's etas.", call. = FALSE)
  }
  if (is.null(data)) data <- .ferx_default_data(model, fit)
  if (is.null(data)) {
    stop(
      "No data supplied. Pass `data`, or add a `[data]` block ",
      "(`path = ...`) to the model file."
    )
  }
  data <- .ferx_data_arg(data)
  stopifnot(
    file.exists(model),
    file.exists(data),
    is.numeric(times),
    length(times) > 0,
    all(is.finite(times))
  )
  times <- as.numeric(times)

  # A refusal is an R error, classed like a refused `ferx_fit()` (#385).
  if (is.null(fit)) {
    res <- .ferx_engine_call(
      ferx_rust_predict_survival(
        model_path = normalizePath(model),
        data_path = normalizePath(data),
        times = times
      ),
      model, data
    )
    return(.ferx_surface_sim_warnings(res, "ferx_predict_survival"))
  }

  fit_pieces <- validate_fit_for_params(fit)
  etas <- .ferx_survival_etas(fit, type)
  # A note on how the fit's occasions were derived (#523) comes back on the
  # simulate paths' `simulation_warnings` channel.
  res <- .ferx_engine_call(
    ferx_rust_predict_survival_from_fit(
      model_path = normalizePath(model),
      data_path = normalizePath(data),
      times = times,
      theta = fit_pieces$theta,
      omega_flat = fit_pieces$omega_flat,
      omega_dim = fit_pieces$omega_dim,
      sigma = fit_pieces$sigma,
      omega_iov_flat = fit_pieces$omega_iov_flat,
      omega_iov_dim = fit_pieces$omega_iov_dim,
      residual_rho = fit_pieces$residual_rho,
      residual_rho_pairs = fit_pieces$residual_rho_pairs,
      mixture_overrides = fit_pieces$mixture_overrides,
      fit_bindings = fit_pieces$fit_bindings,
      iov_occasion = fit_pieces$iov_occasion,
      individual = type == "individual",
      eta_ids = etas$ids,
      eta_names = etas$names,
      eta_flat = etas$flat
    ),
    model, data
  )
  .ferx_surface_sim_warnings(res, "ferx_predict_survival")
}

# Each fitted subject's EBEs for `type = "individual"` (ferx-core #1885): IDs
# verbatim, eta names in the model's order (the engine checks them), values
# row-major. A fit without random effects carries no `ebe_etas`; its subjects
# then get empty rows, and the individual curves are the population ones.
.ferx_survival_etas <- function(fit, type) {
  if (type != "individual") {
    return(list(ids = character(0), names = character(0), flat = numeric(0)))
  }
  ids <- .ferx_fit_subject_id_set(fit, "ferx_predict_survival")$ids
  ebes <- fit$ebe_etas
  if (is.null(ebes) || nrow(ebes) == 0L) {
    return(list(ids = ids, names = character(0), flat = numeric(0)))
  }
  eta_cols <- setdiff(names(ebes), "ID")
  eta_mat <- as.matrix(ebes[, eta_cols, drop = FALSE])
  storage.mode(eta_mat) <- "double"
  list(ids = ids, names = eta_cols, flat = as.numeric(t(eta_mat)))
}
