#' ferx: Nonlinear Mixed Effects Modeling with a Rust Backend
#'
#' ferx provides a user-facing R API for population PK / NLME modeling
#' that delegates to a high-performance Rust engine (ferx-core). Models are
#' written in a compact .ferx DSL (parameters, individual_parameters,
#' structural_model, odes, error_model, initial_values, fit_options,
#' optional scaling).
#'
#' @section Where to start:
#' \tabular{ll}{
#'   \strong{Goal}                                      \tab \strong{Function} \cr
#'   Get a worked example                               \tab \code{\link{ferx_example}} \cr
#'   Inspect dataset columns                            \tab \code{\link{ferx_get_columns}} \cr
#'   Create a model from a template                     \tab \code{\link{ferx_model}(template = ...)} \cr
#'   Inspect / edit an existing model                   \tab \code{\link{ferx_model_show}}, \code{\link{ferx_model_edit}}, \code{\link{ferx_model_validate}}, \code{\link{ferx_model_inspect}} \cr
#'   Sanity-check starting values                       \tab \code{\link{ferx_check_init}}, \code{\link{ferx_inits_from_nca}} \cr
#'   Fit a model (FOCE / FOCEI / GN / SAEM / IMP)       \tab \code{\link{ferx_fit}} \cr
#'   Inspect a fit                                      \tab \code{\link{summary.ferx_fit}}, \code{fit$estimates}, \code{fit$cor_matrix}, \code{\link{check_diagnostics}} \cr
#'   Optimizer trace                                    \tab \code{\link{ferx_trace}}, \code{\link{ferx_plot_trace}} \cr
#'   Covariate screen on EBE ETAs                       \tab \code{fit$eta_cov} \cr
#'   Uncertainty (SIR)                                  \tab \code{\link{ferx_sir}} \cr
#'   Population / individual predictions                \tab \code{\link{ferx_predict}} \cr
#'   Simulate (with or without parameter uncertainty)   \tab \code{\link{ferx_simulate}}, \code{\link{ferx_simulate_with_uncertainty}} \cr
#'   Save / restore a fit                               \tab \code{\link{ferx_save_fit}}, \code{\link{ferx_load_fit}}
#' }
#'
#' Functions are grouped on the reference index under \emph{fitting},
#' \emph{simulation}, \emph{model-editing}, \emph{diagnostics}, and
#' \emph{utilities}; each function's help page links to its siblings
#' via \code{See Also}.
#'
#' @keywords internal
#' @useDynLib ferx, .registration = TRUE
#' @importFrom graphics abline mtext par plot.new
#' @importFrom stats cor.test shapiro.test setNames
#' @importFrom utils read.csv
"_PACKAGE"
