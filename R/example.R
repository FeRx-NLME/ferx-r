#' Get paths to example model and data files
#'
#' Returns file paths to bundled example models and datasets.
#' Called without arguments, lists available example names.
#'
#' @param name Name of the example (e.g. "warfarin"). If NULL, returns
#'   a character vector of available example names.
#'
#' @details
#' Available bundled examples:
#' \describe{
#'   \item{warfarin}{One-compartment oral PK (standard warfarin dataset)}
#'   \item{warfarin_bloq}{One-compartment oral with BLOQ observations (M3 method)}
#'   \item{warfarin_iov}{One-compartment oral with inter-occasion variability (kappa)}
#'   \item{warfarin_block_omega}{One-compartment oral with correlated random effects}
#'   \item{warfarin_saem}{One-compartment oral estimated with SAEM}
#'   \item{warfarin_additive_eta}{One-compartment oral with additive ETA on lag time}
#'   \item{warfarin_logit_f}{One-compartment oral with logit-normal bioavailability}
#'   \item{warfarin_if}{Two-compartment oral with conditional if/else covariate (WT)}
#'   \item{two_cpt_iv}{Two-compartment IV bolus}
#'   \item{two_cpt_oral_cov}{Two-compartment oral with continuous covariates (WT, CRCL)}
#'   \item{three_cpt_iv}{Three-compartment IV bolus}
#'   \item{mm_oral}{One-compartment oral with Michaelis-Menten elimination (ODE)}
#'   \item{warfarin_sde}{One-compartment oral as ODE with SDE process noise
#'     (\code{[diffusion]} block, EKF likelihood; \code{DIFF_CENTRAL} in theta)}
#'   \item{mm_multistart}{One-compartment oral with Michaelis-Menten elimination (ODE);
#'     starting values intentionally far from truth to demonstrate multi-start
#'     (\code{n_starts}, \code{start_sigma}, \code{multi_start_seed} in \code{settings})}
#'   \item{bioavailability}{One-compartment oral with logit-normal bioavailability F
#'     (analytical; THETA_F specified directly on the (0,1) scale)}
#'   \item{bioavailability_ode}{ODE equivalent of \code{bioavailability} — same model
#'     and data via explicit depot/central ODEs; useful for comparing analytical vs ODE paths}
#'   \item{warfarin_dcm}{Two-compartment oral with a Deep Compartment Model
#'     covariate model: a small neural network (\code{[covariate_nn TYPICAL_PK]})
#'     maps WT + CRCL directly to typical CL/V1/Q/V2/KA values, with lognormal
#'     IIV on the final PK params. Requires ferx-r built with the \code{nn}
#'     cargo feature. See \code{inst/examples/ex_warfarin_dcm.R}.}
#' }
#' Call \code{ferx_example()} with no arguments to list all available names.
#'
#' @return If \code{name} is NULL, a character vector of available examples.
#'   Otherwise, a list with components:
#'   \item{model}{Path to the .ferx model file}
#'   \item{data}{Path to the NONMEM-format CSV data file}
#'
#' @examples
#' ferx_example()
#' ex <- ferx_example("warfarin")
#' ex$model
#' ex$data
#'
#' @export
ferx_example <- function(name = NULL) {
  models_dir <- system.file("examples", "models", package = "ferx")
  if (models_dir == "") {
    stop("Example files not found. Is ferx installed?")
  }

  available <- tools::file_path_sans_ext(list.files(models_dir, pattern = "\\.ferx$"))

  if (is.null(name)) {
    return(available)
  }

  if (!name %in% available) {
    stop(
      "Example '", name, "' not found. Available examples: ",
      paste(available, collapse = ", ")
    )
  }

  list(
    model = system.file("examples", "models", paste0(name, ".ferx"), package = "ferx"),
    data = system.file("examples", "data", paste0(name, ".csv"), package = "ferx")
  )
}
