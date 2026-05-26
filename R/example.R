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
#'   \item{bioavailability_ode}{ODE equivalent of \code{bioavailability} - same model
#'     and data via explicit depot/central ODEs; useful for comparing analytical vs ODE paths}
#'   \item{warfarin_dcm}{Two-compartment oral with a Deep Compartment Model
#'     covariate model: a small neural network (\code{[covariate_nn TYPICAL_PK]})
#'     reads WT + CRCL and outputs a multiplicative modulator on the baseline
#'     TVCL/TVV1/TVQ/TVV2/TVKA thetas (i.e. \code{CL = TVCL * TYPICAL_PK.CL *
#'     exp(ETA_CL)}), with lognormal IIV on top. Requires ferx-r built with the
#'     \code{nn} cargo feature.}
#'   \item{warfarin_scaled}{One-compartment oral (warfarin) demonstrating the
#'     \code{[scaling]} block (Form A). The dataset records AMT in micrograms
#'     (100 mg = 100,000 ug) while DV is in mg/L; \code{obs_scale = 1000}
#'     divides model predictions by 1000 so residuals are on the mg/L scale.
#'     Estimates should match the standard \code{warfarin} example.}
#'   \item{warfarin_iov_saem}{One-compartment oral (warfarin) with
#'     inter-occasion variability estimated by SAEM. Identical model to
#'     \code{warfarin_iov} but uses \code{method = saem}; demonstrates that
#'     SAEM fully supports IOV kappa parameters.}
#'   \item{warfarin_ss}{One-compartment oral (warfarin) demonstrating
#'     steady-state dosing. The dose row in the dataset has \code{SS = 1}
#'     and \code{II = 24} (once-daily interval). The engine resolves
#'     steady-state initial conditions automatically; no model changes are
#'     needed.}
#'   \item{transit_2cpt}{Two-compartment ODE model with 3-transit-compartment
#'     absorption, allometric scaling on CL/V1/V2 (reference WT = 70 kg), and
#'     a lag time on the first transit compartment.}
#'   \item{emax_pkpd}{Simultaneous PK/PD: oral 1-compartment PK plus an
#'     effect-compartment Emax PD readout, fit jointly with a SEPARATE
#'     residual error model per endpoint via a per-CMT \code{[error_model]}
#'     block (\code{CMT=2: DV ~ proportional(...)} for plasma,
#'     \code{CMT=3: DV ~ additive(...)} for the PD effect). Demonstrates
#'     multi-endpoint fitting; requires \code{gradient = fd}.}
#' }
#' Call \code{ferx_example()} with no arguments to list all available names.
#'
#' Each example has a matching end-to-end R script installed with the package.
#' List all scripts with:
#' \preformatted{
#' list.files(system.file("examples", package = "ferx"), pattern = "\\.R$")
#' }
#' Open any script directly in RStudio with \code{file.show()} or
#' \code{file.edit()}, for example:
#' \preformatted{
#' file.edit(system.file("examples", "ex1_warfarin.R",     package = "ferx"))
#' file.edit(system.file("examples", "ex_warfarin_dcm.R",  package = "ferx"))
#' }
#'
#' @return If \code{name} is NULL, a character vector of available examples.
#'   Otherwise, a list with components:
#'   \item{model}{Path to the .ferx model file}
#'   \item{data}{Path to the NONMEM-format CSV data file}
#'
#' @examples
#' ferx_example()                        # list all available example names
#' ex <- ferx_example("warfarin")
#' ex$model
#' ex$data
#'
#' # List all bundled example scripts:
#' list.files(system.file("examples", package = "ferx"), pattern = "\\.R$")
#'
#' \dontrun{
#' # Open a script in RStudio to read or run it:
#' file.edit(system.file("examples", "ex1_warfarin.R", package = "ferx"))
#'
#' # Deep Compartment Model (requires nn cargo feature):
#' dcm <- ferx_example("warfarin_dcm")
#' fit <- ferx_fit(dcm$model, dcm$data, method = "focei")
#' fit$neural_networks[[1]]$shape   # NN layer dimensions
#' file.edit(system.file("examples", "ex_warfarin_dcm.R", package = "ferx"))
#' }
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
