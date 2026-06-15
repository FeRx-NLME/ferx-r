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
#'   \item{one_cpt_iv}{One-compartment IV (analytical \code{one_cpt_iv})}
#'   \item{two_cpt_iv}{Two-compartment IV bolus}
#'   \item{two_cpt_oral_cov}{Two-compartment oral with continuous covariates (WT, CRCL)}
#'   \item{three_cpt_iv}{Three-compartment IV bolus}
#'   \item{three_cpt_oral}{Three-compartment oral (analytical \code{three_cpt_oral})}
#'   \item{one_cpt_iv_ode, warfarin_ode, two_cpt_iv_ode, two_cpt_oral_cov_ode,
#'     three_cpt_iv_ode, three_cpt_oral_ode}{ODE-form transcriptions of the
#'     standard analytical models (\code{one_cpt_iv}, \code{warfarin} =
#'     one-cpt oral, \code{two_cpt_iv}, \code{two_cpt_oral_cov},
#'     \code{three_cpt_iv}, \code{three_cpt_oral}). Each shares its analytical
#'     counterpart's dataset and parameters; states are amounts and the
#'     observed concentration is read out with \code{[scaling] obs_scale = V}
#'     (or \code{V1}). The two forms give identical predictions - verified by
#'     \code{test-ode-analytical-equivalence.R}.}
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
#'   \item{warfarin_ltbs}{One-compartment oral (warfarin) with
#'     log-transform-both-sides (LTBS) residual error written as
#'     \code{log(DV) ~ additive(SIGMA_LOG)}. The engine log-transforms
#'     observations and predictions and fits additive (normal) error on
#'     the log scale; the residual is reported as
#'     \code{"additive (log-transformed)"}.}
#'   \item{warfarin_ss}{One-compartment oral (warfarin) demonstrating
#'     steady-state dosing. The dose row in the dataset has \code{SS = 1}
#'     and \code{II = 24} (once-daily interval). The engine resolves
#'     steady-state initial conditions automatically; no model changes are
#'     needed.}
#'   \item{transit_2cpt}{Two-compartment ODE model with 3-transit-compartment
#'     absorption, allometric scaling on CL/V1/V2 (reference WT = 70 kg).
#'     Transit compartments provide the absorption delay without a hard lag time.}
#'   \item{two_cpt_oral_cov_ode_template}{Two-compartment oral PK with covariates
#'     written with \code{ode_template two_cpt_oral(...)} in
#'     \code{[structural_model]}: ferx generates the standard disposition ODE
#'     (states, micro-constant RHS, and \code{obs_scale = V1}) instead of you
#'     hand-writing it. Identical predictions to the analytical
#'     \code{two_cpt_oral_cov} and the hand-ODE \code{two_cpt_oral_cov_ode}
#'     forms. Re-declaring a \code{d/dt(X)} in \code{[odes]} overrides only that
#'     compartment -- the standard way to attach a built-in absorption input
#'     such as \code{transit(...)}. (Requires ferx-core with
#'     FeRx-NLME/ferx-core#363.)}
#'   \item{warfarin_ode_lagtime}{One-compartment oral (warfarin) as an ODE model
#'     demonstrating the \code{LAGTIME} keyword in \code{[individual_parameters]}.
#'     Assigning to the reserved name \code{LAGTIME} delays every dose event for
#'     that subject by the individual lag time; no change to the ODE equations is
#'     needed. IIV on \code{LAGTIME} is log-normal; individual lag times are
#'     recovered in \code{fit$individual_estimates}. Contrast with the
#'     \code{lagtime=} argument used in analytical PK solvers
#'     (\code{warfarin_additive_eta}).}
#'   \item{emax_pkpd}{Simultaneous PK/PD: oral 1-compartment PK plus an
#'     effect-compartment Emax PD readout, fit jointly with a SEPARATE
#'     residual error model per endpoint via a per-CMT \code{[error_model]}
#'     block (\code{CMT=2: DV ~ proportional(...)} for plasma,
#'     \code{CMT=3: DV ~ additive(...)} for the PD effect). Demonstrates
#'     multi-endpoint fitting; requires \code{gradient = fd}.}
#'   \item{warfarin_derived}{One-compartment oral (warfarin) with a
#'     \code{[derived]} block demonstrating per-row expressions (KE, T_HALF),
#'     aggregates (CMAX, TMAX, CMAX_D1), numeric integration (AUC_0_72,
#'     AUC_TAU, AUC_DV_72), and an \code{[output]} block echoing CL, V, KA.
#'     Uses TAFD and TAD as built-in context variables.}
#'   \item{warfarin_derived_pkpd}{One-compartment oral (warfarin) with PD
#'     effect and time-above-MEC derived quantities. Demonstrates TAD-gated
#'     Cmax and \code{integral(1.0, IPRED > MEC, window=24)} for time above
#'     minimum effective concentration.}
#'   \item{two_cpt_oral_derived}{Two-compartment oral PK with micro-rate
#'     constants and day-specific Cmax computed in \code{[derived]}.
#'     Demonstrates \code{integral(IPRED, from=0, to=168)} and periodic AUC.}
#'   \item{warfarin_addl}{One-compartment oral (warfarin) with ADDL-column
#'     dosing. A single EVID=1 row with ADDL=6, II=24 is expanded to seven
#'     daily doses at read time; TAD resets at each expanded dose automatically.
#'     Demonstrates daily trough accumulation over seven doses and CMIN_TAU
#'     using \code{TAD < 1e-10} for floating-point-safe trough detection.}
#'   \item{warfarin_ode_time}{One-compartment oral (warfarin) as an ODE model
#'     demonstrating \code{TIME}, \code{T}, \code{TAFD}, and \code{TAD} inside
#'     \code{[odes]} RHS expressions; accumulator compartments for AUC and
#'     time-above-threshold.}
#'   \item{warfarin_data_selection}{One-compartment oral (warfarin) demonstrating
#'     the \code{[data_selection]} block. Observations below a surrogate LLOQ
#'     (\code{DV < 1.0} mg/L) are excluded at read time without modifying the
#'     CSV, equivalent to NONMEM \code{$DATA IGNORE=}. Exclusion counts are reported
#'     in \code{fit$exclusions} and echoed by \code{print.ferx_fit()} and
#'     \code{ferx_runlog()}.}
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
#' @family utilities
#' @export
ferx_example <- function(name = NULL) {
  models_dir <- system.file("examples", "models", package = "ferx")
  if (models_dir == "") {
    stop("Example files not found. Is ferx installed?")
  }

  # Some example model files share a dataset with an existing example.
  # Maps example name -> data file name (without .csv extension).
  .data_aliases <- list(
    warfarin_derived        = "warfarin",
    warfarin_derived_pkpd   = "warfarin",
    warfarin_ode_time       = "warfarin",

    warfarin_data_selection = "warfarin",
    two_cpt_oral_derived    = "two_cpt_oral_cov",

    # ODE-form siblings of the standard analytical models share their
    # analytical counterpart's dataset, so the two can be verified to give
    # identical predictions (test-ode-analytical-equivalence.R).
    warfarin_ode            = "warfarin",
    one_cpt_iv_ode          = "one_cpt_iv",
    two_cpt_iv_ode          = "two_cpt_iv",
    two_cpt_oral_cov_ode    = "two_cpt_oral_cov",
    three_cpt_iv_ode        = "three_cpt_iv",
    three_cpt_oral_ode      = "three_cpt_oral",

    # `ode_template`-generated disposition; shares the analytical dataset so it
    # can be verified identical to both the analytical and hand-ODE forms.
    two_cpt_oral_cov_ode_template = "two_cpt_oral_cov"
  )

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

  data_name <- .data_aliases[[name]] %||% name
  list(
    model = system.file("examples", "models", paste0(name,      ".ferx"), package = "ferx"),
    data  = system.file("examples", "data",   paste0(data_name, ".csv"),  package = "ferx")
  )
}

#' Inspect column headers of a NONMEM data file
#'
#' Reads only the header row of a NONMEM-format CSV file and prints the column
#' names grouped into required NONMEM columns, optional NONMEM columns, and
#' covariates / user-defined columns. Only the first line is read, so this is
#' fast even on large datasets.
#'
#' @param data One of:
#'   \itemize{
#'     \item A character string path to a NONMEM CSV file.
#'     \item A \code{ferx_fit} object -- uses \code{fit$data_path}.
#'     \item The list returned by \code{\link{ferx_example}} -- uses
#'       \code{x$data}.
#'   }
#'
#' @return A character vector of column names, returned invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_columns(ex)        # pass the ferx_example() list directly
#' ferx_columns(ex$data)   # or pass the path
#'
#' \dontrun{
#' fit <- ferx_fit(ex$model, ex$data)
#' ferx_columns(fit)       # inspect the data used for a fit
#' }
#'
#' @family utilities
#' @export
ferx_columns <- function(data) {
  path <- .ferx_resolve_data_path(data)

  if (!file.exists(path)) {
    stop("Data file not found: '", path, "'")
  }

  con <- tryCatch(
    file(path, open = "r"),
    error = function(e) stop("Cannot open '", path, "': ", conditionMessage(e))
  )
  on.exit(close(con), add = TRUE)

  header_line <- tryCatch(
    readLines(con, n = 1L, warn = FALSE),
    error = function(e) stop("Failed to read '", path, "': ", conditionMessage(e))
  )

  if (length(header_line) == 0L || !nzchar(trimws(header_line))) {
    stop("Data file '", path, "' appears to be empty.")
  }

  cols <- tryCatch(
    names(utils::read.csv(
      text = header_line, nrows = 0L,
      stringsAsFactors = FALSE, check.names = FALSE
    )),
    error = function(e) stop(
      "Could not parse column header of '", path, "': ", conditionMessage(e)
    )
  )

  if (length(cols) == 0L) {
    stop("No columns found in header of '", path, "'.")
  }

  nonmem_required <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT")
  nonmem_optional <- c("RATE", "MDV", "II", "SS", "CENS", "OCC")
  nonmem_all      <- c(nonmem_required, nonmem_optional)

  cols_upper <- toupper(cols)
  idx_req <- which(cols_upper %in% nonmem_required)
  idx_opt <- which(cols_upper %in% nonmem_optional)
  idx_cov <- which(!cols_upper %in% nonmem_all)

  n_total  <- length(cols)
  lbl_w    <- 22L
  indent   <- strrep(" ", lbl_w + 1L)
  con_w    <- getOption("width", 80L)

  fmt_group <- function(label, indices) {
    if (length(indices) == 0L) return(invisible(NULL))
    items   <- sprintf("[%d] %s", indices, cols[indices])
    lbl_str <- formatC(label, width = -lbl_w)
    cat(" ", lbl_str, sep = "")
    cur_pos <- lbl_w + 1L
    for (i in seq_along(items)) {
      sep   <- if (i < length(items)) "  " else ""
      token <- paste0(items[i], sep)
      if (cur_pos > lbl_w + 1L && cur_pos + nchar(items[i]) > con_w) {
        cat("\n", indent, sep = "")
        cur_pos <- lbl_w + 1L
      }
      cat(token)
      cur_pos <- cur_pos + nchar(token)
    }
    cat("\n")
  }

  cat(sprintf("Data file: %s\n", path))
  cat(sprintf("%d column%s\n\n", n_total, if (n_total == 1L) "" else "s"))
  fmt_group("Required NONMEM:", idx_req)
  fmt_group("Optional NONMEM:", idx_opt)
  fmt_group("Covariates / other:", idx_cov)
  cat("\n")

  invisible(cols)
}

.ferx_resolve_data_path <- function(data) {
  if (inherits(data, "ferx_fit")) {
    p <- data$data_path
    if (is.null(p) || length(p) != 1L || is.na(p) || !nzchar(p)) {
      stop("The ferx_fit object does not have a valid data_path.")
    }
    return(p)
  }
  if (is.list(data) && !is.null(names(data)) && "data" %in% names(data)) {
    p <- data[["data"]]
    if (!is.character(p) || length(p) != 1L || is.na(p) || !nzchar(p)) {
      stop("The list's '$data' element is not a valid file path string.")
    }
    return(p)
  }
  if (is.character(data) && length(data) == 1L && !is.na(data) && nzchar(data)) {
    return(data)
  }
  stop(
    "`data` must be a file path (character), a ferx_fit object, ",
    "or a ferx_example() list."
  )
}
