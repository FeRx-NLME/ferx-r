#' Prepare a FREM (Full Random Effects Model) dataset and model
#'
#' Transforms a base model and dataset into a FREM model that treats covariates
#' as additional dependent variables. The extended omega matrix captures
#' covariate-parameter relationships implicitly, enabling covariate screening
#' without stepwise search.
#'
#' The covariates folded into the FREM model - and whether each is continuous or
#' categorical - come from the model's \code{[covariates]} block, which is the
#' single source of truth (and is required).
#'
#' @param model Path to a \code{.ferx} model file, or a \code{\link{ferx_model}}
#'   object. The model must declare its covariates in a \code{[covariates]}
#'   block, each tagged continuous or categorical.
#' @param data Path to a NONMEM-format CSV file containing the covariate
#'   columns.
#' @param covariates Optional character vector used as a \emph{subset filter}
#'   over the covariates declared in the model's \code{[covariates]} block. When
#'   \code{NULL} (the default), \strong{all} declared covariates are included.
#'   When supplied, only the named covariates are included - useful when you do
#'   not want every declared covariate in the FREM model. Each name must be
#'   declared in the block; an undeclared name is an error. This argument cannot
#'   introduce covariates the model has not declared, nor change their
#'   continuous/categorical kind.
#' @param output_dir Directory for the output model and data files. Defaults to
#'   the directory containing \code{model}.
#' @param output_model Optional explicit path for the output \code{.ferx} model
#'   file. When \code{NULL} (default), the file is written to
#'   \code{<output_dir>/<stem>_frem.ferx}.
#' @param output_data Optional explicit path for the output CSV data file. When
#'   \code{NULL} (default), the file is written to
#'   \code{<output_dir>/<stem>_frem_data.csv}.
#' @param fit Optional \code{\link{ferx_fit}} result from fitting \code{model}
#'   (the base model, before FREM conversion). When supplied, its theta and
#'   omega estimates seed the generated FREM model's PK theta inits and PK-PK
#'   omega block, so a subsequent fit of the FREM model warm-starts from
#'   converged parameters instead of \code{model}'s declared inits. Matching is
#'   by name (theta/eta names in \code{fit}); a name not found in \code{model}
#'   keeps its declared init. \code{NULL} (default) leaves every init at the
#'   value declared in \code{model}.
#'
#' @return A \code{\link{ferx_model}} object pointing at the generated FREM
#'   model and dataset, so it composes directly with \code{\link{ferx_fit}} and
#'   the other model helpers. The fixed covariate thetas, covariate omega
#'   initial values, and FREMTYPE mapping are written into the generated model
#'   and data files (and the fit's omega is labelled by eta name), so no
#'   separate metadata object is returned.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Reads the base model and dataset.
#'   \item Adds pseudo-observation rows for each covariate (one per subject),
#'     with \code{DV} set to the covariate value and a \code{FREMTYPE} column
#'     distinguishing covariate rows from PK observations.
#'   \item Generates a new \code{.ferx} model file with an extended omega block
#'     that covers both the original random effects and covariate random effects.
#'   \item Returns a \code{ferx_model} referencing the generated files.
#' }
#'
#' After calling \code{ferx_model_to_frem()}, fit the returned model directly:
#' \code{\link{ferx_fit}(frem)}.
#'
#' @export
ferx_model_to_frem <- function(model,
                         data = NULL,
                         covariates = NULL,
                         output_dir = NULL,
                         output_model = NULL,
                         output_data = NULL,
                         fit = NULL) {

  # --- resolve model path ---
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model_path <- model$model
  } else {
    model_path <- model
  }

  # --- validate inputs ---
  if (!is.character(model_path) || length(model_path) != 1L || !nzchar(model_path)) {
    stop("`model` must be a single non-empty path to a .ferx file or a ferx_model object.")
  }
  if (!file.exists(model_path)) {
    stop("Model file not found: ", model_path)
  }
  if (is.null(data) || !is.character(data) || length(data) != 1L || !nzchar(data)) {
    stop("`data` must be a single non-empty path to a CSV file.")
  }
  if (!file.exists(data)) {
    stop("Data file not found: ", data)
  }
  # `covariates` is an optional subset filter over the model's [covariates]
  # block. NULL (or empty) means "use every declared covariate"; the backend
  # validates that any names given are actually declared.
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("`covariates` must be a character vector (a subset of the model's [covariates] block) or NULL.")
  }
  if (!is.null(fit) && !inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit()) or NULL.")
  }

  # --- extract prior-fit init values (issue #239) ---
  # `fit` is optional; when absent every field below stays empty and the
  # generated FREM model keeps the base model's declared inits.
  if (!is.null(fit)) {
    fit_theta_names <- names(fit$theta) %||% character(0)
    fit_theta_values <- as.numeric(fit$theta)
    fit_eta_names <- if (!is.null(fit$omega)) rownames(fit$omega) else NULL
    fit_eta_names <- fit_eta_names %||% character(0)
    fit_omega_flat <- if (!is.null(fit$omega)) as.numeric(t(fit$omega)) else numeric(0) # row-major
    fit_omega_dim <- if (!is.null(fit$omega)) nrow(fit$omega) else 0L
  } else {
    fit_theta_names <- character(0)
    fit_theta_values <- numeric(0)
    fit_eta_names <- character(0)
    fit_omega_flat <- numeric(0)
    fit_omega_dim <- 0L
  }

  # --- resolve output paths ---
  model_path <- normalizePath(model_path, mustWork = TRUE)
  data <- normalizePath(data, mustWork = TRUE)

  if (is.null(output_dir)) {
    output_dir <- dirname(model_path)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  out_model_path <- if (!is.null(output_model)) {
    as.character(output_model)
  } else {
    ""
  }
  out_data_path <- if (!is.null(output_data)) {
    as.character(output_data)
  } else {
    ""
  }

  # --- call Rust backend ---
  raw <- ferx_rust_prepare_frem(
    model_path  = model_path,
    data_path   = data,
    covariates  = if (!is.null(covariates)) covariates else character(0),
    # Kinds (continuous/categorical) come from the model's [covariates] block;
    # the low-level override is intentionally left empty here.
    categorical_covariates = character(0),
    output_model_path = out_model_path,
    output_data_path  = out_data_path,
    fit_theta_names   = fit_theta_names,
    fit_theta_values  = fit_theta_values,
    fit_eta_names     = fit_eta_names,
    fit_omega_flat    = fit_omega_flat,
    fit_omega_dim     = fit_omega_dim
  )

  # Surface conversion-time advisories (e.g. estimated parameters without a
  # random effect, which IMP/IMPMAP estimate poorly - add an ETA before fitting).
  for (w in raw$warnings) {
    warning(w, call. = FALSE)
  }

  # Return a plain ferx_model referencing the generated FREM files, so the
  # result composes with ferx_fit() and the other model helpers. The covariate
  # statistics / FREMTYPE mapping the backend also computes are baked into the
  # generated model and data files, so they aren't surfaced as a separate object.
  ferx_model(data = raw$data_path, model = raw$model_path)
}
