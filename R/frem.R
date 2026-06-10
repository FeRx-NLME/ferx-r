#' Prepare a FREM (Full Random Effects Model) dataset and model
#'
#' Transforms a base model and dataset into a FREM model that treats covariates
#' as additional dependent variables. The extended omega matrix captures
#' covariate-parameter relationships implicitly, enabling covariate screening
#' without stepwise search.
#'
#' @param model Path to a \code{.ferx} model file, or a \code{\link{ferx_model}}
#'   object.
#' @param data Path to a NONMEM-format CSV file containing the covariate
#'   columns listed in \code{covariates}.
#' @param covariates Character vector of covariate column names to include in
#'   the FREM transformation.
#' @param categorical Character vector of categorical covariate names (must be
#'   a subset of \code{covariates}). These are binarized (one-hot encoded)
#'   before the FREM transformation. Default \code{NULL} (no categoricals).
#' @param output_dir Directory for the output model and data files. Defaults to
#'   the directory containing \code{model}.
#' @param output_model Optional explicit path for the output \code{.ferx} model
#'   file. When \code{NULL} (default), the file is written to
#'   \code{<output_dir>/<stem>_frem.ferx}.
#' @param output_data Optional explicit path for the output CSV data file. When
#'   \code{NULL} (default), the file is written to
#'   \code{<output_dir>/<stem>_frem_data.csv}.
#'
#' @return A list with class \code{"ferx_frem"} containing:
#'   \describe{
#'     \item{model_path}{Path to the generated FREM \code{.ferx} model file.}
#'     \item{data_path}{Path to the generated FREM dataset (CSV).}
#'     \item{covariate_means}{Named numeric vector of population covariate means.}
#'     \item{covariate_variances}{Named numeric vector of population covariate variances.}
#'     \item{fremtype_map}{Named integer vector mapping covariate names to FREMTYPE values.}
#'     \item{n_total_etas}{Total number of etas in the FREM model (base + covariate).}
#'   }
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
#'   \item Returns metadata needed to interpret the FREM fit.
#' }
#'
#' After calling \code{ferx_to_frem()}, fit the returned model with
#' \code{\link{ferx_fit}(result$model_path, result$data_path)}.
#'
#' @export
ferx_to_frem <- function(model,
                         data = NULL,
                         covariates,
                         categorical = NULL,
                         output_dir = NULL,
                         output_model = NULL,
                         output_data = NULL) {

  # --- resolve model path ---
  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model_path <- model$path
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
  if (!is.character(covariates) || length(covariates) < 1L) {
    stop("`covariates` must be a character vector with at least one covariate name.")
  }
  if (!is.null(categorical)) {
    if (!is.character(categorical)) {
      stop("`categorical` must be a character vector or NULL.")
    }
    unknown <- setdiff(categorical, covariates)
    if (length(unknown) > 0L) {
      stop(
        "`categorical` contains names not in `covariates`: ",
        paste(unknown, collapse = ", ")
      )
    }
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
    covariates  = covariates,
    categorical_covariates = if (!is.null(categorical)) categorical else character(0),
    output_model_path = out_model_path,
    output_data_path  = out_data_path
  )

  # --- assemble result ---
  cov_means <- stats::setNames(raw$covariate_mean_values, raw$covariate_mean_names)
  cov_vars  <- stats::setNames(raw$covariate_var_values, raw$covariate_var_names)
  ft_map    <- stats::setNames(raw$fremtype_values, raw$fremtype_names)

  structure(
    list(
      model_path          = raw$model_path,
      data_path           = raw$data_path,
      covariate_means     = cov_means,
      covariate_variances = cov_vars,
      fremtype_map        = ft_map,
      n_total_etas        = raw$n_total_etas
    ),
    class = "ferx_frem"
  )
}
