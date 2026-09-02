#' Inspect the structure of a ferx model file
#'
#' Parses a \code{.ferx} file without fitting and prints a compact summary of
#' its model structure: PK model type, inter-individual variability (IIV),
#' inter-occasion variability (IOV), and residual error type. Useful for
#' verifying that a model file will be interpreted as expected before committing
#' to a potentially long estimation run.
#'
#' Alternatively, pass a \code{ferx_fit} object to display the structure that
#' was auto-derived during fitting (reads \code{fit$model_structure} directly,
#' so no file path is needed post-fit).
#'
#' @param path Path to a \code{.ferx} model file, \emph{or} a
#'   \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#'
#' @return A named list (invisibly) with fields: \code{theta_names}
#'   (character vector of population parameter names), \code{model_type}
#'   (short label such as \code{"1-cpt oral"}, \code{"ODE"}, or
#'   \code{"compartment-free"}, or \code{NULL} when not
#'   unambiguously detectable), \code{iiv} (omega names), \code{iov}
#'   (kappa names), \code{iov_weights} (per-kappa sample-size weight
#'   expressions from \code{kappa K ~ <var> weight = <expr>}, \code{NA} for
#'   an unweighted kappa and \code{character(0)} when no kappa is weighted),
#'   and \code{residual} (error type - one of
#'   \code{"proportional"}, \code{"additive"}, \code{"combined"}, or
#'   \code{"additive (log-transformed)"} for log-transform-both-sides
#'   models written as \code{log(DV) ~ additive(...)} or
#'   \code{DV ~ log_additive(...)}). For multi-endpoint (per-CMT) error
#'   models, \code{residual} is reported as a string of the form
#'   \code{"per-CMT (CMT2=proportional, CMT3=additive)"}.
#'
#' @section Model DSL features detected:
#' \code{ferx_model_inspect()} reflects the parser's view of a
#' \code{.ferx} file. Recently added DSL features that surface here
#' include multi-endpoint residual error (per-CMT \code{[error_model]}),
#' the \code{[scaling]} block for unit conversion / ODE readout (see
#' \code{NEWS.md} for Form A/B/C details), and inter-occasion variability
#' via \code{kappa_*} declarations. Steady-state dosing (\code{SS}/\code{II}
#' columns) is a data-side feature and is documented under
#' \code{\link{ferx_fit}}.
#'
#' A sample-size-weighted kappa - \code{kappa KAPPA_EMAX ~ 2.0 (sd) weight = NARM},
#' meaning \code{kappa_ik ~ N(0, Omega_IOV / N_ik)} - is annotated on the IOV
#' line as \code{KAPPA_EMAX (weight = NARM)} and returned in
#' \code{iov_weights}.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_model_inspect(ex$model)
#'
#' # Programmatic access
#' s <- ferx_model_inspect(ex$model)
#' s$theta_names  # c("TVCL", "TVV", "TVKA")
#' s$model_type   # "1-cpt oral"
#' s$iiv          # c("ETA_CL", "ETA_V", "ETA_KA")
#' s$residual     # "proportional"
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model_edit}},
#'   \code{\link{ferx_fit}}
#' @family model-editing
#' @export
ferx_model_inspect <- function(path) {
  if (inherits(path, "ferx_fit")) {
    s <- path$model_structure
    if (is.null(s)) stop("No model_structure found on this ferx_fit object.")
    label <- if (!is.null(path$model_name) && nzchar(path$model_name))
      paste0(path$model_name, ".ferx") else "ferx_fit"
    cat("Model structure (", label, ")\n", sep = "")
    .ferx_print_structure(s)
    return(invisible(s))
  }
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")

  s <- .ferx_parse_structure(path)

  cat("Model structure (", basename(path), ")\n", sep = "")
  .ferx_print_structure(s)

  invisible(s)
}
