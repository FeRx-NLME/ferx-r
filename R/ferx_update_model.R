#' Write a fitted model with its final estimates as initial values
#'
#' Writes the model a fit was made with as a \code{.ferx} file whose initial
#' values are the fit's final estimates - NONMEM's \code{final.mod}, Pharmpy's
#' \code{update_inits}. Refitting from the written file starts at the fit's
#' optimum: its objective function value at iteration 0 is the fit's OFV.
#'
#' The edit is the engine's own (\code{ModelEdit::SeedInits} in ferx-core),
#' the one every search uses to write its \code{final.ferx}. Only the numbers
#' in the \code{[parameters]} block change:
#' \itemize{
#'   \item each \code{theta} init takes the fitted value; its bounds are kept;
#'   \item each \code{omega} / \code{kappa} takes the fitted variance, or its
#'     square root on a declaration marked \code{(sd)};
#'   \item each \code{sigma} takes the fitted value on the declaration's own
#'     scale (variance, or SD with \code{(sd)});
#'   \item a \code{block_omega} / \code{block_kappa} takes the fitted lower
#'     triangle, keeping its block structure;
#'   \item a \code{FIX}ed parameter is left as declared, as is a vector theta
#'     (\code{theta NAME[...]}) and any estimate that is not finite.
#' }
#' A model with a vector theta therefore refits from the vector's declared
#' values, not its fitted ones, and does not start at the fit's OFV.
#' Whole-line comments, formatting and every other block are kept as written.
#' A comment at the end of a declaration whose value is rewritten is not: the
#' engine writes that line afresh. Residual correlations of a
#' \code{block_sigma} are not rewritten.
#'
#' The model text is the one that was fitted (\code{fit$model_text}, when it
#' still hashes to \code{fit$model_hash}), whatever has happened to the file at
#' \code{fit$model_path} since. A fit that does not carry it falls back to that
#' file, with a warning of class \code{ferx_model_edited_warning} when the file
#' has been edited since the fit.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}} or
#'   \code{\link{ferx_load_fit}}.
#' @param path Where to write the model. \code{NULL} (default) writes nothing
#'   and returns the text instead.
#' @param overwrite Logical. When \code{FALSE} (default), an existing file at
#'   \code{path} is an error.
#' @return With \code{path}, the path written, invisibly. Without, the updated
#'   model text as a single character string.
#'
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' final <- ferx_update_model(fit, tempfile(fileext = ".ferx"))
#' cat(ferx_model_get_section(final, "parameters"), sep = "\n")
#'
#' @seealso \code{\link{ferx_save_fit}} (\code{write_model = TRUE} writes this
#'   file next to the bundle), \code{\link{ferx_fit}}
#' @family model-editing
#' @export
ferx_update_model <- function(fit, path = NULL, overwrite = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (returned by ferx_fit()).")
  }
  if (!is.null(path) &&
      (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))) {
    stop("`path` must be NULL or a single character path.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.")
  }
  if (!is.null(path) && file.exists(path) && !overwrite) {
    stop(path, " already exists. Use overwrite = TRUE to replace it.")
  }

  text <- .ferx_updated_model_text(fit)
  if (is.null(path)) return(text)
  writeBin(charToRaw(text), path)
  invisible(path)
}

# The fitted model's text with the fit's estimates written in as inits. Shared
# by `ferx_update_model()` and `ferx_save_fit(write_model = TRUE)`.
.ferx_updated_model_text <- function(fit) {
  model_text <- .ferx_fit_source_text(fit)
  omega <- fit$omega
  if (!is.matrix(omega) || nrow(omega) != ncol(omega)) {
    stop("`fit$omega` must be a square matrix.")
  }
  theta <- as.numeric(fit$theta)
  theta_names <- as.character(fit$theta_names %||% names(fit$theta) %||%
                                character())
  sigma <- as.numeric(fit$sigma)
  sigma_names <- as.character(fit$sigma_names %||% names(fit$sigma) %||%
                                character())
  eta_names <- as.character(fit$eta_names %||% rownames(omega) %||% character())
  iov <- .ferx_omega_iov_args(fit)
  kappa_names <- if (iov$omega_iov_dim > 0L) {
    as.character(fit$kappa_names %||% rownames(fit$omega_iov) %||% character())
  } else {
    character()
  }
  ferx_rust_update_model(
    model_text     = model_text,
    theta          = theta,
    theta_names    = theta_names,
    omega_flat     = as.numeric(t(omega)),  # row-major
    omega_dim      = as.integer(nrow(omega)),
    eta_names      = eta_names,
    sigma          = sigma,
    sigma_names    = sigma_names,
    omega_iov_flat = iov$omega_iov_flat,
    omega_iov_dim  = iov$omega_iov_dim,
    kappa_names    = kappa_names
  )
}

# The text the fit was made from: `fit$model_text` when it is the fitted one
# (see `.fitrx_fitted_model_text()`), otherwise the file at `fit$model_path`
# as it is now, warning when that file has been edited since the fit - the
# same order `ferx_save_fit()` writes the bundle's model.ferx in.
.ferx_fit_source_text <- function(fit) {
  txt <- .fitrx_fitted_model_text(fit)
  if (!is.null(txt)) return(txt)
  model_path <- fit$model_path
  if (is.character(model_path) && length(model_path) == 1L &&
      !is.na(model_path) && nzchar(model_path) && file.exists(model_path)) {
    .fitrx_warn_model_edited(model_path, fit$model_hash)
    txt <- rawToChar(readBin(model_path, "raw", n = file.size(model_path)))
    Encoding(txt) <- "UTF-8"
    return(txt)
  }
  if (is.character(fit$model_text) && length(fit$model_text) == 1L &&
      !is.na(fit$model_text)) {
    return(fit$model_text)
  }
  stop("The fit carries no model text (fit$model_text) and its model file ",
       "(fit$model_path) cannot be read.")
}
