#' Export a fitted deep compartment model as a plain population PK model
#'
#' Writes a fitted DCM - a model with one or more \code{[covariate_nn]} blocks -
#' out as an ordinary \code{.ferx} model: each network becomes its forward pass,
#' spelled as plain equations at the top of \code{[individual_parameters]} at
#' the fitted weights, and the block itself is removed. The exported model
#' predicts exactly what the fitted network does, bit for bit, has no
#' \code{[covariate_nn]} block, and loads in a ferx build without the
#' \code{nn} cargo feature.
#'
#' The transform is ferx-core's \code{nn::export} (ferx-core #1331, also
#' \code{ferx export-dcm} on the command line). It never evaluates a network,
#' so it runs in every build, including one without \code{nn}.
#'
#' @param model Path to the DCM \code{.ferx} model file, or a
#'   \code{\link{ferx_model}} object.
#' @param fit A \code{\link{ferx_fit}} result of \code{model}: supplies the
#'   fitted weights and every other estimate. The fit must have been made with
#'   this model - the export compares \code{fit$model_text} with
#'   \code{model} and refuses a fit whose networks or equations differ
#'   (comments, layout, \code{[parameters]} values and the \code{[data]} path
#'   may differ). Give \code{fit} or \code{theta}, not both.
#' @param output_model Path the exported \code{.ferx} file is written to. When
#'   \code{NULL} (default), \code{<dir>/<stem>_expanded.ferx} next to
#'   \code{model}. For a \code{\link{ferx_example}()} model that directory is
#'   inside the installed package, so pass a writable path such as
#'   \code{file.path(tempdir(), "dcm_expanded.ferx")}. An existing file is
#'   overwritten; the model file itself is refused.
#' @param weights How the exported model carries the network weights:
#'   \code{"literal"} (default) inlines each weight as a numeric constant, so
#'   the exported network is frozen; \code{"fixed_theta"} declares each one as
#'   \code{theta W_<NN>_<l>_<i>_<j>(value, lower, upper) FIX} in
#'   \code{[parameters]}, so the weights stay named and can be perturbed by
#'   hand. Either way the weights are not estimated: do not refit them.
#' @param network_only Logical. \code{FALSE} (default) writes every other
#'   theta, omega, sigma and kappa at the fit's estimate, so the export is the
#'   fitted model. \code{TRUE} keeps the values \code{model} declares and
#'   rewrites only the networks; needed for a model whose estimates cannot all
#'   be written back (a \code{[mixture]} per-class override, for one), which
#'   the default refuses.
#' @param theta Instead of \code{fit}: a named numeric vector holding at least
#'   the network weights (\code{W_<NN>_*} / \code{B_<NN>_*}), for instance a
#'   saved \code{fit$theta}. Only the networks are rewritten and nothing is
#'   checked against a fit, as with \code{network_only = TRUE}.
#'
#' @return A \code{\link{ferx_model}} pointing at the exported model file, so
#'   it composes with \code{\link{ferx_predict}}, \code{\link{ferx_simulate}}
#'   and the other model helpers. It carries \code{model}'s dataset. The
#'   exported file's own \code{[data]} block is copied verbatim, so a relative
#'   path in it resolves from the exported file's directory.
#'
#' @details
#' The exported equations are wrapped so that they are not mu-referenced; the
#' model is meant for prediction and simulation at the fitted values, not for
#' re-estimation.
#'
#' @seealso \code{\link{ferx_model_to_frem}}, the other model-to-model
#'   transform.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin_dcm")
#' fit <- ferx_fit(ex$model, ex$data)
#' plain <- ferx_model_from_dcm(ex$model, fit,
#'                              output_model = file.path(tempdir(), "dcm.ferx"))
#' ferx_predict(plain$model, ex$data)
#' }
#'
#' @export
ferx_model_from_dcm <- function(model,
                                fit = NULL,
                                output_model = NULL,
                                weights = c("literal", "fixed_theta"),
                                network_only = FALSE,
                                theta = NULL) {
  # --- resolve model path ---
  data <- NULL
  if (inherits(model, "ferx_model")) {
    data <- model$data
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
  weights <- match.arg(weights)
  if (!is.logical(network_only) || length(network_only) != 1L || is.na(network_only)) {
    stop("`network_only` must be TRUE or FALSE.")
  }
  if (is.null(fit) == is.null(theta)) {
    stop("Give exactly one of `fit` (a ferx_fit of `model`) or `theta` ",
         "(a named numeric vector of the network weights).")
  }
  if (!is.null(fit) && !inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit()).")
  }
  if (!is.null(theta) &&
        (!is.numeric(theta) || is.null(names(theta)) || anyNA(names(theta)) ||
           !all(nzchar(names(theta))))) {
    stop("`theta` must be a named numeric vector (e.g. a saved `fit$theta`).")
  }
  if (!is.null(output_model) &&
        (!is.character(output_model) || length(output_model) != 1L ||
           !nzchar(output_model))) {
    stop("`output_model` must be a single non-empty path, or NULL.")
  }

  # --- resolve paths ---
  model_path <- normalizePath(model_path, mustWork = TRUE)
  if (is.null(output_model)) {
    stem <- tools::file_path_sans_ext(basename(model_path))
    output_model <- file.path(dirname(model_path), paste0(stem, "_expanded.ferx"))
  } else if (dir.exists(dirname(output_model))) {
    # Absolute, so the returned ferx_model does not depend on the working
    # directory. An unresolvable directory is left as written; the write
    # then fails against the path the caller gave.
    output_model <- file.path(normalizePath(dirname(output_model), mustWork = FALSE),
                              basename(output_model))
  }
  if (file.exists(output_model) &&
        normalizePath(output_model, mustWork = TRUE) == model_path) {
    stop("`output_model` is the model file itself; write the export to another path.")
  }
  # The dataset the returned ferx_model carries: the ferx_model's own, else the
  # model file's [data] block, resolved against the model file.
  if (is.null(data)) data <- .ferx_model_data_path(model_path, strict = FALSE)

  # --- call Rust backend ---
  # With `fit`, the engine checks the model is the fit's own and writes every
  # estimate; with `theta`, it rewrites the networks alone.
  args <- if (!is.null(fit)) {
    pieces <- validate_fit_for_params(fit)
    list(
      network_only       = network_only,
      from_fit           = TRUE,
      fit_model_text     = as.character(fit$model_text %||% ""),
      theta_names        = names(fit$theta) %||% character(0),
      theta              = pieces$theta,
      omega_flat         = pieces$omega_flat,
      omega_dim          = pieces$omega_dim,
      eta_names          = rownames(fit$omega) %||% character(0),
      sigma              = pieces$sigma,
      sigma_names        = names(fit$sigma) %||% character(0),
      omega_iov_flat     = pieces$omega_iov_flat,
      omega_iov_dim      = pieces$omega_iov_dim,
      kappa_names        = rownames(fit$omega_iov) %||% character(0),
      residual_rho       = pieces$residual_rho,
      residual_rho_pairs = pieces$residual_rho_pairs
    )
  } else {
    list(
      network_only       = TRUE,
      from_fit           = FALSE,
      fit_model_text     = "",
      theta_names        = names(theta),
      theta              = as.numeric(theta),
      omega_flat         = numeric(0),
      omega_dim          = 0L,
      eta_names          = character(0),
      sigma              = numeric(0),
      sigma_names        = character(0),
      omega_iov_flat     = numeric(0),
      omega_iov_dim      = 0L,
      kappa_names        = character(0),
      residual_rho       = numeric(0),
      residual_rho_pairs = integer(0)
    )
  }
  raw <- .ferx_engine_call(
    do.call(ferx_rust_export_dcm, c(list(
      model_path        = model_path,
      output_model_path = output_model,
      weight_form       = weights
    ), args)),
    model_path, NULL
  )

  ferx_model(data = data, model = raw$model_path)
}
