#' Derive NCA-based starting values without fitting
#'
#' \code{ferx_inits_from_nca()} estimates population starting values directly
#' from the data using non-compartmental analysis (NCA), optionally refining
#' parameters NCA cannot estimate with a grid sweep. It runs the same estimation
#' that the \code{inits_from_nca} argument to \code{\link{ferx_fit}} performs,
#' but \emph{returns the suggested values for inspection} instead of running a
#' fit. This is useful for sanity-checking starting values, or for feeding them
#' back into a model file before fitting.
#'
#' All three methods are NCA-based; they differ only in how much refinement runs
#' on top of the NCA estimates:
#' \describe{
#'   \item{\code{"nca"}}{Non-compartmental analysis only (fastest). Leaves
#'     parameters NCA cannot reach (peripheral \code{Q}/\code{V2}, all ODE/PD
#'     thetas) at the model default.}
#'   \item{\code{"nca_sweep"}}{NCA, then a log-space rRMSE grid sweep (population
#'     predictions, etas = 0) over the remaining non-fixed thetas. Model-agnostic:
#'     also covers ODE/PD models. This is the default.}
#'   \item{\code{"nca_ebe"}}{Like \code{"nca_sweep"} but evaluates the grid with
#'     empirical Bayes estimates (etas != 0); more accurate under large
#'     between-subject variability. Falls back to \code{"nca_sweep"} for ODE
#'     models.}
#' }
#'
#' @param model A path to a \code{.ferx} model file, or a \code{ferx_model}
#'   object created with \code{\link{ferx_model}}.
#' @param data Path to a NONMEM-format CSV file. When \code{model} is a
#'   \code{ferx_model} object that already carries a data path, this may be left
#'   \code{NULL}.
#' @param method One of \code{"nca_sweep"} (default), \code{"nca"}, or
#'   \code{"nca_ebe"}. See Details.
#'
#' @return An object of class \code{ferx_inits}: a list with elements
#'   \describe{
#'     \item{theta}{Named numeric vector of suggested theta starting values.}
#'     \item{omega}{Suggested omega matrix (with eta names as dimnames). The CL
#'       eta variance is updated from inter-subject CV-squared; other entries
#'       keep the model default.}
#'     \item{method}{The method actually used (a string).}
#'     \item{warnings}{Character vector of notes about what was and wasn't
#'       estimated.}
#'   }
#'
#' @seealso \code{\link{ferx_fit}}, whose \code{inits_from_nca} argument applies
#'   these values automatically before the optimizer runs.
#'
#' @note The \code{method} argument always wins: a \code{[fit_options]
#'   inits_from_nca} key in the model file is \emph{not} consulted by this
#'   function. To preview the strategy a model file declares, pass that value
#'   explicitly via \code{method}.
#'
#' @examples
#' \donttest{
#' ex <- ferx_example("warfarin")
#'
#' # -- Classic style (path strings) ------------------------------------------
#' inits <- ferx_inits_from_nca(ex$model, ex$data)
#' inits$theta
#' inits$omega
#'
#' # Pick a specific strategy
#' ferx_inits_from_nca(ex$model, ex$data, method = "nca")
#' ferx_inits_from_nca(ex$model, ex$data, method = "nca_ebe")
#'
#' # -- Pipe style (ferx_model object) ----------------------------------------
#'
#' # data |> ferx_model(model) bundles the paths into a ferx_model, which
#' # ferx_inits_from_nca() unpacks automatically (the second arg can stay NULL).
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_inits_from_nca()
#'
#' # Inspect first, then fit using the same strategy in one chain.
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_inits_from_nca(method = "nca_ebe") |>
#'   print()
#'
#' ferx_fit(ex$model, ex$data, inits_from_nca = "nca_ebe")
#'
#' # Peek at the suggested CL/V before committing to a fit.
#' ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_inits_from_nca() |>
#'   (\(x) x$theta[c("TVCL", "TVV")])()
#' }
#'
#' @family fitting
#' @export
ferx_inits_from_nca <- function(model, data = NULL,
                                method = c("nca_sweep", "nca", "nca_ebe")) {
  method <- match.arg(method)

  if (inherits(model, "ferx_model")) {
    if (is.null(data)) {
      data <- model$data
    }
    model <- model$model
  }
  if (is.null(data)) {
    stop("`data` is required. Pass a path to a NONMEM CSV file.")
  }
  stopifnot(file.exists(model), file.exists(data))

  res <- ferx_rust_inits_from_nca(
    normalizePath(model),
    normalizePath(data),
    method
  )
  if (length(res) == 0L) {
    stop("`ferx_inits_from_nca()` failed - see the messages above for details.")
  }

  theta <- as.numeric(res$theta)
  names(theta) <- res$theta_names

  d <- res$omega_dim %||% 0L
  omega <- if (length(res$omega) > 0L && d > 0L) {
    m <- matrix(as.numeric(res$omega), nrow = d, ncol = d)
    eta_nms <- res$eta_names
    if (length(eta_nms) == d) {
      rownames(m) <- colnames(m) <- eta_nms
    }
    m
  } else {
    matrix(numeric(0), 0L, 0L)
  }

  structure(
    list(
      theta = theta,
      omega = omega,
      method = res$method,
      warnings = as.character(res$warnings)
    ),
    class = "ferx_inits"
  )
}

#' @export
print.ferx_inits <- function(x, ...) {
  cat(sprintf(
    "<ferx_inits> NCA-based starting values (method = %s)\n\n", x$method
  ))
  cat("Theta:\n")
  print(x$theta)
  if (!is.null(x$omega) && nrow(x$omega) > 0L) {
    cat("\nOmega:\n")
    print(x$omega)
  }
  if (length(x$warnings) > 0L) {
    cat(sprintf("\n%d note(s):\n", length(x$warnings)))
    for (w in x$warnings) cat(" -", w, "\n")
  }
  invisible(x)
}
