#' Allometric scaling
#'
#' Adds allometric body-size scaling to a model: \code{(WT/70)^0.75} on every
#' clearance the \code{pk} template line binds (\code{cl}, \code{q}, \code{q2},
#' \code{q3}) and \code{(WT/70)^1.0} on every volume (\code{v}, \code{v1},
#' \code{v2}, \code{v3}), written as \code{[covariate_model]} relations. This is
#' Pharmpy's \code{allometry} tool, and it is a convention rather than new
#' machinery - the same thing you would write by hand as
#' \code{power(center = 70, fix = 0.75)}.
#'
#' With \code{fit = FALSE} this is a model transform: you get the scaled model
#' text back and nothing is fitted, so allometry can be one step of a hand-built
#' workflow. With \code{fit = TRUE} (the default) the base and the scaled model
#' are fitted side by side and both outcomes are returned, which is what makes
#' the scaling's cost visible.
#'
#' A parameter that already carries a relation on the size covariate is left
#' alone, as Pharmpy does, with a note saying so.
#'
#' @param model Path to a \code{.ferx} model, or a \code{ferx_model} object.
#'   Omit when \code{config} is given.
#' @param data Path to the dataset. Defaults to the model's \code{[data]}
#'   block. Names the dataset the search runs on, so like \code{model} it
#'   cannot be given beside \code{config} - the file's own \code{data} key
#'   says which dataset that file searches.
#' @param config Path to a \code{.ferxsearch} file carrying an
#'   \code{ALLOMETRY(WT, 70)} statement and an optional \code{[allometry]}
#'   section. Mutually exclusive with the arguments that state the scaling.
#' @param covariate The size covariate (engine default \code{"WT"}).
#' @param reference The reference value it is divided by (engine default 70).
#' @param parameters Parameters to scale. \code{NULL} takes every clearance and
#'   volume the template line binds.
#' @param exponents One exponent per \code{parameters} entry. \code{NULL} uses
#'   the convention: 0.75 for a clearance, 1.0 for a volume.
#' @param estimate Estimate the exponents from those values instead of fixing
#'   them, bounded by \code{lower} and \code{upper}.
#' @param lower,upper Bounds of an estimated exponent (engine defaults 0 and 2).
#' @param fit Fit the base and scaled models, or only build the scaled one.
#' @param threads Total worker threads. \code{NULL} lets the runner choose.
#' @param retries Perturbed restarts per fit on top of the exact one.
#'   \code{NULL} keeps the engine default.
#' @param directory Where the two fits are journalled. \code{NULL} keeps them in
#'   memory.
#'
#' @return An object of class \code{ferx_allometry}:
#'   \describe{
#'     \item{scalings}{One row per scaled parameter: \code{parameter},
#'       \code{exponent}, \code{fixed} and the \code{theta} an estimated
#'       exponent declares.}
#'     \item{model, model_path}{The scaled model as text, and as a file.}
#'     \item{covariate, reference}{The scaling the relations were built on.}
#'     \item{comparison}{With \code{fit = TRUE}: the two fits side by side -
#'       \code{model} (\code{"base"} / \code{"scaled"}), \code{ofv},
#'       \code{converged}, \code{passed} and \code{failures}.}
#'     \item{fit, base_fit}{With \code{fit = TRUE}: the scaled and base fits as
#'       \code{ferx_fit} objects.}
#'     \item{dofv}{\code{OFV(base) - OFV(scaled)}, when both fits exist.}
#'     \item{notes, cancelled}{Parameters left alone and why, and whether the
#'       run was stopped early.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # The transform on its own: the scaled model, nothing fitted
#' scaled <- ferx_allometry(ex$model, ex$data, fit = FALSE)
#' cat(scaled$model)
#'
#' # The tool: both fits, side by side
#' res <- ferx_allometry(ex$model, ex$data, covariate = "WT", reference = 70)
#' res
#' res$comparison
#' }
#'
#' @seealso \code{\link{ferx_covsearch}}, \code{\link{ferx_search_config}}
#' @family search
#' @export
ferx_allometry <- function(model = NULL,
                           data = NULL,
                           config = NULL,
                           covariate = NULL,
                           reference = NULL,
                           parameters = NULL,
                           exponents = NULL,
                           estimate = FALSE,
                           lower = NULL,
                           upper = NULL,
                           fit = TRUE,
                           threads = NULL,
                           retries = NULL,
                           directory = NULL) {
  what <- "ferx_allometry"
  config_path <- .ferx_search_entry_form(
    config, model, data,
    list(
      covariate = covariate,
      reference = reference,
      parameters = parameters,
      exponents = exponents,
      lower = lower,
      upper = upper
    ),
    what
  )

  if (nzchar(config_path)) {
    paths <- list(model = "", data = "")
  } else {
    paths <- .ferx_search_model_data(model, data, what)
  }

  if (!is.null(parameters) && (!is.character(parameters) || anyNA(parameters))) {
    stop(what, ": `parameters` must be a character vector of parameter names")
  }
  if (!is.null(exponents) && (!is.numeric(exponents) || anyNA(exponents))) {
    stop(what, ": `exponents` must be a numeric vector, one per `parameters` entry")
  }
  if (!is.null(covariate) && (!is.character(covariate) || length(covariate) != 1L ||
                              is.na(covariate))) {
    stop(what, ": `covariate` must be a single covariate name")
  }

  dir_arg <- .ferx_search_directory(directory, what)
  raw <- ferx_rust_allometry(
    config_path = config_path,
    model_path  = paths$model,
    data_path   = paths$data,
    covariate   = covariate %||% "",
    reference   = .ferx_search_scalar(reference, "reference", what, positive = TRUE),
    parameters  = as.character(parameters %||% character(0)),
    exponents   = as.numeric(exponents %||% numeric(0)),
    fixed       = !.ferx_search_bool(estimate, "estimate", what),
    lower       = .ferx_search_scalar(lower, "lower", what),
    upper       = .ferx_search_scalar(upper, "upper", what),
    threads     = .ferx_search_count(threads, "threads", what, below = 0L),
    retries     = .ferx_search_count(retries, "retries", what, below = -1L, min = 0L),
    directory   = dir_arg,
    fit         = .ferx_search_bool(fit, "fit", what)
  )

  scalings <- data.frame(
    parameter = as.character(raw$parameter),
    exponent  = as.numeric(raw$exponent),
    fixed     = as.logical(raw$fixed),
    theta     = .ferx_search_chr(raw$theta),
    stringsAsFactors = FALSE
  )

  # The scaled model as a file: written beside the run when it has a directory,
  # otherwise to a temporary file, so it can be fitted or edited as it stands.
  model_path <- if (nzchar(dir_arg)) file.path(dir_arg, "allometric.ferx") else
    tempfile(pattern = "ferx-allometry-", fileext = ".ferx")
  if (nzchar(dir_arg)) dir.create(dir_arg, showWarnings = FALSE, recursive = TRUE)
  writeLines(as.character(raw$model), model_path)

  result <- list(
    scalings   = scalings,
    model      = as.character(raw$model),
    model_path = model_path,
    base_model = as.character(raw$base_model),
    covariate  = as.character(raw$covariate),
    reference  = as.numeric(raw$reference),
    data       = as.character(raw$data),
    directory  = if (nzchar(dir_arg)) dir_arg else NA_character_,
    config     = if (nzchar(config_path)) config_path else NA_character_,
    fitted     = isTRUE(raw$fitted),
    notes      = as.character(raw$notes),
    cancelled  = isTRUE(raw$cancelled %||% FALSE)
  )
  # `fit` is always an element, even when there is nothing in it: assigning
  # NULL would drop it, and `res$fit` would then partially match `res$fitted`
  # and answer FALSE - a fit that reads as a logical.
  result["fit"] <- list(NULL)
  result["base_fit"] <- list(NULL)

  if (isTRUE(raw$fitted)) {
    result$comparison <- data.frame(
      model     = c("base", "scaled"),
      ofv       = .ferx_search_num(c(raw$base_ofv, raw$scaled_ofv)),
      converged = .ferx_search_lgl(c(raw$base_converged, raw$scaled_converged)),
      passed    = as.logical(c(raw$base_passed, raw$scaled_passed)),
      failures  = c(paste(as.character(raw$base_failures), collapse = "; "),
                    paste(as.character(raw$scaled_failures), collapse = "; ")),
      stringsAsFactors = FALSE
    )
    result$comparison$failures <- .ferx_search_chr(result$comparison$failures)
    result$dofv <- .ferx_search_num(raw$dofv)
    if (!is.null(raw$scaled_fit)) {
      result["fit"] <- list(.ferx_fit_from_raw(raw$scaled_fit, model = model_path,
                                               data = as.character(raw$data)))
    }
    if (!is.null(raw$base_fit)) {
      base_path <- tempfile(pattern = "ferx-allometry-base-", fileext = ".ferx")
      writeLines(as.character(raw$base_model), base_path)
      result["base_fit"] <- list(.ferx_fit_from_raw(raw$base_fit, model = base_path,
                                                    data = as.character(raw$data)))
    }
  }

  class(result) <- c("ferx_allometry", "ferx_search_result")
  result
}

#' @param x A \code{ferx_allometry} object.
#' @param digits Significant digits for the printed tables.
#' @param ... Ignored.
#' @rdname ferx_allometry
#' @export
print.ferx_allometry <- function(x, digits = 4, ...) {
  cat("ferx allometric scaling\n")
  cat("  Data:  ", x$data, "\n", sep = "")
  cat("  Model: ", x$model_path, "\n", sep = "")
  if (isTRUE(x$cancelled)) cat("  Cancelled before both fits finished.\n")

  cat(sprintf("\nScaling on %s (reference %s):\n", x$covariate,
              format(signif(x$reference, digits))))
  if (nrow(x$scalings) == 0L) {
    cat("  (nothing scaled)\n")
  } else {
    for (i in seq_len(nrow(x$scalings))) {
      cat(sprintf("  %s ~ %s power(center = %s, %s)\n",
                  x$scalings$parameter[i], x$covariate,
                  format(signif(x$reference, digits)),
                  if (x$scalings$fixed[i]) {
                    sprintf("fix = %s", format(x$scalings$exponent[i]))
                  } else {
                    sprintf("init = %s, estimated", format(x$scalings$exponent[i]))
                  }))
    }
  }

  if (isTRUE(x$fitted)) {
    cat("\nFits:\n")
    .ferx_search_print_table(x$comparison, digits)
    if (!is.na(x$dofv)) {
      cat(sprintf("\n  dOFV (base - scaled): %s\n", format(signif(x$dofv, digits))))
    }
  } else {
    cat("\nNot fitted (`fit = FALSE`); the scaled model is in `$model`.\n")
  }
  .ferx_search_print_notes(x$notes)
  invisible(x)
}
