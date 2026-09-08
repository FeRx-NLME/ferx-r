#' Parse and expand an MFL search space
#'
#' Parses a Model Feature Language (MFL) search space - Pharmpy's grammar, so a
#' space is portable between the two - and, when a base model is supplied,
#' resolves every \code{@}-symbol and wildcard against it. That is what makes an
#' interactive space debuggable: \code{COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])}
#' becomes the explicit parameter x covariate x effect set it stands for on
#' \emph{your} model, before committing to a run.
#'
#' Symbols are resolved against the model's \code{[individual_parameters]}, its
#' template line (\code{pk NAME(...)} / \code{ode_template NAME(...)}), its
#' \code{[covariates]} block and the dataset: \code{@IIV} is every individual
#' parameter carrying an eta, \code{@PK} every parameter bound on the template
#' line, \code{@CONTINUOUS} / \code{@CATEGORICAL} the covariates declared with
#' that type. \code{LET(NAME, [...])} overrides a built-in of the same name.
#' A symbol that resolves to nothing drops its feature with a note (Pharmpy's
#' behaviour) rather than erroring; a name that is neither an individual
#' parameter nor a declared covariate is an error naming it.
#'
#' Resolution needs the dataset (the covariate columns come from it), so
#' \code{data} - or a \code{[data]} block in the model - is required whenever
#' \code{model} is given.
#'
#' @param mfl MFL source text, a \code{ferx_search_config} (its \code{[space]}
#'   is used, together with its \code{base} / \code{data} when \code{model} is
#'   not given), or a \code{ferx_search_space} to re-resolve.
#' @param model Optional path to a \code{.ferx} model, or a \code{ferx_model}
#'   object. \code{NULL} parses the space without resolving it.
#' @param data Optional path to the dataset. Defaults to the model's
#'   \code{[data]} block when the model declares one.
#'
#' @return A data frame of class \code{ferx_search_space} with one row per
#'   feature - \code{feature} (the rendered statement), \code{keyword} and
#'   \code{optional} (\code{TRUE} for an exploratory \code{FEATURE?}) - holding
#'   the resolved (ground) features when a model was supplied and the space as
#'   written otherwise. Attributes: \code{mfl} (the source text),
#'   \code{resolved} (logical), \code{parsed} (the unresolved feature table),
#'   \code{covariate_effects} (a data frame of the ground covariate effects:
#'   \code{parameter}, \code{covariate}, \code{effect}, \code{op},
#'   \code{optional}) and \code{notes} (features the resolver dropped, and
#'   anything else it wants seen once).
#'
#' @examples
#' # Parse only - no model needed
#' ferx_search_space("COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])")
#'
#' # Expanded against a model: what the search would actually explore
#' ex <- ferx_example("two_cpt_oral_cov")
#' sp <- ferx_search_space(
#'   "COVARIATE?(@IIV, @CONTINUOUS, [pow, lin])",
#'   model = ex$model, data = ex$data
#' )
#' attr(sp, "covariate_effects")
#'
#' @seealso \code{\link{ferx_search_coverage}},
#'   \code{\link{ferx_search_config}}
#' @family search
#' @export
ferx_search_space <- function(mfl, model = NULL, data = NULL) {
  if (inherits(mfl, "ferx_search_config")) {
    if (is.null(model)) {
      model <- mfl$base
      if (is.null(data)) data <- mfl$data
    }
    mfl <- mfl$mfl
  } else if (inherits(mfl, "ferx_search_space")) {
    mfl <- attr(mfl, "mfl")
  }
  if (!is.character(mfl) || length(mfl) != 1L || is.na(mfl)) {
    stop("'mfl' must be a single MFL string, a ferx_search_config, or a ",
         "ferx_search_space")
  }

  if (inherits(model, "ferx_model")) {
    if (is.null(data)) data <- model$data
    model <- model$model
  }
  if (is.null(model)) {
    if (!is.null(data)) {
      stop("'data' needs a 'model' - a space is resolved against the two together")
    }
    model_path <- ""
    data_path  <- ""
  } else {
    if (!file.exists(model)) stop("Model file not found: ", model)
    # NULL data defers to the model's [data] block, which the engine resolves.
    if (!is.null(data) && !file.exists(data)) stop("Data file not found: ", data)
    model_path <- normalizePath(model)
    data_path  <- if (is.null(data)) "" else normalizePath(data)
  }

  raw <- ferx_rust_search_space_parse(mfl, model_path, data_path)

  parsed <- .ferx_space_frame(raw$feature, raw$keyword, raw$optional)
  resolved <- isTRUE(raw$resolved)
  space <- if (resolved) {
    .ferx_space_frame(raw$resolved_feature, raw$resolved_keyword,
                      raw$resolved_optional)
  } else {
    parsed
  }

  effects <- data.frame(
    parameter = as.character(raw$effect_parameter),
    covariate = as.character(raw$effect_covariate),
    effect    = as.character(raw$effect_form),
    op        = as.character(raw$effect_op),
    optional  = as.logical(raw$effect_optional),
    stringsAsFactors = FALSE
  )

  structure(
    space,
    mfl               = mfl,
    resolved          = resolved,
    parsed            = parsed,
    covariate_effects = effects,
    notes             = as.character(raw$notes),
    class             = c("ferx_search_space", "data.frame")
  )
}

#' @param x A \code{ferx_search_space}.
#' @param ... Passed to \code{print.data.frame}.
#' @rdname ferx_search_space
#' @export
print.ferx_search_space <- function(x, ...) {
  resolved <- isTRUE(attr(x, "resolved"))
  cat("<ferx search space>",
      if (resolved) " (resolved against the base model)" else " (as written)",
      "\n", sep = "")
  cat("  mfl = ", attr(x, "mfl"), "\n\n", sep = "")
  print(as.data.frame(x), ...)

  effects <- attr(x, "covariate_effects")
  if (!is.null(effects) && nrow(effects) > 0L) {
    cat("\n", nrow(effects), " covariate effect",
        if (nrow(effects) == 1L) "" else "s",
        " (see attr(, \"covariate_effects\"))\n", sep = "")
  }
  notes <- attr(x, "notes")
  if (length(notes)) {
    cat("\nNotes:\n")
    for (n in notes) cat("  * ", n, "\n", sep = "")
  }
  invisible(x)
}

#' Coverage of a search space against what the engine can express
#'
#' The engine refuses to narrow a search space silently: a feature it cannot
#' build a candidate for is an error at load
#' (\code{\link{ferx_search_config}}). This function reports the same
#' information non-fatally, as a table, so an unsupported feature is a row to
#' read rather than a run that aborts - the way to find out that
#' \code{ELIMINATION(MM)} has no candidate before writing the file.
#'
#' Covered features are listed as written. An uncovered feature is listed once
#' per gap, under the name the engine reports - the offending mode
#' (\code{ELIMINATION(MM)}), which for a wildcard such as
#' \code{ELIMINATION(*)} is the more useful identity than the source statement.
#'
#' @param space A \code{ferx_search_space}, a \code{ferx_search_config}, or MFL
#'   source text.
#'
#' @return A data frame with columns \code{feature}, \code{covered} (logical)
#'   and \code{reason} (\code{NA} for a covered feature, otherwise the engine's
#'   explanation).
#'
#' @examples
#' ferx_search_coverage("ABSORPTION([INST, FO]); ELIMINATION(FO)")
#'
#' # A gap is a row, not an error
#' ferx_search_coverage("ELIMINATION([FO, MM])")
#'
#' @seealso \code{\link{ferx_search_space}}, \code{\link{ferx_search_config}}
#' @family search
#' @export
ferx_search_coverage <- function(space) {
  mfl <- if (inherits(space, "ferx_search_space")) {
    attr(space, "mfl")
  } else if (inherits(space, "ferx_search_config")) {
    space$mfl
  } else {
    space
  }
  if (!is.character(mfl) || length(mfl) != 1L || is.na(mfl)) {
    stop("'space' must be a ferx_search_space, a ferx_search_config, or a ",
         "single MFL string")
  }

  raw <- ferx_rust_search_coverage(mfl)
  reason <- as.character(raw$reason)
  data.frame(
    feature = as.character(raw$feature),
    covered = as.logical(raw$covered),
    reason  = ifelse(nzchar(reason), reason, NA_character_),
    stringsAsFactors = FALSE
  )
}
