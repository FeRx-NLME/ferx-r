#' Translate an mrgsolve model file to ferx's native DSL
#'
#' Reads an mrgsolve \code{.cpp} or \code{.mod} model file and emits an
#' equivalent \code{.ferx} source string. The translated text uses the same
#' parser path as hand-written \code{.ferx} models — feeding the result of
#' this function (or its written output file) to \code{\link{ferx_fit}},
#' \code{\link{ferx_simulate}}, or \code{\link{ferx_predict}} will run a
#' standard ferx fit.
#'
#' Supported mrgsolve blocks: \code{$PROB}, \code{$PARAM}, \code{$INPUT},
#' \code{$THETA}, \code{$CMT}, \code{$INIT}, \code{$OMEGA} / \code{$SIGMA}
#' (diagonal + \code{@block} + \code{@labels}), \code{$MAIN} / \code{$PK},
#' \code{$ODE}, \code{$TABLE} / \code{$ERROR}, \code{$CAPTURE},
#' \code{$PKMODEL}, \code{$SET}.
#'
#' Unsupported blocks (\code{$GLOBAL}, \code{$PLUGIN}, \code{$NMXML},
#' \code{$NMEXT}, \code{$PRED}, \code{$EVENT}) and C++ constructs
#' (\code{for}/\code{while}/\code{return}/ternary/C casts/\code{pow(a,b)})
#' raise an error with a \code{file:line:} pointer at the offending line.
#'
#' The \code{$TABLE} block's terminal observation expression must match one
#' of ferx's three error shapes: \code{Y = IPRED * (1 + EPS(p))}
#' (proportional), \code{Y = IPRED + EPS(a)} (additive), or
#' \code{Y = IPRED * (1 + EPS(p)) + EPS(a)} (combined). Other shapes
#' (e.g. log-normal residuals) are rejected with a hint listing the
#' supported forms.
#'
#' @param input Path to a mrgsolve \code{.cpp} or \code{.mod} model file.
#'   The file must exist.
#' @param output Optional path to write the translated \code{.ferx} source
#'   to. When \code{NULL} (default), the text is returned to the caller and
#'   nothing is written to disk.
#' @param overwrite If \code{output} is provided and already exists, set
#'   this to \code{TRUE} to replace it. Default \code{FALSE} (refuse to
#'   overwrite).
#'
#' @return The translated \code{.ferx} source as a single character string.
#'   When \code{output} is provided, the string is also written to that
#'   file and returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Translate an mrgsolve model and print the result
#' ferx_text <- ferx_translate_mrgsolve("warfarin.mod")
#' cat(ferx_text)
#'
#' # Translate and save in one step, then fit
#' ferx_translate_mrgsolve("warfarin.mod", output = "warfarin.ferx")
#' fit <- ferx_fit(ferx_model("data.csv", "warfarin.ferx"))
#' }
#'
#' @seealso \code{\link{ferx_fit}}, \code{\link{ferx_model_validate}}
#' @export
ferx_translate_mrgsolve <- function(input, output = NULL, overwrite = FALSE) {
  if (!is.character(input) || length(input) != 1L) {
    stop("'input' must be a single file path", call. = FALSE)
  }
  if (!file.exists(input)) {
    stop("File not found: ", input, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(input))
  if (!ext %in% c("cpp", "mod")) {
    warning(
      "'", input, "' does not end in .cpp or .mod — proceeding anyway, but ",
      "the translator expects mrgsolve syntax",
      call. = FALSE
    )
  }

  result <- ferx_rust_translate_mrgsolve(normalizePath(input))
  if (!isTRUE(result$ok)) {
    stop("mrgsolve translation failed: ", result$error, call. = FALSE)
  }
  ferx_text <- result$ferx_source

  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L) {
      stop("'output' must be NULL or a single file path", call. = FALSE)
    }
    if (file.exists(output) && !isTRUE(overwrite)) {
      stop(
        "Output file already exists: ", output,
        ". Pass overwrite = TRUE to replace it.",
        call. = FALSE
      )
    }
    writeLines(ferx_text, output)
    return(invisible(ferx_text))
  }

  ferx_text
}
