#' Inspect a section of a ferx model (pipe-friendly)
#'
#' Prints the contents of a named section of a \code{.ferx} model file to the
#' console and returns the input object invisibly so the pipe can continue.
#'
#' When \code{x} is a \code{ferx_model} object the section is read from
#' \code{x$model}. When \code{x} is a plain path string the call is equivalent
#' to \code{\link{ferx_model_section}(x, section)}.
#'
#' @param x A \code{ferx_model} object or a path to a \code{.ferx} file.
#' @param section Name of the section to display, without brackets.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#'
#' # Inspect [parameters] and continue piping: ferx_get_section() prints
#' # the section and passes the ferx_model object through unchanged.
#' ferx_model(ex$data, ex$model) |>
#'   ferx_get_section("parameters")
#'
#' \dontrun{
#' # Mid-pipe peek: print [parameters] then fit without interrupting the chain
#' fit <- ferx_model(ex$data, ex$model) |>
#'   ferx_get_section("parameters") |>
#'   ferx_fit(method = "focei")
#' }
#'
#' @seealso \code{\link{ferx_model_section}}, \code{\link{ferx_set_section}}
#' @family model-editing
#' @export
ferx_get_section <- function(x, section) {
  if (inherits(x, "ferx_model")) {
    ferx_model_section(x$model, section)
  } else {
    ferx_model_section(x, section)
  }
  invisible(x)
}
