#' Extract a section from a ferx model file
#'
#' Returns the lines belonging to a named section of a \code{.ferx} model file,
#' excluding the section header itself. Prints the lines to the console and
#' returns them invisibly.
#'
#' @param path Path to a \code{.ferx} model file.
#' @param section Name of the section to extract, without brackets (e.g.
#'   \code{"parameters"}).
#' @param strip Logical. If \code{TRUE}, leading whitespace is trimmed from each
#'   returned line via \code{\link[base]{trimws}}. Defaults to \code{FALSE} to
#'   preserve the round-trip guarantee with
#'   \code{\link{ferx_model_set_section}}.
#'
#' @return Character vector of lines in the requested section, invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_model_section(ex$model, "parameters")
#' ferx_model_section(ex$model, "parameters", strip = TRUE)
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model_edit}},
#'   \code{\link{ferx_model_set_section}}
#' @family model-editing
#' @export
ferx_model_section <- function(path, section, strip = FALSE) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")

  file_lines <- readLines(path, warn = FALSE)
  hdr        <- ferx_section_headers(file_lines)

  idx <- which(hdr$names == section)
  if (length(idx) == 0L) {
    stop(
      "Section '", section, "' not found. ",
      "Available sections: ", paste(hdr$names, collapse = ", ")
    )
  }

  start <- hdr$positions[idx] + 1L
  end   <- if (idx < length(hdr$positions)) hdr$positions[idx + 1L] - 1L else length(file_lines)
  body  <- if (start <= end) file_lines[start:end] else character(0)
  if (strip) body <- trimws(body, which = "left")

  cat("# [", section, "]\n", sep = "")
  cat(body, sep = "\n")
  cat("\n")
  invisible(body)
}
