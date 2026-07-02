#' Replace a section in a ferx model file
#'
#' Overwrites the body of a named section in a \code{.ferx} file with new
#' lines, leaving all other sections untouched. Use this to modify a model
#' programmatically without opening an editor.
#'
#' @param path Path to a \code{.ferx} model file.
#' @param section Name of the section to replace, without brackets (e.g.
#'   \code{"fit_options"}).
#' @param lines Character vector of replacement lines. These become the new
#'   body of the section (do not include the \code{[section]} header line).
#'
#' @return \code{path}, invisibly.
#'
#' @examples
#' \dontrun{
#' # Switch estimation method without opening an editor
#' ferx_model_set_section("my_model.ferx", "fit_options", c(
#'   "  method     = focei",
#'   "  maxiter    = 500",
#'   "  covariance = false"
#' ))
#'
#' # Read-modify-write a section
#' lines <- ferx_model_section("my_model.ferx", "parameters")
#' lines <- sub("TVCL\\(.*\\)", "TVCL(0.5, 0.001, 10.0)", lines)
#' ferx_model_set_section("my_model.ferx", "parameters", lines)
#' }
#'
#' @seealso \code{\link{ferx_model_section}}, \code{\link{ferx_model_show}},
#'   \code{\link{ferx_model}}
#' @family model-editing
#' @export
ferx_model_set_section <- function(path, section, lines) {
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

  end        <- if (idx < length(hdr$positions)) hdr$positions[idx + 1L] - 1L else length(file_lines)
  tail_lines <- if (end < length(file_lines)) file_lines[seq.int(end + 1L, length(file_lines))] else character(0)

  writeLines(c(file_lines[seq_len(hdr$positions[idx])], lines, tail_lines), path)
  invisible(path)
}
