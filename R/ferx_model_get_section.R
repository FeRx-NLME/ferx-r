#' Get a section from a ferx model
#'
#' Returns the lines belonging to a named section of a \code{.ferx} model,
#' excluding the section header itself. Prints the lines to the console and
#' returns them invisibly.
#'
#' @param x A \code{ferx_model} object or a path to a \code{.ferx} file.
#' @param section Name of the section to extract, without brackets (e.g.
#'   \code{"parameters"}). An instance-named section is addressed as its
#'   header reads, type then instance: \code{"event_model cause_a"} for
#'   \code{[event_model cause_a]}. A bare type such as \code{"event_model"}
#'   addresses the unnamed \code{[event_model]} header when the model has
#'   one, even beside named instances of it. Otherwise it addresses the
#'   type's only instance (\code{"covariate_nn"} for
#'   \code{[covariate_nn TYPICAL_PK]}), and is refused, with the instances to
#'   choose from, when there are several. The
#'   type matches case-insensitively and the instance name exactly, as the
#'   engine reads them.
#' @param strip Logical. If \code{TRUE}, leading whitespace is trimmed from each
#'   returned line via \code{\link[base]{trimws}}. Defaults to \code{FALSE} to
#'   preserve the round-trip guarantee with
#'   \code{\link{ferx_model_set_section}}.
#'
#' @return Character vector of lines in the requested section, invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_model_get_section(ex$model, "parameters")
#' ferx_model_get_section(ex$model, "parameters", strip = TRUE)
#'
#' m <- ferx_model(ex$data, ex$model)
#' ferx_model_get_section(m, "parameters")
#'
#' # An instance-named section
#' tte <- ferx_example("tte_competing_risks")
#' ferx_model_get_section(tte$model, "event_model cause_a")
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model_edit}},
#'   \code{\link{ferx_model_set_section}}
#' @family model-editing
#' @export
ferx_model_get_section <- function(x, section, strip = FALSE) {
  path <- if (inherits(x, "ferx_model")) x$model else x
  body <- .ferx_read_section(path, section, strip = strip)

  cat("# [", attr(body, "section"), "]\n", sep = "")
  attr(body, "section") <- NULL
  cat(body, sep = "\n")
  cat("\n")
  invisible(body)
}

# Low-level section reader backing ferx_model_get_section(). Assumes `path`
# is already resolved from a ferx_model object if needed. The header the
# lines came from (`type` or `type instance`) rides along as attribute
# "section", so a bare-type request prints the header it resolved to.
.ferx_read_section <- function(path, section, strip = FALSE) {
  .ferx_validate_ferx_path(path)

  file_lines <- readLines(path, warn = FALSE)
  hdr        <- ferx_section_headers(file_lines)
  idx        <- .ferx_section_index(hdr, section)

  start <- hdr$positions[idx] + 1L
  end   <- if (idx < length(hdr$positions)) hdr$positions[idx + 1L] - 1L else length(file_lines)
  body  <- if (start <= end) file_lines[start:end] else character(0)
  if (strip) body <- trimws(body, which = "left")
  attr(body, "section") <- hdr$names[idx]
  body
}
