#' Set a section in a ferx model
#'
#' Replaces the body of a named section in a \code{.ferx} model with new
#' lines, leaving all other sections untouched, and returns \code{x} so
#' calls can be chained with \code{|>} or \code{\%>\%}.
#'
#' When \code{x} is a \code{ferx_model} object the section is written to
#' \code{x$model}. \strong{Copy-on-write for bundled examples}: if
#' \code{x$model} points to a file inside the installed \code{ferx} package
#' directory (typically the result of \code{\link{ferx_example}()}), the
#' file is first copied to \code{tempdir()} and \code{x$model} is updated to
#' that copy before the edit. This prevents accidental modification of
#' read-only package examples. When \code{x} is a plain path string the file
#' is edited in place.
#'
#' @param x A \code{ferx_model} object or a path to a \code{.ferx} file.
#' @param section Name of the section to replace, without brackets (e.g.
#'   \code{"fit_options"}).
#' @param lines Character vector of replacement lines (do not include the
#'   \code{[section]} header line).
#'
#' @return \code{x}, invisibly. For a \code{ferx_model} pointing at a
#'   bundled package file, the returned object's \code{$model} field will
#'   point at the temp-directory copy that was actually edited.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Plain path: edited in place
#' ferx_model_set_section(ex$model, "fit_options", c(
#'   "  method     = focei",
#'   "  maxiter    = 500",
#'   "  covariance = true"
#' ))
#'
#' # ferx_model: safe by default. Editing a bundled example transparently
#' # copies the file to tempdir() and edits the copy; ex$model on disk is
#' # untouched.
#' fit <- ferx_model(ex$data, ex$model) |>
#'   ferx_model_set_section("fit_options", c(
#'     "  method     = focei",
#'     "  maxiter    = 500",
#'     "  covariance = true"
#'   )) |>
#'   ferx_fit() |>
#'   summary()
#' }
#'
#' @seealso \code{\link{ferx_model_get_section}}
#' @family model-editing
#' @export
ferx_model_set_section <- function(x, section, lines) {
  if (inherits(x, "ferx_model")) {
    x$model <- .ferx_copy_if_in_pkg(x$model)
    .ferx_write_section(x$model, section, lines)
    invisible(x)
  } else {
    .ferx_write_section(x, section, lines)
    invisible(x)
  }
}

# Low-level section writer shared by ferx_model_set_section()'s two input
# branches. Assumes `path` is already resolved (and copy-on-write applied,
# if applicable) from a ferx_model object if needed.
.ferx_write_section <- function(path, section, lines) {
  .ferx_validate_ferx_path(path)

  file_lines <- readLines(path, warn = FALSE)
  hdr        <- ferx_section_headers(file_lines)
  idx        <- .ferx_section_index(hdr, section)

  end        <- if (idx < length(hdr$positions)) hdr$positions[idx + 1L] - 1L else length(file_lines)
  tail_lines <- if (end < length(file_lines)) file_lines[seq.int(end + 1L, length(file_lines))] else character(0)

  writeLines(c(file_lines[seq_len(hdr$positions[idx])], lines, tail_lines), path)
  invisible(path)
}

# Copy-on-write guard: if `path` is inside the installed ferx package
# directory, copy it to `tempdir()` and return the new path so subsequent
# in-place edits do not mutate the installed package files. Otherwise
# return `path` unchanged.
.ferx_copy_if_in_pkg <- function(path) {
  pkg_dir <- system.file("", package = "ferx")
  if (!nzchar(pkg_dir)) return(path)
  pkg_norm  <- normalizePath(pkg_dir, mustWork = FALSE)
  path_norm <- normalizePath(path,    mustWork = FALSE)
  if (!startsWith(path_norm, pkg_norm)) return(path)

  # A fresh tempfile() per call (rather than a fixed file.path(tempdir(),
  # basename(path))) avoids two independent ferx_model objects wrapping the
  # same bundled file silently clobbering each other's copy.
  dest <- tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(path)), "_"),
    fileext = paste0(".", tools::file_ext(path))
  )
  if (!isTRUE(file.copy(path, dest, overwrite = TRUE))) {
    stop("Failed to copy ", path, " to ", dest)
  }
  message(
    "ferx_model_set_section: copying read-only package model to ", dest,
    " before editing."
  )
  dest
}
