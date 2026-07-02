#' Replace a section in a ferx model (pipe-friendly)
#'
#' Updates the body of a named section in a \code{.ferx} model file and
#' returns the input object so calls can be chained with \code{|>} or
#' \code{\%>\%}.
#'
#' When \code{x} is a \code{ferx_model} object the section is written to
#' \code{x$model} and \code{x} is returned invisibly, allowing further
#' pipe steps. \strong{Copy-on-write for bundled examples}: if \code{x$model}
#' points to a file inside the installed \code{ferx} package directory
#' (typically the result of \code{\link{ferx_example}()}), the file is first
#' copied to \code{tempdir()} and \code{x$model} is updated to that copy
#' before the edit. This prevents accidental modification of read-only
#' package examples. When \code{x} is a plain path string the call is
#' equivalent to \code{\link{ferx_model_set_section}(x, section, lines)}
#' and the file is edited in place.
#'
#' @param x A \code{ferx_model} object or a path to a \code{.ferx} file.
#' @param section Name of the section to replace, without brackets.
#' @param lines Character vector of replacement lines (do not include the
#'   header line).
#'
#' @return \code{x}, invisibly. For a \code{ferx_model} pointing at a
#'   bundled package file, the returned object's \code{$model} field will
#'   point at the temp-directory copy that was actually edited.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Safe by default: editing a bundled example transparently copies the
#' # file to tempdir() and edits the copy. ex$model on disk is untouched.
#' fit <- ferx_model(ex$data, ex$model) |>
#'   ferx_set_section("fit_options", c(
#'     "  method     = focei",
#'     "  maxiter    = 500",
#'     "  covariance = true"
#'   )) |>
#'   ferx_fit() |>
#'   summary()
#' }
#'
#' @seealso \code{\link{ferx_model_set_section}}, \code{\link{ferx_get_section}}
#' @family model-editing
#' @export
ferx_set_section <- function(x, section, lines) {
  if (inherits(x, "ferx_model")) {
    x$model <- .ferx_copy_if_in_pkg(x$model)
    ferx_model_set_section(x$model, section, lines)
    invisible(x)
  } else {
    ferx_model_set_section(x, section, lines)
  }
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

  dest <- file.path(tempdir(), basename(path))
  if (!isTRUE(file.copy(path, dest, overwrite = TRUE))) {
    stop("Failed to copy ", path, " to ", dest)
  }
  message(
    "ferx_set_section: copying read-only package model to ", dest,
    " before editing."
  )
  dest
}
