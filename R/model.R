#' Display a ferx model file in the console
#'
#' Prints the contents of a \code{.ferx} model file to the console.
#'
#' @param path Path to a \code{.ferx} model file.
#'
#' @return \code{path}, invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_model_show(ex$model)
#'
#' @seealso \code{\link{ferx_model_edit}}, \code{\link{ferx_example}}
#' @export
ferx_model_show <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (tools::file_ext(path) != "ferx") stop("'path' must be a .ferx file")
  lines <- readLines(path, warn = FALSE)
  cat("# model:", basename(path), "\n")
  cat(lines, sep = "\n")
  cat("\n")
  invisible(path)
}

#' Open a ferx model file in an editor
#'
#' Opens a \code{.ferx} model file for editing. If the file lives inside the
#' installed ferx package directory (i.e. a bundled read-only example), a copy
#' is written to \code{dest} first and that copy is opened instead.
#'
#' @param path Path to a \code{.ferx} model file.
#' @param dest Directory to copy read-only package files into before editing.
#'   Defaults to the current working directory. Ignored when \code{path} is
#'   already a writable user-owned file.
#'
#' @return The path of the file that was opened (i.e. \code{path} for
#'   user-owned files, or the copied path for package examples), invisibly.
#'
#' @examples
#' \dontrun{
#' ex <- ferx_example("warfarin")
#'
#' # Edit a copy of the bundled example in a temp directory
#' my_model <- ferx_model_edit(ex$model, dest = tempdir())
#'
#' # Edit a user-owned model directly
#' ferx_model_edit("my_model.ferx")
#' }
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_example}}
#' @export
ferx_model_edit <- function(path, dest = ".") {
  if (!file.exists(path)) stop("File not found: ", path)

  pkg_dir <- system.file("", package = "ferx")
  in_pkg  <- nzchar(pkg_dir) && startsWith(normalizePath(path), normalizePath(pkg_dir))

  if (in_pkg) {
    dest_path <- file.path(dest, basename(path))
    file.copy(path, dest_path, overwrite = FALSE)
    message("Copied to ", dest_path, "; editing your copy.")
    path <- dest_path
  }

  utils::file.edit(path)
  invisible(path)
}
