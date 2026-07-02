#' Open a ferx model file in an editor
#'
#' Opens a \code{.ferx} model file for editing. If the file lives inside the
#' installed ferx package directory (i.e. a bundled read-only example), a copy
#' is written to \code{dest} first and that copy is opened instead.
#'
#' After the editor closes, an optional \code{save_as} step lets you copy the
#' edited file to a new path - useful when iterating on a model to keep
#' versioned copies.
#'
#' @param path Path to a \code{.ferx} model file.
#' @param dest Directory to copy read-only package files into before editing.
#'   Defaults to the current working directory. Ignored when \code{path} is
#'   already a writable user-owned file.
#' @param overwrite Logical. If \code{TRUE}, overwrite an existing file in
#'   \code{dest} when copying a package example. If \code{FALSE} (default) and
#'   the destination file already exists, an error is raised.
#' @param save_as Controls post-edit save-as behaviour:
#'   \itemize{
#'     \item \code{NULL} or \code{FALSE} (default) - no extra action after
#'       editing. Accepting \code{FALSE} lets you pass expressions like
#'       \code{save_as = interactive()}.
#'     \item \code{TRUE} - interactively prompt the user for a destination path.
#'     \item A character string - silently copy the edited file to that path.
#'   }
#'   When a copy is made the \emph{copy} path is returned; otherwise the edited
#'   file path is returned.
#' @param .editor Function used to open the file. Defaults to
#'   \code{utils::file.edit}. Override in tests or non-interactive contexts to
#'   replace or suppress the editor call.
#'
#' @return The path of the file in its final location, invisibly. When
#'   \code{save_as} produces a copy, that copy's path is returned; otherwise
#'   the path of the edited file is returned.
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
#'
#' # After editing, save a copy to a new versioned path
#' ferx_model_edit("run1.ferx", save_as = "run2.ferx")
#'
#' # After editing, interactively ask for the destination
#' ferx_model_edit("run1.ferx", save_as = TRUE)
#' }
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model}},
#'   \code{\link{ferx_example}}
#' @family model-editing
#' @export
ferx_model_edit <- function(path, dest = ".", overwrite = FALSE, save_as = NULL,
                            .editor = utils::file.edit) {
  if (!file.exists(path)) stop("File not found: ", path)

  # Validate save_as up front so invalid input fails fast - before opening an
  # editor or copying any files. FALSE collapses to NULL so callers can pass
  # `save_as = interactive()` and have it no-op in non-interactive sessions.
  if (isFALSE(save_as)) save_as <- NULL
  if (!is.null(save_as) && !isTRUE(save_as) &&
      !(is.character(save_as) && length(save_as) == 1L && !is.na(save_as))) {
    stop("'save_as' must be NULL, TRUE, FALSE, or a single character string.")
  }

  pkg_dir <- system.file("", package = "ferx")
  in_pkg  <- nzchar(pkg_dir) && startsWith(normalizePath(path), normalizePath(pkg_dir))

  if (in_pkg) {
    dest_path <- file.path(dest, basename(path))
    if (file.exists(dest_path) && !overwrite) {
      stop(
        dest_path, " already exists. ",
        "Use overwrite = TRUE to replace it, or choose a different dest."
      )
    }
    if (!isTRUE(file.copy(path, dest_path, overwrite = overwrite))) {
      stop("Failed to copy ", path, " to ", dest_path)
    }
    message("Copied to ", dest_path, "; editing your copy.")
    path <- dest_path
  }

  .editor(path)

  if (is.null(save_as)) return(invisible(path))

  if (isTRUE(save_as)) {
    save_as <- trimws(readline(
      prompt = paste0("Save a copy of '", basename(path), "' to: ")
    ))
    if (!nzchar(save_as)) {
      message("No path entered; keeping edited file at ", path)
      return(invisible(path))
    }
  }

  if (file.exists(save_as) && !overwrite) {
    stop(
      save_as, " already exists. ",
      "Use overwrite = TRUE to replace it."
    )
  }
  if (!isTRUE(file.copy(path, save_as, overwrite = overwrite))) {
    stop("Failed to copy ", path, " to ", save_as)
  }
  message("Saved copy to ", save_as)
  invisible(save_as)
}
