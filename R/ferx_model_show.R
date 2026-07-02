#' Display a ferx model file in the console
#'
#' Prints the contents of a \code{.ferx} model file to the console. When the
#' console supports ANSI colour (and the \pkg{cli} package is installed), the
#' output is syntax-highlighted: section headers (\code{[parameters]}, ...) are
#' bold yellow, declaration keywords (\code{theta}, \code{omega}, \code{sigma},
#' \code{kappa}, ...) are cyan, and comments are dimmed. In a non-colour context
#' (file, pipe, \code{NO_COLOR}, or no \pkg{cli}) the raw text is printed
#' unchanged.
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
#' @family model-editing
#' @export
ferx_model_show <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")
  lines <- readLines(path, warn = FALSE)
  cat("# model:", basename(path), "\n")
  if (.ferx_use_color()) {
    lines <- unname(vapply(lines, .ferx_highlight_line, character(1L)))
  }
  cat(lines, sep = "\n")
  cat("\n")
  invisible(path)
}

# TRUE when the console can render ANSI colour and cli is available. Honours
# NO_COLOR, non-tty output, etc. via cli::num_ansi_colors().
.ferx_use_color <- function() {
  requireNamespace("cli", quietly = TRUE) && cli::num_ansi_colors() > 1L
}

# Syntax-highlight one line of a .ferx model file. Only called when colour is
# active (so cli is guaranteed present). Plain text is always recoverable via
# cli::ansi_strip().
.ferx_highlight_line <- function(line) {
  # Section header on its own line, e.g. "[parameters]".
  if (grepl("^\\s*\\[[^]]+\\]\\s*$", line)) {
    return(cli::style_bold(cli::col_yellow(line)))
  }
  # Peel off a trailing comment (from the first '#') and dim it.
  comment <- ""
  hash <- regexpr("#", line, fixed = TRUE)
  if (hash > 0L) {
    comment <- cli::style_dim(substr(line, hash, nchar(line)))
    line <- substr(line, 1L, hash - 1L)
  }
  # Colour a leading declaration keyword.
  keywords <- c("theta", "omega", "sigma", "kappa", "block_omega",
                "pk", "ode", "odes")
  g <- regmatches(line, regexec("^(\\s*)(\\S+)(.*)$", line))[[1L]]
  if (length(g) == 4L && g[[3L]] %in% keywords) {
    line <- paste0(g[[2L]], cli::col_cyan(g[[3L]]), g[[4L]])
  }
  paste0(line, comment)
}
