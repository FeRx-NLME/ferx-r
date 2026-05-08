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
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")
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
#' @param overwrite Logical. If \code{TRUE}, overwrite an existing file in
#'   \code{dest} when copying a package example. If \code{FALSE} (default) and
#'   the destination file already exists, an error is raised.
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
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model_new}},
#'   \code{\link{ferx_example}}
#' @export
ferx_model_edit <- function(path, dest = ".", overwrite = FALSE) {
  if (!file.exists(path)) stop("File not found: ", path)

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
    file.copy(path, dest_path, overwrite = overwrite)
    message("Copied to ", dest_path, "; editing your copy.")
    path <- dest_path
  }

  utils::file.edit(path)
  invisible(path)
}

# Parse section header positions and names from a character vector of file
# lines. Returns a list with $positions (integer indices) and $names (strings).
ferx_section_headers <- function(lines) {
  pos   <- grep("^\\s*\\[.+\\]\\s*$", lines)
  names <- gsub("^\\s*\\[|\\]\\s*$", "", lines[pos])
  list(positions = pos, names = names)
}


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
#'   \code{\link{ferx_model_new}}
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

  start      <- hdr$positions[idx] + 1L
  end        <- if (idx < length(hdr$positions)) hdr$positions[idx + 1L] - 1L else length(file_lines)
  tail_lines <- if (end < length(file_lines)) file_lines[seq.int(end + 1L, length(file_lines))] else character(0)

  writeLines(c(file_lines[seq_len(hdr$positions[idx])], lines, tail_lines), path)
  invisible(path)
}

#' Create a new ferx model file from a skeleton template
#'
#' Writes a new \code{.ferx} file pre-filled with a skeleton for the chosen
#' model type, then opens it in an editor. Pass \code{print = TRUE} instead of
#' supplying a \code{path} to print the skeleton to the console without writing
#' any file — useful for copy-pasting or piping in a scripted workflow.
#'
#' @param path Path for the new \code{.ferx} file. Must not already exist
#'   unless \code{overwrite = TRUE}. Ignored when \code{print = TRUE}.
#' @param template One of \code{"1cpt_oral"} (default), \code{"1cpt_iv"},
#'   \code{"2cpt_oral"}, \code{"2cpt_iv"}, or \code{"ode"}.
#' @param overwrite Logical. Overwrite \code{path} if it already exists?
#'   Defaults to \code{FALSE}.
#' @param edit Logical. Open the file in an editor after writing? Defaults to
#'   \code{TRUE}. Ignored when \code{print = TRUE}.
#' @param print Logical. If \code{TRUE}, print the skeleton to the console
#'   instead of writing a file. No file is created and no editor is opened.
#'   Defaults to \code{FALSE}.
#'
#' @return \code{path} invisibly, or \code{NULL} invisibly when
#'   \code{print = TRUE}.
#'
#' @examples
#' # Print a skeleton to the console (no file written, no editor opened)
#' ferx_model_new(print = TRUE)
#' ferx_model_new(template = "2cpt_iv", print = TRUE)
#'
#' \dontrun{
#' # Write a file and open it for editing
#' ferx_model_new("my_model.ferx")
#'
#' # Write a file without opening an editor
#' ferx_model_new("my_model.ferx", edit = FALSE)
#' }
#'
#' @seealso \code{\link{ferx_model_edit}}, \code{\link{ferx_model_show}},
#'   \code{\link{ferx_model_set_section}}
#' @export
ferx_model_new <- function(path = NULL, template = "1cpt_oral",
                            overwrite = FALSE, edit = TRUE, print = FALSE) {
  if (!print) {
    if (is.null(path)) stop("'path' is required when print = FALSE")
    if (tools::file_ext(path) != "ferx") stop("'path' must end in .ferx")
    if (file.exists(path) && !overwrite) {
      stop(path, " already exists. Use overwrite = TRUE to replace it.")
    }
  }

  templates <- list(
    `1cpt_oral` = c(
      "# One-compartment oral PK model",
      "",
      "[parameters]",
      "  theta TVCL(1.0, 0.001, 100.0)",
      "  theta TVV(10.0, 0.1, 1000.0)",
      "  theta TVKA(1.0, 0.01, 50.0)",
      "",
      "  omega ETA_CL ~ 0.09",
      "  omega ETA_V  ~ 0.09",
      "  omega ETA_KA ~ 0.25",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "  KA = TVKA * exp(ETA_KA)",
      "",
      "[structural_model]",
      "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[initial_values]",
      "  theta = [1.0, 10.0, 1.0]",
      "  omega = [0.09, 0.09, 0.25]",
      "  sigma = [0.01]",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 300",
      "  covariance = true"
    ),
    `1cpt_iv` = c(
      "# One-compartment IV bolus PK model",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV(20.0, 1.0, 500.0)",
      "",
      "  omega ETA_CL ~ 0.09",
      "  omega ETA_V  ~ 0.09",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "",
      "[structural_model]",
      "  pk one_cpt_iv_bolus(cl=CL, v=V)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[initial_values]",
      "  theta = [5.0, 20.0]",
      "  omega = [0.09, 0.09]",
      "  sigma = [0.01]",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 300",
      "  covariance = true"
    ),
    `2cpt_oral` = c(
      "# Two-compartment oral PK model",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV1(50.0, 1.0, 500.0)",
      "  theta TVQ(10.0, 0.1, 100.0)",
      "  theta TVV2(100.0, 1.0, 1000.0)",
      "  theta TVKA(1.2, 0.01, 10.0)",
      "",
      "  omega ETA_CL ~ 0.10",
      "  omega ETA_V1 ~ 0.10",
      "  omega ETA_Q  ~ 0.05",
      "  omega ETA_V2 ~ 0.05",
      "  omega ETA_KA ~ 0.15",
      "",
      "  sigma PROP_ERR ~ 0.02",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V1 = TVV1 * exp(ETA_V1)",
      "  Q  = TVQ  * exp(ETA_Q)",
      "  V2 = TVV2 * exp(ETA_V2)",
      "  KA = TVKA * exp(ETA_KA)",
      "",
      "[structural_model]",
      "  pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2, ka=KA)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[initial_values]",
      "  theta = [5.0, 50.0, 10.0, 100.0, 1.2]",
      "  omega = [0.10, 0.10, 0.05, 0.05, 0.15]",
      "  sigma = [0.02]",
      "",
      "[fit_options]",
      "  method     = focei",
      "  maxiter    = 500",
      "  covariance = true"
    ),
    `2cpt_iv` = c(
      "# Two-compartment IV bolus PK model",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV1(15.0, 1.0, 500.0)",
      "  theta TVQ(3.0, 0.01, 100.0)",
      "  theta TVV2(30.0, 1.0, 500.0)",
      "",
      "  omega ETA_CL ~ 0.10",
      "  omega ETA_V1 ~ 0.10",
      "  omega ETA_Q  ~ 0.10",
      "  omega ETA_V2 ~ 0.10",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V1 = TVV1 * exp(ETA_V1)",
      "  Q  = TVQ  * exp(ETA_Q)",
      "  V2 = TVV2 * exp(ETA_V2)",
      "",
      "[structural_model]",
      "  pk two_cpt_iv_bolus(cl=CL, v1=V1, q=Q, v2=V2)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[initial_values]",
      "  theta = [5.0, 15.0, 3.0, 30.0]",
      "  omega = [0.10, 0.10, 0.10, 0.10]",
      "  sigma = [0.01]",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 500",
      "  covariance = true"
    ),
    ode = c(
      "# ODE-based PK model",
      "",
      "[parameters]",
      "  theta TVPARAM(1.0, 0.001, 1000.0)",
      "",
      "  omega ETA_PARAM ~ 0.09",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  PARAM = TVPARAM * exp(ETA_PARAM)",
      "",
      "[structural_model]",
      "  ode(obs_cmt=central, states=[depot, central])",
      "",
      "[odes]",
      "  d/dt(depot)   = 0  # replace with your equations",
      "  d/dt(central) = 0",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[initial_values]",
      "  theta = [1.0]",
      "  omega = [0.09]",
      "  sigma = [0.01]",
      "",
      "[fit_options]",
      "  method     = focei",
      "  maxiter    = 500",
      "  covariance = true"
    )
  )

  valid <- names(templates)
  if (!template %in% valid) {
    stop(
      "Unknown template '", template, "'. ",
      "Choose one of: ", paste(valid, collapse = ", ")
    )
  }

  skeleton <- templates[[template]]

  if (print) {
    cat(skeleton, sep = "\n")
    cat("\n")
    return(invisible(NULL))
  }

  writeLines(skeleton, path)
  message("Created ", path)

  if (edit) utils::file.edit(path)
  invisible(path)
}

# Extract all named [section] blocks from a .ferx file.
# Returns a named list: section name → character vector of (comment-stripped, trimmed) lines.
.ferx_extract_blocks <- function(path) {
  raw <- readLines(path, warn = FALSE)
  blocks  <- list()
  current <- NULL
  for (line in raw) {
    stripped <- trimws(sub("(#|//).*$", "", line))
    if (!nzchar(stripped)) next
    m <- regmatches(stripped, regexpr("^\\[(\\w+)\\]$", stripped, perl = TRUE))
    if (length(m) > 0L) {
      current <- tolower(gsub("^\\[|\\]$", "", m))
      if (is.null(blocks[[current]])) blocks[[current]] <- character(0)
      next
    }
    if (!is.null(current)) blocks[[current]] <- c(blocks[[current]], stripped)
  }
  blocks
}

# Detect an unambiguous model-type label from [structural_model] lines.
# Returns a short string ("1-cpt oral", "ODE", etc.) or NULL when unrecognised.
# Labels and the recognised function names are kept in sync with the Rust
# parser (ferx-nlme src/parser/model_parser.rs `pk_func_name` match arms) so
# pre-fit `ferx_model_inspect(path)` reports the same string the engine would
# attach to a fitted result.
.ferx_model_type <- function(lines) {
  s <- paste(lines, collapse = " ")
  if (grepl("\\bode\\(", s, perl = TRUE)) return("ODE")

  m <- regmatches(s, regexpr("\\bpk\\s+([A-Za-z_][A-Za-z_0-9]*)\\s*\\(",
                             s, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  fn <- sub("^pk\\s+([A-Za-z_][A-Za-z_0-9]*)\\s*\\(.*", "\\1", m, perl = TRUE)

  # Long `*_compartment_*` aliases collapse to their `*_cpt_*` equivalents.
  fn <- sub("_compartment_", "_cpt_", fn, fixed = TRUE)

  switch(fn,
    one_cpt_iv_bolus    = "1-cpt IV bolus",
    one_cpt_infusion    = "1-cpt IV infusion",
    one_cpt_oral        = "1-cpt oral",
    two_cpt_iv_bolus    = "2-cpt IV bolus",
    two_cpt_infusion    = "2-cpt IV infusion",
    two_cpt_oral        = "2-cpt oral",
    three_cpt_iv_bolus  = "3-cpt IV bolus",
    three_cpt_infusion  = "3-cpt IV infusion",
    three_cpt_oral      = "3-cpt oral",
    NULL
  )
}

# Format the structural display string from a model_structure list.
# Combines the optional type label with the population parameter (theta) names.
.ferx_format_structural <- function(ms) {
  theta_str <- if (length(ms$theta_names) > 0L)
    paste(ms$theta_names, collapse = ", ")
  else
    NULL
  if (!is.null(ms$model_type) && !is.null(theta_str))
    sprintf("%s  (%s)", ms$model_type, theta_str)
  else if (!is.null(ms$model_type))
    ms$model_type
  else if (!is.null(theta_str))
    theta_str
  else
    "unknown"
}

# Print the four structural detail lines (Structural / IIV / IOV / Residual).
# The caller is responsible for any preceding header line.
.ferx_print_structure <- function(ms) {
  cat(sprintf("  Structural:  %s\n", .ferx_format_structural(ms)))
  cat(sprintf("  IIV:         %s\n",
    if (length(ms$iiv) > 0L) paste(ms$iiv, collapse = ", ") else "none"))
  cat(sprintf("  IOV:         %s\n",
    if (length(ms$iov) > 0L) paste(ms$iov, collapse = ", ") else "none"))
  cat(sprintf("  Residual:    %s\n", ms$residual))
  invisible(NULL)
}

# Parse a .ferx file and return a named list describing model structure.
# Fields: theta_names (pop param names), model_type (label or NULL),
#         iiv, iov, residual.
# Used by ferx_model_inspect() (pre-fit) and attached to ferx_fit() results.
.ferx_parse_structure <- function(path) {
  b <- .ferx_extract_blocks(path)

  # Population (theta) parameter names from [parameters]
  params      <- b[["parameters"]] %||% character(0)
  theta_lines <- grep("^theta\\s", params, value = TRUE)
  thetas <- if (length(theta_lines) > 0L)
    sub("^theta\\s+(\\w+).*", "\\1", theta_lines)
  else
    character(0)

  # Optional model-type label from [structural_model]
  struct_lines <- b[["structural_model"]] %||% character(0)
  model_type   <- if (length(struct_lines) > 0L) .ferx_model_type(struct_lines) else NULL

  # IIV: omega lines
  omega_lines <- grep("^omega\\s", params, value = TRUE)
  iiv <- if (length(omega_lines) > 0L)
    sub("^omega\\s+(\\w+).*", "\\1", omega_lines)
  else
    character(0)

  # IOV: kappa lines
  kappa_lines <- grep("^kappa\\s", params, value = TRUE)
  iov <- if (length(kappa_lines) > 0L)
    sub("^kappa\\s+(\\w+).*", "\\1", kappa_lines)
  else
    character(0)

  # Residual error type from [error_model]
  err_line <- (b[["error_model"]] %||% character(0))[1L]
  residual <- if (is.na(err_line) || !nzchar(err_line %||% "")) {
    "unknown"
  } else if (grepl("proportional", err_line, ignore.case = TRUE)) {
    "proportional"
  } else if (grepl("additive",     err_line, ignore.case = TRUE)) {
    "additive"
  } else if (grepl("combined",     err_line, ignore.case = TRUE)) {
    "combined"
  } else {
    warning("Unrecognised residual error type; reporting as \"unknown\". Line: ",
            err_line, call. = FALSE)
    "unknown"
  }

  list(
    theta_names = thetas,
    model_type  = model_type,
    iiv         = iiv,
    iov         = iov,
    residual    = residual
  )
}

#' Inspect the structure of a ferx model file
#'
#' Parses a \code{.ferx} file without fitting and prints a compact summary of
#' its model structure: PK model type, inter-individual variability (IIV),
#' inter-occasion variability (IOV), and residual error type. Useful for
#' verifying that a model file will be interpreted as expected before committing
#' to a potentially long estimation run.
#'
#' Alternatively, pass a \code{ferx_fit} object to display the structure that
#' was auto-derived during fitting (reads \code{fit$model_structure} directly,
#' so no file path is needed post-fit).
#'
#' @param path Path to a \code{.ferx} model file, \emph{or} a
#'   \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#'
#' @return A named list (invisibly) with fields: \code{theta_names}
#'   (character vector of population parameter names), \code{model_type}
#'   (short label such as \code{"1-cpt oral"} or \code{NULL} when not
#'   unambiguously detectable), \code{iiv} (omega names), \code{iov}
#'   (kappa names), and \code{residual} (error type).
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_model_inspect(ex$model)
#'
#' # Programmatic access
#' s <- ferx_model_inspect(ex$model)
#' s$theta_names  # c("TVCL", "TVV", "TVKA")
#' s$model_type   # "1-cpt oral"
#' s$iiv          # c("ETA_CL", "ETA_V", "ETA_KA")
#' s$residual     # "proportional"
#'
#' @seealso \code{\link{ferx_model_show}}, \code{\link{ferx_model_edit}},
#'   \code{\link{ferx_fit}}
#' @export
ferx_model_inspect <- function(path) {
  if (inherits(path, "ferx_fit")) {
    s <- path$model_structure
    if (is.null(s)) stop("No model_structure found on this ferx_fit object.")
    label <- if (!is.null(path$model_name) && nzchar(path$model_name))
      paste0(path$model_name, ".ferx") else "ferx_fit"
    cat("Model structure (", label, ")\n", sep = "")
    .ferx_print_structure(s)
    return(invisible(s))
  }
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")

  s <- .ferx_parse_structure(path)

  cat("Model structure (", basename(path), ")\n", sep = "")
  .ferx_print_structure(s)

  invisible(s)
}
