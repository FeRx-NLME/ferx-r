#' Validate a ferx model file for syntax and structural errors
#'
#' Parses a \code{.ferx} model file using the Rust engine and checks for
#' required sections, without running the optimizer. Useful for catching
#' syntax errors and missing sections before committing to a long estimation
#' run. This is the in-R counterpart to the \code{ferx check} CLI shipped
#' with ferx-core and shares the same parser and structural diagnostics.
#'
#' @param path Path to a \code{.ferx} model file. The file must exist and have
#'   a \code{.ferx} extension; otherwise an error is raised (these are caller
#'   errors, not validation failures).
#' @param data Optional path to a NONMEM-format CSV. When supplied, the engine
#'   additionally runs data-dependent checks (covariate columns present,
#'   per-CMT scaling/error-model coverage, steady-state II sanity, lagtime
#'   signs). \code{NULL} runs the model-only checks.
#'
#' @return Invisibly returns a list with \code{ok} (logical), \code{model},
#'   \code{data}, and a \code{diagnostics} data frame with one row per finding
#'   (\code{severity}, \code{code}, \code{message}, \code{block}, \code{line},
#'   \code{suggestion}). The function always prints a report to the console.
#'   Codes are stable identifiers (\code{E_*} for errors, \code{W_*} for
#'   warnings) suitable for programmatic handling - the registry lives in
#'   ferx-core's \code{docs/src/file-formats/check-report.md}. A missing file
#'   or non-\code{.ferx} extension raises an error rather than returning a
#'   failed diagnostic.
#'
#' @examples
#' # Valid model
#' ex <- ferx_example("warfarin")
#' ferx_model_validate(ex$model)
#'
#' # With data-dependent checks
#' ferx_model_validate(ex$model, data = ex$data)
#'
#' # Inspect findings programmatically
#' res <- ferx_model_validate(ex$model)
#' res$ok
#' res$diagnostics
#'
#' \dontrun{
#' # Invalid model (missing required sections)
#' bad <- tempfile(fileext = ".ferx")
#' writeLines(c(
#'   "[parameters]",
#'   "  theta TVCL(1.0, 0.001, 100.0)",
#'   "[structural_model]",
#'   "  pk one_cpt_oral(cl=CL, v=V, ka=KA)"
#' ), bad)
#' ferx_model_validate(bad)
#' # Validating: <file>.ferx
#' #
#' # Sections present:
#' #   parameters                     [ok]
#' #   individual_parameters          [MISSING]
#' #   structural_model               [ok]
#' #   error_model                    [MISSING]
#' #
#' # Result: INVALID
#' #   * Missing required section: [individual_parameters]
#' #   * Missing required section: [error_model]
#' }
#'
#' @seealso \code{\link{ferx_model_inspect}}, \code{\link{ferx_model_show}}
#' @family model-editing
#' @export
ferx_model_validate <- function(path, data = NULL) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")
  if (!is.null(data)) {
    if (!is.character(data) || length(data) != 1L) {
      stop("'data' must be a single file path or NULL")
    }
    if (!file.exists(data)) stop("Data file not found: ", data)
  }

  required_sections <- c(
    "parameters", "individual_parameters", "structural_model",
    "error_model"
  )
  optional_sections <- c("odes", "fit_options", "scaling", "initial_values",
                         "covariate_nn", "diffusion", "derived", "output")

  blocks   <- .ferx_extract_blocks(path)
  present  <- names(blocks)
  missing  <- setdiff(required_sections, present)
  unknown  <- setdiff(present, c(required_sections, optional_sections))

  data_arg <- if (is.null(data)) "" else normalizePath(data)
  rust_result <- ferx_rust_validate_model(normalizePath(path), data_arg)

  diag <- data.frame(
    severity   = as.character(rust_result$severity),
    code       = as.character(rust_result$code),
    message    = as.character(rust_result$message),
    block      = ifelse(nzchar(rust_result$block), rust_result$block, NA_character_),
    line       = ifelse(rust_result$line == 0L, NA_integer_, as.integer(rust_result$line)),
    suggestion = ifelse(nzchar(rust_result$suggestion), rust_result$suggestion, NA_character_),
    stringsAsFactors = FALSE
  )

  ok <- isTRUE(rust_result$ok) && length(missing) == 0L

  cat("Validating:", basename(path), "\n")
  if (!is.null(data)) cat("       data:", basename(data), "\n")
  cat("\n")

  cat("Sections present:\n")
  for (s in required_sections) {
    status <- if (s %in% present) "[ok]" else "[MISSING]"
    cat(sprintf("  %-30s %s\n", s, status))
  }
  for (s in optional_sections) {
    if (s %in% present) cat(sprintf("  %-30s [ok] (optional)\n", s))
  }
  if (length(unknown) > 0L) {
    for (s in unknown) cat(sprintf("  %-30s [unknown section]\n", s))
  }
  cat("\n")

  if (ok && nrow(diag) == 0L) {
    cat("Result: VALID\n")
  } else if (ok) {
    cat("Result: VALID (with warnings)\n")
  } else {
    cat("Result: INVALID\n")
  }
  if (length(missing) > 0L) {
    for (s in missing) cat("  * Missing required section: [", s, "]\n", sep = "")
  }
  if (nrow(diag) > 0L) {
    for (i in seq_len(nrow(diag))) {
      tag <- if (diag$severity[i] == "error") "ERROR" else "warning"
      loc <- if (!is.na(diag$block[i])) {
        if (!is.na(diag$line[i])) sprintf(" [%s:%d]", diag$block[i], diag$line[i])
        else sprintf(" [%s]", diag$block[i])
      } else ""
      cat(sprintf("  * %s %s%s: %s\n", tag, diag$code[i], loc, diag$message[i]))
      if (!is.na(diag$suggestion[i])) cat("      hint:", diag$suggestion[i], "\n")
    }
  }

  invisible(list(
    ok = ok,
    model = as.character(rust_result$model),
    data = if (is.null(data)) NULL else as.character(rust_result$data),
    diagnostics = diag
  ))
}
