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
#' @return Invisibly returns a list with \code{ok} (logical - FALSE when the
#'   engine reports an error, a required section is missing, or the file
#'   carries a section name this build of the engine does not accept - one it
#'   never knew, one it has retired, or one behind a disabled feature),
#'   \code{model},
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
  # Optional sections come from the engine (ferx-core's `known_block_names()`),
  # never a list maintained here. The duplicated vector this replaces had
  # drifted: it omitted covariates, event_model, binary_model, markov_model,
  # data_selection, adaptive_dosing, mixture, data and simulation - all of them
  # blocks the parser reads and all of them used by bundled examples, so
  # `ferx_model_validate(ferx_example("two_cpt_oral_cov")$model)` reported
  # `covariates [unknown section]` - and it still listed `initial_values`,
  # which the engine dropped in ferx-core e5e934d. See ferx-core #1040.
  optional_sections <- setdiff(ferx_rust_known_blocks(), required_sections)

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

  # An unrecognised section counts against `ok`. It used to be printed as
  # `[unknown section]` and then left out of the returned status, so `res$ok`
  # was TRUE for a model carrying a block the engine would ignore - the exact
  # silent-drop ferx-core #1040 closes. A current engine already errors on it
  # (`E_UNKNOWN_BLOCK`), which is what `rust_result$ok` carries; folding it in
  # here keeps the status honest against an older pinned engine too.
  ok <- isTRUE(rust_result$ok) && length(missing) == 0L && length(unknown) == 0L

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
  # `ferx_rust_known_blocks()` is build-dependent and omits names the engine
  # still recognises: a retired block (`E_DEPRECATED_BLOCK`) and one gated
  # behind a cargo feature this binary lacks (`E_BLOCK_FEATURE_DISABLED`) both
  # land in `unknown`. Printing `[unknown section]` for those contradicts the
  # engine's own diagnostic two lines further down, which correctly says the
  # name was retired or needs a feature flag - so take the label from that
  # diagnostic when the engine has already named the block. `[unknown section]`
  # is left for a header nothing explains.
  if (length(unknown) > 0L) {
    for (s in unknown) {
      codes <- diag$code[!is.na(diag$block) & diag$block == s]
      label <- if ("E_DEPRECATED_BLOCK" %in% codes) {
        "[retired section]"
      } else if ("E_BLOCK_FEATURE_DISABLED" %in% codes) {
        "[feature not enabled]"
      } else {
        "[unknown section]"
      }
      cat(sprintf("  %-30s %s\n", s, label))
    }
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
