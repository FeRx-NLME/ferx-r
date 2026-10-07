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
#'   engine reports an error, including a section this model's family requires
#'   (\code{E_MISSING_BLOCK}; a binary or time-to-event model needs no
#'   \code{[individual_parameters]}, \code{[structural_model]} or
#'   \code{[error_model]}), or the file
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
#' # Invalid model (no [error_model]: the engine names the first missing section)
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
#' #   structural_model               [ok]
#' #   error_model                    [MISSING]
#' #
#' # Result: INVALID
#' #   * ERROR E_MISSING_BLOCK [error_model]: Missing [error_model] block
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

  # Which sections exist, and which a model needs, both come from the engine.
  # The known names are ferx-core's `known_block_names()`; an R copy of that
  # list drifted (ferx-core #1040). Which ones are *required* depends on the
  # model family - a binary or time-to-event model has no
  # [individual_parameters], [structural_model] or [error_model] - and the
  # engine's parser says so itself with `E_MISSING_BLOCK`, naming the block.
  # A hard-coded required list here reported four valid bundled examples
  # INVALID (#306).
  known_sections <- ferx_rust_known_blocks()

  blocks   <- .ferx_extract_blocks(path)
  present  <- names(blocks)
  unknown  <- setdiff(present, known_sections)

  data_arg <- if (is.null(data)) "" else normalizePath(data)
  rust_result <- ferx_rust_validate_model(normalizePath(path), data_arg)

  diag <- .ferx_diagnostics_frame(rust_result)
  missing <- unique(diag$block[diag$code == "E_MISSING_BLOCK" & !is.na(diag$block)])

  # An unrecognised section counts against `ok`. It used to be printed as
  # `[unknown section]` and then left out of the returned status, so `res$ok`
  # was TRUE for a model carrying a block the engine would ignore - the exact
  # silent-drop ferx-core #1040 closes. A current engine already errors on it
  # (`E_UNKNOWN_BLOCK`), which is what `rust_result$ok` carries; folding it in
  # here keeps the status honest against an older pinned engine too.
  ok <- isTRUE(rust_result$ok) && length(unknown) == 0L

  cat("Validating:", basename(path), "\n")
  if (!is.null(data)) cat("       data:", basename(data), "\n")
  cat("\n")

  cat("Sections present:\n")
  for (s in intersect(present, known_sections)) {
    cat(sprintf("  %-30s [ok]\n", s))
  }
  for (s in setdiff(missing, present)) {
    cat(sprintf("  %-30s [MISSING]\n", s))
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

# -- Engine diagnostics as a data frame --------------------------------------

# Shape what `ferx_rust_validate_model()` returns (parallel character / integer
# vectors) into the `diagnostics` data frame `ferx_model_validate()` documents.
# Shared with `.ferx_engine_error()` so the code attached to a refused fit is
# the same identifier `ferx_model_validate()` prints for the same file.
.ferx_diagnostics_frame <- function(rust_result) {
  data.frame(
    severity   = as.character(rust_result$severity),
    code       = as.character(rust_result$code),
    message    = as.character(rust_result$message),
    block      = ifelse(nzchar(rust_result$block), rust_result$block, NA_character_),
    line       = ifelse(rust_result$line == 0L, NA_integer_, as.integer(rust_result$line)),
    suggestion = ifelse(nzchar(rust_result$suggestion), rust_result$suggestion, NA_character_),
    stringsAsFactors = FALSE
  )
}

# Re-raise an engine error from a fit with the stable diagnostic code attached.
#
# A refused fit used to carry only the engine's prose, while
# `ferx_model_validate()` on the same file returned `E_UNKNOWN_BLOCK` (or
# whichever code applied). A script could therefore branch on the code only on
# the path most users reach second (#367). This runs the engine's own
# validation over the same model / data, finds the diagnostic behind the
# failure, and returns a condition of class `ferx_engine_error` carrying
# `code`, `block`, `line` and `suggestion`, with the code appended to the
# message for whoever reads the console. The engine's prose is preserved
# verbatim at the front of the message, so existing `tryCatch()` /
# `expect_error()` matches on it keep working. Returns the original condition
# unchanged when no diagnostic can be tied to the failure.
#
# Only the failure path pays for the extra validation pass; a fit that runs
# never calls this. Nor does a refusal that already carries ferx-core's code
# (#498, `.ferx_engine_coded_error()`): the validation pass is the fallback for
# the ones that do not - `ferx_fit()`'s, and an uncoded one from elsewhere.
#
# `fallback_stages`: the single-error fallback below labels a failure with the
# one error validation found, without a text match. `NULL` allows it for any
# message, which is sound where every failure is a failure of the model or the
# data (`ferx_fit()`). Otherwise a regular expression naming the messages it
# may apply to; see `.ferx_engine_call()`.
.ferx_engine_error <- function(e, model, data, fallback_stages = NULL) {
  msg  <- conditionMessage(e)
  coded <- .ferx_engine_coded_error(e, msg)
  if (!is.null(coded)) return(coded)
  diag <- tryCatch(
    .ferx_diagnostics_frame(ferx_rust_validate_model(
      normalizePath(model),
      if (is.null(data)) "" else normalizePath(data)
    )),
    error = function(...) NULL
  )
  if (is.null(diag) || nrow(diag) == 0L) return(e)
  errs <- diag[diag$severity == "error", , drop = FALSE]
  if (nrow(errs) == 0L) return(e)

  # Prefer the diagnostic whose message the engine actually raised: both paths
  # render the same prose, so a substring match identifies the finding exactly.
  # Only when nothing matches do we fall back, and then only if validation
  # found a single error - with several, guessing which one stopped the fit
  # would risk labelling the failure with an unrelated code.
  #
  # Byte-wise: the diagnostic arrives marked UTF-8 and the condition message as
  # unmarked bytes, so outside a UTF-8 locale a message with a non-ASCII
  # character in it (the engine writes em-dashes) never matched itself.
  hit <- which(vapply(
    errs$message,
    function(m) nzchar(m) && grepl(m, msg, fixed = TRUE, useBytes = TRUE),
    logical(1)
  ))
  fallback_ok <- is.null(fallback_stages) || grepl(fallback_stages, msg, useBytes = TRUE)
  i <- if (length(hit) > 0L) {
    hit[1L]
  } else if (fallback_ok && nrow(errs) == 1L) {
    1L
  } else {
    return(e)
  }

  structure(
    class = c("ferx_engine_error", "error", "condition"),
    list(
      message    = sprintf("%s [%s]", msg, errs$code[i]),
      call       = conditionCall(e),
      code       = errs$code[i],
      block      = errs$block[i],
      line       = errs$line[i],
      suggestion = errs$suggestion[i]
    )
  )
}

# The condition for a refusal ferx-core raised with its own diagnostic code
# (ferx-r #498), or NULL when it carried none and `.ferx_engine_error()` has to
# re-validate to find one. Since ferx-core #1746 the prediction, simulation,
# NPDE, SIR and covariance entry points return the `Diagnostic` `ferx check`
# reports for the same refusal; the glue keeps it beside the text it raised.
# The record is used only for the condition whose message is that text, byte
# for byte, so a refusal the glue raised without a code never borrows the
# previous one's. `message` is the engine's message without its suggestion,
# which is the `suggestion` field's alone, so the advice is never shown twice.
.ferx_engine_coded_error <- function(e, msg) {
  d <- tryCatch(ferx_rust_take_engine_diagnostic(), error = function(...) NULL)
  if (is.null(d) || !identical(charToRaw(d$text), charToRaw(msg))) return(NULL)
  structure(
    class = c("ferx_engine_error", "error", "condition"),
    list(
      message    = sprintf("%s [%s]", d$message, d$code),
      call       = conditionCall(e),
      code       = d$code,
      block      = if (nzchar(d$block)) d$block else NA_character_,
      line       = if (d$line == 0L) NA_integer_ else d$line,
      suggestion = if (nzchar(d$suggestion)) d$suggestion else NA_character_
    )
  )
}

# Evaluate a call into the Rust glue so that a refusal reaches the caller the
# way a refused `ferx_fit()` does: as an R error, classed `ferx_engine_error`
# when `.ferx_engine_error()` can tie a diagnostic code to it (#385). `call` is
# evaluated lazily inside the `tryCatch()`, so it covers a failure the glue
# raises itself and one extendr raises for it (a panic in the engine) alike.
#
# Unlike `ferx_fit()`, these entry points also fail for reasons the validation
# pass never sees - a `fit` whose theta does not fit the model, a TTE simulation
# with no horizon - and with one unrelated error in the data the single-error
# fallback labelled those with that error's code ("theta length 2 does not match
# model [E_DOSE_CMT_NOT_INFUSABLE]", #386 review). So the fallback is confined
# to the two stages validation re-runs with the engine's own parser and reader,
# which the glue names in its prefix; there a text match is not to be had (the
# parser adds "(line N)", validation wraps the reader's message in the file
# name) and the single error is the one that stopped the call. A message from
# any other stage gets a code only on a text match: a missing code, never a
# wrong one.
#
# Most of the glue spells those two prefixes "Error parsing model: " /
# "Error reading data: ". `ferx_rust_simulate_adaptive()` puts its own name in
# front and lower-cases the word ("ferx_simulate_adaptive: error parsing
# model: ", #390), so the optional `ferx_*: ` prefix is what lets those two
# stages reach the fallback there too. Nothing else is widened: that function's
# other refusals - a model with no `[adaptive_dosing]` block, and whatever the
# controller or the regimen rejects - carry the same prefix but neither stage
# name, so they still get a code only on a text match.
.ferx_engine_fallback_stages <-
  "^(ferx_[a-z_]+: )?[Ee]rror (parsing model|reading data): "

.ferx_engine_call <- function(call, model, data) {
  tryCatch(
    call,
    error = function(e) {
      stop(.ferx_engine_error(e, model, data,
                              fallback_stages = .ferx_engine_fallback_stages))
    }
  )
}
