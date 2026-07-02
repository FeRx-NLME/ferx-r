# Parse section header positions and names from a character vector of file
# lines. Returns a list with $positions (integer indices) and $names (strings).
ferx_section_headers <- function(lines) {
  pos   <- grep("^\\s*\\[.+\\]\\s*$", lines)
  names <- gsub("^\\s*\\[|\\]\\s*$", "", lines[pos])
  list(positions = pos, names = names)
}


# Extract all named [section] blocks from a .ferx file.
# Returns a named list: section name ? character vector of (comment-stripped, trimmed) lines.
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
# parser (ferx-core src/parser/model_parser.rs `pk_func_name` match arms) so
# pre-fit `ferx_model_inspect(path)` reports the same string the engine would
# attach to a fitted result.
# Format a pk function name (snake_case) into a readable label.
# e.g. "one_cpt_oral" -> "1-cpt oral", "two_cpt_iv" -> "2-cpt IV".
#
# ferx-core #176 retired the split `*_iv_bolus` / `*_infusion` names -
# the bolus-vs-infusion choice is now per dose from RATE, so `*_iv`
# covers both. Retired names are passed through here so pre-fit
# `ferx_model_inspect()` on a stale `.ferx` file can still produce a
# best-effort label; the engine itself will then reject the parse with
# a migration error.
.ferx_fmt_pk_name <- function(fn) {
  label <- fn
  label <- sub("^one_",   "1_",   label)
  label <- sub("^two_",   "2_",   label)
  label <- sub("^three_", "3_",   label)
  label <- sub("_cpt_",   "-cpt ", label, fixed = TRUE)
  label <- sub("_cpt$",   "-cpt", label)
  label <- gsub("_", " ", label, fixed = TRUE)
  label <- gsub("(?<![a-z])iv(?![a-z])", "IV", label, perl = TRUE)
  # `*_iv_bolus` (retired) - the lone "IV" is already capitalised above;
  # nothing extra to inject here.
  # `*_infusion` (retired) - the bare function name has no `iv` token,
  # so inject "IV " to keep the R-side label intelligible. The Rust
  # label for these retired forms is unreachable (parser errors first).
  label <- sub("\\binfusion\\b", "IV infusion", label, perl = TRUE)
  label
}

.ferx_model_type <- function(lines) {
  s <- paste(lines, collapse = " ")
  if (grepl("\\bode\\(", s, perl = TRUE)) return("ODE")

  m <- regmatches(s, regexpr("\\bpk\\s+([A-Za-z_][A-Za-z_0-9]*)\\s*\\(",
                             s, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  fn <- sub("^pk\\s+([A-Za-z_][A-Za-z_0-9]*)\\s*\\(.*", "\\1", m, perl = TRUE)

  # Long `*_compartment_*` aliases collapse to their `*_cpt_*` equivalents.
  fn <- sub("_compartment_", "_cpt_", fn, fixed = TRUE)

  .ferx_fmt_pk_name(fn)
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

  # Residual error type from [error_model]. Multi-endpoint blocks use a
  # `CMT=N:` prefix per line; report them as "per-CMT (CMT2=proportional,
  # CMT3=additive)" to match the label the Rust engine attaches post-fit.
  err_lines <- b[["error_model"]] %||% character(0)
  err_lines <- err_lines[nzchar(trimws(err_lines)) & !grepl("^\\s*#", err_lines)]
  err_type <- function(line) {
    # Log-transform-both-sides (LTBS): `log(DV) ~ additive(...)` (engine logs DV)
    # or `DV ~ log_additive(...)` (DV already log). Both are additive error on
    # the log scale - detect them BEFORE the plain `additive` branch, since both
    # forms contain the substring "additive".
    if (grepl("log_additive", line, ignore.case = TRUE) ||
        grepl("log\\s*\\(\\s*DV\\s*\\)", line, ignore.case = TRUE)) {
      "additive (log-transformed)"
    } else if (grepl("proportional", line, ignore.case = TRUE)) {
      "proportional"
    } else if (grepl("additive", line, ignore.case = TRUE)) {
      "additive"
    } else if (grepl("combined", line, ignore.case = TRUE)) {
      "combined"
    } else {
      NA_character_
    }
  }
  cmt_lines <- grep("^\\s*CMT\\s*=", err_lines, value = TRUE)
  residual <- if (length(err_lines) == 0L) {
    "unknown"
  } else if (length(cmt_lines) > 0L) {
    cmts  <- suppressWarnings(as.integer(sub("^\\s*CMT\\s*=\\s*([0-9]+).*", "\\1", cmt_lines)))
    types <- vapply(cmt_lines, err_type, character(1L))
    if (anyNA(cmts) || anyNA(types)) {
      # An unparseable CMT index or unrecognised error type would otherwise
      # render as "CMTNA=..." / "CMT2=NA". Warn and fall back, mirroring the
      # single-endpoint path, so the residual label never contains NA.
      warning(
        "Unrecognised per-CMT residual error specification; reporting as ",
        "\"unknown\". Lines: ", paste(cmt_lines, collapse = " | "),
        call. = FALSE
      )
      "unknown"
    } else {
      ord <- order(cmts)
      paste0(
        "per-CMT (",
        paste(sprintf("CMT%d=%s", cmts[ord], types[ord]), collapse = ", "),
        ")"
      )
    }
  } else {
    t1 <- err_type(err_lines[1L])
    if (is.na(t1)) {
      warning("Unrecognised residual error type; reporting as \"unknown\". Line: ",
              err_lines[1L], call. = FALSE)
      "unknown"
    } else {
      t1
    }
  }

  list(
    theta_names = thetas,
    model_type  = model_type,
    iiv         = iiv,
    iov         = iov,
    residual    = residual
  )
}
