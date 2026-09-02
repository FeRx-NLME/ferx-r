# Parse section header positions and names from a character vector of file
# lines. Returns a list with $positions (integer indices) and $names (strings).
ferx_section_headers <- function(lines) {
  pos   <- grep("^\\s*\\[.+\\]\\s*$", lines)
  names <- gsub("^\\s*\\[|\\]\\s*$", "", lines[pos])
  list(positions = pos, names = names)
}

# Shared validation for ferx_model_get_section()/ferx_model_set_section():
# `path` must exist and have a .ferx extension. Used by both so the two
# functions can't drift on what counts as a valid model file.
.ferx_validate_ferx_path <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (tolower(tools::file_ext(path)) != "ferx") stop("'path' must be a .ferx file")
}

# Locate a named section's header index within already-parsed `hdr`
# (from ferx_section_headers()). Errors with the available section names
# when `section` isn't present. Shared by ferx_model_get_section() and
# ferx_model_set_section() so the "not found" message can't drift between
# the two.
.ferx_section_index <- function(hdr, section) {
  idx <- which(hdr$names == section)
  if (length(idx) == 0L) {
    stop(
      "Section '", section, "' not found. ",
      "Available sections: ", paste(hdr$names, collapse = ", ")
    )
  }
  idx
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
# Returns a short string ("1-cpt oral", "ODE", "compartment-free", etc.) or
# NULL when unrecognised.
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
  if (length(m) == 0L) {
    # A compartment-free structural model declares the prediction directly.
    # The engine requires a y assignment and rejects mixing this form with a
    # pk/ode declaration, so y is the unambiguous marker on this pre-fit path.
    y_re <- "^\\s*y(?:\\s*\\[\\s*CMT\\s*=\\s*\\d+\\s*\\])?\\s*="
    if (any(grepl(y_re, lines, ignore.case = TRUE, perl = TRUE))) {
      return("compartment-free")
    }
    return(NULL)
  }
  fn <- sub("^pk\\s+([A-Za-z_][A-Za-z_0-9]*)\\s*\\(.*", "\\1", m, perl = TRUE)

  # Long `*_compartment_*` aliases collapse to their `*_cpt_*` equivalents.
  fn <- sub("_compartment_", "_cpt_", fn, fixed = TRUE)

  .ferx_fmt_pk_name(fn)
}

# Peel an optional trailing `weight = <expr>` modifier off a `[parameters]`
# kappa declaration (ferx-core #1031), returning the weight expression or
# NA_character_ when the line carries none.
#
# This deliberately mirrors ferx-core's `split_weight_modifier()`
# (parser/model_parser.rs) rather than approximating it with a regex, because
# the two are read side by side: this one feeds ferx_model_inspect() pre-fit
# while the engine's feeds model_structure post-fit, and any divergence shows
# up as inspect() calling a weighted model unweighted. The engine's matching
# rules, all reproduced here:
#   * case-insensitive, so `WEIGHT =` is the modifier;
#   * whole word, so `WEIGHTED` and `X_weight` are not;
#   * only at bracket depth 0, so a covariate named `weight` inside a
#     magnitude expression is not mistaken for the modifier;
#   * followed by a single `=`, so a `==` comparison is skipped;
#   * first match wins.
# The engine rejects an empty right-hand side and a modifier with no statement
# in front of it; here both report "no weight", since such a file fails to
# parse at fit time anyway and inspect() must not invent a label for it.
.ferx_split_weight_modifier <- function(line) {
  none <- NA_character_
  ch <- strsplit(line, "", fixed = TRUE)[[1L]]
  n  <- length(ch)
  kw <- c("w", "e", "i", "g", "h", "t")
  k  <- length(kw)
  if (n < k) return(none)
  is_ident <- function(c) grepl("^[A-Za-z0-9_]$", c)
  depth <- 0L
  for (i in seq_len(n)) {
    c_i <- ch[i]
    if (c_i == "(" || c_i == "[") {
      depth <- depth + 1L
      next
    }
    if (c_i == ")" || c_i == "]") {
      depth <- depth - 1L
      next
    }
    if (depth != 0L || i + k - 1L > n) next
    if (!identical(tolower(ch[i:(i + k - 1L)]), kw)) next
    # Whole-word match on both sides.
    if (i > 1L && is_ident(ch[i - 1L])) next
    j <- i + k
    if (j <= n && is_ident(ch[j])) next
    # Skip whitespace, then require a single `=` (not `==`).
    while (j <= n && grepl("^[ \t]$", ch[j])) j <- j + 1L
    if (j > n || ch[j] != "=") next
    if (j + 1L <= n && ch[j + 1L] == "=") next
    expr <- trimws(paste(ch[seq_len(n) > j], collapse = ""))
    stmt <- trimws(paste(ch[seq_len(i - 1L)], collapse = ""))
    if (!nzchar(expr) || !nzchar(stmt)) return(none)
    return(expr)
  }
  none
}

# Parse a .ferx file and return a named list describing model structure.
# Fields: theta_names (pop param names), model_type (label or NULL),
#         iiv, iov, residual.
# Used by ferx_model_inspect() (pre-fit) and attached to ferx_fit() results.
.ferx_parse_structure <- function(path) {
  b <- .ferx_extract_blocks(path)

  # Population (theta) parameter names from [parameters].
  #
  # Every declaration keyword is matched case-insensitively because the
  # engine's declaration regexes all carry `(?i)` (model_parser.rs: theta_re,
  # omega_re, sigma_re, kappa_re and the three block_* forms). A model written
  # `THETA TVCL(1.0, 0.001, 100.0)` fits exactly like the lowercase spelling,
  # but matching it case-sensitively here reported *no* thetas, no IIV and no
  # IOV - an entirely blank structure for a model the engine parses fine.
  params      <- b[["parameters"]] %||% character(0)
  theta_lines <- grep("^theta\\s", params, value = TRUE, ignore.case = TRUE)
  thetas <- if (length(theta_lines) > 0L)
    sub("^theta\\s+(\\w+).*", "\\1", theta_lines, ignore.case = TRUE)
  else
    character(0)

  # Optional model-type label from [structural_model]
  struct_lines <- b[["structural_model"]] %||% character(0)
  model_type   <- if (length(struct_lines) > 0L) .ferx_model_type(struct_lines) else NULL

  # IIV: omega lines
  omega_lines <- grep("^omega\\s", params, value = TRUE, ignore.case = TRUE)
  iiv <- if (length(omega_lines) > 0L)
    sub("^omega\\s+(\\w+).*", "\\1", omega_lines, ignore.case = TRUE)
  else
    character(0)

  # IOV: kappa lines
  kappa_lines <- grep("^kappa\\s", params, value = TRUE, ignore.case = TRUE)
  iov <- if (length(kappa_lines) > 0L)
    sub("^kappa\\s+(\\w+).*", "\\1", kappa_lines, ignore.case = TRUE)
  else
    character(0)

  # Sample-size-weighted IOV (ferx-core #1031): `kappa K ~ 2.0 (sd) weight = NARM`
  # declares `kappa_ik ~ N(0, Omega_IOV / N_ik)`. Capture the weight expression
  # so ferx_model_inspect() shows it pre-fit, mirroring the `iov_weights` the
  # engine attaches to model_structure post-fit. Left as character(0) when no
  # kappa is weighted, so an ordinary IOV model's structure list is unchanged.
  iov_weights <- if (length(kappa_lines) > 0L) {
    w <- vapply(kappa_lines, .ferx_split_weight_modifier, character(1L),
                USE.NAMES = FALSE)
    if (all(is.na(w))) character(0) else w
  } else {
    character(0)
  }

  # Safety net for the next time this parser drifts from the engine's. Every
  # field above is a best-effort regex read of a grammar that lives in Rust,
  # so a form the engine accepts and these patterns miss is reported as an
  # absence rather than an error - the failure mode that let a whole
  # `[parameters]` block read as empty (see the case-insensitivity note
  # above). A block with content but not one recognised declaration in it is
  # not a plausible model; say so rather than hand back a blank structure.
  # Deliberately keyword-agnostic: it fires on whatever the next divergence
  # turns out to be, not just on the ones already known.
  decl_re <- "^(block_)?(theta|omega|sigma|kappa)\\s*[[:space:](]"
  if (length(params) > 0L &&
        !any(grepl(decl_re, params, ignore.case = TRUE))) {
    warning(
      "No parameter declarations recognised in the [parameters] block of '",
      basename(path), "' (", length(params), " non-empty lines). Reporting an ",
      "empty model structure; the engine may still parse this file. If it ",
      "fits, this is a bug in ferx's model inspector - please report it.",
      call. = FALSE
    )
  }

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
    iov_weights = iov_weights,
    residual    = residual
  )
}

# Resolve the dataset path declared in a model file's `[data]` block (#254).
# Returns the resolved path (character(1)) or NULL when the model file declares
# no `[data] path`. Relative paths are resolved by the engine relative to the
# model file's directory. `model_path` must be a path to a .ferx file.
.ferx_model_data_path <- function(model_path) {
  if (is.null(model_path) || !is.character(model_path) ||
      length(model_path) != 1L || !file.exists(model_path)) {
    return(NULL)
  }
  p <- tryCatch(
    ferx_rust_model_data_path(normalizePath(model_path)),
    error = function(e) ""
  )
  if (length(p) != 1L || is.na(p) || !nzchar(p)) return(NULL)
  p
}
