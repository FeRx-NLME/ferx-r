# rust-code.awk - read Rust source one line at a time and keep only its code.
#
# Shared by the glue guards (tools/check-glue-raise.sh,
# tools/check-glue-r-lookups.sh). Each guard prepends this file to its own awk
# program; run it with LC_ALL=C so `length` and `substr` count bytes.
#
#   code = rust_code($0)
#
# returns the line with every comment removed and every literal replaced by a
# placeholder: a string, byte string or raw string becomes `"N"`, where
# rust_lit[N] holds its contents (escapes left as written, lines joined by
# "\n"); a char literal becomes `' '`. A literal or block comment that spans
# lines leaves its placeholder on the line where it closes.
#
# One scan, left to right, so what opens first wins: a `/*` or `//` inside a
# string is text, a `"` inside a comment is text, and block comments nest as
# Rust's do. Char literals are told apart from lifetimes and labels (`'a`),
# which stay code. Stripping comments with a regex first - `/* .. */`, then
# `//` - is what this replaced: a `/*` in a line comment or a string switched
# that stripper into block-comment mode for the rest of the file, and a `//`
# in a string hid the rest of its line (ferx-r #781 review findings 1, 2).
#
# At the end of a file, rust_scan_open() names what the scan is still inside
# ("" when it is back in code). Anything else means the scan lost track and
# every line after that point went unchecked, so a guard reports it as a
# failure rather than passing the file.

function rust_ident(ch) { return ch ~ /^[A-Za-z0-9_]$/ }

function rust_scan_open() {
  if (rust_state == "block") return "a block comment"
  if (rust_state == "str") return "a string literal"
  if (rust_state == "raw") return "a raw string"
  return ""
}

function rust_close_lit() {
  rust_nlit++
  rust_lit[rust_nlit] = rust_buf
  rust_buf = ""
  rust_state = ""
  return "\"" rust_nlit "\""
}

function rust_code(line,    n, i, c, c2, out, prev, rest, skip, nx, close_at, j) {
  n = length(line); out = ""; i = 1
  # A literal that continues from the previous line keeps its line break.
  if (rust_state == "str" || rust_state == "raw") rust_buf = rust_buf "\n"
  while (i <= n) {
    c = substr(line, i, 1); c2 = substr(line, i, 2)
    if (rust_state == "block") {
      if (c2 == "/*") { rust_depth++; i += 2 }
      else if (c2 == "*/") { rust_depth--; i += 2; if (rust_depth == 0) { rust_state = ""; out = out " " } }
      else i++
      continue
    }
    if (rust_state == "str") {
      if (c == "\\") { rust_buf = rust_buf substr(line, i, 2); i += 2 }
      else if (c == "\"") { out = out rust_close_lit(); i++ }
      else { rust_buf = rust_buf c; i++ }
      continue
    }
    if (rust_state == "raw") {
      if (c == "\"" && substr(line, i + 1, rust_nhash) == rust_hashes) { out = out rust_close_lit(); i += 1 + rust_nhash }
      else { rust_buf = rust_buf c; i++ }
      continue
    }
    # -- code --
    if (c2 == "//") break
    if (c2 == "/*") { rust_state = "block"; rust_depth = 1; i += 2; continue }
    if (c == "\"") { rust_state = "str"; rust_buf = ""; i++; continue }
    prev = (i > 1) ? substr(line, i - 1, 1) : ""
    if (!rust_ident(prev) && (c == "r" || c2 == "br")) {
      skip = (c == "r") ? 1 : 2
      rest = substr(line, i + skip)
      if (match(rest, /^#*"/)) {
        rust_nhash = RLENGTH - 1; rust_hashes = substr(rest, 1, rust_nhash)
        rust_state = "raw"; rust_buf = ""; i += skip + RLENGTH
        continue
      }
    }
    if (c == "\x27") {
      # A char literal, or a lifetime / label, which stays code.
      nx = substr(line, i + 1, 1)
      close_at = 0
      if (nx == "\\") {
        for (j = i + 3; j <= n; j++) if (substr(line, j, 1) == "\x27") { close_at = j; break }
      } else if (substr(line, i + 2, 1) == "\x27") {
        close_at = i + 2
      } else if (nx >= "\200") {
        # A multi-byte UTF-8 char: its closing quote is within 4 bytes.
        for (j = i + 2; j <= i + 5 && j <= n; j++) if (substr(line, j, 1) == "\x27") { close_at = j; break }
      }
      if (close_at) { out = out "\x27 \x27"; i = close_at + 1; continue }
    }
    out = out c; i++
  }
  return out
}
