#!/usr/bin/env bash
# Check that the extendr glue never builds a table with extendr's `data_frame!`.
#
# Why: `data_frame!` builds its result by evaluating `data.frame(..)` in R, and
# that lookup starts in the global environment. A user with a function of that
# name there - by accident or on purpose - makes every ferx call whose glue uses
# the macro panic on an `unwrap()` and print the whole table (#780; the VPC
# glue first, #778 review finding 5). The glue builds every table with
# `r_data_frame(columns, nrow)` instead, which lays the list out exactly as
# `data.frame()` would (names, class, compact row names `c(NA, -n)`,
# `integer(0)` for 0 rows) without calling R. A `data_frame!` call compiles and
# passes every test that does not shadow `data.frame`, so only a text rule
# keeps it out.
#
# The rule: the identifier `data_frame` does not occur in code. That covers a
# call with any delimiter (`data_frame!(..)`, `![..]`, `!{..}`), a qualified
# one (`extendr_api::data_frame!(..)`), and an import or rename in any layout
# (`use extendr_api::{\n List,\n data_frame as df,\n};` is flagged on the line
# that names it, so a later `df!(..)` cannot hide). `r_data_frame` and other
# longer identifiers are not it.
#
# "In code" is decided by one left-to-right scan per file that tracks `//` and
# (nested) `/* .. */` comments, string and byte-string literals with their
# escapes, raw strings `r#".."#`, and char literals (told apart from
# lifetimes). Comments and literal contents are dropped before the rule runs,
# so naming the macro in prose or in a message is free, and a `/*` or `//`
# inside a string or a line comment cannot blind the rest of the file (#781
# review findings 1, 2, 8). A file the scan ends inside a comment or literal
# is an error, not a pass: the scan lost track, and everything after that
# point went unchecked.
#
# Every `.rs` file under the directory is scanned, nested modules included.
# Exercised by tools/test-check-glue-data-frame.sh.
#
# Usage: tools/check-glue-data-frame.sh [dir]   (default src/rust/src)
# Exits 0 when no glue source names data_frame in code, 1 otherwise.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="${1:-$ROOT/src/rust/src}"
status=0

fail() { # file, line, message
  if [[ "${GITHUB_ACTIONS:-}" == true ]]; then
    echo "::error file=$1,line=$2::$3"
  else
    echo "error: $1:$2: $3" >&2
  fi
  status=1
}

sources=()
if [[ -d "$DIR" ]]; then
  while IFS= read -r f; do sources+=("$f"); done < <(find "$DIR" -type f -name '*.rs' | LC_ALL=C sort)
fi
if [[ ${#sources[@]} -eq 0 ]]; then
  echo "error: no .rs sources under $DIR" >&2
  exit 1
fi

MSG_CALL="data_frame! builds the table by evaluating data.frame() from the global environment, so a user's own data.frame there makes this call panic (#780). Build it with r_data_frame(columns, nrow), which gives data.frame()'s exact layout without calling R."
MSG_NAME="data_frame is named here (an import or a rename). The glue builds tables with r_data_frame(columns, nrow), never extendr's data_frame! (#780)."

# `file<TAB>line<TAB>kind`, kind one of call / name / lost:<what>.
scan() {
  for f in "${sources[@]}"; do
    LC_ALL=C awk -v file="${f#"$ROOT"/}" '
      function ident(ch) { return ch ~ /^[A-Za-z0-9_]$/ }
      BEGIN { state = "code" }   # code | block | str | raw
      {
        line = $0; n = length(line); out = ""; i = 1
        while (i <= n) {
          c = substr(line, i, 1); c2 = substr(line, i, 2)
          if (state == "block") {
            if (c2 == "/*") { depth++; i += 2 }
            else if (c2 == "*/") { depth--; i += 2; if (depth == 0) { state = "code"; out = out " " } }
            else i++
            continue
          }
          if (state == "str") {
            if (c == "\\") i += 2
            else if (c == "\"") { state = "code"; out = out " "; i++ }
            else i++
            continue
          }
          if (state == "raw") {
            if (c == "\"" && substr(line, i + 1, nhash) == hashes) { state = "code"; out = out " "; i += 1 + nhash }
            else i++
            continue
          }
          # -- code --
          if (c2 == "//") break
          if (c2 == "/*") { state = "block"; depth = 1; i += 2; continue }
          if (c == "\"") { state = "str"; i++; continue }
          prev = (i > 1) ? substr(line, i - 1, 1) : ""
          if (!ident(prev) && (c == "r" || c2 == "br")) {
            rest = substr(line, i + (c == "r" ? 1 : 2))
            if (match(rest, /^#*"/)) {
              nhash = RLENGTH - 1; hashes = substr(rest, 1, nhash)
              state = "raw"; i += (c == "r" ? 1 : 2) + RLENGTH
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
            if (close_at) { out = out " "; i = close_at + 1; continue }
          }
          out = out c; i++
        }
        code = " " out " "
        while (match(code, /[^A-Za-z0-9_]data_frame[^A-Za-z0-9_]/)) {
          after = substr(code, RSTART + RLENGTH - 1)
          kind = (after ~ /^[[:space:]]*!/) ? "call" : "name"
          printf "%s\t%d\t%s\n", file, NR, kind
          code = substr(code, RSTART + RLENGTH - 1)
        }
      }
      END {
        if (state != "code") {
          what = (state == "block") ? "a block comment" : (state == "raw") ? "a raw string" : "a string literal"
          printf "%s\t%d\tlost:%s\n", file, NR, what
        }
      }
    ' "$f"
  done
}

n=0
while IFS=$'\t' read -r file line kind; do
  case "$kind" in
    call) n=$((n + 1)); fail "$file" "$line" "$MSG_CALL" ;;
    name) n=$((n + 1)); fail "$file" "$line" "$MSG_NAME" ;;
    lost:*) fail "$file" "$line" "the scan reached the end of the file inside ${kind#lost:}, so it lost track of what is code and checked nothing after that point. Fix the source, or teach this script the syntax it misread." ;;
  esac
done < <(scan)

if [[ $status -eq 0 ]]; then
  echo "glue data frames OK: data_frame named nowhere in code in ${#sources[@]} source(s) under ${DIR#"$ROOT"/}; tables go through r_data_frame()."
elif [[ $n -gt 0 ]]; then
  echo "$n data_frame use(s) in the glue." >&2
fi
exit $status
