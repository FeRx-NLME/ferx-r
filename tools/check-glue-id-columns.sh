#!/usr/bin/env bash
# Check that every `ID` column the extendr glue returns is built by
# `subject_id_column()`, so it is the subjects' ID text, verbatim (#786).
#
# Why: `fit$sdtab$ID` used to be the engine's number for the subject - the ID
# parsed as a number, else its 1-based position (`sdtab_id()`). "S1" came back
# as 1 and "001" as 1, subjects "1" and "1.0" merged, and the column's type
# depended on the data, so a join of `ferx_simulate(fit = )` (text IDs) onto
# `fit$sdtab` matched nothing. A new numeric ID column compiles and passes
# every test that uses numeric IDs, so a text rule keeps them out:
#
#   tuple      `("ID", <value>)` - a column of a table - whose value is not a
#              `subject_id_column(..)` call. The value may start on the next
#              line. A table built from parallel name/value vectors names
#              "ID" in another shape; those sites call the helper too, and the
#              R sweep in tests/testthat/test-subject-id-type.R holds them.
#   sdtab_id   the deleted position fallback, named anywhere in code.
#   parse      `.id.parse::<..>()` - a subject ID read as a number.
#
# Code is read with tools/rust-code.awk, so comments and string contents never
# trigger a rule, and a file the scan ends inside a comment or literal is an
# error, not a pass. Every `.rs` file under the directory is scanned.
# Exercised by tools/test-check-glue-id-columns.sh.
#
# Usage: tools/check-glue-id-columns.sh [dir]   (default src/rust/src)
# Exits 0 when every returned ID column is the subjects' text, 1 otherwise.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="${1:-$ROOT/src/rust/src}"
SCAN="$(cat "$ROOT/tools/rust-code.awk")"
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

# `file<TAB>line<TAB>kind`; kind one of tuple / sdtab_id / parse / lost.
scan() {
  for f in "${sources[@]}"; do
    LC_ALL=C awk -v file="${f#"$ROOT"/}" "$SCAN"'
      function report(kind, at) { printf "%s\t%d\t%s\n", file, at, kind }
      # Whether `rest` (what follows `("ID",`) opens with the helper call.
      function helper(rest) { return rest ~ /^[[:space:]]*subject_id_column[[:space:]]*\(/ }
      {
        code = rust_code($0)
        if (pending) {
          if (code ~ /^[[:space:]]*$/) next
          if (!helper(code)) report("tuple", pending)
          pending = 0
        }
        # A tuple may open at the end of the line before its `"ID",`.
        rest = (opened ? "(" : "") code
        if (code !~ /^[[:space:]]*$/) opened = (code ~ /\([[:space:]]*$/)
        while (match(rest, /\([[:space:]]*"[0-9]+"[[:space:]]*,/)) {
          hit = substr(rest, RSTART, RLENGTH)
          rest = substr(rest, RSTART + RLENGTH)
          n = hit; gsub(/[^0-9]/, "", n)
          if (rust_lit[n + 0] != "ID") continue
          if (rest ~ /^[[:space:]]*$/) { pending = NR; break }
          if (!helper(rest)) report("tuple", NR)
        }
        rest = " " code " "
        while (match(rest, /[^A-Za-z0-9_]sdtab_id[^A-Za-z0-9_]/)) {
          report("sdtab_id", NR)
          rest = substr(rest, RSTART + RLENGTH - 1)
        }
        if (code ~ /\.id[[:space:]]*\.[[:space:]]*parse[[:space:]]*::/) report("parse", NR)
      }
      END {
        if (pending) report("tuple", pending)
        if ((what = rust_scan_open()) != "") report("lost", 0)
      }
    ' "$f"
  done
}

msg() { # kind
  case "$1" in
    tuple) echo "an ID column built without subject_id_column(). Every ID column ferx returns is the subjects' ID text, verbatim, as fit\$subject_ids holds it (#786): build it with (\"ID\", subject_id_column(ids))." ;;
    sdtab_id) echo "sdtab_id() numbered a subject by its ID parsed as a number, else its position, so \"S1\" and \"001\" both came back as 1 (#786). Return the subject's text with subject_id_column()." ;;
    parse) echo "a subject ID parsed as a number. An ID is a key, not a quantity: \"001\" would become 1 and \"1\" and \"1.0\" one subject (#786). Return the text with subject_id_column()." ;;
    lost) echo "the scan reached the end of the file inside a comment or literal, so it lost track of what is code and checked nothing after that point. Fix the source, or teach tools/rust-code.awk the syntax it misread." ;;
  esac
}

n=0
while IFS=$'\t' read -r file line kind; do
  n=$((n + 1))
  fail "$file" "$line" "$(msg "$kind")"
done < <(scan)

if [[ $status -eq 0 ]]; then
  echo "glue ID columns OK: every ID column in ${#sources[@]} source(s) under ${DIR#"$ROOT"/} is built by subject_id_column()."
else
  echo "$n ID column problem(s) in the glue." >&2
fi
exit $status
