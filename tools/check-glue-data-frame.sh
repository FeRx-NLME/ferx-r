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
# Matches `data_frame!` followed by any macro delimiter - `(`, `[` or `{` -
# with or without a path in front (`extendr_api::data_frame!(..)`). `//` and
# `/* .. */` comments are stripped first, so naming the macro in prose is free.
# A `use .. data_frame as other` rename is not caught; none exists, and the
# import itself is flagged (`data_frame` named in a `use`).
# Exercised by tools/test-check-glue-data-frame.sh.
#
# Usage: tools/check-glue-data-frame.sh [dir]   (default src/rust/src)
# Exits 0 when no glue source uses the macro, 1 otherwise.

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

shopt -s nullglob
sources=("$DIR"/*.rs)
if [[ ${#sources[@]} -eq 0 ]]; then
  echo "error: no .rs sources under $DIR" >&2
  exit 1
fi

MSG_CALL="data_frame! builds the table by evaluating data.frame() from the global environment, so a user's own data.frame there makes this call panic (#780). Build it with r_data_frame(columns, nrow), which gives data.frame()'s exact layout without calling R."
MSG_USE="data_frame is imported here. The glue builds tables with r_data_frame(columns, nrow), never extendr's data_frame! (#780)."

scan() {
  for f in "${sources[@]}"; do
    awk -v file="${f#"$ROOT"/}" '
      {
        code = $0
        if (in_block) {
          if (match(code, /\*\//)) { code = substr(code, RSTART + 2); in_block = 0 }
          else code = ""
        }
        while (match(code, /\/\*/)) {
          head = substr(code, 1, RSTART - 1)
          tail = substr(code, RSTART + 2)
          if (match(tail, /\*\//)) code = head substr(tail, RSTART + 2)
          else { code = head; in_block = 1; break }
        }
        sub(/[[:space:]]*\/\/.*$/, "", code)
      }
      code ~ /(^|[^A-Za-z0-9_])data_frame[[:space:]]*![[:space:]]*(\(|\[|\{)/ {
        printf "%s\t%d\tcall\n", file, NR
      }
      code ~ /^[[:space:]]*(pub(\([^)]*\))?[[:space:]]+)?use[[:space:]]/ && code ~ /(^|[^A-Za-z0-9_])data_frame([^A-Za-z0-9_]|$)/ {
        printf "%s\t%d\tuse\n", file, NR
      }
    ' "$f"
  done
}

n=0
while IFS=$'\t' read -r file line kind; do
  n=$((n + 1))
  if [[ "$kind" == use ]]; then fail "$file" "$line" "$MSG_USE"; else fail "$file" "$line" "$MSG_CALL"; fi
done < <(scan)

if [[ $status -eq 0 ]]; then
  echo "glue data frames OK: no data_frame! in ${#sources[@]} source(s) under ${DIR#"$ROOT"/}; tables go through r_data_frame()."
else
  echo "$n data_frame! use(s) in the glue." >&2
fi
exit $status
