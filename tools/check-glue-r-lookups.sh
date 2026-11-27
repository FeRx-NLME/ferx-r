#!/usr/bin/env bash
# Check that the extendr glue never looks an R function or variable up from
# the global environment.
#
# Why: extendr's `data_frame!` built its table by evaluating `data.frame(..)`
# in the global environment, so a user with their own `data.frame` there - by
# accident or on purpose - made every ferx call whose glue used the macro
# panic and print the whole table (#780; the VPC glue first, #778 review
# finding 5). The glue now builds tables with `r_data_frame(columns, nrow)`,
# which lays the list out exactly as `data.frame()` would without calling R.
# The same lookup sits behind other extendr-api 0.9.0 entry points, and any
# of them compiles and passes every test that does not shadow the name, so
# only a text rule keeps them out (#781 review finding 5):
#
#   macros     call!  data_frame!  factor!  R!  Rraw!   (parse and evaluate R
#              text there)   var!  global!   (read variables there)
#              lang!   (a call whose function `.eval()` resolves there)
#   functions  global_var  local_var  global_function  global_env
#   methods    .eval(  .eval_blind(   (evaluate in the global environment)
#   eval_string / eval_string_with_params - allowed only when the first
#              argument is a string literal that starts with a
#              namespace-qualified name, as `eval_string("base::tempdir()")`.
#              That still resolves `::` itself from the global environment;
#              a user who redefines `::` is past what this guard defends.
#
# A name is flagged as a macro call (`name!`, any delimiter, path-qualified
# too), as a function or method call, or when a `use` statement names it in
# any layout - so a rename (`use extendr_api::{\n data_frame as df,\n};`)
# cannot hide a later `df!(..)` (#781 review finding 3). `data_frame` is
# flagged wherever it occurs in code.
#
# Code is read with tools/rust-code.awk, one left-to-right scan per file:
# comments and literal contents never trigger a rule, a `/*` or `//` inside a
# string or a line comment cannot blind the lines after it (#781 review
# findings 1, 2, 8), and a file the scan ends inside a comment or literal is
# an error, not a pass. Every `.rs` file under the directory is scanned,
# nested modules included (finding 4).
# Exercised by tools/test-check-glue-r-lookups.sh.
#
# Usage: tools/check-glue-r-lookups.sh [dir]   (default src/rust/src)
# Exits 0 when no glue source looks a name up globally, 1 otherwise.

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

# `file<TAB>line<TAB>kind<TAB>name`; kind one of
#   macro / fn / method / use / eval / lost.
scan() {
  for f in "${sources[@]}"; do
    LC_ALL=C awk -v file="${f#"$ROOT"/}" "$SCAN"'
      BEGIN {
        split("call data_frame factor R Rraw var global lang", m, " ")
        for (k in m) MACRO[m[k]] = 1
        split("global_var local_var global_function global_env", g, " ")
        for (k in g) FN[g[k]] = 1
        for (k in MACRO) NAMED[k] = 1
        for (k in FN) NAMED[k] = 1
        NAMED["eval_string"] = 1; NAMED["eval_string_with_params"] = 1
      }
      function report(kind, name) { printf "%s\t%d\t%s\t%s\n", file, NR, kind, name }
      {
        code = rust_code($0)
        trimmed = code
        sub(/^[[:space:]]+/, "", trimmed)
        # A `use` statement runs to its `;`, over as many lines as it takes.
        if (trimmed ~ /^(pub(\([^)]*\))?[[:space:]]+)?use[[:space:]]/) in_use = 1
        rest = " " code " "
        while (match(rest, /[A-Za-z_][A-Za-z0-9_]*/)) {
          name = substr(rest, RSTART, RLENGTH)
          before = substr(rest, RSTART - 1, 1)
          after = substr(rest, RSTART + RLENGTH)
          rest = after
          if (before ~ /[A-Za-z0-9_]/) continue
          bang = (after ~ /^[[:space:]]*!/)
          paren = (after ~ /^[[:space:]]*\(/)
          if (name == "data_frame") { report(bang ? "macro" : (in_use ? "use" : "name"), name); continue }
          if (in_use && (name in NAMED)) { report("use", name); continue }
          if (bang && (name in MACRO)) { report("macro", name); continue }
          if (!paren) continue
          if (before == "." && (name == "eval" || name == "eval_blind")) { report("method", name); continue }
          if (before == ".") continue
          if (name in FN) { report("fn", name); continue }
          if (name == "eval_string" || name == "eval_string_with_params") {
            arg = after
            sub(/^[[:space:]]*\([[:space:]]*/, "", arg)
            ok = 0
            if (match(arg, /^"[0-9]+"/)) {
              lit = rust_lit[substr(arg, 2, RLENGTH - 2) + 0]
              ok = (lit ~ /^[[:space:]]*[A-Za-z][A-Za-z0-9._]*:::?[A-Za-z.`]/)
            }
            if (!ok) report("eval", name)
          }
        }
        if (in_use && code ~ /;/) in_use = 0
      }
      END { if ((what = rust_scan_open()) != "") report("lost", what) }
    ' "$f"
  done
}

msg() { # kind, name
  case "$1" in
    macro)
      if [[ "$2" == data_frame ]]; then
        echo "data_frame! builds the table by evaluating data.frame() from the global environment, so a user's own data.frame there makes this call panic (#780). Build it with r_data_frame(columns, nrow), which gives data.frame()'s exact layout without calling R."
      else
        echo "$2! looks R names up from the global environment, where a user's own definition replaces base's (#780). Build the object in Rust (r_data_frame() for a table), or reach the function through its namespace."
      fi ;;
    name) echo "data_frame is named here. The glue builds tables with r_data_frame(columns, nrow), never extendr's data_frame! (#780)." ;;
    use) echo "$2 is imported or renamed here. It looks R names up from the global environment, where a user's own definition replaces base's (#780); a rename would hide every later use from this guard." ;;
    fn) echo "$2() reads the global environment, where a user's own definition replaces base's (#780). Reach the object through its namespace instead." ;;
    method) echo ".$2() evaluates in the global environment, where a user's own definitions replace base's (#780). Use eval_with_env() with base_env() or the namespace the code needs." ;;
    eval) echo "$2() evaluates its text in the global environment (#780). Only a string literal that starts with a namespace-qualified name, as eval_string(\"base::tempdir()\"), is allowed here - write it on the call's line." ;;
    lost) echo "the scan reached the end of the file inside $2, so it lost track of what is code and checked nothing after that point. Fix the source, or teach tools/rust-code.awk the syntax it misread." ;;
  esac
}

n=0
while IFS=$'\t' read -r file line kind name; do
  [[ "$kind" == lost ]] || n=$((n + 1))
  fail "$file" "$line" "$(msg "$kind" "$name")"
done < <(scan)

if [[ $status -eq 0 ]]; then
  echo "glue R lookups OK: no global-environment lookup in ${#sources[@]} source(s) under ${DIR#"$ROOT"/}; tables go through r_data_frame()."
elif [[ $n -gt 0 ]]; then
  echo "$n global-environment lookup(s) in the glue." >&2
fi
exit $status
