#!/usr/bin/env bash
# Check that the extendr glue raises R errors in exactly one place.
#
# Why: `Rf_error()` reads its first argument as a printf format string and
# longjmps out of Rust without unwinding. Handing it the engine's text garbles
# a message that contains `%` and can abort the R session (ferx-r #388), and
# raising from inside an entry point body leaks every local that body is
# holding (ferx-r #389). Both are closed by the same shape: an `#[extendr]`
# body is a closure returning `Result<_, String>`, and `entry()` raises once -
# after that closure frame has returned - through a `"%s"` format the glue owns.
#
# Most of that shape the compiler holds on its own: lib.rs shadows the prelude
# `throw_r_error` with a local `fn throw_r_error(_: Infallible) -> !`, so a bare
# `throw_r_error(format!(..))` is `error[E0308]` on every machine, not a CI
# finding. This script covers only what the compiler cannot see:
#
#   1. the shadow itself is there, and no other `throw_r_error(` is in the code
#      (a fully qualified `extendr_api::throw_r_error(..)` walks past the shadow);
#   2. `raise_verbatim(` is called from exactly one place, so no body can raise
#      while its locals are still alive;
#   3. every `#[extendr]` fn body opens with `entry(`.
#
# `//` comments are stripped before every test, so naming a helper in prose is
# free. Exercised by tools/test-check-glue-raise.sh.
#
# Usage: tools/check-glue-raise.sh [dir]   (default src/rust/src)
# Exits 0 when the glue is in shape, 1 otherwise.

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

# `file<TAB>line<TAB>kind<TAB>message`, one line per finding, over every source
# at once so the raise_verbatim count is a count over the whole glue and not
# per file. `@shadow` and `@raise` are tallies rather than findings.
scan() {
  for f in "${sources[@]}"; do
    awk -v file="${f#"$ROOT"/}" '
      {
        code = $0
        sub(/[[:space:]]*\/\/.*$/, "", code)
        trimmed = code
        sub(/^[[:space:]]+/, "", trimmed)
      }

      code ~ /throw_r_error[[:space:]]*\(/ {
        if (trimmed ~ /^fn[[:space:]]+throw_r_error[[:space:]]*\(/)
          printf "%s\t%d\t@shadow\t-\n", file, NR
        else
          printf "%s\t%d\tfinding\t%s\n", file, NR, "throw_r_error() is reachable here. Raising goes through entry(); the local fn of that name exists only so that a bare call does not compile."
      }

      code ~ /raise_verbatim[[:space:]]*\(/ \
        && trimmed !~ /^fn[[:space:]]+raise_verbatim[[:space:]]*\(/ {
        printf "%s\t%d\t@raise\t-\n", file, NR
      }

      trimmed == "#[extendr]" { pending = 1; in_sig = 0; expect_entry = 0; next }
      pending && !in_sig && trimmed ~ /^fn[[:space:]]/ {
        in_sig = 1; fn_line = NR
        name = trimmed
        sub(/^fn[[:space:]]+/, "", name)
        sub(/[^A-Za-z0-9_].*$/, "", name)
        fn_name = name
      }
      in_sig {
        if (code ~ /\{[[:space:]]*$/) { in_sig = 0; expect_entry = 1 }
        next
      }
      expect_entry && trimmed == "" { next }
      expect_entry {
        if (trimmed !~ /^entry[[:space:]]*\(/)
          printf "%s\t%d\tfinding\t%s\n", file, fn_line, fn_name " does not open its body with entry(). Every #[extendr] body is entry(move || { ..; Ok(value) }), so the raise happens once, after the body has returned."
        expect_entry = 0; pending = 0
      }
    ' "$f"
  done
}

shadows=0
raises=0
while IFS=$'\t' read -r file line kind msg; do
  case "$kind" in
    @shadow) shadows=$((shadows + 1)) ;;
    @raise)
      raises=$((raises + 1))
      if [[ $raises -gt 1 ]]; then
        fail "$file" "$line" "raise_verbatim() is called from more than one place. Only entry() may raise - a body that raises directly still holds its locals when Rf_error() longjmps over them (#389)."
      fi
      ;;
    *) fail "$file" "$line" "$msg" ;;
  esac
done < <(scan)

if [[ $shadows -eq 0 ]]; then
  fail "${sources[0]#"$ROOT"/}" 1 \
    "the local fn throw_r_error(_: Infallible) shadow is gone. Without it a bare throw_r_error(format!(..)) compiles again and hands the message to Rf_error() as its format string (#388)."
fi

if [[ $raises -eq 0 ]]; then
  fail "${sources[0]#"$ROOT"/}" 1 \
    "nothing calls raise_verbatim(). entry() is the one place an R error is raised from; a glue that raises nowhere is a glue whose errors went somewhere else."
fi

if [[ $status -eq 0 ]]; then
  echo "glue raise shape OK: entry() opens every #[extendr] body, one raise_verbatim() call site, throw_r_error() shadowed."
fi
exit $status
