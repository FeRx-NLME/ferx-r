#!/usr/bin/env bash
# Check that the extendr glue raises R errors in exactly one place.
#
# Why: `Rf_error()` reads its first argument as a printf format string and
# longjmps out of Rust without unwinding. Handing it the engine's text garbles
# a message that contains `%` and can abort the R session (ferx-r #388), and
# raising from inside an entry point body leaks every local that body is
# holding (ferx-r #389). Both are closed by the same shape: an `#[extendr]`
# body is a closure returning `Result<_, String>`, and `entry()` raises once -
# after that closure frame has returned. It raises by unwinding into extendr's
# wrapper with the text `%`-escaped for extendr-api 0.9.0's `throw_r_error`,
# so that extendr's own frame - which holds every argument `Robj` - unwinds too
# instead of being longjmp'd over (ferx-r #394).
#
# Part of that shape the compiler holds on its own. lib.rs shadows the prelude
# `throw_r_error` with a local `fn throw_r_error(_: Infallible) -> !`, so a bare
# `throw_r_error(format!(..))` is `error[E0308]`; and `Rf_error` is declared
# inside `mod raise`, so `Rf_error(..)` from elsewhere in the file is
# `error[E0425]` and `raise::Rf_error(..)` is `error[E0603]`.
#
# The module is not a substitute for check 4 below and check 4 is not a
# substitute for the module - they are complements. Measured on 2026-09-20: a
# second `extern "C" { fn Rf_error(..) -> !; }` block anywhere in the file
# re-declares the symbol at file scope and compiles with no warning, which puts
# #388 back with both the compiler and checks 1-3 green. The module closes the
# accidental call; only a text rule closes the deliberate re-declaration.
#
# So, what the compiler cannot see:
#
#   1. the shadow itself is there, and no other `throw_r_error(` is in the code
#      (a fully qualified `extendr_api::throw_r_error(..)` walks past the shadow);
#   2. `raise_verbatim(` is called from exactly one place, so no body can raise
#      while its locals are still alive;
#   3. every `#[extendr]` fn body opens with `entry(`;
#   4. `Rf_error` is named nowhere outside `mod raise` - not called, and not
#      declared a second time;
#   5. Cargo.lock resolves extendr-api and extendr-macros to 0.9.0 *from
#      crates.io* - a git or path replacement can report 0.9.0 and still carry
#      a different wrapper or `throw_r_error`. The `%`
#      doubling in `raise::format_escaped` is right only while `throw_r_error`
#      passes its text to `Rf_error` as the format, which extendr `main`
#      (b0cb8a81, extendr/extendr#1058) changes to a `"%s"` argument. A lock
#      that moves extendr has to come past that function: drop the doubling
#      for a release carrying b0cb8a81, then move the pin here and the `=`
#      pins in src/rust/Cargo.toml. Those keep `cargo update` - and the
#      ferx-core bump script's whole-graph update - from moving extendr on
#      their own; this check explains why when someone moves them by hand.
#
# `//` and `/* .. */` comments are stripped before every test, so naming a
# helper in prose is free. Neither strip knows about string literals, so the
# first `https://` inside a message would blind every check on that line - keep
# URLs in comments, where they cost nothing.
# Exercised by tools/test-check-glue-raise.sh.
#
# Usage: tools/check-glue-raise.sh [dir [lock]]
#   (defaults src/rust/src and src/rust/Cargo.lock)
# Exits 0 when the glue is in shape, 1 otherwise.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="${1:-$ROOT/src/rust/src}"
LOCK="${2:-$ROOT/src/rust/Cargo.lock}"
# The extendr release whose `throw_r_error` raise_verbatim escapes for.
EXTENDR_RAISE_VERSION="0.9.0"
EXTENDR_RAISE_SOURCE="registry+https://github.com/rust-lang/crates.io-index"
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
        trimmed = code
        sub(/^[[:space:]]+/, "", trimmed)
      }

      # The extent of `mod raise`, by brace depth. Not by a plain `}` test: the
      # extern block and raise_verbatim each close with one, and the module
      # would end three lines in.
      {
        if (in_raise) {
          raise_depth += gsub(/\{/, "{", code) - gsub(/\}/, "}", code)
          if (raise_depth <= 0) in_raise = 0
        } else if (trimmed ~ /^mod[[:space:]]+raise([[:space:]]|\{|$)/) {
          in_raise = 1
          raise_depth = gsub(/\{/, "{", code) - gsub(/\}/, "}", code)
        }
      }

      !in_raise && code ~ /Rf_error/ {
        printf "%s\t%d\tfinding\t%s\n", file, NR, "Rf_error is named outside mod raise. It is declared in there, privately, so that nothing but raise_verbatim can reach it - but a second extern \"C\" block re-declaring the symbol compiles, and hands the engine text to Rf_error as its format string again (#388)."
      }

      code ~ /throw_r_error[[:space:]]*\(/ {
        if (trimmed ~ /^fn[[:space:]]+throw_r_error[[:space:]]*\(/)
          printf "%s\t%d\t@shadow\t-\n", file, NR
        else
          printf "%s\t%d\tfinding\t%s\n", file, NR, "throw_r_error() is reachable here. Raising goes through entry(); the local fn of that name exists only so that a bare call does not compile."
      }

      code ~ /raise_verbatim[[:space:]]*\(/ \
        && trimmed !~ /^(pub(\([^)]*\))?[[:space:]]+)?fn[[:space:]]+raise_verbatim[[:space:]]*\(/ {
        printf "%s\t%d\t@raise\t-\n", file, NR
      }

      trimmed ~ /^#\[extendr(\]|\(|[[:space:]])/ {
        pending = 1; in_sig = 0; expect_entry = 0
        # `#[extendr(` may wrap over several lines; run to the one that closes
        # the attribute, or the signature below is never seen.
        if (trimmed !~ /\][[:space:]]*$/) in_attr = 1
        next
      }
      in_attr { if (trimmed ~ /\][[:space:]]*$/) in_attr = 0; next }

      pending && !in_sig && trimmed ~ /^impl([[:space:]]|<)/ {
        printf "%s\t%d\tfinding\t%s\n", file, NR, "#[extendr] on an impl block. This guard only understands #[extendr] fn - it would check only the first method in that block and skip every other one. Teach it that shape before using it here."
        pending = 0
        next
      }

      pending && !in_sig && trimmed ~ /^fn[[:space:]]/ {
        in_sig = 1; fn_line = NR
        name = trimmed
        sub(/^fn[[:space:]]+/, "", name)
        sub(/[^A-Za-z0-9_].*$/, "", name)
        fn_name = name
      }
      in_sig {
        # The first `{` on a signature line opens the body - a Rust signature
        # has no other brace. What follows it is either nothing (the body opens
        # on the next line) or the whole body, on this one.
        if (match(code, /\{/)) {
          rest = substr(code, RSTART + 1)
          sub(/^[[:space:]]+/, "", rest)
          in_sig = 0
          if (rest == "") {
            expect_entry = 1
          } else {
            if (rest !~ /^entry[[:space:]]*\(/)
              printf "%s\t%d\tfinding\t%s\n", file, fn_line, fn_name " does not open its body with entry(). Every #[extendr] body is entry(move || { ..; Ok(value) }), so the raise happens once, after the body has returned."
            pending = 0
          }
        }
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

# Check 5. `name = ".."` is followed by its `version = ".."` in every
# Cargo.lock package table, then by `source = ".."` - absent for a path
# package, which is what a `[patch.crates-io]` path replacement leaves.
lock_rel="${LOCK#"$ROOT"/}"
if [[ ! -f "$LOCK" ]]; then
  fail "$lock_rel" 1 "no Cargo.lock at $LOCK. raise_verbatim's % escape is right for one extendr release only, so the lock has to say which one."
else
  for crate in extendr-api extendr-macros; do
    # `version<TAB>source` of the crate's table, source "-" when there is none.
    entry=$(awk -v want="name = \"$crate\"" '
      $0 == want { found = 1; v = ""; s = "-"; next }
      found && $1 == "version" { v = $3; gsub(/"/, "", v); next }
      found && $1 == "source" { s = $3; gsub(/"/, "", s); next }
      found && ($0 == "" || $0 ~ /^\[/) { print v "\t" s; found = 0 }
      END { if (found) print v "\t" s }
    ' "$LOCK")
    got="${entry%%$'\t'*}"
    src="${entry#*$'\t'}"
    if [[ -z "$entry" ]]; then
      fail "$lock_rel" 1 "$crate is not in the lock. raise_verbatim unwinds into extendr's wrapper and escapes % for extendr-api $EXTENDR_RAISE_VERSION's throw_r_error; without extendr that is unchecked."
    elif [[ "$got" != "$EXTENDR_RAISE_VERSION" ]]; then
      fail "$lock_rel" 1 "$crate is $got in the lock, but raise::format_escaped doubles % for $EXTENDR_RAISE_VERSION, whose throw_r_error hands the text to Rf_error as its format. A release carrying extendr b0cb8a81 passes it as a \"%s\" argument, and every % would print as %%. Update format_escaped for $got, then the = pins in src/rust/Cargo.toml and EXTENDR_RAISE_VERSION here (#394)."
    elif [[ "$src" != "$EXTENDR_RAISE_SOURCE" ]]; then
      fail "$lock_rel" 1 "$crate $got comes from '$src' in the lock, not crates.io. raise::format_escaped and the unwind into extendr's wrapper are written against the published $EXTENDR_RAISE_VERSION; a git or path replacement at the same version can differ in both. Point it back at crates.io, or re-check format_escaped against that source and teach this check its source (#394)."
    fi
  done
fi

if [[ $status -eq 0 ]]; then
  echo "glue raise shape OK: entry() opens every #[extendr] body, one raise_verbatim() call site, throw_r_error() shadowed, Rf_error named only inside mod raise, extendr at $EXTENDR_RAISE_VERSION."
fi
exit $status
