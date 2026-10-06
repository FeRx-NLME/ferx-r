#!/usr/bin/env bash
# Drive tools/check-agents-md.sh over stub repositories, one per shape it is
# supposed to catch and one per shape it must not.
#
# The real checkout is always in shape, so nothing there would fail if one of
# the checks were deleted. Each fixture is a throwaway `git init` whose index
# holds a guidance file in good shape plus one change, so every check has a
# fixture that reddens only it - and every bad fixture must exit 1 with its own
# message, so a check firing for the wrong reason does not pass. Each case runs
# twice - GITHUB_ACTIONS unset and set - because CI takes the `::error` branch.
# Output is captured and printed only for a failing case, indented, so no
# `::error` line becomes an annotation.
#
# Usage: tools/test-check-agents-md.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="$ROOT/tools/check-agents-md.sh"
# An explicit template, so TMPDIR is honoured on BSD mktemp (macOS) as on GNU.
TMP="$(mktemp -d "${TMPDIR:-/tmp}/test-check-agents-md.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
# A fixture that is not a repository must stay one even when TMPDIR sits inside
# a git checkout: git's discovery stops at $TMP instead of climbing into it.
export GIT_CEILING_DIRECTORIES="$TMP"
# Built at runtime, so this harness does not trip the check it drives.
LEGACY="CLAUDE"".md"
LEGACY_ESCAPED="CLAUDE""\\.md"
LOWER="claude"".md"
failures=0

report() { # headline, captured output
  echo "$1"
  failures=$((failures + 1))
  printf '%s\n' "$2" | sed 's/^/     | /'
}

# A repository in shape: AGENTS.md with its title, a NEWS.md entry recording
# the rename (history, exempt), and an ordinary source file. Everything in the
# index; the caller then applies its one change.
fixture() { # name
  local dir="$TMP/$1"
  mkdir -p "$dir/R"
  git -C "$dir" init -q
  printf '# AGENTS.md\n\nGuidance for agents.\n' > "$dir/AGENTS.md"
  printf '# ferx\n\n- Renamed %s to AGENTS.md (#489).\n' "$LEGACY" > "$dir/NEWS.md"
  printf 'f <- function() 1\n' > "$dir/R/f.R"
  git -C "$dir" add -A
  echo "$dir"
}

run_check() { # mode, dir -> sets out, rc
  rc=0
  if [[ "$1" == ci ]]; then
    out=$(GITHUB_ACTIONS=true bash "$CHECK" "$2" 2>&1) || rc=$?
  else
    out=$(env -u GITHUB_ACTIONS bash "$CHECK" "$2" 2>&1) || rc=$?
  fi
}

# file and line, when given, must appear as each mode prints a location:
# `file:line: ` locally, `::error file=...,line=...::` in CI, where the file is
# escaped (ci_file, defaulting to file).
expect() { # name, dir, expected exit, message fragment[, file, line[, ci_file]]
  local name=$1 dir=$2 want=$3 fragment=$4 file=${5:-} line=${6:-} ci_file=${7:-${5:-}}
  local mode out rc where
  for mode in local ci; do
    run_check "$mode" "$dir"
    if [[ "$rc" != "$want" ]]; then
      report "FAIL [$mode] $name: expected exit $want, got $rc" "$out"
      continue
    fi
    if ! printf '%s' "$out" | grep -qF -- "$fragment"; then
      report "FAIL [$mode] $name: output does not mention \"$fragment\"" "$out"
      continue
    fi
    if [[ -n "$file" ]]; then
      if [[ "$mode" == ci ]]; then where="::error file=$ci_file,line=$line::"; else where="$file:$line: "; fi
      if ! printf '%s' "$out" | grep -qF -- "$where"; then
        report "FAIL [$mode] $name: output does not locate the hit as \"$where\"" "$out"
        continue
      fi
    fi
    echo "ok   [$mode] $name"
  done
}

expect_count() { # name, dir, fragment, number of output lines holding it
  local name=$1 dir=$2 fragment=$3 want=$4
  local mode out rc n
  for mode in local ci; do
    run_check "$mode" "$dir"
    n=$(printf '%s\n' "$out" | grep -cF -- "$fragment" || true)
    if [[ "$n" != "$want" ]]; then
      report "FAIL [$mode] $name: $n lines mention \"$fragment\", expected $want" "$out"
      continue
    fi
    echo "ok   [$mode] $name"
  done
}

OK="AGENTS.md is the one agent-guidance file: OK"
RENAMED="which was renamed to AGENTS.md in #489"
BACK="is back at the repository root: keep one guidance file, AGENTS.md (#489)"
NESTED="in a subdirectory is loaded as guidance for it: keep one guidance file, AGENTS.md (#489)"
UNTRACKED="AGENTS.md is not tracked at the repository root: it is the one agent-guidance file (#489)"

# -- green --------------------------------------------------------------------

expect "a repository in shape passes" "$(fixture in-shape)" 0 "$OK"

expect "this checkout passes" "$ROOT" 0 "$OK"

# Untracked files exist on one machine only; the scan must not see them.
d=$(fixture untracked-mention)
printf 'see %s\n' "$LEGACY" > "$d/notes.txt"
mkdir -p "$d/scratch"
printf 'x\n' > "$d/scratch/$LEGACY"
expect "an untracked file naming the old file, or under it below the root, is not scanned" "$d" 0 "$OK"

# -- red: contents ------------------------------------------------------------

d=$(fixture tracked-mention)
printf 'f <- function() 1\n# The label convention (%s).\n' "$LEGACY" > "$d/R/f.R"
expect "a tracked file naming the old file fails, with file and line" "$d" 1 \
  "names $LEGACY, $RENAMED" R/f.R 2

# The spelling the one reference .Rbuildignore held used: a regex, dot escaped.
d=$(fixture escaped-mention)
printf '^\\.github$\n^%s$\n' "$LEGACY_ESCAPED" > "$d/.Rbuildignore"
git -C "$d" add -A
expect "a regex-escaped mention of the old file fails" "$d" 1 \
  "names $LEGACY_ESCAPED, $RENAMED" .Rbuildignore 2

d=$(fixture lowercase-mention)
printf 'f <- function() 1\n# see %s\n' "$LOWER" > "$d/R/f.R"
expect "a lowercase mention of the old file fails" "$d" 1 \
  "names $LOWER, $RENAMED" R/f.R 2

# A history exemption keyed on a basename would let this through.
d=$(fixture nested-news)
mkdir -p "$d/docs"
printf 'see %s\n' "$LEGACY" > "$d/docs/NEWS.md"
git -C "$d" add -A
expect "only the root NEWS.md is exempt" "$d" 1 \
  "names $LEGACY, $RENAMED" docs/NEWS.md 1

# One invalid UTF-8 byte must not hide the rest of the file from the scan.
d=$(fixture binary-mention)
printf '\377\376\000bytes\nsee %s\n' "$LEGACY" > "$d/blob.bin"
git -C "$d" add -A
expect "a binary file naming the old file fails" "$d" 1 \
  "names $LEGACY, $RENAMED" blob.bin 2

# Without -z git quotes such a path ("sub:dir/we\"ird,na%me.R"), and splitting
# at the first `:` cuts it in two; in CI, `:` and `,` would end the file=
# property, and an unescaped `%` would be read as the start of an escape.
d=$(fixture awkward-path)
mkdir -p "$d/sub:dir"
printf 'x\nsee %s\n' "$LEGACY" > "$d/sub:dir/we\"ird,na%me.R"
git -C "$d" add -A
expect "a path holding \`:\`, \`,\`, \`%\` and a quote is reported verbatim, escaped in CI" "$d" 1 \
  "names $LEGACY, $RENAMED" 'sub:dir/we"ird,na%me.R' 2 'sub%3Adir/we"ird%2Cna%25me.R'

d=$(fixture two-on-a-line)
printf 'f <- function() 1\n# %s, or %s\n' "$LEGACY" "$LEGACY_ESCAPED" > "$d/R/f.R"
expect_count "two mentions on one line are one finding" "$d" "$RENAMED" 1

# -- red: the file itself -----------------------------------------------------

d=$(fixture legacy-back)
printf '# %s\n' "$LEGACY" > "$d/$LEGACY"
expect "the old file recreated at the root fails, even untracked" "$d" 1 \
  "$LEGACY $BACK" "$LEGACY" 1

# Dangling, so `-e` (which follows the link) is false and only `-L` sees it.
d=$(fixture legacy-symlink)
ln -s nowhere "$d/$LEGACY"
expect "a dangling symlink under the old name fails" "$d" 1 "$LEGACY $BACK"

# Matched by name, not by `-e`: on a case-insensitive filesystem `-e` would find
# it under the old spelling too, on Linux it would not.
d=$(fixture legacy-lowercase)
printf 'x\n' > "$d/$LOWER"
expect "the old file at the root in lowercase fails" "$d" 1 "$LOWER $BACK" "$LOWER" 1

d=$(fixture legacy-nested)
mkdir -p "$d/inst"
printf 'x\n' > "$d/inst/$LEGACY"
git -C "$d" add -A
expect "a tracked file under the old name in a subdirectory fails" "$d" 1 \
  "$LEGACY $NESTED" "inst/$LEGACY" 1

d=$(fixture legacy-nested-lowercase)
mkdir -p "$d/a/b"
printf 'x\n' > "$d/a/b/$LOWER"
git -C "$d" add -A
expect "a tracked lowercase file under the old name, two levels down, fails" "$d" 1 \
  "$LOWER $NESTED" "a/b/$LOWER" 1

# -- red: AGENTS.md -----------------------------------------------------------

d=$(fixture wrong-heading)
printf '# Guidance\n\nGuidance for agents.\n' > "$d/AGENTS.md"
expect "AGENTS.md without its title fails" "$d" 1 \
  'AGENTS.md must open with its own name as the title: `# AGENTS.md`' AGENTS.md 1

d=$(fixture agents-missing)
git -C "$d" rm -q --cached AGENTS.md
rm "$d/AGENTS.md"
expect "a repository without AGENTS.md fails" "$d" 1 "$UNTRACKED"

d=$(fixture agents-untracked)
git -C "$d" rm -q --cached AGENTS.md
expect "an AGENTS.md on disk but not in git fails" "$d" 1 "$UNTRACKED"

# Not a repository at all: git errors, and that must not read as "no match".
d="$TMP/not-a-repo"
mkdir -p "$d"
printf '# AGENTS.md\n' > "$d/AGENTS.md"
expect "a directory outside git fails rather than passing vacuously" "$d" 1 \
  "git ls-files failed (exit 128) in $d"

if [[ "$failures" -gt 0 ]]; then
  echo "$failures case(s) failed"
  exit 1
fi
echo "all cases passed"
