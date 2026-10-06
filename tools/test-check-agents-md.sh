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
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
# Built at runtime, so this harness does not trip the check it drives.
LEGACY="CLAUDE"".md"
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

expect() { # name, dir, expected exit, message fragment[, file:line]
  local name=$1 dir=$2 want=$3 fragment=$4 at=${5:-}
  local mode out rc where
  for mode in local ci; do
    rc=0
    if [[ "$mode" == ci ]]; then
      out=$(GITHUB_ACTIONS=true bash "$CHECK" "$dir" 2>&1) || rc=$?
    else
      out=$(env -u GITHUB_ACTIONS bash "$CHECK" "$dir" 2>&1) || rc=$?
    fi
    if [[ "$rc" != "$want" ]]; then
      report "FAIL [$mode] $name: expected exit $want, got $rc" "$out"
      continue
    fi
    if ! printf '%s' "$out" | grep -qF -- "$fragment"; then
      report "FAIL [$mode] $name: output does not mention \"$fragment\"" "$out"
      continue
    fi
    # The location, in the form each mode prints it.
    if [[ -n "$at" ]]; then
      if [[ "$mode" == ci ]]; then where="::error file=${at%%:*},line=${at#*:}::"; else where="$at: "; fi
      if ! printf '%s' "$out" | grep -qF -- "$where"; then
        report "FAIL [$mode] $name: output does not locate the hit as \"$where\"" "$out"
        continue
      fi
    fi
    echo "ok   [$mode] $name"
  done
}

OK="AGENTS.md is the one agent-guidance file: OK"

# -- green --------------------------------------------------------------------

expect "a repository in shape passes" "$(fixture in-shape)" 0 "$OK"

expect "this checkout passes" "$ROOT" 0 "$OK"

# Untracked files exist on one machine only; the scan must not see them.
d=$(fixture untracked-mention)
printf 'see %s\n' "$LEGACY" > "$d/notes.txt"
expect "an untracked file naming the old file is not scanned" "$d" 0 "$OK"

# -- red ----------------------------------------------------------------------

d=$(fixture tracked-mention)
printf 'f <- function() 1\n# The label convention (%s).\n' "$LEGACY" > "$d/R/f.R"
expect "a tracked file naming the old file fails, with file and line" "$d" 1 \
  "names $LEGACY, which was renamed to AGENTS.md in #489" R/f.R:2

# The spelling the one reference .Rbuildignore held used: a regex, dot escaped.
d=$(fixture escaped-mention)
printf '^\\.github$\n^%s\\.md$\n' "CLAUDE" > "$d/.Rbuildignore"
git -C "$d" add -A
expect "a regex-escaped mention of the old file fails" "$d" 1 \
  "names $LEGACY, which was renamed to AGENTS.md in #489" .Rbuildignore:2

# A history exemption keyed on a basename would let this through.
d=$(fixture nested-news)
mkdir -p "$d/docs"
printf 'see %s\n' "$LEGACY" > "$d/docs/NEWS.md"
git -C "$d" add -A
expect "only the root NEWS.md is exempt" "$d" 1 \
  "names $LEGACY, which was renamed to AGENTS.md in #489" docs/NEWS.md:1

# One invalid UTF-8 byte must not hide the rest of the file from the scan.
d=$(fixture binary-mention)
printf '\377\376\000bytes\nsee %s\n' "$LEGACY" > "$d/blob.bin"
git -C "$d" add -A
expect "a binary file naming the old file fails" "$d" 1 \
  "names $LEGACY, which was renamed to AGENTS.md in #489" blob.bin:2

BACK="$LEGACY is back at the repository root: keep one guidance file, AGENTS.md (#489)"
UNTRACKED="AGENTS.md is not tracked at the repository root: it is the one agent-guidance file (#489)"

d=$(fixture legacy-back)
printf '# %s\n' "$LEGACY" > "$d/$LEGACY"
expect "the old file recreated at the root fails, even untracked" "$d" 1 "$BACK" "$LEGACY:1"

# Dangling, so `-e` (which follows the link) is false and only `-L` sees it.
d=$(fixture legacy-symlink)
ln -s nowhere "$d/$LEGACY"
expect "a dangling symlink under the old name fails" "$d" 1 "$BACK"

d=$(fixture wrong-heading)
printf '# Guidance\n\nGuidance for agents.\n' > "$d/AGENTS.md"
expect "AGENTS.md without its title fails" "$d" 1 \
  'AGENTS.md must open with its own name as the title: `# AGENTS.md`' AGENTS.md:1

d=$(fixture agents-missing)
git -C "$d" rm -q --cached AGENTS.md
rm "$d/AGENTS.md"
expect "a repository without AGENTS.md fails" "$d" 1 "$UNTRACKED"

d=$(fixture agents-untracked)
git -C "$d" rm -q --cached AGENTS.md
expect "an AGENTS.md on disk but not in git fails" "$d" 1 "$UNTRACKED"

# Not a repository at all: git grep errors, and that must not read as "no match".
d="$TMP/not-a-repo"
mkdir -p "$d"
printf '# AGENTS.md\n' > "$d/AGENTS.md"
expect "a directory outside git fails rather than passing vacuously" "$d" 1 \
  "git grep failed"

if [[ "$failures" -gt 0 ]]; then
  echo "$failures case(s) failed"
  exit 1
fi
echo "all cases passed"
