#!/usr/bin/env bash
# Check that AGENTS.md is the repository's one agent-guidance file (#489,
# mirroring ferx-core #1528 / tests/agents_md_canonical.rs).
#
# Why: two copies of a rulebook drift the way two copies of a formula do. ferx-core
# carried a Claude-specific guidance file next to an AGENTS.md copied from it once,
# and three weeks later the copy described a source layout that no longer existed
# while every code comment pointed at the other file. So this checks that:
#
#   - AGENTS.md is tracked at the root and opens with `# AGENTS.md`;
#   - the old Claude-specific file is not back at the root, tracked or not, in any
#     letter case (a case-insensitive filesystem opens it under either), and not
#     even as a dangling symlink;
#   - no tracked file in a subdirectory carries the old name either: agents load
#     those as guidance for that directory, so each is a second rulebook;
#   - no tracked file names the old file, in any letter case, plain or
#     regex-escaped (the spelling .Rbuildignore used). NEWS.md is exempt: it is
#     history, and the entry recording the rename names the old file by necessity.
#
# Tracked files, not a directory walk: a walk would also see untracked and
# gitignored files that exist on one machine only (.claude/, local notes).
# The scan reads bytes (`git grep -a`), so a binary or non-UTF-8 file is scanned too.
# Paths travel NUL-separated (`-z`), so a path holding `:`, a quote or a tab is
# reported as it is, and in CI the annotation's properties are escaped.
#
# Runs no R and no cargo. Used by the R-CMD-check workflow and exercised by
# tools/test-check-agents-md.sh. Exits 0 when the guidance file is in shape, 1 otherwise.
#
# Usage: tools/check-agents-md.sh [repo root]

set -euo pipefail

ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
# Built at runtime, so this script does not match its own scan.
LEGACY="CLAUDE"".md"
# The same name as a regex spells it, dot escaped: that is how .Rbuildignore
# held it, and a fixed-string match on the plain name does not see it.
LEGACY_ESCAPED="CLAUDE""\\.md"
HISTORY=NEWS.md
status=0
# Command substitution drops NUL bytes, so the -z listings go through files.
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/check-agents-md.XXXXXX")"
trap 'rm -rf "$SCRATCH"' EXIT

lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }
LEGACY_LOWER="$(lower "$LEGACY")"

# GitHub workflow commands: a property value escapes %, CR, LF, `:` and `,`; the
# message escapes the first three. Unescaped, a path holding `,` or `:` ends the
# file= property early and the annotation lands on the wrong file.
gh_data() { local s=${1//%/%25}; s=${s//$'\r'/%0D}; printf '%s' "${s//$'\n'/%0A}"; }
gh_prop() { local s; s=$(gh_data "$1"); s=${s//:/%3A}; printf '%s' "${s//,/%2C}"; }

# Runs `git -C ROOT <args>` into a file. Exit 1 is "no match" for grep and never
# happens for ls-files; anything above it is an error, which must not read as clean.
git_to() { # outfile, git args...
  local out=$1 rc=0
  shift
  git -C "$ROOT" "$@" > "$out" || rc=$?
  if [[ "$rc" -gt 1 ]]; then
    echo "git $1 failed (exit $rc) in $ROOT" >&2
    exit 1
  fi
}

fail() { # file, line, message
  if [[ "${GITHUB_ACTIONS:-}" == true ]]; then
    echo "::error file=$(gh_prop "$1"),line=$(gh_prop "$2")::$(gh_data "$3")"
  else
    echo "$1:$2: $3"
  fi
  status=1
}

if ! git -C "$ROOT" ls-files --error-unmatch -- AGENTS.md >/dev/null 2>&1; then
  fail AGENTS.md 1 "AGENTS.md is not tracked at the repository root: it is the one agent-guidance file (#489)"
elif [[ "$(head -n 1 "$ROOT/AGENTS.md")" != "# AGENTS.md" ]]; then
  fail AGENTS.md 1 "AGENTS.md must open with its own name as the title: \`# AGENTS.md\`"
fi

# The root, on disk: tracked or not, any letter case, a dangling symlink included
# (the glob lists the link itself, whatever it points at).
for entry in "$ROOT"/*; do
  name=${entry##*/}
  if [[ "$(lower "$name")" == "$LEGACY_LOWER" ]] && [[ -e "$entry" || -L "$entry" ]]; then
    fail "$name" 1 "$name is back at the repository root: keep one guidance file, AGENTS.md (#489)"
  fi
done

# Below the root: every tracked path, by basename, any letter case.
git_to "$SCRATCH/tracked" ls-files -z
while IFS= read -r -d '' path; do
  [[ "$path" == */* ]] || continue
  name=${path##*/}
  if [[ "$(lower "$name")" == "$LEGACY_LOWER" ]]; then
    fail "$path" 1 "$name in a subdirectory is loaded as guidance for it: keep one guidance file, AGENTS.md (#489)"
  fi
done < "$SCRATCH/tracked"

# Contents. `-z -n -o` prints `path NUL line NUL match LF` per match, the path
# verbatim; `-o` keeps a binary line's own NULs out of the record.
git_to "$SCRATCH/hits" grep -z -n -o -i -a -F --no-color -e "$LEGACY" -e "$LEGACY_ESCAPED" \
  -- . ":(exclude)$HISTORY"
seen=""
while IFS= read -r -d '' file && IFS= read -r -d '' line && IFS= read -r match; do
  # Two matches on one line are one finding.
  [[ "$file:$line" == "$seen" ]] && continue
  seen="$file:$line"
  fail "$file" "$line" "names $match, which was renamed to AGENTS.md in #489"
done < "$SCRATCH/hits"

if [[ "$status" -eq 0 ]]; then
  echo "AGENTS.md is the one agent-guidance file: OK"
fi
exit "$status"
