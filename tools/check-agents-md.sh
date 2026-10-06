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
#   - the old Claude-specific file is not back at the root (not even as a symlink);
#   - no tracked file names the old file. NEWS.md is exempt: it is history, and
#     the entry recording the rename names the old file by necessity.
#
# Tracked files, not a directory walk: a walk would also see untracked and
# gitignored files that exist on one machine only (.claude/, local notes).
# The scan reads bytes (`git grep -a`), so a binary or non-UTF-8 file is scanned too.
#
# Runs no R and no cargo. Used by the R-CMD-check workflow and exercised by
# tools/test-check-agents-md.sh. Exits 0 when the guidance file is in shape, 1 otherwise.
#
# Usage: tools/check-agents-md.sh [repo root]

set -euo pipefail

ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
# Built at runtime, so this script does not match its own scan.
LEGACY="CLAUDE"".md"
HISTORY=NEWS.md
status=0

fail() { # file, line, message
  if [[ "${GITHUB_ACTIONS:-}" == true ]]; then
    echo "::error file=$1,line=$2::$3"
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

if [[ -e "$ROOT/$LEGACY" || -L "$ROOT/$LEGACY" ]]; then
  fail "$LEGACY" 1 "$LEGACY is back at the repository root: keep one guidance file, AGENTS.md (#489)"
fi

# git grep exits 1 for "no match" and >1 for an error; an error must not read as clean.
rc=0
hits=$(git -C "$ROOT" grep -n -a -F --no-color -e "$LEGACY" -- . ":(exclude)$HISTORY") || rc=$?
if [[ "$rc" -gt 1 ]]; then
  echo "git grep failed (exit $rc) in $ROOT" >&2
  exit 1
fi
while IFS= read -r hit; do
  [[ -z "$hit" ]] && continue
  file=${hit%%:*}
  rest=${hit#*:}
  line=${rest%%:*}
  fail "$file" "$line" "names $LEGACY, which was renamed to AGENTS.md in #489"
done <<< "$hits"

if [[ "$status" -eq 0 ]]; then
  echo "AGENTS.md is the one agent-guidance file: OK"
fi
exit "$status"
