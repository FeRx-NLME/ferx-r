#!/usr/bin/env bash
# Feed tools/check-ferx-core-pin.sh the lock states a local sibling build can
# leave behind, derived from the committed src/rust/Cargo.lock, and assert each
# verdict. CI never produces a damaged lock itself, so without this nothing would
# fail if one of the script's checks were deleted.
#
# Every bad fixture must exit 1 *with its own message*, so a check that fires for
# the wrong reason does not pass. The checks run with GITHUB_ACTIONS unset and
# their output is captured and indented, so the `::error` lines they print for bad
# fixtures never become annotations on a green run.
#
# Usage: tools/test-check-ferx-core-pin.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="$ROOT/tools/check-ferx-core-pin.sh"
LOCK="$ROOT/src/rust/Cargo.lock"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
failures=0

# A literal prefix, matched with index(): `awk -v` would process backslash
# escapes in a regex, and `+` / `.` in it would stop meaning themselves.
PIN_PREFIX='source = "git+https://github.com/FeRx-NLME/ferx-core'

expect() { # name, lockfile, expected exit, expected message fragment ("" = none)
  local name=$1 lock=$2 want=$3 fragment=$4 out rc=0
  out=$(env -u GITHUB_ACTIONS bash "$CHECK" "$lock" 2>&1) || rc=$?
  if [[ "$rc" -ne "$want" ]]; then
    echo "FAIL $name: exit $rc, expected $want"; failures=$((failures + 1))
  elif [[ -n "$fragment" && "$out" != *"$fragment"* ]]; then
    echo "FAIL $name: exit $rc as expected, but without \"$fragment\""; failures=$((failures + 1))
  else
    echo "ok   $name"; return 0
  fi
  printf '%s\n' "$out" | sed 's/^/     | /'
}

# The committed lock must itself be pinned, or every fixture below is meaningless.
expect "committed lock" "$LOCK" 0 "ferx-core and ferx-tools pinned to"

# An applied [patch]: cargo drops both source lines.
awk -v pfx="$PIN_PREFIX" 'index($0, pfx) != 1' "$LOCK" > "$TMP/stripped.lock"
expect "both pins stripped" "$TMP/stripped.lock" 1 "ferx-core lacks a git source pin"

# Only ferx-tools stripped: must not be excused by the revision check skipping it.
awk -v pfx="$PIN_PREFIX" '
  $0 == "name = \"ferx-tools\"" { tools = 1 }
  tools && index($0, pfx) == 1 { tools = 0; next }
  { print }
' "$LOCK" > "$TMP/tools-stripped.lock"
expect "ferx-tools pin stripped" "$TMP/tools-stripped.lock" 1 "ferx-tools lacks a git source pin"

# A git pin, but to another repository (a fork): not the pin CI must build.
sed 's#github.com/FeRx-NLME/ferx-core#github.com/someone-else/ferx-core#' "$LOCK" > "$TMP/fork.lock"
expect "pinned to a fork" "$TMP/fork.lock" 1 "ferx-core lacks a git source pin"

# An unused [patch]: pins intact, [[patch.unused]] appended.
{ cat "$LOCK"; printf '\n[[patch.unused]]\nname = "ferx-core"\nversion = "0.0.0"\n'; } > "$TMP/unused.lock"
expect "[[patch.unused]] appended" "$TMP/unused.lock" 1 "[[patch.unused]]"

# The two crates on different revisions of the one repository.
awk -v pfx="$PIN_PREFIX" '
  $0 == "name = \"ferx-tools\"" { tools = 1 }
  tools && index($0, pfx) == 1 { sub(/#[0-9a-f]+/, "#0000000000000000000000000000000000000000"); tools = 0 }
  { print }
' "$LOCK" > "$TMP/split.lock"
expect "revisions split" "$TMP/split.lock" 1 "pinned to different revisions"

expect "missing lock" "$TMP/does-not-exist.lock" 1 "not found"

if [[ "$failures" -ne 0 ]]; then
  echo "$failures check-ferx-core-pin.sh case(s) failed"
  exit 1
fi
echo "all check-ferx-core-pin.sh cases passed"
