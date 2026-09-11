#!/usr/bin/env bash
# Feed tools/check-ferx-core-pin.sh the lock states a local sibling build can
# leave behind, derived from the committed src/rust/Cargo.lock, and assert each
# verdict. CI never produces a damaged lock itself, so without this nothing would
# fail if one of the script's checks were deleted.
#
# Every bad fixture must exit 1 *with its own message* and with repair advice, so
# a check that fires for the wrong reason does not pass. Each case runs twice:
# with GITHUB_ACTIONS unset (plain `error:` lines, local advice) and with
# GITHUB_ACTIONS=true (`::error` lines, CI advice), because CI's own guard step
# takes the second branch. Output is captured and printed only for a failing
# case, indented, so no `::error` line becomes an annotation.
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

report() { # headline, captured output
  echo "$1"
  failures=$((failures + 1))
  printf '%s\n' "$2" | sed 's/^/     | /'
}

expect() { # name, lockfile, expected exit, message fragment, advice|no-advice
  local name=$1 lock=$2 want=$3 fragment=$4 advice_wanted=$5 mode out rc advice
  for mode in local ci; do
    rc=0
    if [[ "$mode" == ci ]]; then
      out=$(GITHUB_ACTIONS=true bash "$CHECK" "$lock" 2>&1) || rc=$?
      advice="The committed Cargo.lock is broken"
    else
      out=$(env -u GITHUB_ACTIONS bash "$CHECK" "$lock" 2>&1) || rc=$?
      advice="Read \`git diff src/rust/Cargo.lock\` before restoring"
    fi
    if [[ "$rc" -ne "$want" ]]; then
      report "FAIL $name ($mode): exit $rc, expected $want" "$out"
    elif [[ "$out" != *"$fragment"* ]]; then
      report "FAIL $name ($mode): exit $rc as expected, but without \"$fragment\"" "$out"
    elif [[ "$advice_wanted" == advice && "$out" != *"$advice"* ]]; then
      report "FAIL $name ($mode): no repair advice (\"$advice\")" "$out"
    else
      echo "ok   $name ($mode)"
    fi
  done
}

# The committed lock must itself be pinned, or every fixture below is meaningless.
expect "committed lock" "$LOCK" 0 "ferx-core and ferx-tools pinned to" no-advice

# An applied [patch]: cargo drops both source lines.
awk -v pfx="$PIN_PREFIX" 'index($0, pfx) != 1' "$LOCK" > "$TMP/stripped.lock"
expect "both pins stripped" "$TMP/stripped.lock" 1 "ferx-core lacks a git source pin" advice

# One pin stripped at a time. ferx-tools' table follows ferx-core's, so a scan that
# ran past the end of ferx-core's table would borrow ferx-tools' pin; and a
# stripped ferx-tools must not be excused by the revision check skipping it.
for crate in ferx-core ferx-tools; do
  awk -v pfx="$PIN_PREFIX" -v crate="$crate" '
    $0 == "name = \"" crate "\"" { hit = 1 }
    hit && index($0, pfx) == 1 { hit = 0; next }
    { print }
  ' "$LOCK" > "$TMP/$crate-stripped.lock"
  expect "$crate pin stripped" "$TMP/$crate-stripped.lock" 1 "$crate lacks a git source pin" advice
done

# A git pin, but to another repository (a fork): not the pin CI must build.
sed 's#github.com/FeRx-NLME/ferx-core#github.com/someone-else/ferx-core#' "$LOCK" > "$TMP/fork.lock"
expect "pinned to a fork" "$TMP/fork.lock" 1 "ferx-core lacks a git source pin" advice

# An unused [patch]: pins intact, [[patch.unused]] appended.
{ cat "$LOCK"; printf '\n[[patch.unused]]\nname = "ferx-core"\nversion = "0.0.0"\n'; } > "$TMP/unused.lock"
expect "[[patch.unused]] appended" "$TMP/unused.lock" 1 "[[patch.unused]]" advice

# The two crates on different revisions of the one repository.
awk -v pfx="$PIN_PREFIX" '
  $0 == "name = \"ferx-tools\"" { tools = 1 }
  tools && index($0, pfx) == 1 { sub(/#[0-9a-f]+/, "#0000000000000000000000000000000000000000"); tools = 0 }
  { print }
' "$LOCK" > "$TMP/split.lock"
expect "revisions split" "$TMP/split.lock" 1 "pinned to different revisions" advice

expect "missing lock" "$TMP/does-not-exist.lock" 1 "not found" no-advice

if [[ "$failures" -ne 0 ]]; then
  echo "$failures check-ferx-core-pin.sh case(s) failed"
  exit 1
fi
echo "all check-ferx-core-pin.sh cases passed"
