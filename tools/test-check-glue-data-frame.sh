#!/usr/bin/env bash
# Drive tools/check-glue-data-frame.sh over stub glue sources, one per shape it
# must catch and one per shape it must not. The real lib.rs is clean, so it
# would stay green if the check's pattern were broken; the stubs are what
# redden. Each case runs with GITHUB_ACTIONS unset and set (CI takes the
# `::error` branch). Output is printed only for a failing case, indented, so no
# `::error` line becomes an annotation.
#
# Usage: tools/test-check-glue-data-frame.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="$ROOT/tools/check-glue-data-frame.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
failures=0
cases=0

# case <name> <expected exit> <expected text or ""> <source>
case_() {
  local name="$1" want="$2" text="$3" src="$4"
  local dir="$TMP/$name"
  mkdir -p "$dir"
  printf '%s\n' "$src" > "$dir/lib.rs"
  for gha in "" true; do
    cases=$((cases + 1))
    local out rc=0
    out="$(GITHUB_ACTIONS="$gha" bash "$CHECK" "$dir" 2>&1)" || rc=$?
    if [[ $rc -ne $want ]] || { [[ -n "$text" ]] && ! grep -qF -- "$text" <<<"$out"; }; then
      echo "FAIL: $name (GITHUB_ACTIONS='$gha'): exit $rc, wanted $want${text:+ with '$text'}"
      printf '%s\n' "$out" | sed 's/^/     | /'
      failures=$((failures + 1))
    fi
  done
}

CLEAN='use extendr_api::prelude::*;
fn r_data_frame(columns: Vec<(&str, Robj)>, nrow: usize) -> Robj { todo!() }
fn table() -> Robj { r_data_frame(vec![("ID", 1.into())], 1) }'

case_ clean 0 "glue data frames OK" "$CLEAN"
case_ prose_line_comment 0 "" "$CLEAN
// extendr's data_frame!(ID = id) evaluates data.frame in the global env.
/// Unlike \`data_frame!(..)\`, this does not call R."
case_ prose_block_comment 0 "" "$CLEAN
/* data_frame!(
   ID = id) */
fn other() {}"
case_ a_different_macro 0 "" "$CLEAN
fn f() { my_data_frame!(x = 1); data_frames(); }"

case_ paren_call 1 "data_frame! builds the table" "$CLEAN
fn f() -> Robj { data_frame!(ID = id, TIME = time).into() }"
case_ multiline_call 1 "data_frame! builds the table" "$CLEAN
fn f() -> Robj {
    let df = data_frame!(
        ID = id
    );
}"
case_ spaced_call 1 "data_frame! builds the table" "$CLEAN
fn f() { let d = data_frame! (a = 1); }"
case_ bracket_call 1 "data_frame! builds the table" "$CLEAN
fn f() { let d = data_frame![a = 1]; }"
case_ qualified_call 1 "data_frame! builds the table" "$CLEAN
fn f() { let d = extendr_api::data_frame!(a = 1); }"
case_ after_block_comment 1 "data_frame! builds the table" "$CLEAN
fn f() { /* note */ let d = data_frame!(a = 1); }"
case_ import 1 "data_frame is imported here" "$CLEAN
use extendr_api::data_frame;"
case_ import_as 1 "data_frame is imported here" "$CLEAN
use extendr_api::{data_frame as df, List};"

# The finding names the line of the call, not of the statement it opens -
# locally as `file:line:` and in CI as an `::error` annotation.
mkdir -p "$TMP/where"
printf '%s\n' "$CLEAN" 'fn f() {' '    let d = data_frame!(' '        a = 1' '    );' '}' > "$TMP/where/lib.rs"
for want in "lib.rs:5: data_frame!" "::error file="; do
  cases=$((cases + 1))
  if [[ "$want" == ::* ]]; then
    out="$(GITHUB_ACTIONS=true bash "$CHECK" "$TMP/where" 2>/dev/null || true)"
    want_re='^::error file=[^,]*lib\.rs,line=5::data_frame!'
  else
    out="$(GITHUB_ACTIONS='' bash "$CHECK" "$TMP/where" 2>&1 || true)"
    want_re='lib\.rs:5: data_frame!'
  fi
  if ! grep -qE -- "$want_re" <<<"$out"; then
    echo "FAIL: location ($want): no line matching /$want_re/"
    printf '%s\n' "$out" | sed 's/^/     | /'
    failures=$((failures + 1))
  fi
done

# No sources at all is an error, not a pass.
cases=$((cases + 1))
mkdir -p "$TMP/empty"
if bash "$CHECK" "$TMP/empty" >/dev/null 2>&1; then
  echo "FAIL: empty dir passed"
  failures=$((failures + 1))
fi

if [[ $failures -gt 0 ]]; then
  echo "$failures of $cases glue data-frame guard case(s) failed."
  exit 1
fi
echo "glue data-frame guard: all $cases cases behave."
