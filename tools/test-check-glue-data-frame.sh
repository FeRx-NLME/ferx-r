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
    local out rc=0 expect="$text"
    # A location is written `lib.rs:N: ` locally and `lib.rs,line=N::` in CI.
    if [[ "$gha" == true && "$expect" =~ ^lib\.rs:([0-9]+):\ (.*)$ ]]; then
      expect="lib.rs,line=${BASH_REMATCH[1]}::${BASH_REMATCH[2]}"
    fi
    out="$(GITHUB_ACTIONS="$gha" bash "$CHECK" "$dir" 2>&1)" || rc=$?
    if [[ $rc -ne $want ]] || { [[ -n "$expect" ]] && ! grep -qF -- "$expect" <<<"$out"; }; then
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
case_ import 1 "data_frame is named here" "$CLEAN
use extendr_api::data_frame;"
case_ import_as 1 "data_frame is named here" "$CLEAN
use extendr_api::{data_frame as df, List};"
# #781 review finding 3: the rename split over lines, then used under its
# new name, which no `data_frame!` rule can see.
case_ import_as_multiline 1 "lib.rs:6: data_frame is named here" "$CLEAN
use extendr_api::{
    List,
    data_frame as df,
};
fn f() { df!(a = 1); }"

# -- What the scan must treat as not-code, and what it must not lose track
#    in (#781 review findings 1, 2, 8). Each "hit" case puts the call on a
#    line after the tricky syntax, so a scan that went blind there passes it.
NL=$'\n'
case_ slashstar_in_line_comment 1 "lib.rs:5: data_frame!" "$CLEAN$NL"'// scans tools/*.sh'"$NL"'fn f() { data_frame!(a = 1); }'
case_ slashstar_in_doc_comment 1 "lib.rs:5: data_frame!" "$CLEAN$NL"'/// reads inst/examples/*.R'"$NL"'fn f() { data_frame!(a = 1); }'
case_ slashstar_in_string 1 "lib.rs:5: data_frame!" "$CLEAN$NL"'fn g() { let p = "inst/*.R"; }'"$NL"'fn f() { data_frame!(a = 1); }'
case_ slashslash_in_string 1 "lib.rs:4: data_frame!" "$CLEAN$NL"'fn f() { let u = "http://x"; let d = data_frame!(a = 1); }'
case_ escaped_quote_in_string 1 "lib.rs:4: data_frame!" "$CLEAN$NL"'fn f() { let s = "a \"data_frame!(\" b"; let d = data_frame!(a = 1); }'
case_ quote_char_literal 1 "lib.rs:4: data_frame!" "$CLEAN$NL"$'fn f() { let q = \'"\'; let d = data_frame!(a = 1); }'
case_ escaped_quote_in_string_only 0 "glue data frames OK" "$CLEAN$NL"'fn f() { let s = "a \" data_frame!(x) \" b"; }'
case_ escaped_char_literal 1 "lib.rs:4: data_frame!" "$CLEAN$NL"$'fn f() { let q = \'\\"\'; let d = data_frame!(a = 1); }'
case_ lifetime_then_call 1 "lib.rs:4: data_frame!" "$CLEAN$NL"$'fn f<\'a>(x: &\'a str) { let d = data_frame!(a = 1); }'
case_ raw_string_with_quote 1 "lib.rs:5: data_frame!" "$CLEAN$NL"'fn g() { let r = r#"say "hi" /* "#; }'"$NL"'fn f() { data_frame!(a = 1); }'
# Rust block comments nest: the first `*/` closes only the inner one.
case_ nested_block_comment 0 "glue data frames OK" "$CLEAN$NL"'/* outer /* inner */'"$NL"'   still a comment: data_frame!(a = 1) */'"$NL"'fn f() {}'
case_ after_nested_block_comment 1 "lib.rs:6: data_frame!" "$CLEAN$NL"'/* outer /* inner */'"$NL"'   still a comment */'"$NL"'fn f() { data_frame!(a = 1); }'
case_ call_in_string 0 "glue data frames OK" "$CLEAN$NL"'fn f() { let m = "use data_frame!(..) never"; }'
case_ call_in_byte_string 0 "glue data frames OK" "$CLEAN$NL"'fn f() { let m = b"data_frame!(a)"; }'
case_ call_in_raw_string 0 "glue data frames OK" "$CLEAN$NL"'fn f() { let m = r##"a "# data_frame!(x) "##; }'
case_ call_in_multiline_string 0 "glue data frames OK" "$CLEAN$NL"'fn f() { let m = "first \'"$NL"'    data_frame!(a = 1) \'"$NL"'    last"; }'
case_ unterminated_string 1 "inside a string literal" "$CLEAN$NL"'fn f() { let m = "never closed; }'
case_ unterminated_block 1 "inside a block comment" "$CLEAN$NL"'/* never closed'
case_ unterminated_raw 1 "inside a raw string" "$CLEAN$NL"'fn f() { let m = r#"closed with the wrong count"; }'

# #781 review finding 4: a nested module is glue too.
mkdir -p "$TMP/nested/inner"
printf '%s\n' "$CLEAN" > "$TMP/nested/lib.rs"
printf '%s\n' 'fn f() { data_frame!(a = 1); }' > "$TMP/nested/inner/mod.rs"
cases=$((cases + 1))
out="$(GITHUB_ACTIONS='' bash "$CHECK" "$TMP/nested" 2>&1 || true)"
if ! grep -qE 'inner/mod\.rs:1: data_frame!' <<<"$out"; then
  echo "FAIL: nested module: no finding at inner/mod.rs:1"
  printf '%s\n' "$out" | sed 's/^/     | /'
  failures=$((failures + 1))
fi

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
