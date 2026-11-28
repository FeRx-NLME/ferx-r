#!/usr/bin/env bash
# Drive tools/check-glue-id-columns.sh over stub glue sources, one per shape it
# must catch and one per shape it must not. The real lib.rs is clean, so it
# would stay green if the check's pattern were broken; the stubs are what
# redden. Each case runs with GITHUB_ACTIONS unset and set (CI takes the
# `::error` branch). Output is printed only for a failing case, indented, so no
# `::error` line becomes an annotation.
#
# Usage: tools/test-check-glue-id-columns.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="$ROOT/tools/check-glue-id-columns.sh"
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

CLEAN='fn subject_id_column<I>(ids: I) -> Robj { todo!() }
fn table(ids: Vec<String>, time: Vec<f64>) -> Robj {
    r_data_frame(vec![("ID", subject_id_column(ids)), ("TIME", time.into())], 1)
}'

case_ clean 0 "glue ID columns OK" "$CLEAN"
case_ helper_on_next_line 0 "" "$CLEAN
fn f() { pairs.push((
    \"ID\",
    subject_id_column(rows.iter().map(|r| &r.id)),
)); }"
case_ helper_after_blank_line 0 "" "$CLEAN
fn f() { let c = (\"ID\",

    subject_id_column(ids)); }"
case_ prose_line_comment 0 "" "$CLEAN
// (\"ID\", ids.into()) gave a number column; sdtab_id(id, 0) was the fallback.
/// Unlike \`id.parse::<f64>()\`, this keeps the text."
case_ prose_block_comment 0 "" "$CLEAN
/* (\"ID\", ids.into())
   sdtab_id(&s.id, si) */
fn other() {}"
case_ in_a_string 0 "" "$CLEAN
fn f() -> &'static str { \"(\\\"ID\\\", ids.into()) and sdtab_id(x)\" }"
case_ other_columns 0 "" "$CLEAN
fn f() { cols.push((\"IDX\", ids.into())); cols.push((\"SUBJECT_ID\", n.into())); }"
case_ similar_names 0 "" "$CLEAN
fn f() { let sdtab_ids = 1; my_sdtab_id(2); let n = s.idx.parse::<f64>(); }"
case_ id_match_arm 0 "" "$CLEAN
fn f(name: &str) -> f64 { match name { \"ID\" => f64::NAN, _ => 0.0 } }"

case_ planted_tuple 1 "lib.rs:5: an ID column built without subject_id_column()" "$CLEAN
fn f() { cols.push((\"ID\", ids.into())); }"
case_ planted_tuple_vec 1 "an ID column built without subject_id_column()" "$CLEAN
fn f() { r_data_frame(vec![(\"TIME\", t.into()), ( \"ID\" , id.into())], n) }"
case_ planted_value_next_line 1 "lib.rs:6: an ID column built without subject_id_column()" "$CLEAN
fn f() { pairs.push((
    \"ID\",
    ids.into(),
)); }"
case_ planted_collect 1 "an ID column built without subject_id_column()" "$CLEAN
fn f() { (\"ID\", l.iter().map(|e| e.subject.clone()).collect::<Vec<String>>().into()) }"
case_ planted_at_end_of_file 1 "an ID column built without subject_id_column()" "$CLEAN
fn f() { (\"ID\","
case_ helper_later_on_the_line 1 "an ID column built without subject_id_column()" "$CLEAN
fn f() { (\"ID\", ids.into()); subject_id_column(ids); }"
case_ planted_sdtab_id 1 "lib.rs:5: sdtab_id() numbered a subject" "$CLEAN
fn f() { id.push(sdtab_id(&subject.id, si)); }"
case_ planted_sdtab_id_definition 1 "sdtab_id() numbered a subject" "$CLEAN
fn sdtab_id(id: &str, position: usize) -> f64 { 0.0 }"
case_ planted_parse 1 "lib.rs:5: a subject ID parsed as a number" "$CLEAN
fn f() { let n = subj.id.parse::<f64>().unwrap_or(1.0); }"
case_ lost_in_comment 1 "lost track" "$CLEAN
/* never closed"

if [[ $failures -gt 0 ]]; then
  echo "$failures of $cases glue ID-column guard case(s) failed." >&2
  exit 1
fi
echo "glue ID-column guard: all $cases cases pass."
