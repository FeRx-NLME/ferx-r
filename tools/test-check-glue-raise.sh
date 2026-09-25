#!/usr/bin/env bash
# Drive tools/check-glue-raise.sh over stub glue sources, one per shape it is
# supposed to catch and one per shape it must not.
#
# The real src/rust/src/lib.rs is always in shape, so nothing there would fail
# if one of the checks were deleted - and every check exists only because the
# compiler cannot see the shape it holds (a fully qualified
# `extendr_api::throw_r_error(..)` compiles; a second `raise_verbatim(` call
# site compiles; a `#[extendr]` body that never returns through `entry()`
# compiles; a second extern "C" block re-declaring Rf_error compiles; a
# Cargo.lock moving extendr past the release raise_verbatim escapes for
# resolves). The stubs are the only thing that reddens when a check goes away,
# so each of the five has a fixture of its own - including "nothing raises at all", whose
# absence let that check be deleted with the harness still green.
#
# Every bad fixture must exit 1 with its own message, so a check firing for the
# wrong reason does not pass. Each case runs twice - GITHUB_ACTIONS unset and
# set - because CI takes the `::error` branch. Output is captured and printed
# only for a failing case, indented, so no `::error` line becomes an annotation.
#
# Usage: tools/test-check-glue-raise.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="$ROOT/tools/check-glue-raise.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
failures=0

report() { # headline, captured output
  echo "$1"
  failures=$((failures + 1))
  printf '%s\n' "$2" | sed 's/^/     | /'
}

# The shape every fixture starts from: the shadow, the single raise_verbatim()
# call inside entry(), and one wrapped entry point.
preamble() {
  cat <<'RS'
use extendr_api::prelude::*;

mod raise {
    extern "C" {
        fn Rf_error(fmt: *const std::ffi::c_char, ...) -> !;
    }

    pub(super) fn raise_verbatim(msg: String) -> ! {
        std::panic::resume_unwind(Box::new(msg.replace('%', "%%")))
    }
}

use raise::raise_verbatim;

fn entry<T>(f: impl FnOnce() -> Result<T, String>) -> T {
    let msg = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(Ok(value)) => return value,
        Ok(Err(msg)) => msg,
        Err(payload) => panic_message(payload, std::panic::Location::caller()),
    };
    raise_verbatim(msg)
}

fn throw_r_error(never: std::convert::Infallible) -> ! {
    match never {}
}

#[extendr]
fn ferx_rust_ok(path: &str) -> Robj {
    entry(move || {
        if path.is_empty() {
            return Err("empty".to_string());
        }
        Ok(NULL.into())
    })
}
RS
}

fixture() { # name, extra source appended to the preamble
  local dir="$TMP/$1"
  mkdir -p "$dir"
  { preamble; printf '%s\n' "$2"; } > "$dir/lib.rs"
  echo "$dir"
}

expect() { # name, dir, expected exit, message fragment[, Cargo.lock]
  local name=$1 dir=$2 want=$3 fragment=$4 lock=${5:-$ROOT/src/rust/Cargo.lock}
  local mode out rc
  for mode in local ci; do
    rc=0
    if [[ "$mode" == ci ]]; then
      out=$(GITHUB_ACTIONS=true bash "$CHECK" "$dir" "$lock" 2>&1) || rc=$?
    else
      out=$(env -u GITHUB_ACTIONS bash "$CHECK" "$dir" "$lock" 2>&1) || rc=$?
    fi
    if [[ "$rc" != "$want" ]]; then
      report "FAIL [$mode] $name: expected exit $want, got $rc" "$out"
      continue
    fi
    if [[ -n "$fragment" ]] && ! printf '%s' "$out" | grep -qF -- "$fragment"; then
      report "FAIL [$mode] $name: output does not mention \"$fragment\"" "$out"
      continue
    fi
    echo "ok   [$mode] $name"
  done
}

# -- green --------------------------------------------------------------------

expect "a glue in shape passes" \
  "$(fixture in-shape '')" 0 "glue raise shape OK"

expect "naming throw_r_error in a comment is free" \
  "$(fixture commented '
// Never call throw_r_error(msg) from a body: see #388.
/// raise_verbatim(msg) belongs to entry() alone.
#[extendr]
fn ferx_rust_second(x: &str) -> Robj {
    entry(move || Ok(x.into()))
}')" 0 "glue raise shape OK"

# ferx_rust_fit is exactly this shape: an #[allow] between the attribute and a
# signature that spans many lines.
expect "an #[allow] and a multi-line signature still count as wrapped" \
  "$(fixture multiline '
#[extendr]
#[allow(clippy::too_many_arguments)]
fn ferx_rust_wide(
    a: &str,
    b: &str,
) -> Robj {
    entry(move || {
        Ok([a, b].concat().into())
    })
}')" 0 "glue raise shape OK"

# `#[extendr(r_name = ..)]`, `#[extendr(use_try_from = true)]`: the attribute
# takes arguments, and the check used to compare it for exact equality.
expect "an #[extendr] with arguments counts as an entry point" \
  "$(fixture attr-args '
#[extendr(r_name = "ferx_rust_renamed")]
fn ferx_rust_named(x: &str) -> Robj {
    entry(move || Ok(x.into()))
}')" 0 "glue raise shape OK"

# No #[extendr] body in lib.rs is a one-liner today, but the check used to wait
# for a line that ends in `{`, so any body written on its signature line was
# swallowed whole.
expect "a wrapped one-line body passes" \
  "$(fixture one-line-ok '
#[extendr]
fn ferx_rust_blocks() -> Robj { entry(move || Ok(NULL.into())) }')" 0 "glue raise shape OK"

# `#[extendr(` may wrap. The attribute used to have to close on its own line
# for the signature below it to be seen at all.
expect "a wrapped #[extendr(..) attribute still finds the body" \
  "$(fixture attr-wrapped '
#[extendr(
    use_try_from = true
)]
fn ferx_rust_wrapped_attr(x: &str) -> Robj {
    entry(move || Ok(x.into()))
}')" 0 "glue raise shape OK"

# The `//` strip knows nothing of /* .. */, so a body opening after a block
# comment used to be read as the body itself.
expect "a block comment before the body is not the body" \
  "$(fixture block-comment '
#[extendr]
fn ferx_rust_blocky(x: &str) -> Robj {
    /* the engine refuses an empty path, and says so itself */
    entry(move || Ok(x.into()))
}')" 0 "glue raise shape OK"

# -- red ----------------------------------------------------------------------

# #388 rebuilt by hand: `mod raise` keeps Rf_error out of scope, but a second
# extern block re-declares the symbol at file scope and compiles clean.
expect "a second extern block re-declaring Rf_error is refused" \
  "$(fixture redeclared '
extern "C" {
    fn Rf_error(fmt: *const std::ffi::c_char, ...) -> !;
}

#[extendr]
fn ferx_rust_redeclares(x: &str) -> Robj {
    entry(move || {
        let c = std::ffi::CString::new(x).unwrap();
        unsafe { Rf_error(c.as_ptr()) }
    })
}')" 1 "Rf_error is named outside mod raise"

# #[extendr] impl: the guard reads the first method and nothing else, so it
# says so rather than passing the block.
expect "#[extendr] on an impl block is refused" \
  "$(fixture extendr-impl '
#[extendr]
impl Fit {
    fn ofv(&self) -> f64 {
        entry(move || Ok(0.0))
    }
    fn theta(&self) -> f64 {
        0.0
    }
}')" 1 "on an impl block"


expect "a bare throw_r_error() call is refused" \
  "$(fixture bare '
#[extendr]
fn ferx_rust_bare(x: &str) -> Robj {
    entry(move || {
        throw_r_error(format!("no: {x}"));
    })
}')" 1 "throw_r_error() is reachable here"

expect "a fully qualified extendr_api::throw_r_error() is refused" \
  "$(fixture qualified '
#[extendr]
fn ferx_rust_qualified(x: &str) -> Robj {
    entry(move || {
        extendr_api::throw_r_error(format!("no: {x}"));
    })
}')" 1 "throw_r_error() is reachable here"

expect "a second raise_verbatim() call site is refused" \
  "$(fixture second-raise '
#[extendr]
fn ferx_rust_raises(x: &str) -> Robj {
    entry(move || {
        raise_verbatim(format!("no: {x}"));
    })
}')" 1 "called from more than one place"

expect "an #[extendr] with arguments whose body is unwrapped is refused" \
  "$(fixture attr-args-unwrapped '
#[extendr(use_try_from = true)]
fn ferx_rust_try_from(x: &str) -> Robj {
    let y = x.to_string();
    Ok(y.into()).unwrap()
}')" 1 "does not open its body with entry()"

expect "a wrapped #[extendr(..) attribute over an unwrapped body is refused" \
  "$(fixture attr-wrapped-unwrapped '
#[extendr(
    use_try_from = true
)]
fn ferx_rust_wrapped_bare(x: &str) -> Robj {
    let y = x.to_string();
    Ok(y.into()).unwrap()
}')" 1 "does not open its body with entry()"

expect "an unwrapped one-line body is refused" \
  "$(fixture one-line-unwrapped '
#[extendr]
fn ferx_rust_inline() -> Robj { NULL.into() }')" 1 "does not open its body with entry()"

expect "an #[extendr] body that does not open with entry() is refused" \
  "$(fixture unwrapped '
#[extendr]
fn ferx_rust_unwrapped(x: &str) -> Robj {
    let y = x.to_string();
    Ok(y.into()).unwrap()
}')" 1 "does not open its body with entry()"

# The shadow is what makes a bare call a compile error, so its deletion has to
# be a finding of its own - otherwise the red case above silently stops firing.
shadowless="$TMP/shadowless"
mkdir -p "$shadowless"
preamble | grep -v 'fn throw_r_error' | grep -v 'match never' > "$shadowless/lib.rs"
expect "deleting the shadow is refused" "$shadowless" 1 "shadow is gone"

# entry() is the one raise site, so a glue with no raise site at all is a glue
# whose errors leave by some other door. Without this fixture the `raises -eq 0`
# check could be deleted with the harness still printing "all cases behaved".
raiseless="$TMP/raiseless"
mkdir -p "$raiseless"
preamble | grep -v '^    raise_verbatim(msg)$' > "$raiseless/lib.rs"
expect "a glue that calls raise_verbatim() nowhere is refused" "$raiseless" 1 \
  "nothing calls raise_verbatim()"

# `mod raise` deleted outright: Rf_error back at file scope, which is where it
# was before the module and is what check 4 is ultimately holding.
moduleless="$TMP/moduleless"
mkdir -p "$moduleless"
preamble | grep -v '^mod raise {$' > "$moduleless/lib.rs"
expect "deleting mod raise is refused" "$moduleless" 1 \
  "Rf_error is named outside mod raise"

# Check 5: the % escape in raise_verbatim is right for extendr-api 0.9.0's
# throw_r_error only (#394), so a lock that moves extendr is a finding.
lock_fixture() { # name, extendr-api version, extendr-macros version ("" = absent)
  local f="$TMP/$1.lock"
  {
    echo 'version = 4'
    echo
    if [[ -n "$2" ]]; then
      printf '[[package]]\nname = "extendr-api"\nversion = "%s"\n\n' "$2"
    fi
    if [[ -n "$3" ]]; then
      printf '[[package]]\nname = "extendr-macros"\nversion = "%s"\n\n' "$3"
    fi
  } > "$f"
  echo "$f"
}
in_shape="$(fixture lock-shape '')"

expect "a lock at the extendr the escape was written for passes" \
  "$in_shape" 0 "extendr at 0.9.0" "$(lock_fixture lock-ok 0.9.0 0.9.0)"

expect "a lock that moves extendr-api is refused" \
  "$in_shape" 1 "extendr-api is 0.9.1 in the lock" \
  "$(lock_fixture lock-api 0.9.1 0.9.0)"

expect "a lock that moves extendr-macros alone is refused" \
  "$in_shape" 1 "extendr-macros is 0.10.0 in the lock" \
  "$(lock_fixture lock-macros 0.9.0 0.10.0)"

expect "a lock without extendr is refused" \
  "$in_shape" 1 "extendr-api is not in the lock" \
  "$(lock_fixture lock-none '' 0.9.0)"

expect "a missing lock is refused" \
  "$in_shape" 1 "no Cargo.lock at" "$TMP/no-such.lock"

empty="$TMP/empty"
mkdir -p "$empty"
expect "a directory with no sources is refused" "$empty" 1 "no .rs sources"

if [[ $failures -gt 0 ]]; then
  echo "$failures check(s) failed" >&2
  exit 1
fi
echo "all glue-raise guard cases behaved"
