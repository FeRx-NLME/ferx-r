#!/bin/sh
# ad-preflight.sh — verify the `enzyme` toolchain actually differentiates
# before committing to the ~15-30 min fat-LTO build.
#
# A toolchain assembled by copying a standalone Enzyme plugin into a prebuilt
# nightly's sysroot PASSES `rustup toolchain list | grep enzyme` but then hangs
# forever in llvm::Constant::getNullValue the first time it differentiates
# anything. The trivial `E0601` check (`rustc - </dev/null`) never runs
# differentiation, so it cannot see this. The only reliable test is to compile
# and run a real #[autodiff_forward] function — under a timeout, because the
# failure mode is an infinite hang rather than an error.
#
# Exit codes:
#   0  toolchain verified good (or cached good for this toolchain)
#   1  toolchain present but broken (hang/timeout or wrong result)
#   2  could not run the check (e.g. `rustc +enzyme` unavailable)
#
# Environment:
#   FERX_AD_PREFLIGHT_TIMEOUT  seconds allowed for each probe step (default 90)
#   FERX_AD_PREFLIGHT_SKIP=1   skip the check entirely (already verified)

set -u

TOOLCHAIN=enzyme
TIMEOUT="${FERX_AD_PREFLIGHT_TIMEOUT:-90}"

[ "${FERX_AD_PREFLIGHT_SKIP:-0}" = "1" ] && exit 0

# Resolve the toolchain the build will actually use and fingerprint it (commit
# hash + LLVM line) so the probe re-runs when the toolchain changes but is
# skipped on identical reinstalls.
RUSTC_VV="$(rustc "+$TOOLCHAIN" -vV 2>/dev/null)" || {
  echo "ferx: AD preflight could not invoke 'rustc +$TOOLCHAIN'." >&2
  exit 2
}
FINGERPRINT="$(printf '%s' "$RUSTC_VV" | grep -E 'commit-hash|LLVM' \
                | tr -d ' ' | tr '\n' '-')"
CACHE_DIR="${TMPDIR:-/tmp}/ferx-ad-preflight"
CACHE_FILE="$CACHE_DIR/${FINGERPRINT:-unknown}.ok"
[ -f "$CACHE_FILE" ] && exit 0

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ferx-adprobe.XXXXXX")" || exit 2
trap 'rm -rf "$WORK"' EXIT

# The same real differentiation the install docs verify with: d/dx of
# f(x,y) = exp(x*y) + x*x is y*exp(x*y) + 2x, so df(1.5, 1, 2) = 2*e^3 + 3.
cat > "$WORK/ad_check.rs" <<'RS'
#![feature(autodiff)]
use std::autodiff::autodiff_forward;
#[autodiff_forward(df, Dual, Const, Dual)]
#[inline(never)]
fn f(x: f64, y: f64) -> f64 { (x * y).exp() + x * x }
fn main() {
    let (val, dfdx) = df(1.5, 1.0, 2.0);
    let expect = 2.0 * (1.5 * 2.0_f64).exp() + 2.0 * 1.5;
    assert!((dfdx - expect).abs() < 1e-9);
    println!("autodiff OK: f={val:.4}, df/dx={dfdx:.4}");
}
RS

# run_guarded CMD... : run CMD in the background under a sleep/kill watchdog so
# an infinite hang becomes a clean failure. Portable — does not depend on
# `timeout`/`gtimeout` (absent on stock macOS). Returns CMD's exit status, or
# 137 (128 + SIGKILL) when the watchdog fired. rustc's Enzyme hang is in-process
# during the codegen pass, so killing its PID is sufficient — no process-group
# handling (which would need non-portable setsid/`kill -- -pgid`).
run_guarded() {
  "$@" &
  cmd_pid=$!
  ( sleep "$TIMEOUT"; kill -9 "$cmd_pid" 2>/dev/null ) &
  wd_pid=$!
  wait "$cmd_pid" 2>/dev/null
  status=$?
  kill "$wd_pid" 2>/dev/null
  wait "$wd_pid" 2>/dev/null
  return "$status"
}

# 1) Compile — this is where the LLVM/Enzyme mismatch hangs.
run_guarded rustc "+$TOOLCHAIN" -Zautodiff=Enable -Clto=fat -Copt-level=2 \
  "$WORK/ad_check.rs" -o "$WORK/ad_check" >"$WORK/rustc.log" 2>&1
status=$?
if [ "$status" -ne 0 ]; then
  if [ "$status" -eq 137 ]; then
    echo "ferx: Enzyme AD preflight TIMED OUT after ${TIMEOUT}s." >&2
    echo "      The '$TOOLCHAIN' toolchain compiles but HANGS the first time it" >&2
    echo "      differentiates — the standalone-plugin / CI-LLVM mismatch." >&2
    echo "      Rebuild the toolchain from source with --set llvm.download-ci-llvm=false:" >&2
    echo "      https://ferx-nlme.github.io/ferx-core/installation.html#do-not-build-a-standalone-enzyme-plugin-against-a-separately-installed-llvm" >&2
  else
    echo "ferx: Enzyme AD preflight failed to compile (rustc exit $status):" >&2
    sed 's/^/      /' "$WORK/rustc.log" >&2
  fi
  exit 1
fi

# 2) Run it — a subtler mismatch can hang at run time rather than compile time.
run_guarded "$WORK/ad_check" >"$WORK/run.log" 2>&1
status=$?
if [ "$status" -ne 0 ] || ! grep -q 'autodiff OK' "$WORK/run.log"; then
  echo "ferx: Enzyme AD preflight produced no valid result (exit $status)." >&2
  sed 's/^/      /' "$WORK/run.log" >&2
  exit 1
fi

mkdir -p "$CACHE_DIR" && : > "$CACHE_FILE"
exit 0
