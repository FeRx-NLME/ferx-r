#!/usr/bin/env bash
# Drive tools/sibling-cargo-build.sh through every sibling-build outcome, with a
# stub `cargo` that rewrites Cargo.lock the way a real resolve does and then
# skips the 15-minute compile. Nothing in CI has a sibling ferx-core checkout,
# so without this the whole patched-build path - the one that unpins the lock -
# would ship untested.
#
# Each case asserts the exit status, the verdict, and that src/rust/Cargo.lock
# came out byte-identical to what went in. The interesting cases are the ones
# the round-5 review of #349 found: a lock that was already stripped before the
# build, a mixed-revision build, and a shell whose stdout reader goes away (dash
# does not run an EXIT trap on SIGPIPE). The pipe and concurrency cases run
# under every shell found on the box, because that is where the shells differ.
#
# Usage: tools/test-sibling-cargo-build.sh   (from anywhere; exits 1 on failure)

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/sibling-cargo-build.sh"
REAL_LOCK="$ROOT/src/rust/Cargo.lock"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
failures=0

PIN_PREFIX='source = "git+https://github.com/FeRx-NLME/ferx-core'

# The lock version both crates carry in the committed lock; the fixtures below
# match or mismatch it deliberately.
LOCK_VERSION="$(awk '
  $0 == "name = \"ferx-core\"" { in_pkg = 1; next }
  in_pkg && /^version = "/ { gsub(/^version = "|".*$/, ""); print; exit }
' "$REAL_LOCK")"

report() { # headline, captured output
  echo "FAIL $1"
  failures=$((failures + 1))
  printf '%s\n' "$2" | sed 's/^/     | /'
}

# -- fixtures ----------------------------------------------------------------

# A package tree with a stub cargo on PATH and a sibling checkout beside it.
#   make_case <name> <core version> <tools version> [lock variant]
# The lock variant is `pinned` (the committed lock) or `stripped` (what a direct
# cargo run with the old persistent [patch] left behind).
make_case() { # name, core version, tools version, lock variant
  local dir="$TMP/$1" core=$2 tools=$3 variant=${4:-pinned}
  mkdir -p "$dir/pkg/src/rust" "$dir/ferx-core/crates/ferx-tools" "$dir/bin"

  if [[ "$variant" == stripped ]]; then
    awk -v pfx="$PIN_PREFIX" 'index($0, pfx) != 1' "$REAL_LOCK" > "$dir/pkg/src/rust/Cargo.lock"
  else
    cp "$REAL_LOCK" "$dir/pkg/src/rust/Cargo.lock"
  fi
  cp "$dir/pkg/src/rust/Cargo.lock" "$dir/lock.before"

  printf '[package]\nname = "ferx-core"\nversion = "%s"\n' "$core" > "$dir/ferx-core/Cargo.toml"
  printf '[package]\nname = "ferx-tools"\nversion = "%s"\n' "$tools" > "$dir/ferx-core/crates/ferx-tools/Cargo.toml"

  # The stub records how it was called, then rewrites ./Cargo.lock (cargo's cwd
  # is src/rust) the way the resolve it stands in for would have.
  cat > "$dir/bin/cargo" <<'STUB'
#!/bin/sh
echo "$@" > "$STUB_DIR/args"
: > "$STUB_DIR/ran"
[ -n "${STUB_SLEEP:-}" ] && sleep "$STUB_SLEEP"
# Stand in for whatever removes the snapshot mid-build (a tmp reaper, a stray
# rm). Scoped to the case's own TMPDIR so no real build's snapshot is touched.
[ -n "${EAT_SNAPSHOT:-}" ] && rm -f "$EAT_SNAPSHOT"/ferx-Cargo.lock.*
strip_pin() { # crate
  awk -v pfx='source = "git+https://github.com/FeRx-NLME/ferx-core' -v crate="$1" '
    $0 == "name = \"" crate "\"" { hit = 1 }
    hit && index($0, pfx) == 1 { hit = 0; next }
    { print }
  ' Cargo.lock > Cargo.lock.new && mv Cargo.lock.new Cargo.lock
}
mark_unused() { # crate
  printf '\n[[patch.unused]]\nname = "%s"\nversion = "0.0.0"\n' "$1" >> Cargo.lock
}
case "${STUB_MODE:-apply}" in
  apply)  strip_pin ferx-core; strip_pin ferx-tools ;;
  unused) mark_unused ferx-core; mark_unused ferx-tools ;;
  mixed)  strip_pin ferx-core; mark_unused ferx-tools ;;
  noop)   ;;
esac
# Sit on the rewritten lock until released, so a second build can be started
# with a resolve already in flight - the ordering a `sleep` before the rewrite
# cannot produce.
if [ -n "${STUB_HOLD:-}" ]; then
  : > "$STUB_DIR/resolved"
  while [ ! -f "$STUB_HOLD" ]; do sleep 0.2; done
fi
echo "stub cargo: built (mode ${STUB_MODE:-apply})"
exit "${STUB_EXIT:-0}"
STUB
  chmod +x "$dir/bin/cargo"
  echo "$dir"
}

lock_unchanged() { # case dir
  cmp -s "$1/lock.before" "$1/pkg/src/rust/Cargo.lock"
}

# Run the script under one shell, from the case's src/ directory.
#   run_case <case dir> <shell> [extra env assignments...]
run_case() {
  local dir=$1 shell=$2
  shift 2
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" FERX_CORE_SIBLING="$dir/ferx-core" "$@" \
      "$shell" "$SCRIPT" build --release 2>&1 )
}

# -- assertions --------------------------------------------------------------

# check <name> <case dir> <expected exit> <output> <actual exit> <fragment>...
# Always also asserts the lock is byte-identical to the one the case started
# with: that is the acceptance criterion, in every outcome.
check() {
  local name=$1 dir=$2 want=$3 out=$4 rc=$5
  shift 5
  if [[ "$rc" -ne "$want" ]]; then
    report "$name: exit $rc, expected $want" "$out"
    return
  fi
  if ! lock_unchanged "$dir"; then
    report "$name: Cargo.lock was left rewritten" "$(diff "$dir/lock.before" "$dir/pkg/src/rust/Cargo.lock" | head -20)"
    return
  fi
  local fragment
  for fragment in "$@"; do
    if [[ "$fragment" == !* ]]; then
      if [[ "$out" == *"${fragment#!}"* ]]; then
        report "$name: output must not contain \"${fragment#!}\"" "$out"
        return
      fi
    elif [[ "$out" != *"$fragment"* ]]; then
      report "$name: output lacks \"$fragment\"" "$out"
      return
    fi
  done
  echo "ok   $name"
}

# -- the cases ---------------------------------------------------------------

# 1. The patch applies to both crates: the happy path.
dir=$(make_case applied "$LOCK_VERSION" "$LOCK_VERSION")
out=$(STUB_MODE=apply run_case "$dir" /bin/sh); rc=$?
check "patch applies to both" "$dir" 0 "$out" "$rc" \
  "supplied BOTH ferx-core and ferx-tools" "Cargo.lock restored to its pin"
# and the patch really was passed per crate, to that one cargo run
args=$(cat "$dir/args")
for want in 'patch."https://github.com/FeRx-NLME/ferx-core".ferx-core.path' \
            'patch."https://github.com/FeRx-NLME/ferx-core".ferx-tools.path'; do
  if [[ "$args" != *"$want"* ]]; then
    report "patch applies to both: cargo was not given $want" "$args"
  else
    echo "ok   cargo received --config for ${want##*core\".}"
  fi
done

# 2. Both versions differ from the lock, so cargo ignores the patch entirely.
dir=$(make_case unused 9.9.9 9.9.9)
out=$(STUB_MODE=unused run_case "$dir" /bin/sh); rc=$?
check "patch unused for both" "$dir" 0 "$out" "$rc" \
  "the sibling was NOT used" "patch entries unused: ferx-core ferx-tools" \
  "sibling 9.9.9" "!supplied BOTH"

# 3. A mixed-revision build that cargo's own resolve produced. The version
#    pre-check cannot see it (both manifests match the lock), so only the
#    per-crate verdict read back from the lock catches it.
dir=$(make_case mixed "$LOCK_VERSION" "$LOCK_VERSION")
out=$(STUB_MODE=mixed run_case "$dir" /bin/sh); rc=$?
check "mixed revisions refused" "$dir" 1 "$out" "$rc" \
  "MIXED revisions" "ferx-core came from the sibling" \
  "ferx-tools from the pinned revision" "!the sibling was NOT used"

# 4. A mixed build the versions predict: refused before cargo runs at all.
dir=$(make_case mixed-preflight "$LOCK_VERSION" 9.9.9)
out=$(run_case "$dir" /bin/sh); rc=$?
check "mixed versions refused before building" "$dir" 1 "$out" "$rc" \
  "would supply only one of the two crates" "Refusing to build"
if [[ -f "$dir/ran" ]]; then
  report "mixed versions refused before building: cargo ran anyway" "$out"
else
  echo "ok   mixed versions refused before cargo ran"
fi

# 5. A lock some earlier direct cargo run already stripped: reported as found,
#    never announced as restored to a pin.
dir=$(make_case pre-stripped "$LOCK_VERSION" "$LOCK_VERSION" stripped)
out=$(STUB_MODE=apply run_case "$dir" /bin/sh); rc=$?
check "pre-stripped lock reported" "$dir" 0 "$out" "$rc" \
  "ALREADY unpinned before this build" "!restored to its pin"

# 6. cargo fails: its status is passed on, and the lock still goes back.
dir=$(make_case cargo-fails "$LOCK_VERSION" "$LOCK_VERSION")
out=$(STUB_MODE=apply STUB_EXIT=101 run_case "$dir" /bin/sh); rc=$?
check "cargo failure propagated" "$dir" 101 "$out" "$rc" "cargo exited 101"

# 7. The snapshot disappears mid-build. The lock cannot be put back, so the
#    only acceptable behaviour is to fail and say which file it is about -
#    never to exit 0 on a lock left as cargo wrote it.
dir=$(make_case snapshot-eaten "$LOCK_VERSION" "$LOCK_VERSION")
mkdir -p "$dir/tmp"
out=$( ( cd "$dir/pkg/src" &&
  env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply TMPDIR="$dir/tmp" \
    EAT_SNAPSHOT="$dir/tmp" FERX_CORE_SIBLING="$dir/ferx-core" \
    /bin/sh "$SCRIPT" build 2>&1 ) ) && rc=0 || rc=$?
if [[ "$rc" -ne 1 ]]; then
  report "a vanished snapshot fails loudly: exit $rc, expected 1" "$out"
elif [[ "$out" != *"vanished before it could be put back"* ]]; then
  report "a vanished snapshot fails loudly: without saying so" "$out"
elif lock_unchanged "$dir"; then
  report "a vanished snapshot fails loudly: the lock is intact, so the case tested nothing" "$out"
else
  echo "ok   a vanished snapshot fails loudly"
fi

# 7b. A guard left behind by a build that was killed outright. The checkout must
#     not be blocked forever by a directory whose owner is gone.
dir=$(make_case stale-guard "$LOCK_VERSION" "$LOCK_VERSION")
mkdir -p "$dir/pkg/src/rust/target/.ferx-lock-guard"
# A pid that cannot be running: allocated, then reaped.
sh -c 'exit 0' & dead_pid=$!; wait "$dead_pid" 2>/dev/null
echo "$dead_pid" > "$dir/pkg/src/rust/target/.ferx-lock-guard/pid"
out=$(STUB_MODE=apply run_case "$dir" /bin/sh); rc=$?
check "a guard from a dead build is taken over" "$dir" 0 "$out" "$rc" \
  "taking over the build guard left behind by process $dead_pid" \
  "Cargo.lock restored to its pin"
if [[ -d "$dir/pkg/src/rust/target/.ferx-lock-guard" ]]; then
  report "a guard from a dead build is taken over: the guard was not released" "$out"
else
  echo "ok   the guard is released on the way out"
fi

# 8. A sibling predating the ferx-tools workspace split.
dir=$(make_case no-ferx-tools "$LOCK_VERSION" "$LOCK_VERSION")
rm -rf "$dir/ferx-core/crates"
out=$(run_case "$dir" /bin/sh); rc=$?
check "sibling without crates/ferx-tools" "$dir" 2 "$out" "$rc" "has no crates/ferx-tools"

# 9. A POSIX path with a backslash in it: a legal directory name, but not
#    something a TOML string can carry, and not something to rewrite either -
#    the build would then point at a different directory.
dir=$(make_case backslash-path "$LOCK_VERSION" "$LOCK_VERSION")
if mv "$dir/ferx-core" "$dir/ferx\\core" 2>/dev/null; then
  out=$( ( cd "$dir/pkg/src" && env PATH="$dir/bin:$PATH" STUB_DIR="$dir" \
      FERX_CORE_SIBLING="$dir/ferx\\core" /bin/sh "$SCRIPT" build 2>&1 ) ) && rc=0 || rc=$?
  check "backslash in the sibling path refused" "$dir" 2 "$out" "$rc" \
    "contains a backslash"
else
  echo "skip backslash in the sibling path refused (this filesystem will not take the name)"
fi

# 9b. A Windows path, though, is translated: cargo reads forward slashes there.
#     The directory does not exist here, and the error naming it forward-slashed
#     is what shows the translation happened before anything else looked at it.
dir=$(make_case windows-path "$LOCK_VERSION" "$LOCK_VERSION")
out=$( ( cd "$dir/pkg/src" && env PATH="$dir/bin:$PATH" STUB_DIR="$dir" \
    FERX_CORE_SIBLING='C:\Users\Someone\ferx-core' /bin/sh "$SCRIPT" build 2>&1 ) ) && rc=0 || rc=$?
check "a drive-letter path is translated, not refused" "$dir" 2 "$out" "$rc" \
  "no sibling ferx-core checkout at C:/Users/Someone/ferx-core"

# 10. The shell-dependent cases: a reader that goes away mid-build, and two
#    builds racing in one checkout. dash is /bin/sh on Debian and Ubuntu and
#    kills itself on SIGPIPE instead of running the EXIT trap, which is exactly
#    what the script's `trap '' PIPE` is for.
shells=(/bin/sh)
seen=" $(readlink -f /bin/sh 2>/dev/null || echo /bin/sh) "
for candidate in /bin/dash /bin/bash /usr/local/bin/dash /opt/homebrew/bin/dash; do
  [[ -x "$candidate" ]] || continue
  # On Debian and Ubuntu /bin/sh *is* dash; running it twice buys nothing.
  resolved="$(readlink -f "$candidate" 2>/dev/null || echo "$candidate")"
  [[ "$seen" == *" $resolved "* ]] && continue
  seen+="$resolved "
  shells+=("$candidate")
done

for shell in "${shells[@]}"; do
  label="${shell##*/}"

  # The output reader exits while the build is still writing.
  dir=$(make_case "pipe-$label" "$LOCK_VERSION" "$LOCK_VERSION")
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply STUB_SLEEP=2 \
      FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build 2>&1 ) | (sleep 0.2; exit 0)
  # Give the orphaned build time to finish writing and restoring.
  for _ in 1 2 3 4 5 6 7 8 9 10; do
    [[ -f "$dir/ran" ]] && lock_unchanged "$dir" && break
    sleep 0.5
  done
  if lock_unchanged "$dir"; then
    echo "ok   closed pipe restores the lock ($label)"
  else
    report "closed pipe restores the lock ($label)" "Cargo.lock left rewritten"
  fi

  # The reader is already gone when the script writes its first line. Without
  # `trap '' PIPE` the shell is killed by SIGPIPE right there and cargo never
  # runs at all; with it the write fails harmlessly and the build goes ahead.
  # The reader exits at once and the writer starts half a second later, so the
  # order does not depend on scheduling.
  dir=$(make_case "dead-reader-$label" "$LOCK_VERSION" "$LOCK_VERSION")
  ( sleep 0.5
    cd "$dir/pkg/src" &&
      env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply \
        FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build 2>&1 ) | (exit 0)
  for _ in 1 2 3 4 5 6 7 8 9 10; do
    [[ -f "$dir/ran" ]] && lock_unchanged "$dir" && break
    sleep 0.5
  done
  if [[ ! -f "$dir/ran" ]]; then
    report "a reader gone before the first line does not kill the build ($label)" \
      "cargo never ran: SIGPIPE killed the shell"
  elif ! lock_unchanged "$dir"; then
    report "a reader gone before the first line does not kill the build ($label)" \
      "Cargo.lock left rewritten"
  else
    echo "ok   a reader gone before the first line does not kill the build ($label)"
  fi

  # Two builds started moments apart must not share, or delete, one snapshot.
  dir=$(make_case "race-$label" "$LOCK_VERSION" "$LOCK_VERSION")
  out_a="$dir/out-a"; out_b="$dir/out-b"
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply STUB_SLEEP=2 \
      FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build >"$out_a" 2>&1 ) &
  pid_a=$!
  sleep 0.3
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply STUB_SLEEP=1 \
      FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build >"$out_b" 2>&1 ) &
  pid_b=$!
  wait "$pid_a"; wait "$pid_b"
  both=$(cat "$out_a" "$out_b")
  if ! lock_unchanged "$dir"; then
    report "two concurrent builds restore the lock ($label)" "$both"
  elif [[ "$both" == *"No such file"* ]]; then
    report "two concurrent builds restore the lock ($label): a snapshot went missing" "$both"
  else
    echo "ok   two concurrent builds restore the lock ($label)"
  fi

  # The ordering that breaks a snapshot on its own: A resolves first and its
  # rewritten lock is on disk *before* B starts. Unguarded, B snapshots A's
  # stripped lock, A restores the pin, and B then puts the strip back - both
  # exiting 0 on an unpinned checkout. B must instead wait for A.
  dir=$(make_case "ordered-race-$label" "$LOCK_VERSION" "$LOCK_VERSION")
  out_a="$dir/out-a"; out_b="$dir/out-b"
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply STUB_HOLD="$dir/release" \
      FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build >"$out_a" 2>&1 ) &
  pid_a=$!
  # A has rewritten the lock and is holding it there.
  for _ in $(seq 1 100); do
    [[ -f "$dir/resolved" ]] && break
    sleep 0.2
  done
  if [[ ! -f "$dir/resolved" ]] || lock_unchanged "$dir"; then
    report "a build started mid-resolve waits its turn ($label)" \
      "the stub never got the lock into its rewritten state, so the case tested nothing"
    kill "$pid_a" 2>/dev/null
    wait "$pid_a" 2>/dev/null
    continue
  fi
  # B's stub is slow, so B is the one that finishes last and gets the last word
  # on the lock - the order in which an unguarded B puts A's strip back for good.
  ( cd "$dir/pkg/src" &&
    env PATH="$dir/bin:$PATH" STUB_DIR="$dir" STUB_MODE=apply STUB_SLEEP=4 \
      FERX_CORE_SIBLING="$dir/ferx-core" "$shell" "$SCRIPT" build >"$out_b" 2>&1 ) &
  pid_b=$!
  sleep 1                 # long enough for B to reach the point it snapshots at
  : > "$dir/release"      # let A finish and restore, while B is still going
  wait "$pid_a"; rc_a=$?
  # B blocks until A is done; a wrapper that never releases must not hang CI.
  b_done=no
  for _ in $(seq 1 120); do
    kill -0 "$pid_b" 2>/dev/null || { b_done=yes; break; }
    sleep 0.5
  done
  if [[ "$b_done" == no ]]; then
    kill -9 "$pid_b" 2>/dev/null
    report "a build started mid-resolve waits its turn ($label)" "B never finished (60s)"
  fi
  wait "$pid_b"; rc_b=$?
  both=$(cat "$out_a" "$out_b")
  if [[ "$b_done" == no ]]; then
    : # already reported
  elif ! lock_unchanged "$dir"; then
    report "a build started mid-resolve waits its turn ($label): Cargo.lock left rewritten" "$both"
  elif [[ "$rc_a" -ne 0 || "$rc_b" -ne 0 ]]; then
    report "a build started mid-resolve waits its turn ($label): exits $rc_a / $rc_b" "$both"
  elif ! grep -q "waiting for another sibling build" "$out_b"; then
    report "a build started mid-resolve waits its turn ($label): B never said it waited" "$both"
  else
    echo "ok   a build started mid-resolve waits its turn ($label)"
  fi
done

if [[ "$failures" -ne 0 ]]; then
  echo "$failures sibling-cargo-build.sh case(s) failed"
  exit 1
fi
echo "all sibling-cargo-build.sh cases passed"
