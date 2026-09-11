#!/usr/bin/env bash
# Bump the ferx-core / ferx-tools pin in src/rust/Cargo.lock to current main HEAD.
#
# Both crates live in the one ferx-core repository and are patched to the local
# sibling checkout together, so they are bumped together and must end up on the
# same revision.
#
# Why this script exists: src/rust/.cargo/config.toml carries a [patch] that
# redirects ferx-core to a sibling ../ferx-core checkout when present. If you
# run any `cargo` command that resolves while that patch is applied, cargo
# writes a *path*-style lock entry with no `source = "git+..."` line — which
# silently unpins ferx-core for everyone who builds without the sibling (CI,
# downstream users). That happened twice, both times in a lock bump
# (commits 1ce7f59, b96c867).
#
# This script temporarily removes the patch so cargo resolves ferx-core from
# GitHub and writes the correct git+https pin. It runs a whole-graph
# `cargo update`, which also moves every registry crate to its latest
# compatible version. With the patch absent, `cargo update -p ferx-core
# --precise <sha>` also works (measured 2026-09-11): it moves ferx-core and
# ferx-tools together, since they share one git source, plus whatever registry
# crates that revision needs beyond the current lock - nothing else, when the
# lock already satisfies it. Use that to pin a specific revision.
#
# Run from the repo root.

set -euo pipefail

cd "$(dirname "$0")/.."
RUST_DIR="src/rust"
CONFIG="$RUST_DIR/.cargo/config.toml"
LOCK="$RUST_DIR/Cargo.lock"

# The lock as it was before this run, so a failed update never leaves a broken
# pin behind, and the check's captured messages.
LOCK_BEFORE="$(mktemp)"
CHECK_ERR="$(mktemp)"
cp "$LOCK" "$LOCK_BEFORE"
BACKUP=""
cleanup() {
  rm -f "$LOCK_BEFORE" "$CHECK_ERR"
  if [[ -n "$BACKUP" ]]; then mv "$BACKUP" "$CONFIG"; fi
}
trap cleanup EXIT

# No config.toml (a fresh clone or worktree that has not been built yet) means no
# [patch] to take out of the way, so there is nothing to back up.
if [[ -f "$CONFIG" ]]; then
  cp "$CONFIG" "$CONFIG.bumplock.bak"
  BACKUP="$CONFIG.bumplock.bak"
  rm "$CONFIG"
fi

( cd "$RUST_DIR" && cargo update )

# Verify the lock pins both crates via a git source line, on one revision - the
# same check the R-CMD-check workflow runs. Only its error lines are shown: its
# generic repair advice would point back at this script.
if ! bash tools/check-ferx-core-pin.sh "$LOCK" >/dev/null 2>"$CHECK_ERR"; then
  grep '^error:' "$CHECK_ERR" >&2 || true
  cp "$LOCK_BEFORE" "$LOCK"
  echo "error: cargo update left a broken pin, so $LOCK is restored to what it was" >&2
  echo "error: before this run. Re-running will not help; investigate the errors above." >&2
  exit 2
fi

SHA=$(grep -m1 '^source = "git+https://github.com/FeRx-NLME/ferx-core' "$RUST_DIR/Cargo.lock" \
  | sed -E 's/.*#([0-9a-f]+).*/\1/')

echo
echo "ferx-core and ferx-tools now pinned to: $SHA"
echo "Commit suggestion:"
echo "  chore(deps): update Cargo.lock to ferx-core main (${SHA:0:7})"
