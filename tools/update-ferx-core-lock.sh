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
# --precise <sha>` also works (measured 2026-09-11) and moves only ferx-core and
# ferx-tools, which share one git source; use that to pin a specific revision.
#
# Run from the repo root.

set -euo pipefail

cd "$(dirname "$0")/.."
RUST_DIR="src/rust"
CONFIG="$RUST_DIR/.cargo/config.toml"

# No config.toml (a fresh clone or worktree that has not been built yet) means no
# [patch] to take out of the way, so there is nothing to back up.
if [[ -f "$CONFIG" ]]; then
  BACKUP="$CONFIG.bumplock.bak"
  cp "$CONFIG" "$BACKUP"
  trap 'mv "$BACKUP" "$CONFIG"' EXIT
  rm "$CONFIG"
fi

( cd "$RUST_DIR" && cargo update )

# Verify the lock pins both crates via a git source line, on one revision - the
# same check the R-CMD-check workflow runs.
if ! bash tools/check-ferx-core-pin.sh "$RUST_DIR/Cargo.lock" >/dev/null; then
  echo "error: refusing to leave a broken pin after update" >&2
  exit 2
fi

SHA=$(grep -m1 '^source = "git+https://github.com/FeRx-NLME/ferx-core' "$RUST_DIR/Cargo.lock" \
  | sed -E 's/.*#([0-9a-f]+).*/\1/')

echo
echo "ferx-core and ferx-tools now pinned to: $SHA"
echo "Commit suggestion:"
echo "  chore(deps): update Cargo.lock to ferx-core main (${SHA:0:7})"
