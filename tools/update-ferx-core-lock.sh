#!/usr/bin/env bash
# Bump the ferx-core pin in src/rust/Cargo.lock to current main HEAD.
#
# Why this script exists: src/rust/.cargo/config.toml carries a [patch] that
# redirects ferx-core to a sibling ../ferx-core checkout when present. If you
# run any `cargo` command that rewrites the lock while that patch is active,
# cargo writes a *path*-style lock entry with no `source = "git+..."` line —
# which silently unpins ferx-core for everyone who builds without the sibling
# (CI, downstream users). That happened twice (commits 1ce7f59, b96c867).
#
# This script temporarily removes the patch so cargo resolves ferx-core from
# GitHub and writes the correct git+https pin. `cargo update -p ferx-core`
# does NOT work here because with the patch absent the package spec is the
# full git URL, not a bare name; plain `cargo update` is used instead.
#
# Run from the repo root.

set -euo pipefail

cd "$(dirname "$0")/.."
RUST_DIR="src/rust"
CONFIG="$RUST_DIR/.cargo/config.toml"

if [[ ! -f "$CONFIG" ]]; then
  echo "error: $CONFIG not found — run from a fresh checkout where Makevars has populated it" >&2
  exit 1
fi

BACKUP="$CONFIG.bumplock.bak"
cp "$CONFIG" "$BACKUP"
trap 'mv "$BACKUP" "$CONFIG"' EXIT
rm "$CONFIG"

( cd "$RUST_DIR" && cargo update )

# Verify the lock still pins ferx-core via a git source line.
if ! grep -q '^source = "git+https://github.com/FeRx-NLME/ferx-core' "$RUST_DIR/Cargo.lock"; then
  echo "error: Cargo.lock has no git source for ferx-core after update — refusing to commit a broken pin" >&2
  exit 2
fi

SHA=$(grep -A2 '^name = "ferx-core"$' "$RUST_DIR/Cargo.lock" \
      | grep '^source = ' \
      | sed -E 's/.*#([0-9a-f]+).*/\1/')

echo
echo "ferx-core now pinned to: $SHA"
echo "Commit suggestion:"
echo "  chore(deps): update Cargo.lock to ferx-core main (${SHA:0:7})"
