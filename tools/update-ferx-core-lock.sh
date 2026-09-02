#!/usr/bin/env bash
# Bump the ferx-core / ferx-tools pin in src/rust/Cargo.lock to current main HEAD.
#
# Both crates live in the one ferx-core repository and are patched to the local
# sibling checkout together, so they are bumped together and must end up on the
# same revision.
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

# Verify the lock still pins both crates via a git source line, on one revision.
source_of() {
  awk -v pkg="$1" '
    $0 == "name = \"" pkg "\"" { in_pkg = 1; next }
    in_pkg && /^source = / { print; exit }
    in_pkg && /^\[\[package\]\]/ { exit }
  ' "$RUST_DIR/Cargo.lock"
}

for pkg in ferx-core ferx-tools; do
  case "$(source_of "$pkg")" in
    'source = "git+https://github.com/FeRx-NLME/ferx-core'*) ;;
    *)
      echo "error: Cargo.lock has no git source for $pkg after update - refusing to commit a broken pin" >&2
      exit 2
      ;;
  esac
done

SHA=$(source_of ferx-core | sed -E 's/.*#([0-9a-f]+).*/\1/')
TOOLS_SHA=$(source_of ferx-tools | sed -E 's/.*#([0-9a-f]+).*/\1/')
if [[ "$SHA" != "$TOOLS_SHA" ]]; then
  echo "error: ferx-core ($SHA) and ferx-tools ($TOOLS_SHA) landed on different revisions" >&2
  exit 2
fi

echo
echo "ferx-core and ferx-tools now pinned to: $SHA"
echo "Commit suggestion:"
echo "  chore(deps): update Cargo.lock to ferx-core main (${SHA:0:7})"
