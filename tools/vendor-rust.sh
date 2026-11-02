#!/usr/bin/env bash
# Vendor the Rust dependencies into src/rust/vendor.tar.xz, so a source tarball
# built afterwards (`R CMD build .`) installs with no network and a read-only
# CARGO_HOME (ferx-r #666).
#
# The tarball is a release artifact, not a tracked file: it is ~18 MB and would
# grow git history by that much on every ferx-core pin bump. It is gitignored
# but not Rbuildignored, so `R CMD build` picks it up when it is present and a
# GitHub install (which builds from the repository) never sees it. src/Makevars
# takes the offline path when the tarball is there and the online one - cargo
# fetching the revision Cargo.lock pins - when it is not. Either way the build
# gets exactly the ferx-core commit the lock pins.
#
# The tarball holds three things:
#   vendor/              every crate in Cargo.lock, crates.io and the ferx-core
#                        git source alike (cargo vendor --locked)
#   vendor-Cargo.lock    the lock it was vendored from (Makevars checks it)
#   vendor-config.toml   the [source] replacement cargo vendor printed for this
#                        lock, with `directory` made relative. Captured rather
#                        than hand-written because the git source's key embeds
#                        the dependency spec (`?branch=main`), which a
#                        hand-written copy would silently fall out of step with.
#
# A tarball, not a vendor/ directory, because R CMD INSTALL drops hidden files
# and every vendored crate carries a .cargo-checksum.json.
#
# Vendors what Cargo.lock pins, never a sibling ../ferx-core: cargo vendor
# takes no [patch] from the wrapper, and this script refuses a config.toml that
# still carries one (a checkout not rebuilt since #353).
#
# Run from anywhere. Usage:
#   tools/vendor-rust.sh          write src/rust/vendor.tar.xz
#   tools/vendor-rust.sh --clean  remove it (and any unpacked vendor/)

set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"
RUST_DIR="$TOOLS_DIR/../src/rust"
cd "$RUST_DIR"

if [[ "${1:-}" == "--clean" ]]; then
  rm -rf vendor vendor.tar.xz vendor-config.toml vendor-Cargo.lock
  # An offline build appended the [source] replacement to the generated
  # config.toml; left there, it points plain cargo at the vendor/ just removed.
  # Reset it to what an online build writes.
  if [[ -f .cargo/config.toml ]] && grep -q 'vendored-sources' .cargo/config.toml; then
    echo '[build]' > .cargo/config.toml
  fi
  echo "vendor: removed src/rust/vendor.tar.xz"
  exit 0
fi

if [[ -f .cargo/config.toml ]] && grep -q '^\[patch' .cargo/config.toml; then
  echo "vendor: error: src/rust/.cargo/config.toml carries a [patch] (a checkout" >&2
  echo "vendor:   not rebuilt since #353); vendoring would pick up the sibling." >&2
  echo "vendor:   Delete that file and re-run." >&2
  exit 1
fi

"$TOOLS_DIR/check-ferx-core-pin.sh" Cargo.lock

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# --locked: vendor exactly what the lock pins, and fail rather than re-resolve.
# --versioned-dirs: stable directory names, independent of how many versions
# of one crate the graph holds.
cargo vendor --locked --versioned-dirs "$STAGE/vendor" > "$STAGE/raw-config.toml"

# cargo vendor prints `directory = "<absolute stage path>"`; Makevars needs it
# relative to src/rust, where the tarball is unpacked. Cargo resolves a relative
# `directory` in .cargo/config.toml against the directory holding .cargo/.
sed -E 's|^directory = ".*"$|directory = "vendor"|' "$STAGE/raw-config.toml" \
  | grep -v '^To use vendored sources' > "$STAGE/vendor-config.toml"
grep -q '^directory = "vendor"$' "$STAGE/vendor-config.toml" || {
  echo "vendor: error: cargo vendor printed no [source.vendored-sources] directory" >&2
  exit 1
}

# The lock must come out of this byte-for-byte unchanged.
"$TOOLS_DIR/check-ferx-core-pin.sh" Cargo.lock

# The lock this set was vendored from, so src/Makevars can refuse a tarball
# that a later pin bump has left stale, by name, instead of letting cargo fail
# on some missing crate.
cp Cargo.lock "$STAGE/vendor-Cargo.lock"

rm -f vendor.tar.xz
# COPYFILE_DISABLE keeps macOS tar from adding ._* AppleDouble members.
COPYFILE_DISABLE=1 tar -C "$STAGE" -cJf vendor.tar.xz vendor vendor-config.toml vendor-Cargo.lock

n=$(find "$STAGE/vendor" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
rev=$(grep -A2 '^name = "ferx-core"$' Cargo.lock | sed -nE 's/.*#([0-9a-f]+)".*/\1/p')
echo "vendor: wrote src/rust/vendor.tar.xz ($(du -h vendor.tar.xz | awk '{print $1}'), $n crates, ferx-core $rev)"
echo "vendor: now run R CMD build . - the source tarball will build offline."
echo "vendor: tools/vendor-rust.sh --clean removes it again."
