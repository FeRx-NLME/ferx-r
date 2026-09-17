#!/usr/bin/env bash
# Check that src/rust/Cargo.lock still pins ferx-core and ferx-tools to the
# ferx-core GitHub repository, on one revision, with no [[patch.unused]] tables.
#
# Why: any cargo command that resolves while a [patch] to a sibling ../ferx-core
# checkout is in place rewrites the lock. An applied patch deletes both
# `source = "git+..."` lines (unpinning the crates for CI and everyone who
# builds without the sibling); an unused one appends [[patch.unused]] tables
# instead. Since #353 only tools/sibling-cargo-build.sh passes such a patch, to
# its own cargo run, and it restores the lock afterwards - but a checkout that
# has not been rebuilt since then still has the old [patch] in
# src/rust/.cargo/config.toml, and `cargo update` and a `kill -9` reach the lock
# either way.
#
# Runs no cargo, so it is always safe. Used by the R-CMD-check workflow and by
# tools/update-ferx-core-lock.sh, and exercised by tools/test-check-ferx-core-pin.sh.
# Exits 0 when the pin is intact, 1 otherwise.
#
# Usage: tools/check-ferx-core-pin.sh [path/to/Cargo.lock]

set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"
# The lock parsing is shared with tools/sibling-cargo-build.sh, which reports
# per crate what a local sibling build resolved.
# shellcheck source=tools/ferx-core-lock-lib.sh
. "$TOOLS_DIR/ferx-core-lock-lib.sh"

LOCK="${1:-$TOOLS_DIR/../src/rust/Cargo.lock}"
status=0

fail() {
  if [[ "${GITHUB_ACTIONS:-}" == true ]]; then
    echo "::error file=src/rust/Cargo.lock::$1"
  else
    echo "error: $1" >&2
  fi
  status=1
}

if [[ ! -f "$LOCK" ]]; then
  fail "$LOCK not found"
  exit 1
fi

for pkg in ferx-core ferx-tools; do
  case "$(ferx_lock_source_of "$pkg" "$LOCK")" in
    "$FERX_GIT_SOURCE_PREFIX"*) ;;
    *) fail "$pkg lacks a git source pin - a cargo run with the sibling ../ferx-core [patch] applied strips it." ;;
  esac
done

core_rev=$(ferx_lock_source_of ferx-core "$LOCK" | sed -E 's/.*#([0-9a-f]+).*/\1/')
tools_rev=$(ferx_lock_source_of ferx-tools "$LOCK" | sed -E 's/.*#([0-9a-f]+).*/\1/')
if [[ -n "$core_rev" && -n "$tools_rev" && "$core_rev" != "$tools_rev" ]]; then
  fail "ferx-core ($core_rev) and ferx-tools ($tools_rev) are pinned to different revisions of the same repository."
fi

if grep -q '^\[\[patch\.unused\]\]' "$LOCK"; then
  fail "Cargo.lock carries [[patch.unused]] tables - written by a cargo run whose sibling [patch] went unused; they belong to no change."
fi

if [[ "$status" -ne 0 ]]; then
  if [[ "${GITHUB_ACTIONS:-}" == true ]]; then
    cat <<'EOF'
The committed Cargo.lock is broken, so there is nothing local to check out. Locally,
restore it from main (`git checkout origin/main -- src/rust/Cargo.lock`) and redo any
wanted lock change the way CLAUDE.md's "ferx-core dependency" section describes, or, if
this PR bumps the pin, re-run tools/update-ferx-core-lock.sh. Then commit the lock.
EOF
  else
    cat >&2 <<'EOF'
Read `git diff src/rust/Cargo.lock` before restoring:
- only the damage above: `git checkout -- src/rust/Cargo.lock`.
- a pin bump you meant: re-run tools/update-ferx-core-lock.sh.
- any other change you meant (e.g. a new dependency): `git checkout -- src/rust/Cargo.lock`,
  run cargo in src/rust (e.g. `cargo metadata --format-version 1 >/dev/null`), and run this
  check again. A plain cargo run there carries no [patch], so it records that change and
  nothing else. Not the bump script: its whole-graph `cargo update` also moves the pin to
  ferx-core main.
EOF
  fi
  exit 1
fi

echo "ferx-core and ferx-tools pinned to $core_rev"
