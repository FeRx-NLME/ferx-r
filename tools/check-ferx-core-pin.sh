#!/usr/bin/env bash
# Check that src/rust/Cargo.lock still pins ferx-core and ferx-tools to the
# ferx-core GitHub repository, on one revision, with no [[patch.unused]] tables.
#
# Why: src/Makevars patches both crates to a sibling ../ferx-core checkout when
# one exists, and any cargo command that resolves while that [patch] is in place
# rewrites the lock - R CMD INSTALL, roxygenize() and load_all() included, as
# well as cargo run by hand or by an editor's rust-analyzer. An applied patch
# deletes both `source = "git+..."` lines (unpinning the crates for CI and
# everyone who builds without the sibling); an unused one appends
# [[patch.unused]] tables instead.
#
# Runs no cargo, so it is always safe. Used by the R-CMD-check workflow and by
# tools/update-ferx-core-lock.sh, and exercised by tools/test-check-ferx-core-pin.sh.
# Exits 0 when the pin is intact, 1 otherwise.
#
# Usage: tools/check-ferx-core-pin.sh [path/to/Cargo.lock]

set -euo pipefail

LOCK="${1:-$(dirname "$0")/../src/rust/Cargo.lock}"
GIT_SOURCE='source = "git+https://github.com/FeRx-NLME/ferx-core'
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

# The `source = ` line of a package's [[package]] table; empty for a path
# package, which is what an applied [patch] leaves behind.
source_of() {
  awk -v pkg="$1" '
    $0 == "name = \"" pkg "\"" { in_pkg = 1; next }
    in_pkg && /^source = / { print; exit }
    in_pkg && /^\[/ { exit }
  ' "$LOCK"
}

for pkg in ferx-core ferx-tools; do
  case "$(source_of "$pkg")" in
    "$GIT_SOURCE"*) ;;
    *) fail "$pkg lacks a git source pin - a cargo run with the sibling ../ferx-core [patch] applied strips it." ;;
  esac
done

core_rev=$(source_of ferx-core | sed -E 's/.*#([0-9a-f]+).*/\1/')
tools_rev=$(source_of ferx-tools | sed -E 's/.*#([0-9a-f]+).*/\1/')
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
  move src/rust/.cargo/config.toml aside, run cargo in src/rust (e.g. `cargo metadata
  --format-version 1 >/dev/null`), move config.toml back, and run this check again.
  Not the bump script: its whole-graph `cargo update` also moves the pin to ferx-core main.
EOF
  fi
  exit 1
fi

echo "ferx-core and ferx-tools pinned to $core_rev"
