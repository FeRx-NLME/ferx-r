# POSIX sh helpers for reading src/rust/Cargo.lock.
#
# Sourced by tools/check-ferx-core-pin.sh (which verifies a committed lock) and
# by tools/sibling-cargo-build.sh (which reports, per crate, what a local
# sibling build actually resolved). Both need the same answer to "where did this
# package come from", so the parsing lives here once.
#
# Not executable, and it runs nothing on its own: source it, then call the
# functions with an explicit lock path.

# The prefix every ferx-core / ferx-tools `source =` line must carry. Matched as
# a literal, never as a regex: it contains `+`, `.` and `/`.
FERX_GIT_SOURCE_PREFIX='source = "git+https://github.com/FeRx-NLME/ferx-core'

# The `source = ` line of a package's [[package]] table, verbatim, or empty for
# a path package - which is what an applied [patch] leaves behind.
#
# The scan stops at the end of the package's own table, so a stripped crate
# cannot borrow the `source =` of the table that follows it. `[[patch.unused]]`
# tables repeat the crate names, but cargo writes them after every
# [[package]] table, so the first `name = ` hit is always the real one.
ferx_lock_source_of() { # package, lockfile
  awk -v pkg="$1" '
    $0 == "name = \"" pkg "\"" { in_pkg = 1; next }
    in_pkg && /^source = / { print; exit }
    in_pkg && /^\[/ { exit }
  ' "$2"
}

# The version a package's [[package]] table declares, e.g. 0.4.0. Empty when the
# lock has no such package.
ferx_lock_version_of() { # package, lockfile
  awk -v pkg="$1" '
    $0 == "name = \"" pkg "\"" { in_pkg = 1; next }
    in_pkg && /^version = "/ { gsub(/^version = "|".*$/, ""); print; exit }
    in_pkg && /^\[/ { exit }
  ' "$2"
}

# The crate names under the lock`s [[patch.unused]] tables, one per line. Cargo
# writes one such table per patch entry it was given and did not use, so this is
# the per-crate record of which half of the sibling was ignored.
ferx_lock_patch_unused_names() { # lockfile
  awk '
    /^\[\[patch\.unused\]\]/ { in_patch = 1; next }
    in_patch && /^name = "/ { gsub(/^name = "|".*$/, ""); print; in_patch = 0; next }
    in_patch && /^\[/ { in_patch = 0 }
  ' "$1"
}

# "git" when the package is pinned to the ferx-core repository, "path" when an
# applied [patch] replaced it with a local directory, "other" for a git pin to
# some other repository, and "absent" when the lock has no such package.
ferx_lock_origin_of() { # package, lockfile
  ferx_lock_source_of_line=$(ferx_lock_source_of "$1" "$2")
  if [ -z "$ferx_lock_source_of_line" ]; then
    if [ -z "$(ferx_lock_version_of "$1" "$2")" ]; then
      echo absent
    else
      echo path
    fi
  else
    case "$ferx_lock_source_of_line" in
      "$FERX_GIT_SOURCE_PREFIX"*) echo git ;;
      *) echo other ;;
    esac
  fi
}
