#!/bin/sh
# Run cargo in src/rust against a sibling ../ferx-core checkout, and leave
# src/rust/Cargo.lock exactly as it was found.
#
# Why this exists (ferx-r #353). The [patch] that redirects ferx-core and
# ferx-tools to the sibling used to be written into the persistent
# src/rust/.cargo/config.toml, so *every* later cargo resolve rewrote the lock:
# R CMD INSTALL, roxygen2::roxygenize(), pkgload::load_all(), cargo run by hand,
# an editor's rust-analyzer. An applied patch deletes both `source = "git+..."`
# pins - unpinning the crates for CI and for everyone who builds without the
# sibling - while an unused one appends [[patch.unused]] tables instead. Here
# the patch is handed to this one cargo invocation with --config, and the lock
# is snapshotted and put back around it, so nothing else in the checkout ever
# sees a patch.
#
# Building against the sibling is therefore an explicit opt-in: src/Makevars
# calls this script, and by hand it is
#
#   cd src && sh ../tools/sibling-cargo-build.sh check
#
# The sibling defaults to ../../ferx-core relative to the src/ directory this is
# run from; FERX_CORE_SIBLING overrides it with an absolute path.
#
# Usage (from the package's src/ directory, which is where Makevars runs):
#   sh ../tools/sibling-cargo-build.sh <cargo args...>
#
# Exits with cargo's status, or 1 when the build mixed revisions or the lock
# could not be put back.

set -eu

if [ "$#" -lt 1 ]; then
  echo "usage: $0 <cargo args...>   (run from the package's src/ directory)" >&2
  exit 2
fi

TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=tools/ferx-core-lock-lib.sh
. "$TOOLS_DIR/ferx-core-lock-lib.sh"

GIT_URL="https://github.com/FeRx-NLME/ferx-core"
RUST_DIR="rust"
LOCK="$RUST_DIR/Cargo.lock"
CHECK="$TOOLS_DIR/check-ferx-core-pin.sh"

# `|| true`: with SIGPIPE ignored (see the PIPE trap below) a write to a pipe
# nobody reads fails instead of killing the shell, and `set -e` would then end
# the build over an unread progress line.
say() { echo "ferx: $*" || true; }
err() { echo "ferx: $*" >&2 || true; }

if [ ! -f "$LOCK" ]; then
  err "ERROR $PWD/$LOCK not found - run this from the package's src/ directory."
  exit 2
fi

if [ -n "${FERX_CORE_SIBLING:-}" ]; then
  SIBLING="$FERX_CORE_SIBLING"
elif [ -d ../../ferx-core ]; then
  SIBLING="$(cd ../../ferx-core && pwd)"
else
  SIBLING="$(cd .. && pwd)/../ferx-core" # only ever reached by the error below
fi
# The path goes into a TOML string handed to `cargo --config`, where a quote
# ends the string and a backslash starts an escape. A Windows path is nothing
# but backslashes, and cargo reads forward slashes there just as well, so
# translate a drive-letter path and refuse anything else that carries one - on
# a POSIX filesystem a backslash is an ordinary character in a name, and
# rewriting it would point the build at some other directory.
case "$SIBLING" in
  *\"*)
    err "ERROR the sibling path contains a quote, which cannot be passed through"
    err "  cargo --config as a TOML path: $SIBLING"
    exit 2
    ;;
esac
case "$SIBLING" in
  *\\*)
    case "$SIBLING" in
      [A-Za-z]:[/\\]*) SIBLING=$(printf '%s' "$SIBLING" | tr '\\' '/') ;;
      *)
        err "ERROR the sibling path contains a backslash, which cannot be passed"
        err "  through cargo --config as a TOML path: $SIBLING"
        exit 2
        ;;
    esac
    ;;
esac

if [ ! -f "$SIBLING/Cargo.toml" ]; then
  err "ERROR no sibling ferx-core checkout at $SIBLING."
  err "  Without one there is nothing to patch: run cargo directly, and it will"
  err "  build the revision Cargo.lock pins."
  exit 2
fi
if [ ! -f "$SIBLING/crates/ferx-tools/Cargo.toml" ]; then
  err "ERROR the sibling ferx-core checkout has no crates/ferx-tools:"
  err "  $SIBLING"
  err "  (a checkout predating the ferx-tools workspace split). Patching ferx-core"
  err "  alone would build that checkout's engine against ferx-tools from GitHub"
  err "  main - two revisions of one workspace, with no compile error to say so."
  exit 2
fi

# -- what cargo will do with the patch ---------------------------------------
#
# Cargo applies a [patch] only when the patched crate's version *equals* the one
# Cargo.lock holds; anything else (an older sibling, a patch-level difference)
# leaves the entry unused, with no more than a "was not used in the crate graph"
# warning. ferx-core and ferx-tools carry their own version numbers, so the two
# can disagree - and a build that patches one and pins the other is the
# mixed-revision build this script exists to keep out. Predict it here, before
# spending a compile on it; the post-build verdict below re-checks what cargo
# actually resolved, which is what decides.

# The version a manifest declares, resolving `version.workspace = true` against
# the sibling's [workspace.package]. Empty when neither spells it out.
manifest_version() { # manifest
  mv_value=$(awk '
    /^\[/ { in_pkg = ($0 == "[package]"); next }
    in_pkg && /^version *= *"/ { gsub(/^version *= *"|".*$/, ""); print; exit }
    in_pkg && /^version\.workspace *= *true/ { print "@workspace"; exit }
  ' "$1")
  if [ "$mv_value" = "@workspace" ]; then
    mv_value=$(awk '
      /^\[/ { in_ws = ($0 == "[workspace.package]"); next }
      in_ws && /^version *= *"/ { gsub(/^version *= *"|".*$/, ""); print; exit }
    ' "$SIBLING/Cargo.toml")
  fi
  echo "$mv_value"
}

# yes / no / unknown: whether cargo will take this crate from the sibling.
would_patch() { # sibling version, locked version
  if [ -z "$1" ] || [ -z "$2" ]; then
    echo unknown
  elif [ "$1" = "$2" ]; then
    echo yes
  else
    echo no
  fi
}

core_locked_version=$(ferx_lock_version_of ferx-core "$LOCK")
tools_locked_version=$(ferx_lock_version_of ferx-tools "$LOCK")
core_sibling_version=$(manifest_version "$SIBLING/Cargo.toml")
tools_sibling_version=$(manifest_version "$SIBLING/crates/ferx-tools/Cargo.toml")
core_expected=$(would_patch "$core_sibling_version" "$core_locked_version")
tools_expected=$(would_patch "$tools_sibling_version" "$tools_locked_version")

versions_line="ferx-core: lock ${core_locked_version:-?}, sibling ${core_sibling_version:-?}; ferx-tools: lock ${tools_locked_version:-?}, sibling ${tools_sibling_version:-?}"

if { [ "$core_expected" = yes ] && [ "$tools_expected" = no ]; } ||
   { [ "$core_expected" = no ] && [ "$tools_expected" = yes ]; }; then
  err "ERROR the sibling checkout would supply only one of the two crates:"
  err "  $versions_line"
  err "  cargo applies a [patch] only when the sibling's version equals the locked"
  err "  one, so this build would take one crate from $SIBLING"
  err "  and the other from the revision Cargo.lock pins - two revisions of one"
  err "  workspace, with no compile error to say so. Refusing to build."
  err "  Fix it by bumping the lock (tools/update-ferx-core-lock.sh), by moving the"
  err "  sibling to a matching revision, or by building wholly from the pin with"
  err "  MAKEFLAGS=\"LOCAL_FERX_CORE=\" R CMD INSTALL ."
  exit 1
fi

# -- snapshot the lock -------------------------------------------------------

# Whether the lock was pinned *before* this build. A lock that some earlier
# direct cargo run already stripped is put back exactly as found, and saying
# "restored to its pin" about it would be a lie.
pin_before=unknown
if [ -f "$CHECK" ] && command -v bash >/dev/null 2>&1; then
  if bash "$CHECK" "$LOCK" >/dev/null 2>&1; then
    pin_before=intact
  else
    pin_before=broken
  fi
fi

# mktemp, not a fixed name under target/: two builds in one checkout would
# otherwise share one snapshot, and the second to finish would restore - or fail
# to find - the first one's file.
SNAPSHOT=$(mktemp "${TMPDIR:-/tmp}/ferx-Cargo.lock.XXXXXX")
cp "$LOCK" "$SNAPSHOT"
restored=unknown

restore_lock() {
  [ -n "${SNAPSHOT:-}" ] && [ -f "$SNAPSHOT" ] || return 0
  if cmp -s "$SNAPSHOT" "$LOCK" 2>/dev/null; then
    restored=untouched
  elif cp "$SNAPSHOT" "$LOCK" 2>/dev/null; then
    restored=rewritten
  else
    restored=failed
    err "ERROR could not restore $PWD/$LOCK from $SNAPSHOT - the snapshot is kept."
    return 0
  fi
  rm -f "$SNAPSHOT"
}

on_signal() {
  restore_lock
  trap - EXIT
  exit 130
}

trap restore_lock EXIT
trap on_signal HUP INT TERM
# dash - /bin/sh on Debian and Ubuntu - does not run an EXIT trap when a write
# hits a closed pipe: SIGPIPE kills the shell outright. That is not exotic here;
# `R CMD INSTALL . 2>&1 | head` and an interrupted load_all() (whose SIGKILL
# reaches only the R process, leaving make, sh and cargo writing into a pipe
# nobody reads) both produce it, and the lock would stay stripped. Ignoring
# SIGPIPE turns the write into an ordinary EIO failure, which does run the trap.
trap '' PIPE

# -- build -------------------------------------------------------------------

say "patching ferx-core + ferx-tools to the sibling checkout, for this cargo run only"
say "  $SIBLING"
if [ "$core_expected" = no ] && [ "$tools_expected" = no ]; then
  say "  (both versions differ from the lock, so cargo will most likely ignore it: $versions_line)"
fi

status=0
(
  cd "$RUST_DIR" &&
    cargo "$@" \
      --config "patch.\"$GIT_URL\".ferx-core.path=\"$SIBLING\"" \
      --config "patch.\"$GIT_URL\".ferx-tools.path=\"$SIBLING/crates/ferx-tools\""
) || status=$?

# -- per-crate verdict, read from the lock cargo just wrote ------------------

core_origin=$(ferx_lock_origin_of ferx-core "$LOCK")
tools_origin=$(ferx_lock_origin_of ferx-tools "$LOCK")
unused=$(ferx_lock_patch_unused_names "$LOCK" | tr '\n' ' ' | sed 's/ *$//')

# Put the lock back before saying anything about the build: a verdict printed
# first is a verdict that can be lost with the shell that was about to print it.
restore_lock

case "$pin_before,$restored" in
  *,failed)
    err "ERROR Cargo.lock is left as cargo wrote it. Run tools/check-ferx-core-pin.sh."
    exit 1
    ;;
  *,unknown)
    err "ERROR the lock snapshot $SNAPSHOT vanished before it could be put back, so"
    err "  Cargo.lock is as cargo left it. Run tools/check-ferx-core-pin.sh."
    exit 1
    ;;
  intact,untouched) say "Cargo.lock untouched by this build" ;;
  intact,rewritten) say "Cargo.lock restored to its pin" ;;
  broken,*)
    say "WARNING Cargo.lock was ALREADY unpinned before this build, so it has been"
    say "  put back as it was found, not to a pin. Some earlier cargo run wrote it."
    say "  Run tools/check-ferx-core-pin.sh for what is wrong, and restore the lock"
    say "  as CLAUDE.md's \"ferx-core dependency\" section describes."
    ;;
  *) say "Cargo.lock put back as found (its pin was not checked: no bash on PATH)" ;;
esac

if [ "$status" -ne 0 ]; then
  err "cargo exited $status"
  exit "$status"
fi

# `path` means an applied patch: cargo replaced the git pin with the sibling
# directory. `git` means the patch went unused and the pinned revision was
# built. Deciding per crate is the whole point - a verdict taken from one of
# them calls the mixed build "the sibling was not used".
origin_word() { # git | path | other | absent
  case "$1" in
    path) echo "the sibling" ;;
    git) echo "the pinned revision" ;;
    *) echo "somewhere else ($1)" ;;
  esac
}

case "$core_origin,$tools_origin" in
  path,path)
    say "the sibling supplied BOTH ferx-core and ferx-tools"
    ;;
  git,git)
    say "WARNING the sibling was NOT used: ferx-core and ferx-tools both come from"
    say "  the revision Cargo.lock pins."
    if [ "$core_expected" = yes ] && [ "$tools_expected" = yes ]; then
      say "  Both versions match the lock ($versions_line), so cargo was expected to"
      say "  take the patch. It did not, which means this run did not resolve - an"
      say "  unchanged lock and a warm cache, or a cargo subcommand that reads the"
      say "  lock without writing it."
    else
      say "  cargo applies a [patch] only when the sibling's version equals the"
      say "  locked one - $versions_line"
    fi
    if [ -n "$unused" ]; then say "  cargo reported these patch entries unused: $unused"; fi
    ;;
  *)
    err "ERROR this build MIXED revisions: ferx-core came from $(origin_word "$core_origin"),"
    err "  ferx-tools from $(origin_word "$tools_origin")."
    err "  $versions_line"
    if [ -n "$unused" ]; then err "  cargo reported these patch entries unused: $unused"; fi
    err "  The two crates are halves of one workspace and move together, so what was"
    err "  just built is a combination that exists nowhere. Bump the lock"
    err "  (tools/update-ferx-core-lock.sh), move the sibling to a matching revision,"
    err "  or build wholly from the pin with MAKEFLAGS=\"LOCAL_FERX_CORE=\"."
    exit 1
    ;;
esac
