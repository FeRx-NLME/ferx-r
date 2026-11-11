#!/usr/bin/env bash
# Drive tools/docker-release-tags.sh with a stub `docker` and a stub `git`,
# one case per outcome, and assert each exit code and message. The script only
# runs for real on a v* tag push, so without this nothing in CI would fail if
# one of its rules were deleted.
#
# Usage: tools/test-docker-release-tags.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/docker-release-tags.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
failures=0
cases=0

# Stub docker: `docker buildx imagetools inspect <ref> --format <fmt>`, answered
# from $STUB/<ref with / : @ turned into _>.{err,mt,kids}.
mkdir -p "$TMP/bin"
cat > "$TMP/bin/docker" <<'EOF'
#!/bin/sh
[ "$1 $2 $3" = "buildx imagetools inspect" ] || { echo "stub docker: unexpected $*" >&2; exit 99; }
key="$(printf '%s' "$4" | tr '/:@' '___')"
f="$STUB/$key"
if [ -f "$f.err" ]; then cat "$f.err" >&2; exit 1; fi
[ -f "$f.mt" ] || { echo "ERROR: $4: not found" >&2; exit 1; }
case "$6" in
  *MediaType*) cat "$f.mt" ;;
  *Manifests*) cat "$f.kids" 2>/dev/null || true ;;
  *) echo "stub docker: unexpected format $6" >&2; exit 99 ;;
esac
EOF
# Stub git: `git ls-remote ...` prints $STUB/ls-remote, or fails when absent.
cat > "$TMP/bin/git" <<'EOF'
#!/bin/sh
[ "$1" = "ls-remote" ] || { echo "stub git: unexpected $*" >&2; exit 99; }
[ -f "$STUB/ls-remote" ] || { echo "fatal: unable to access 'https://github.com/x/y/': Could not resolve host" >&2; exit 128; }
cat "$STUB/ls-remote"
EOF
chmod +x "$TMP/bin/docker" "$TMP/bin/git"

IMG=ghcr.io/ferx-nlme/ferx-r
INDEX=application/vnd.oci.image.index.v1+json
MANIFEST=application/vnd.oci.image.manifest.v1+json
d() { printf 'sha256:%064d' "$1"; }   # a fake digest

new_case() { STUB="$TMP/case$cases"; rm -rf "$STUB"; mkdir -p "$STUB/digests"; export STUB; }
key() { printf '%s' "$1" | tr '/:@' '___'; }
index() { # ref, child digests...
  local k; k="$(key "$1")"; shift
  echo "$INDEX" > "$STUB/$k.mt"; printf '%s\n' "$@" > "$STUB/$k.kids"
}
single() { echo "$MANIFEST" > "$STUB/$(key "$1").mt"; }
error() { echo "$2" > "$STUB/$(key "$1").err"; }
# The run's two per-arch images: index A (amd64 + attestation), index B.
run_images() {
  index "$IMG@$(d 1)" "$(d 11)" "$(d 12)"; touch "$STUB/digests/$(printf '%064d' 1)"
  index "$IMG@$(d 2)" "$(d 21)" "$(d 22)"; touch "$STUB/digests/$(printf '%064d' 2)"
}
tags() { # tag lines for the ls-remote stub
  local t; : > "$STUB/ls-remote"
  for t in "$@"; do printf '%040d\trefs/tags/%s\n' 0 "$t" >> "$STUB/ls-remote"; done
}

# Every POSIX shell on the box: CI's /bin/sh is dash, a Mac's is bash. (Not
# zsh: it is not POSIX by default, and the script is #!/bin/sh.)
SHELLS=""
for s in sh dash bash ksh; do command -v "$s" >/dev/null 2>&1 && SHELLS="$SHELLS $s"; done

expect() { # name, want-exit, want-fragment, args...
  local name=$1 want=$2 frag=$3 out rc sh; shift 3
  for sh in $SHELLS; do
    cases=$((cases + 1)); rc=0
    out="$(PATH="$TMP/bin:$PATH" "$sh" "$SCRIPT" "$@" 2>&1)" || rc=$?
    if [ "$rc" -ne "$want" ] || ! printf '%s' "$out" | grep -qF -- "$frag"; then
      echo "FAIL [$sh] $name: exit $rc (want $want), want text: $frag"
      printf '%s\n' "$out" | sed 's/^/     | /'
      failures=$((failures + 1))
    else
      echo "ok   [$sh] $name"
    fi
  done
}

# ---- guard ---------------------------------------------------------------
new_case
expect "guard: neither tag published -> free" 0 "ferx-r:0.5.0 is not published yet" \
  guard "$IMG" v0.5.0

new_case; index "$IMG:v0.5.0" "$(d 11)" "$(d 12)" "$(d 21)" "$(d 22)"
expect "guard (early): v-tag published -> refused" 1 "Re-run failed jobs" \
  guard "$IMG" v0.5.0

new_case; index "$IMG:0.5.0" "$(d 11)"
expect "guard (early): plain tag published -> refused" 1 "ferx-r:0.5.0 is already published" \
  guard "$IMG" v0.5.0

new_case; run_images; index "$IMG:v0.5.0" "$(d 22)" "$(d 11)" "$(d 21)" "$(d 12)"
expect "guard: v-tag holds exactly this run's images -> resume" 0 "completing an earlier publish" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; run_images; index "$IMG:v0.5.0" "$(d 11)" "$(d 12)" "$(d 21)" "$(d 99)"
expect "guard: v-tag holds other images -> refused" 1 "already published with different images" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; run_images; index "$IMG:v0.5.0" "$(d 11)" "$(d 12)"
expect "guard: v-tag holds a subset of this run -> refused" 1 "already published with different images" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; run_images; index "$IMG:v0.5.0" "$(d 11)" "$(d 12)" "$(d 21)" "$(d 22)"
index "$IMG:0.5.0" "$(d 11)" "$(d 12)" "$(d 21)" "$(d 99)"
expect "guard: v-tag same, plain tag different -> refused" 1 "ferx-r:0.5.0 is already published with different images" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; run_images; single "$IMG:v0.5.0"
expect "guard: tag is a single manifest -> refused" 1 "single manifest" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; error "$IMG:v0.5.0" "ERROR: failed to authorize: 401 Unauthorized"
expect "guard: registry says 401 -> could not check" 2 "could not tell whether" \
  guard "$IMG" v0.5.0

new_case; run_images; index "$IMG:v0.5.0" "$(d 11)"; error "$IMG@$(d 2)" "ERROR: 503 Service Unavailable"
expect "guard: this run's image unreadable -> could not check" 2 "could not read this run's image" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case; index "$IMG:v0.5.0" "$(d 11)"
expect "guard: empty digests dir -> could not check" 2 "no digests in" \
  guard "$IMG" v0.5.0 "$STUB/digests"

new_case
expect "guard: pre-release tags are guarded too" 0 "ferx-r:0.5.0-rc1 is not published yet" \
  guard "$IMG" v0.5.0-rc1

# ---- is-newest -----------------------------------------------------------
new_case; tags v0.3.0 v0.4.0 v0.5.0 v0.5.0-rc1
expect "is-newest: highest final release -> true" 0 "true" is-newest https://github.com/x/y v0.5.0

new_case; tags v0.4.0 v0.4.1 v0.5.0
expect "is-newest: backport to an older line -> false" 0 "false" is-newest https://github.com/x/y v0.4.1

new_case; tags v0.9.0 v0.10.0
expect "is-newest: version order, not text order -> true" 0 "newest final release tag: v0.10.0" \
  is-newest https://github.com/x/y v0.10.0

new_case; tags v0.5.0 v0.6.0-rc1
expect "is-newest: pre-release -> false, no lookup" 0 "not a final release tag" \
  is-newest https://github.com/x/y v0.6.0-rc1

new_case
expect "is-newest: tag lookup fails -> could not check" 2 "could not list the tags" \
  is-newest https://github.com/x/y v0.5.0

new_case; tags v0.3.0 v0.4.0
expect "is-newest: own tag missing from the list -> could not check" 2 "the lookup is broken" \
  is-newest https://github.com/x/y v0.5.0

# ---- usage ---------------------------------------------------------------
new_case
expect "usage: unknown mode -> 2" 2 "usage" bogus

echo "$((cases - failures))/$cases passed"
[ "$failures" -eq 0 ]
