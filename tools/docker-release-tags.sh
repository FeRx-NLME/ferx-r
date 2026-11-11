#!/bin/sh
# Release-tag rules for the published image (.github/workflows/docker.yml).
#
#   tools/docker-release-tags.sh guard <image> <git-tag> [<digests-dir>]
#   tools/docker-release-tags.sh is-newest <repo-url> <git-tag>
#
# guard: release tags (:vX.Y.Z and :X.Y.Z) are permanent - validated work
# cites them - so a release image is published once and never replaced. A
# rebuild is not byte-identical (apt and CRAN move), so a re-pushed git tag or
# a full re-run must not repoint them. For each of the two tags:
#   - not in the registry                              -> free
#   - there, and <digests-dir> given, holding exactly the images this run
#     would publish (same platform and attestation manifests)  -> free: this
#     completes an earlier publish of the same bytes that stopped part-way
#     ("Re-run failed jobs" reuses the build job's digests)
#   - there, otherwise                                 -> refused
#   - registry gives any other answer                  -> could not check
# <digests-dir> holds one empty file per pushed per-architecture image, named
# by its sha256 digest (the build job's "Export digest" step). Without it, any
# existing tag is refused: the early check in the build job, before an hour
# of builds and before anything is pushed.
#
# is-newest: prints `true` when <git-tag> is a final release (plain vX.Y.Z,
# no -rc1 etc.) and the highest such tag in <repo-url>, else `false`. That
# decides the floating :release tag, so a backport (v0.4.1 after v0.5.0) or a
# pre-release leaves it alone. The tag that triggered the run must be in the
# remote's list: if the lookup fails or it is missing, that is "could not
# check" (exit 2), never a quiet `false`.
#
# Exit codes: 0 ok, 1 refused (guard only), 2 could not check.
# tools/test-docker-release-tags.sh drives both modes with stub docker/git.
set -eu

err() { echo "docker-release-tags: $*" >&2; }

TMPD="$(mktemp -d)"
trap 'rm -rf "$TMPD"' EXIT

# inspect <ref> <format>: 0 = found (output on stdout), 1 = not found,
# 2 = any other answer (printed on stderr).
inspect() {
  if docker buildx imagetools inspect "$1" --format "$2" >"$TMPD/out" 2>"$TMPD/err"; then
    cat "$TMPD/out"; return 0
  fi
  case "$(cat "$TMPD/err")" in
    *"not found"*|*"manifest unknown"*|*NOT_FOUND*|*MANIFEST_UNKNOWN*) return 1 ;;
  esac
  sed 's/^/  /' "$TMPD/err" >&2
  return 2
}

# children <ref> [<own digest>]: the sorted digests of the manifests an index
# holds (platform images and attestations); a single manifest stands for
# itself when its digest is known. 0 = listed, 1 = not found, 2 = could not
# read, 3 = a single manifest whose digest is unknown (a tag).
children() {
  rc=0; mt="$(inspect "$1" '{{.Manifest.MediaType}}')" || rc=$?
  [ "$rc" -eq 0 ] || return "$rc"
  case "$mt" in
    *image.index*|*manifest.list*)
      rc=0; list="$(inspect "$1" '{{range .Manifest.Manifests}}{{.Digest}}{{println}}{{end}}')" || rc=$?
      [ "$rc" -eq 0 ] || { [ "$rc" -eq 1 ] && return 2; return "$rc"; }
      printf '%s\n' "$list" | grep -E '^sha256:[0-9a-f]{64}$' | sort -u ;;
    *)
      [ -n "${2:-}" ] || return 3
      printf '%s\n' "$2" ;;
  esac
}

guard() {
  image="$1"; tag="$2"; digests="${3:-}"
  expected=""
  for t in "$tag" "${tag#v}"; do
    rc=0; have="$(children "$image:$t")" || rc=$?
    case "$rc" in
      0|3) ;;
      1) echo "$image:$t is not published yet."; continue ;;
      *) err "could not tell whether $image:$t exists (see above)"; return 2 ;;
    esac
    if [ -z "$digests" ]; then
      err "$image:$t is already published. Release images are never replaced: tag a new version. (If an earlier publish of this release stopped part-way, use 'Re-run failed jobs' on that run instead.)"
      return 1
    fi
    if [ "$rc" -eq 3 ]; then
      err "$image:$t is already published as a single manifest, not this run's multi-arch image. Release images are never replaced."
      return 1
    fi
    if [ -z "$expected" ]; then
      set -- "$digests"/*
      [ -e "$1" ] || { err "no digests in $digests"; return 2; }
      for f in "$@"; do
        d="sha256:$(basename "$f")"
        c=0; kids="$(children "$image@$d" "$d")" || c=$?
        [ "$c" -eq 0 ] || { err "could not read this run's image $image@$d"; return 2; }
        expected="$expected$kids
"
      done
      expected="$(printf '%s' "$expected" | sort -u)"
    fi
    if [ "$have" = "$expected" ]; then
      echo "$image:$t already points at exactly this run's images: completing an earlier publish of the same bytes."
    else
      err "$image:$t is already published with different images. Release images are never replaced: tag a new version."
      printf '  published:\n%s\n  this run:\n%s\n' "$(printf '%s\n' "$have" | sed 's/^/    /')" "$(printf '%s\n' "$expected" | sed 's/^/    /')" >&2
      return 1
    fi
  done
  return 0
}

PLAIN='^v[0-9]+\.[0-9]+\.[0-9]+$'

is_newest() {
  url="$1"; tag="$2"
  if ! printf '%s\n' "$tag" | grep -Eq "$PLAIN"; then
    echo "$tag is not a final release tag: :release stays where it is." >&2
    echo false; return 0
  fi
  list="$(git ls-remote --tags --refs "$url" 'v*')" \
    || { err "could not list the tags of $url"; return 2; }
  plain="$(printf '%s\n' "$list" | sed 's|.*refs/tags/||' | grep -E "$PLAIN" || true)"
  if ! printf '%s\n' "$plain" | grep -qxF "$tag"; then
    err "$tag, the tag this run is for, is not in the tag list of $url: the lookup is broken"
    return 2
  fi
  newest="$(printf '%s\n' "$plain" | sort -V | tail -n 1)"
  echo "this tag: $tag; newest final release tag: $newest" >&2
  if [ "$newest" = "$tag" ]; then echo true; else echo false; fi
}

mode="${1:-}"
[ $# -gt 0 ] && shift
case "$mode" in
  guard)
    if [ $# -lt 2 ] || [ $# -gt 3 ]; then err "usage: guard <image> <git-tag> [<digests-dir>]"; exit 2; fi
    guard "$@" ;;
  is-newest)
    [ $# -eq 2 ] || { err "usage: is-newest <repo-url> <git-tag>"; exit 2; }
    is_newest "$@" ;;
  *)
    err "usage: tools/docker-release-tags.sh guard|is-newest ..."; exit 2 ;;
esac
