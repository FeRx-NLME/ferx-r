#!/bin/sh
# Scan a built Docker image for credentials and license files that must never
# ship: GitHub/AWS/Slack tokens, private keys, .netrc/.git-credentials/.env
# files, license files (*.lic: NONMEM's nonmem.lic and friends) and NONMEM
# license banners. Exits 1 with every hit listed, 0 when clean.
#
#   tools/scan-image-secrets.sh <image>
#   tools/scan-image-secrets.sh --attestations <registry ref>
#
# .github/workflows/docker.yml runs it on every image it builds (both
# Dockerfiles) before anything is pushed, and again on the exact image it
# pushed (pulled back by digest) before any tag points at it. An image scan
# looks at three places:
#   1. the image config: environment variables and labels;
#   2. the build history (every RUN/ENV/ARG line baked into the layers);
#   3. the filesystem, as root, outside /proc, /sys and /dev.
# The CA bundles every base image carries (/etc/ssl/certs,
# /usr/share/ca-certificates) are public certificates, not keys: they are not
# flagged, because the key check reads file content (PRIVATE KEY), not names.
#
# --attestations scans what a registry publishes next to the image: the
# provenance (build args, the Dockerfile, build metadata) and SBOM
# attestations of a pushed index. `docker pull` never fetches those, so an
# image scan alone would not see them. Needs buildx and registry access.
set -eu

mode=image
if [ "${1:-}" = "--attestations" ]; then mode=attestations; shift; fi
image="${1:?usage: tools/scan-image-secrets.sh [--attestations] <image>}"
hits=0
report() { echo "SECRET-SCAN: $*" >&2; hits=1; }

# Token shapes, as whole words: a token is never glued to more base64, which
# keeps random runs inside embedded WebAssembly/base64 blobs from matching.
TOKEN_RE='\b(gh[pousr]_[A-Za-z0-9]{36}|github_pat_[A-Za-z0-9_]{40,}|AKIA[0-9A-Z]{16}|xox[baprs]-[A-Za-z0-9-]{10,})\b'
# A PEM private key header. A file is only reported when a base64 body line
# follows it (checked in the image), so PEM-parsing code that merely contains
# the header string does not match.
PEM_RE='-----BEGIN [A-Z ]*PRIVATE KEY-----'
# Assignments of a non-empty secret-looking variable in config or history.
ASSIGN_RE='(GITHUB_TOKEN|GH_TOKEN|GITHUB_PAT|NPM_TOKEN|AWS_SECRET_ACCESS_KEY|AWS_SESSION_TOKEN|PASSWORD|PASSWD|SECRET|API_KEY|NONMEM_LICENSE|NMLICENSE)=[^[:space:]"'"'"']'
# What a report prints: token values and assigned values are redacted, so the
# scanner never copies a live secret into a CI log.
redact() {
  sed -E 's/(gh[pousr]_|github_pat_|AKIA|xox[baprs]-)[A-Za-z0-9_-]+/\1<redacted>/g;
          s/((TOKEN|PAT|KEY|PASSWORD|PASSWD|SECRET|LICENSE)=)[^[:space:]"'"'"']+/\1<redacted>/g'
}

if [ "$mode" = attestations ]; then
  # A secret-named key with a non-empty value, as JSON (provenance records
  # build args as "build-arg:NAME": "value", not NAME=value).
  JSON_ASSIGN_RE='"(build-arg:)?(GITHUB_TOKEN|GH_TOKEN|GITHUB_PAT|NPM_TOKEN|AWS_SECRET_ACCESS_KEY|AWS_SESSION_TOKEN|PASSWORD|PASSWD|SECRET|API_KEY|NONMEM_LICENSE|NMLICENSE)"[[:space:]]*:[[:space:]]*"[^"]'
  att="$(docker buildx imagetools inspect "$image" --format '{{json .Provenance}}{{println}}{{json .SBOM}}')" \
    || { echo "SECRET-SCAN: could not read the attestations of $image" >&2; exit 2; }
  if printf '%s\n' "$att" | grep -Eq "$TOKEN_RE|$PEM_RE|$ASSIGN_RE|$JSON_ASSIGN_RE"; then
    report "attestations (provenance/SBOM) of $image:"
    printf '%s\n' "$att" | grep -Eo ".{0,40}($TOKEN_RE|$PEM_RE|$ASSIGN_RE|$JSON_ASSIGN_RE).{0,20}" \
      | redact | sed -E 's/("(build-arg:)?[A-Z_]+"[[:space:]]*:[[:space:]]*")[^"]+/\1<redacted>/g; s/^/  /' >&2
  fi
  [ "$hits" -eq 0 ] && echo "secret scan: attestations of $image are clean (provenance, SBOM)"
  exit "$hits"
fi

# 1. Config: environment and labels.
config="$(docker image inspect --format '{{range .Config.Env}}{{println .}}{{end}}{{range $k, $v := .Config.Labels}}{{$k}}={{$v}}{{println}}{{end}}' "$image")"
if printf '%s\n' "$config" | grep -Eq "$TOKEN_RE|$PEM_RE|$ASSIGN_RE"; then
  report "image config (ENV/LABEL) of $image:"
  printf '%s\n' "$config" | grep -E "$TOKEN_RE|$PEM_RE|$ASSIGN_RE" | redact | sed 's/^/  /' >&2
fi

# 2. Build history.
history="$(docker image history --no-trunc --format '{{.CreatedBy}}' "$image")"
if printf '%s\n' "$history" | grep -Eq "$TOKEN_RE|$PEM_RE|$ASSIGN_RE"; then
  report "build history of $image:"
  printf '%s\n' "$history" | grep -Eo ".{0,40}($TOKEN_RE|$PEM_RE|$ASSIGN_RE).{0,20}" | redact | sed 's/^/  /' >&2
fi

# 3. Filesystem, inside the image as root (the image's default user cannot
#    read everything). The scan program is passed on stdin, so nothing is
#    mounted and the image's own tools (find, grep from the base) do the work.
fs_hits="$(docker run --rm -i --user 0 --entrypoint sh --network none "$image" -s -- "$TOKEN_RE" "$PEM_RE" <<'SCAN'
TOKEN_RE="$1"
PEM_RE="$2"
prune='-path /proc -prune -o -path /sys -prune -o -path /dev -prune -o'
# Files that are secrets by name, wherever they are.
eval "find / -xdev $prune -type f \\( \
  -name '*.lic' -o -name 'nmlicense*' -o -name '.netrc' -o -name '.git-credentials' \
  -o -name '.env' -o -name '.env.*' -o -name '.npmrc' -o -name '*.p12' -o -name '*.pfx' \
  -o -name 'id_rsa' -o -name 'id_ecdsa' -o -name 'id_ed25519' \\) -print" 2>/dev/null \
  | sed 's/^/name: /'
# Token or private-key content in text files, outside the public CA bundles.
files() {
  eval "find / -xdev $prune -path /etc/ssl/certs -prune -o -path /usr/share/ca-certificates -prune -o \
    -type f -size -5M -print0" 2>/dev/null
}
files | xargs -0 grep -IlE "$TOKEN_RE" 2>/dev/null | sed 's/^/token: /'
# A private key only when a base64 body follows the header (after optional
# Proc-Type/DEK-Info headers): PEM-handling code that just names the header
# is not a key. Two layouts: real newlines (a .pem/.key file, docs) and the
# escaped `\n` of a key embedded in JSON (cloud service-account files).
files | xargs -0 grep -IlE -- "$PEM_RE" 2>/dev/null \
  | while IFS= read -r f; do
      if grep -IPzq -- "$PEM_RE\\r?\\n(?:[A-Za-z-]+:[^\\n]*\\n)*\\r?\\n?[A-Za-z0-9+/]{40,}" "$f" 2>/dev/null \
         || grep -IEq -- "$PEM_RE"'(\\r)?\\n[A-Za-z0-9+/]{40,}' "$f" 2>/dev/null; then
        echo "private-key: $f"
      fi
    done
# A NONMEM license banner (what nmfe prints into .lst files).
for d in /usr/local /opt /home /root /work /srv; do
  [ -d "$d" ] && grep -rIl 'License Registered to' "$d" 2>/dev/null | sed 's/^/nonmem-banner: /'
done
exit 0
SCAN
)"
# Known public third-party files that look like secrets but are not. Each
# entry is a whole hit line (kind: path glob) with its reason; anything else
# fails. Only add an entry after reading the file, and only for a file that is
# public upstream (a distro package or a CRAN package), never for anything
# from this repository.
known_public() {
  case "$1" in
    # Ubuntu cmake-data: CMake's dummy Windows app-signing template.
    # (Globs stay unquoted - a quoted * is literal in `case` - so only the
    # space after the kind is escaped.)
    name:\ /usr/share/cmake-*/Templates/Windows/Windows_TemporaryKey.pfx) return 0 ;;
    # CRAN openssl: example keys in the package's own documentation.
    private-key:\ /usr/local/lib/R/site-library/openssl/doc/keys.html) return 0 ;;
    # CRAN xpose (a ferx Suggests; RStudio image only): example NONMEM output
    # whose banner names Uppsala University's 2018 license. A registration
    # line, not a license key; public on CRAN.
    nonmem-banner:\ /usr/local/lib/R/site-library/xpose/extdata/*.lst) return 0 ;;
  esac
  return 1
}
fs_real=""
fs_known=""
if [ -n "$fs_hits" ]; then
  while IFS= read -r line; do
    if known_public "$line"; then fs_known="$fs_known  $line
"; else fs_real="$fs_real  $line
"; fi
  done <<EOF
$fs_hits
EOF
fi
if [ -n "$fs_known" ]; then
  echo "secret scan: ignored known public third-party files in $image:"
  printf '%s' "$fs_known"
fi
if [ -n "$fs_real" ]; then
  report "filesystem of $image:"
  printf '%s' "$fs_real" >&2
fi

if [ "$hits" -eq 0 ]; then
  echo "secret scan: $image is clean (config, history, filesystem)"
fi
exit "$hits"
