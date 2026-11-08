#!/bin/bash
# Regression test for tools/scan-image-secrets.sh: plant one fake secret per
# detection route in a throwaway image built on <base image>, and require the
# scanner to fail, to name every plant, to leave a decoy alone, and to print no
# secret value. Run by .github/workflows/docker.yml after the slim build.
#
#   tools/test-scan-image-secrets.sh [--with-attestations] <base image>
#
# --with-attestations also tests `scan-image-secrets.sh --attestations`
# against a throwaway local registry (registry:2 on 127.0.0.1:5055, and a
# host-network buildx builder that can push to it): a planted build-arg token
# in provenance is found, clean provenance passes, a missing provenance blob is
# "could not read" (exit 2, never "clean"), and an image pushed without
# attestations passes. Needs docker with buildx and python3.
#
# Fake values are generated per run (never committed), so nothing in this
# repository looks like a live token.
set -u
with_att=no
if [ "${1:-}" = "--with-attestations" ]; then with_att=yes; shift; fi
base="${1:?usage: tools/test-scan-image-secrets.sh [--with-attestations] <base image>}"
here="$(cd "$(dirname "$0")" && pwd)"
work="$(mktemp -d)"
cleanup() {
  rm -rf "$work"
  docker image rm -f ferx-scan-test:planted >/dev/null 2>&1
  docker rm -f ferx-scan-test-registry >/dev/null 2>&1
  docker buildx rm ferx-scan-test-builder >/dev/null 2>&1
}
trap cleanup EXIT

gh="ghp_$(LC_ALL=C tr -dc 'A-Za-z0-9' < /dev/urandom | head -c 36)"
aws="AKIA$(LC_ALL=C tr -dc 'A-Z0-9' < /dev/urandom | head -c 16)"
b64="$(LC_ALL=C tr -dc 'A-Za-z0-9+/' < /dev/urandom | head -c 64)"

cat > "$work/Dockerfile" <<EOF
FROM $base
ARG GHTOK
ARG AWSKEY
ARG B64
ENV GH_TOKEN=\${GHTOK}
USER root
RUN set -e; mkdir -p /opt/nonmem /srv /home/scan; \\
    printf 'my token %s here\\n' "\$GHTOK" > /home/scan/notes.txt; \\
    echo 'FAKE LICENSE' > /opt/nonmem/nonmem.lic; \\
    echo 'machine example.org login x password y' > /home/scan/.netrc; \\
    printf -- '-----BEGIN PRIVATE KEY-----\\n%s\\n-----END PRIVATE KEY-----\\n' "\$B64" > /srv/key.pem; \\
    printf '{"private_key": "-----BEGIN PRIVATE KEY-----\\\\n%s\\\\n-----END PRIVATE KEY-----\\\\n"}\\n' "\$B64" > /srv/sa.json; \\
    printf 'License Registered to: Test Org\\nExpiration Date: 1 JAN 2099\\n' > /opt/run.lst; \\
    printf 'aws %s\\n' "\$AWSKEY" > /srv/aws.txt; \\
    printf 'const HDR = "-----BEGIN PRIVATE KEY-----"; // a parser constant\\n' > /srv/pemcode.js
EOF

docker build -q -t ferx-scan-test:planted \
  --build-arg GHTOK="$gh" --build-arg AWSKEY="$aws" --build-arg B64="$b64" "$work" >/dev/null \
  || { echo "test-scan-image-secrets: FAIL could not build the planted image"; exit 1; }

out="$("$here/scan-image-secrets.sh" ferx-scan-test:planted 2>&1)"
code=$?
fail=0
check() { if printf '%s\n' "$out" | grep -qF -- "$1"; then echo "  ok   $2"; else echo "  FAIL $2 (expected: $1)"; fail=1; fi; }

echo "test-scan-image-secrets: planted image on $base"
[ "$code" -eq 1 ] && echo "  ok   scanner exits 1" || { echo "  FAIL scanner exited $code, expected 1"; fail=1; }
check "image config (ENV/LABEL)"          "token in ENV (config)"
check "build history"                     "token in ENV/ARG (history)"
check "token: /home/scan/notes.txt"       "GitHub token in a file"
check "token: /srv/aws.txt"               "AWS key in a file"
check "name: /opt/nonmem/nonmem.lic"      "license file by name"
check "name: /home/scan/.netrc"           ".netrc by name"
check "private-key: /srv/key.pem"         "PEM private key"
check "private-key: /srv/sa.json"         "private key embedded in JSON"
check "nonmem-banner: /opt/run.lst"       "NONMEM license banner"
if printf '%s\n' "$out" | grep -qF "pemcode.js"; then echo "  FAIL decoy (PEM header string in code) was flagged"; fail=1; else echo "  ok   decoy not flagged"; fi
if printf '%s\n' "$out" | grep -qF -e "${gh#ghp_}" -e "${aws#AKIA}" -e "$b64"; then
  echo "  FAIL a secret value appears in the scanner output"; fail=1
else
  echo "  ok   no secret value in the scanner output"
fi
if [ "$with_att" = yes ]; then
  echo "test-scan-image-secrets: attestations, against a throwaway registry"
  reg=127.0.0.1:5055
  docker run -d --rm --name ferx-scan-test-registry -p "$reg:5000" registry:2 >/dev/null \
    && docker buildx create --name ferx-scan-test-builder --driver docker-container \
         --driver-opt network=host >/dev/null \
    || { echo "  FAIL could not start the test registry/builder"; exit 1; }
  # push <name> <dockerfile dir> <provenance setting> [build args...] -> prints the pushed digest
  push() {
    name="$1"; dir="$2"; prov="$3"; shift 3
    docker buildx build --builder ferx-scan-test-builder --provenance="$prov" "$@" \
      --output "type=image,name=$reg/$name,push-by-digest=true,name-canonical=true,push=true,registry.insecure=true" \
      --metadata-file "$work/meta.json" -q "$dir" >/dev/null || return 1
    sed -n 's/.*"containerimage.digest": *"\([^"]*\)".*/\1/p' "$work/meta.json"
  }
  att_check() {  # att_check <expected exit> <label> <ref> [text that must not appear]
    out="$("$here/scan-image-secrets.sh" --attestations "$3" 2>&1)"; code=$?
    if [ "$code" -eq "$1" ]; then echo "  ok   $2 (exit $code)"; else echo "  FAIL $2: exit $code, expected $1"; printf '%s\n' "$out" | sed 's/^/       /'; fail=1; fi
    if [ -n "${4:-}" ] && printf '%s\n' "$out" | grep -qF -- "$4"; then echo "  FAIL $2: a secret value appears in the output"; fail=1; fi
  }
  mkdir -p "$work/att"
  printf 'FROM busybox:1.36\nARG GH_TOKEN\nRUN echo build\n' > "$work/att/Dockerfile"

  d="$(push planted "$work/att" mode=max --build-arg GH_TOKEN="$gh")" && [ -n "$d" ] \
    && att_check 1 "token passed as a build arg is found in provenance" "$reg/planted@$d" "${gh#ghp_}" \
    || { echo "  FAIL could not push the planted image"; fail=1; }

  d="$(push clean "$work/att" mode=max)" && [ -n "$d" ] \
    && att_check 0 "clean provenance passes" "$reg/clean@$d" \
    || { echo "  FAIL could not push the clean image"; fail=1; }

  # The defect this guards: buildx prints {} (exit 0) when the provenance blob
  # cannot be fetched. Delete that blob from the registry's storage.
  if [ -n "${d:-}" ]; then
    att="$(docker buildx imagetools inspect "$reg/clean@$d" --raw \
             | python3 -I -c 'import json,sys; print([m["digest"] for m in json.load(sys.stdin)["manifests"] if m.get("platform",{}).get("os")=="unknown"][0])')"
    blob="$(docker buildx imagetools inspect "$reg/clean@$att" --raw \
             | python3 -I -c 'import json,sys; print([l["digest"] for l in json.load(sys.stdin)["layers"] if "in-toto" in l["mediaType"]][0])')"
    h="${blob#sha256:}"
    docker exec ferx-scan-test-registry rm -rf "/var/lib/registry/docker/registry/v2/blobs/sha256/${h:0:2}/$h"
    att_check 2 "unreadable provenance is 'could not read', not clean" "$reg/clean@$d"
  fi

  d="$(push noatt "$work/att" false)" && [ -n "$d" ] \
    && att_check 0 "an image pushed without attestations passes" "$reg/noatt@$d" \
    || { echo "  FAIL could not push the image without attestations"; fail=1; }
fi

[ "$fail" -eq 0 ] && echo "test-scan-image-secrets: PASS" || echo "test-scan-image-secrets: FAIL"
exit "$fail"
