#!/bin/bash
# Regression test for tools/scan-image-secrets.sh: plant one fake secret per
# detection route in a throwaway image built on <base image>, and require the
# scanner to fail, to name every plant, to leave a decoy alone, and to print no
# secret value. Run by .github/workflows/docker.yml after the slim build.
#
#   tools/test-scan-image-secrets.sh <base image>
#
# Fake values are generated per run (never committed), so nothing in this
# repository looks like a live token.
set -u
base="${1:?usage: tools/test-scan-image-secrets.sh <base image>}"
here="$(cd "$(dirname "$0")" && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"; docker image rm -f ferx-scan-test:planted >/dev/null 2>&1' EXIT

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
[ "$fail" -eq 0 ] && echo "test-scan-image-secrets: PASS" || echo "test-scan-image-secrets: FAIL"
exit "$fail"
