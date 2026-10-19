#!/usr/bin/env bash
# Drive tools/changelog.sh over scratch trees, one per shape it is supposed to
# catch and one per shape it must not.
#
# The real checkout is always in shape, so nothing there would fail if one of
# the checks were deleted. Each fixture is a throwaway `git init` holding a
# NEWS.md, a changelog.d/ in good shape and one change, pointed at through
# CHANGELOG_ROOT; every bad fixture must exit 1 with its own message, so a check
# firing for the wrong reason does not pass. Output is captured and printed only
# for a failing case, indented.
#
# Usage: tools/test-changelog.sh   (from anywhere; exits 1 on a failure)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$ROOT/tools/changelog.sh"
# An explicit template, so TMPDIR is honoured on BSD mktemp (macOS) as on GNU.
TMP="$(mktemp -d "${TMPDIR:-/tmp}/test-changelog.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
# A fixture must stay its own repository even when TMPDIR sits inside a git
# checkout: git's discovery stops at $TMP instead of climbing into it.
export GIT_CEILING_DIRECTORIES="$TMP"
export GIT_AUTHOR_NAME=test GIT_AUTHOR_EMAIL=test@example.com
export GIT_COMMITTER_NAME=test GIT_COMMITTER_EMAIL=test@example.com
failures=0
cases=0

report() { # headline, captured output
  echo "$1"
  failures=$((failures + 1))
  printf '%s\n' "$2" | sed 's/^/     | /'
}

# A tree in shape: an empty development section above a release, the README,
# one fragment, DESCRIPTION at the next release's version, and a source file;
# committed on `main`. The caller then applies its one change.
fixture() { # name
  local dir="$TMP/$1"
  mkdir -p "$dir/changelog.d" "$dir/R" "$dir/tests/testthat"
  cat >"$dir/NEWS.md" <<'EOF'
# ferx 0.4.0.9000 (development version)

# ferx 0.4.0

## Bug fixes

- **An old fix** (#1).

  ```r
  # a comment inside a code block is no heading
  ```
EOF
  printf '# NEWS fragments\n' >"$dir/changelog.d/README.md"
  printf -- '- **Pending fix** (#10).\n' >"$dir/changelog.d/10.fix.md"
  printf 'Package: ferx\nVersion: 0.5.0\n' >"$dir/DESCRIPTION"
  printf 'f <- function() 1\n' >"$dir/R/f.R"
  printf 'test_that("f", expect_equal(f(), 1))\n' >"$dir/tests/testthat/test-f.R"
  git -C "$dir" init -q -b main
  git -C "$dir" add -A
  git -C "$dir" commit -q -m base
  echo "$dir"
}

# A branch off `main` in fixture $1, for `require`.
branch() { git -C "$1" checkout -q -b pr; }
commit() { git -C "$1" add -A && git -C "$1" commit -q -m change; }

run() { # dir, args... -> sets out, rc
  local dir=$1
  shift
  rc=0
  out=$(cd "$dir" && CHANGELOG_ROOT="$dir" bash "$SCRIPT" "$@" 2>&1) || rc=$?
}

expect() { # name, dir, expected exit, message fragment, args...
  local name=$1 dir=$2 want=$3 fragment=$4
  shift 4
  cases=$((cases + 1))
  run "$dir" "$@"
  if [[ "$rc" != "$want" ]]; then
    report "FAIL $name: expected exit $want, got $rc" "$out"
  elif ! printf '%s' "$out" | grep -qF -- "$fragment"; then
    report "FAIL $name: output does not mention \"$fragment\"" "$out"
  fi
}

# -- check ---------------------------------------------------------------------

d=$(fixture good)
expect "a tree in shape passes" "$d" 0 "1 fragment(s) OK" check

for good in 10-2.fix.md 11.breaking.md 12.feature.md 13.docs.md 14.internal.md; do
  d=$(fixture "ok-$good")
  printf -- '- **Entry** (#99).\n\n  continued, with a nested list:\n\n  - item\n' >"$d/changelog.d/$good"
  expect "accepts $good" "$d" 0 "2 fragment(s) OK" check
done

bad_name() { # file name, message
  local d
  d=$(fixture "name-$1")
  printf -- '- **Entry** (#99).\n' >"$d/changelog.d/$1"
  expect "rejects the name $1" "$d" 1 "$2" check
}
bad_name 0.fix.md "name must be <N>.<category>.md"
bad_name 10-1.fix.md "name must be <N>.<category>.md"
bad_name 12.Fix.md "name must be <N>.<category>.md"
bad_name fix.md "name must be <N>.<category>.md"
bad_name 10.fixed.md "unknown category 'fixed'"
bad_name 10.added.md "unknown category 'added'"

bad_body() { # name, body, message
  local d
  d=$(fixture "body-$1")
  printf '%b' "$2" >"$d/changelog.d/20.fix.md"
  expect "rejects $1" "$d" 1 "$3" check
}
bad_body empty '\n\n' "empty"
bad_body prose 'Fixed a thing (#20).\n' "line 1 must start a bullet"
bad_body two-bullets '- one (#20).\n- two.\n' "must be a single bullet"
bad_body one-space '- one (#20).\n - nested by one space.\n' "must be a single bullet"
bad_body heading '- one (#20).\n## Bug fixes\n' "must be a single bullet"
bad_body crlf '- one (#20).\r\n  more.\r\n' "CR line endings"
bad_body no-ref '- one, with no reference.\n' "no issue/PR reference"
bad_body zero-ref '- one (#0).\n' "no issue/PR reference"

d=$(fixture dir-fragment)
mkdir "$d/changelog.d/30.fix.md"
expect "rejects a directory named like a fragment" "$d" 1 "not a regular file" check

d=$(fixture symlink-fragment)
ln -s 10.fix.md "$d/changelog.d/31.fix.md"
expect "rejects a symlinked fragment" "$d" 1 "not a regular file" check

d=$(fixture no-readme)
rm "$d/changelog.d/README.md"
expect "requires changelog.d/README.md" "$d" 1 "README.md: missing" check

d=$(fixture no-dir)
rm -r "$d/changelog.d"
expect "requires changelog.d/" "$d" 1 "no changelog.d/ directory" check

d=$(fixture dev-entry)
awk 'NR == 1 { print; print ""; print "- **Written straight into NEWS.md** (#40)."; next } { print }' \
  "$d/NEWS.md" >"$d/NEWS.tmp" && mv "$d/NEWS.tmp" "$d/NEWS.md"
expect "rejects an entry under the development heading" "$d" 1 \
  "development section must stay empty" check

d=$(fixture two-devs)
printf '# ferx 0.4.0.9001 (development version)\n\n' | cat - "$d/NEWS.md" >"$d/NEWS.tmp" &&
  mv "$d/NEWS.tmp" "$d/NEWS.md"
expect "rejects two development headings" "$d" 1 "at most one" check

d=$(fixture dev-below-release)
printf '# ferx 0.4.0\n\n# ferx 0.4.0.9000 (development version)\n' >"$d/NEWS.md"
expect "rejects a development heading below a release" "$d" 1 \
  "must come before every release heading" check

d=$(fixture no-dev)
sed -i.bak '1,2d' "$d/NEWS.md" && rm "$d/NEWS.md.bak"
expect "accepts NEWS.md without a development heading" "$d" 0 "fragment(s) OK" check

# -- preview -------------------------------------------------------------------

d=$(fixture preview)
printf -- '- **A feature** (#50).\n' >"$d/changelog.d/50.feature.md"
printf -- '- **A newer fix** (#60).\n' >"$d/changelog.d/60.fix.md"
printf -- '- **Second on 60** (#60).\n' >"$d/changelog.d/60-2.fix.md"
printf -- '- **Tenth on 60** (#60).\n' >"$d/changelog.d/60-10.fix.md"
cases=$((cases + 1))
run "$d" preview
want='# ferx (development version)

## New features

- **A feature** (#50).

## Bug fixes

- **A newer fix** (#60).

- **Second on 60** (#60).

- **Tenth on 60** (#60).

- **Pending fix** (#10).'
[[ "$rc" == 0 && "$out" == "$want" ]] ||
  report "FAIL preview: sections in category order, newest first" "$out"

# -- assemble ------------------------------------------------------------------

d=$(fixture assemble)
printf -- '- **A breaking change** (#70).\n\n  With a second paragraph.\n' >"$d/changelog.d/70.breaking.md"
cases=$((cases + 1))
run "$d" assemble 0.5.0
want='# ferx 0.5.0

## Breaking changes

- **A breaking change** (#70).

  With a second paragraph.

## Bug fixes

- **Pending fix** (#10).

# ferx 0.4.0

## Bug fixes

- **An old fix** (#1).

  ```r
  # a comment inside a code block is no heading
  ```'
if [[ "$rc" != 0 ]]; then
  report "FAIL assemble: expected exit 0, got $rc" "$out"
elif [[ "$(cat "$d/NEWS.md")" != "$want" ]]; then
  report "FAIL assemble: the development heading becomes the release section" "$(cat "$d/NEWS.md")"
elif [[ -n "$(find "$d/changelog.d" -name '*.md' ! -name README.md)" ]]; then
  report "FAIL assemble: fragments left behind" "$(ls "$d/changelog.d")"
elif [[ ! -f "$d/changelog.d/README.md" ]]; then
  report "FAIL assemble: removed the README" "$(ls "$d/changelog.d")"
fi
expect "assembled tree passes check" "$d" 0 "0 fragment(s) OK" check

d=$(fixture assemble-no-dev)
sed -i.bak '1,2d' "$d/NEWS.md" && rm "$d/NEWS.md.bak"
cases=$((cases + 1))
run "$d" assemble 0.5.0
[[ "$rc" == 0 && "$(head -n 1 "$d/NEWS.md")" == "# ferx 0.5.0" &&
  "$(grep -c '^# ferx 0.4.0$' "$d/NEWS.md")" == 1 ]] ||
  report "FAIL assemble without a development heading: goes above the last release" "$out"

refuse_assemble() { # name, message, version, setup command...
  local name=$1 msg=$2 version=$3 d before frags
  shift 3
  d=$(fixture "refuse-$name")
  (cd "$d" && "$@")
  before="$(cat "$d/NEWS.md")"
  frags="$(ls -a "$d/changelog.d")"
  expect "assemble refuses $name" "$d" 1 "$msg" assemble "$version"
  [[ "$(cat "$d/NEWS.md")" == "$before" && "$(ls -a "$d/changelog.d")" == "$frags" ]] ||
    report "FAIL assemble refusing $name changed NEWS.md or the fragments" ""
}
refuse_assemble bad-version "not a release version" 0.5 true
refuse_assemble dev-version "not a release version" 0.5.0.9000 true
refuse_assemble leading-zero "not a release version" 0.05.0 true
refuse_assemble description-mismatch "DESCRIPTION says Version: 0.5.0" 0.6.0 true
refuse_assemble no-fragments "no fragments" 0.5.0 rm changelog.d/10.fix.md
refuse_assemble existing-version "already has a '# ferx 0.5.0' section" 0.5.0 \
  sed -i.bak '1s/.*/# ferx 0.5.0/' NEWS.md
refuse_assemble invalid-fragment "problem(s) found" 0.5.0 \
  sh -c 'printf "no bullet\n" > changelog.d/11.fix.md'
expect "assemble takes a leading v" "$(fixture v-prefix)" 0 "assembled 1 fragment(s)" assemble v0.5.0

# -- require -------------------------------------------------------------------

d=$(fixture req-missing)
branch "$d"
printf 'f <- function() 2\n' >"$d/R/f.R"
commit "$d"
expect "require fails an R/ change without a fragment" "$d" 1 \
  "Missing: changelog.d/77.<category>.md" require main --pr 77
expect "require names the user-facing path" "$d" 1 "  R/f.R" require main --pr 77
expect "require honours the opt-out" "$d" 0 "label set" require main --pr 77 --opt-out true
expect "require rejects a non-boolean opt-out" "$d" 1 "must be true or false" \
  require main --opt-out yes

for path in src/rust/Cargo.lock inst/examples/x.ferx NAMESPACE DESCRIPTION; do
  d=$(fixture "req-path-${path//\//-}")
  branch "$d"
  mkdir -p "$d/$(dirname "$path")"
  printf 'changed\n' >>"$d/$path"
  commit "$d"
  expect "require counts $path as user-facing" "$d" 1 "Missing:" require main
done

d=$(fixture req-fragment)
branch "$d"
printf 'f <- function() 2\n' >"$d/R/f.R"
printf -- '- **f returns 2** (#77).\n' >"$d/changelog.d/77.fix.md"
commit "$d"
expect "require passes an R/ change with a fragment" "$d" 0 "changelog.d/77.fix.md" require main

d=$(fixture req-reword-only)
branch "$d"
printf 'f <- function() 2\n' >"$d/R/f.R"
printf -- '- **Pending fix, reworded** (#10).\n' >"$d/changelog.d/10.fix.md"
commit "$d"
expect "rewording another fragment does not count as adding one" "$d" 1 "Missing:" require main

d=$(fixture req-tests-only)
branch "$d"
printf 'test_that("f", expect_true(TRUE))\n' >>"$d/tests/testthat/test-f.R"
mkdir -p "$d/tools" && printf 'echo\n' >"$d/tools/x.sh"
commit "$d"
expect "require passes tests/ and tools/ changes" "$d" 0 "no user-facing paths" require main

d=$(fixture req-delete)
branch "$d"
git -C "$d" rm -q changelog.d/10.fix.md
commit "$d"
expect "require fails deleting another PR's fragment" "$d" 1 "deletes pending" require main
expect "the opt-out does not excuse deleting one" "$d" 1 "deletes pending" \
  require main --opt-out true

d=$(fixture req-rename)
branch "$d"
git -C "$d" mv changelog.d/10.fix.md changelog.d/11.fix.md
commit "$d"
expect "require fails renaming another PR's fragment" "$d" 1 "renames pending" require main

d=$(fixture req-release)
branch "$d"
(cd "$d" && CHANGELOG_ROOT="$d" bash "$SCRIPT" assemble 0.5.0 >/dev/null)
commit "$d"
expect "require passes a release PR that assembled every fragment" "$d" 0 "release PR" require main

d=$(fixture req-release-leftover)
branch "$d"
(cd "$d" && CHANGELOG_ROOT="$d" bash "$SCRIPT" assemble 0.5.0 >/dev/null)
printf -- '- **Late** (#88).\n' >"$d/changelog.d/88.fix.md"
commit "$d"
expect "require fails a release PR that leaves a fragment pending" "$d" 1 \
  "leaves pending ones" require main

d=$(fixture req-release-dropped)
branch "$d"
(cd "$d" && CHANGELOG_ROOT="$d" bash "$SCRIPT" assemble 0.5.0 >/dev/null)
sed -i.bak '/Pending fix/d' "$d/NEWS.md" && rm "$d/NEWS.md.bak"
commit "$d"
expect "require fails a release PR that drops an assembled entry" "$d" 1 \
  "deletes pending" require main

if [[ "$failures" -gt 0 ]]; then
  echo "test-changelog: $failures of $cases case(s) failed"
  exit 1
fi
echo "test-changelog: all $cases case(s) passed"
