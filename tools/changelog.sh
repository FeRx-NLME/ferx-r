#!/usr/bin/env bash
#
# NEWS fragments (#600, after ferx-core #1545). `--help` prints usage; this
# header is the why.
#
# Every user-facing PR used to add a bullet under the development heading of
# NEWS.md, at the same insertion point. With several PRs open and green at once,
# merging the first made every other one conflict on NEWS.md, so each needed a
# local rebase, a conflict resolution and a fresh CI run before it could merge -
# a new merge and a new CI run for nothing but a changelog line. That is
# structural: the convention forced N concurrent PRs to edit the same lines of
# the same file. `merge=union` does not fix it - GitHub's merge button and
# conflict detection ignore custom merge drivers.
#
# So each PR adds a NEW file instead, `changelog.d/<N>.<category>.md`, and the
# release step assembles them into NEWS.md (the towncrier / scriv / changie
# pattern, and the same convention as ferx-core's changelog.d/). Two PRs adding
# two different files cannot conflict.
#
# Plain bash + awk + sed on purpose: no new dependency, and it runs unchanged
# under macOS's bash 3.2 and CI's bash 5.
#
# `CHANGELOG_ROOT` (default: the repo root) points every subcommand at another
# tree; tools/test-changelog.sh uses it to run against scratch copies instead of
# the real NEWS.md.
set -euo pipefail

script_path="${BASH_SOURCE[0]}"
repo_root="$(cd "$(dirname "$script_path")/.." && pwd)"
root="${CHANGELOG_ROOT:-$repo_root}"
frag_dir="$root/changelog.d"
news="$root/NEWS.md"
description="$root/DESCRIPTION"

# The `##` sections NEWS.md has used since 0.4.0, in the order a release lists
# them. `heading_of` maps each to its heading.
CATEGORIES=(breaking feature fix docs internal)

# Files in changelog.d/ that are not fragments.
NON_FRAGMENTS=(README.md .gitkeep)

# Paths whose change is presumed user-facing: the R API, the Rust glue and its
# pinned engine (src/rust/Cargo.lock), the shipped examples, and the package
# metadata. tests/, tools/, .github/ and the docs are not.
USER_FACING_PREFIXES=(R/ src/ inst/ NAMESPACE DESCRIPTION)

OPT_OUT_LABEL="no-changelog"

# A release heading (`# ferx 0.4.0`) and the development one
# (`# ferx 0.4.0.9000 (development version)`). Anything else starting `# ` -
# a comment inside a fenced code block of an old entry - is not a heading.
RELEASE_RE='^# ferx [0-9]'
DEV_RE='^# ferx [0-9][0-9.]* [(]development version[)][[:space:]]*$'

# An awk function, prepended to the scripts that need it: `fence_step(line)`
# returns 1 if the line opens or closes a fenced code block, and keeps
# `fence_n` > 0 while inside one. CommonMark's rule: a fence is 3+ backticks or
# tildes, and only a run of the SAME character, at least as long, with nothing
# after it, closes it - so a ``` line inside a ~~~ block is content.
# shellcheck disable=SC2016 # awk source, not shell
FENCE_AWK='
function fence_step(line,   s, c, n, rest) {
  s = line
  sub(/^[ \t]*/, "", s)
  c = substr(s, 1, 1)
  if (c != "`" && c != "~") return 0
  n = 0
  while (substr(s, n + 1, 1) == c) n++
  if (n < 3) return 0
  rest = substr(s, n + 1)
  if (!fence_n) { fence_c = c; fence_n = n; return 1 }
  if (c == fence_c && n >= fence_n && rest ~ /^[ \t]*$/) { fence_n = 0; return 1 }
  return 0
}
'

# Every `# ferx ...` heading of NEWS.md ($1, default the real one) as
# `<line>\t<text>`, skipping fenced code blocks: an old entry's example may
# well show a NEWS heading, and counting it would misplace `assemble`'s splice.
# Every consumer below locates headings through this one function.
news_headings() {
  awk -v rel="$RELEASE_RE" "$FENCE_AWK"'
    fence_step($0)        { next }
    !fence_n && $0 ~ rel  { print FNR "\t" $0 }
  ' "${1:-$news}"
}

# The line of the development heading, or nothing.
dev_line() {
  news_headings "$@" | awk -F '\t' -v dev="$DEV_RE" '$2 ~ dev { print $1; exit }'
}

# The line of the first heading after line $1 (0: the first heading), or nothing.
heading_after() {
  local after=$1
  shift
  news_headings "$@" | awk -F '\t' -v after="$after" '$1 + 0 > after + 0 { print $1; exit }'
}

usage() {
  cat <<EOF
NEWS fragments (#600): each user-facing PR adds a file under changelog.d/
instead of editing NEWS.md, so concurrent PRs never conflict.

  tools/changelog.sh check                     validate every fragment, and that
                                               NEWS.md's development section is empty
  tools/changelog.sh preview                   print the section the fragments
                                               would assemble into
  tools/changelog.sh assemble <version>        release time only: write the
                                               fragments into NEWS.md as
                                               '# ferx <version>' (replacing the
                                               development heading) and delete them
  tools/changelog.sh require <base-ref> [--pr N] [--opt-out true|false]
                                               fail if the diff base...HEAD touches
                                               user-facing code and adds no fragment,
                                               unless --opt-out true (CI passes
                                               whether the '$OPT_OUT_LABEL' label is set)

Fragment: changelog.d/<N>.<category>.md, where <N> is the issue or PR number
(append -2, -3, ... for a second entry on the same number, e.g. 541-2.fix.md)
and <category> is one of:
  ${CATEGORIES[*]}
Its body is exactly one Markdown bullet, as it will appear in NEWS.md: the
first line starts with '- ', every further line is blank or indented.
EOF
}

die() {
  echo "changelog: $*" >&2
  exit 1
}

is_category() {
  local c
  for c in "${CATEGORIES[@]}"; do
    [ "$c" = "$1" ] && return 0
  done
  return 1
}

is_non_fragment() {
  local n
  for n in "${NON_FRAGMENTS[@]}"; do
    [ "$n" = "$1" ] && return 0
  done
  return 1
}

heading_of() {
  case "$1" in
    breaking) echo "Breaking changes" ;;
    feature) echo "New features" ;;
    fix) echo "Bug fixes" ;;
    docs) echo "Documentation" ;;
    internal) echo "Internal" ;;
    *) die "internal error: no heading for category '$1'" ;;
  esac
}

# Every fragment file name, newest first as NEWS.md lists its entries: by the
# leading number descending, then a number's own entries in order (541, 541-2,
# 541-10). A name without that shape sorts last, for `check` to reject.
fragment_names() {
  [ -d "$frag_dir" ] || return 0
  # `find`, not a glob: every glob spelling misses some dotfile shape (`..x`),
  # and a file that `check` never sees but `require` counts is a hole in both.
  #
  # A reserved name (README.md, .gitkeep) is skipped only when it is a plain
  # file: a `.gitkeep/` DIRECTORY is listed, so `check` rejects it rather than
  # letting `require` count a path nested inside it.
  local name
  find "$frag_dir" -mindepth 1 -maxdepth 1 | sed 's#.*/##' | while IFS= read -r name; do
    if is_non_fragment "$name" && [ -f "$frag_dir/$name" ] && [ ! -L "$frag_dir/$name" ]; then
      continue
    fi
    printf '%s\n' "$name"
  done | awk '{
    n = -1; k = 1
    if (match($0, /^[0-9]+/)) {
      n = substr($0, 1, RLENGTH) + 0
      rest = substr($0, RLENGTH + 1)
      if (match(rest, /^-[0-9]+/)) k = substr(rest, 2, RLENGTH - 1) + 0
    }
    printf "%d\t%d\t%s\n", n, k, $0
  }' | LC_ALL=C sort -t "$(printf '\t')" -k1,1nr -k2,2n -k3,3 | cut -f3-
}

# The fragment manifest, read ONCE per run into FRAGMENTS. Every consumer -
# check, render, verify, delete - walks this same list, so a file appearing
# mid-assemble is neither rendered nor deleted.
#
# Split on newlines with globbing off, not a `while read ... <<EOF $(...)` loop:
# a here-document needs a temp file, and when bash cannot create one the loop
# silently runs zero times - `check` then reports every fragment OK having read
# none, and exits 0.
FRAGMENTS=()
load_fragments() {
  local names
  names="$(fragment_names)" || die "could not list changelog.d/"
  FRAGMENTS=()
  [ -n "$names" ] || return 0
  local IFS=$'\n'
  set -f
  # shellcheck disable=SC2206 # word-splitting on newlines is the point
  FRAGMENTS=($names)
  set +f
}

# Fragment body with trailing blank lines removed and a final newline ensured.
fragment_body() {
  awk '{ lines[NR] = $0 } NF { last = NR }
       END { for (i = 1; i <= last; i++) print lines[i] }' "$1"
}

# The non-blank lines of NEWS.md's development section (excluding its heading),
# each prefixed with its line number.
dev_content() {
  local from to
  from="$(dev_line)"
  [ -n "$from" ] || return 0
  to="$(heading_after "$from")"
  awk -v from="$from" -v to="${to:-0}" '
    to + 0 > 0 && FNR >= to + 0 { exit }
    FNR > from + 0 && NF        { print FNR ": " $0 }
  ' "$news"
}

cmd_check() {
  local errors=0 name stem cat body first bad devs
  local ref_re='#[1-9][0-9]*'
  [ -f "$news" ] || die "no NEWS.md at $news"
  # The directory and its README are part of the convention, not optional: an
  # accidental `git rm changelog.d/*` after assembly would otherwise delete the
  # guide and, with it, the tracked directory.
  [ -d "$frag_dir" ] && [ ! -L "$frag_dir" ] ||
    die "no changelog.d/ directory at $frag_dir - restore it (with its README.md)"
  if [ ! -f "$frag_dir/README.md" ] || [ -L "$frag_dir/README.md" ]; then
    echo "changelog.d/README.md: missing or not a regular file - restore it" >&2
    errors=$((errors + 1))
  fi

  # At most one development heading, and if there is one it is the file's first
  # `# ferx` heading. Between `assemble` and the post-release version bump there
  # is none, which is fine; two, or one below a release, means a hand edit that
  # `assemble` would splice into the wrong place.
  devs="$(news_headings | cut -f2- | grep -Ec "$DEV_RE" || true)"
  if [ "$devs" -gt 1 ]; then
    echo "NEWS.md: expected at most one '# ferx X.Y.Z.9000 (development version)' heading, found $devs" >&2
    errors=$((errors + 1))
  elif [ "$devs" -eq 1 ] && [ "$(dev_line)" != "$(heading_after 0)" ]; then
    echo "NEWS.md: the development heading must come before every release heading" >&2
    errors=$((errors + 1))
  fi

  load_fragments
  for name in ${FRAGMENTS[@]+"${FRAGMENTS[@]}"}; do
    local path="$frag_dir/$name"
    if [ -L "$path" ] || [ ! -f "$path" ]; then
      echo "changelog.d/$name: not a regular file; fragments are plain files directly" \
        "under changelog.d/ (no directories, no symlinks)" >&2
      errors=$((errors + 1))
      continue
    fi
    # <N>[-<k>].<category>.md, k >= 2: the second entry on one number is `-2`.
    if ! printf '%s' "$name" | grep -Eq '^[1-9][0-9]*(-([2-9]|[1-9][0-9]+))?\.[a-z]+\.md$'; then
      echo "changelog.d/$name: name must be <N>.<category>.md (e.g. 541.fix.md;" \
        "541-2.fix.md for a second entry)" >&2
      errors=$((errors + 1))
      continue
    fi
    stem="${name%.md}"
    cat="${stem##*.}"
    if ! is_category "$cat"; then
      echo "changelog.d/$name: unknown category '$cat' - expected one of: ${CATEGORIES[*]}" >&2
      errors=$((errors + 1))
      continue
    fi
    body="$(fragment_body "$path")"
    if [ -z "$(printf '%s' "$body" | tr -d '[:space:]')" ]; then
      echo "changelog.d/$name: empty - write the bullet as it should read in NEWS.md" >&2
      errors=$((errors + 1))
      continue
    fi
    # No reader of $body below stops early: under pipefail, one that does (head,
    # grep -q, awk's exit) can close the pipe while printf is still writing a
    # long fragment, and the SIGPIPE fails the check on a fragment in shape.
    first="${body%%$'\n'*}"
    if ! printf '%s' "$first" | grep -Eq '^- [^[:space:]]'; then
      echo "changelog.d/$name: line 1 must start a bullet ('- ...'), found: $first" >&2
      errors=$((errors + 1))
      continue
    fi
    # Every later non-blank line must be indented by at least two spaces (or a
    # tab) - the content column of `- `, so Markdown keeps it inside the one
    # bullet. A column-0 line is a second bullet, a heading, or prose outside
    # the list; and ` - x` with ONE space is still a sibling bullet, not a
    # nested one.
    bad="$(printf '%s\n' "$body" | awk '!found && NR > 1 && NF && !/^(  |\t)/ { print NR ": " $0; found = 1 }')"
    if [ -n "$bad" ]; then
      echo "changelog.d/$name: must be a single bullet; line $bad is not indented" \
        "by two spaces (continuation lines need at least two; one entry per file)" >&2
      errors=$((errors + 1))
      continue
    fi
    if [[ $body == *$'\r'* ]]; then
      echo "changelog.d/$name: contains CR line endings; save it with LF" >&2
      errors=$((errors + 1))
      continue
    fi
    # A real number: `#0` / `#00` is a placeholder left in, not a reference.
    if ! [[ $body =~ $ref_re ]]; then
      echo "changelog.d/$name: no issue/PR reference - cite it as (#NN)" >&2
      errors=$((errors + 1))
      continue
    fi
  done

  local stray
  stray="$(dev_content)"
  if [ -n "$stray" ]; then
    {
      echo "NEWS.md: the development section must stay empty between releases (#600) -"
      echo "it is assembled from changelog.d/ at release time. Move each entry into its"
      echo "own changelog.d/<N>.<category>.md. Offending lines:"
      printf '%s\n' "$stray" | sed -n '1,5s/^/  /p'
    } >&2
    errors=$((errors + 1))
  fi

  if [ "$errors" -gt 0 ]; then
    echo "changelog: $errors problem(s) found" >&2
    exit 1
  fi
  echo "changelog: ${#FRAGMENTS[@]} fragment(s) OK"
}

# The `## Category` blocks for every fragment, in category order, each fragment
# followed by a blank line as NEWS.md spaces its entries.
render_sections() {
  local cat name first=1 matched
  for cat in "${CATEGORIES[@]}"; do
    matched=()
    for name in ${FRAGMENTS[@]+"${FRAGMENTS[@]}"}; do
      case "$name" in *."$cat".md) matched+=("$name") ;; esac
    done
    [ ${#matched[@]} -gt 0 ] || continue
    [ "$first" -eq 1 ] || echo
    first=0
    echo "## $(heading_of "$cat")"
    for name in "${matched[@]}"; do
      echo
      fragment_body "$frag_dir/$name"
    done
  done
}

cmd_preview() {
  cmd_check >/dev/null # also loads FRAGMENTS
  echo "# ferx (development version)"
  echo
  render_sections
}

# 0 if $1 is a release version: X.Y.Z, no leading zeros - not a development
# X.Y.Z.9000. Single-line input only (the caller checks); grep matches per line.
is_release_version() {
  local num='(0|[1-9][0-9]*)'
  printf '%s' "$1" | grep -Eq "^$num\\.$num\\.$num\$"
}

cmd_assemble() {
  [ $# -eq 1 ] || die "assemble: usage: tools/changelog.sh assemble <version>"
  local version="${1#v}"
  # grep matches per LINE, so `$'junk\n1.2.3'` would pass the regex below and
  # write a split heading. Refuse any line break first.
  case "$version" in
    *$'\n'* | *$'\r'*) die "assemble: the version must be single-line" ;;
  esac
  is_release_version "$version" ||
    die "assemble: '$version' is not a release version (X.Y.Z, no leading zeros)"

  # The release PR sets DESCRIPTION's Version and assembles NEWS.md together; a
  # mismatch is a typo in one of the two.
  if [ -f "$description" ]; then
    local described
    described="$(sed -n 's/^Version:[[:space:]]*//p' "$description" | tr -d '[:space:]')"
    [ "$described" = "$version" ] ||
      die "assemble: DESCRIPTION says Version: $described - set it to $version first"
  fi

  cmd_check >/dev/null # also loads FRAGMENTS: the one manifest used below
  [ ${#FRAGMENTS[@]} -gt 0 ] || die "assemble: no fragments in changelog.d/ - nothing to release"
  if news_headings | cut -f2- | grep -Fqx "# ferx $version"; then
    die "assemble: NEWS.md already has a '# ferx $version' section"
  fi

  # The rewrite goes to a copy beside NEWS.md (`cp -p`, so the `mv` keeps its
  # mode - a `mktemp` file is 0600) and replaces it only once verified.
  # Globals, not locals: the EXIT trap runs after this function has returned.
  section="$(mktemp)"
  tmp="$news.assemble.$$"
  cp -p "$news" "$tmp"
  trap 'rm -f "$section" "$tmp"' EXIT
  {
    echo "# ferx $version"
    echo
    render_sections
  } >"$section"

  # The new section takes the development heading's place (its body is empty,
  # `check` saw to that); without one, it goes before the first release heading,
  # or closes a file that has none. Lines come from `news_headings`, so a
  # heading quoted inside a code block is never the splice point.
  local from to
  from="$(dev_line)"
  if [ -n "$from" ]; then
    to="$(heading_after "$from")"
  else
    from="$(heading_after 0)"
    to="$from"
  fi
  awk -v section="$section" -v from="${from:-0}" -v to="${to:-0}" '
    function emit_section(   line) {
      while ((getline line < section) > 0) print line
      close(section)
      print ""
      done = 1
    }
    # [from, to): the development heading and its empty body are dropped; with
    # no development heading from == to and nothing is.
    !done && from + 0 > 0 && FNR == from + 0 { emit_section() }
    from + 0 > 0 && FNR >= from + 0 && (to + 0 == 0 || FNR < to + 0) { next }
    { print }
    END { if (!done) emit_section() }
  ' "$news" >"$tmp"

  # Never delete a fragment that did not make it into the file: the new heading
  # and every fragment's first line must be present before anything is removed.
  grep -Fqx "# ferx $version" "$tmp" ||
    die "assemble: internal error: the new section was not written; nothing changed"
  local name
  for name in ${FRAGMENTS[@]+"${FRAGMENTS[@]}"}; do
    grep -Fqx -- "$(sed -n 1p "$frag_dir/$name")" "$tmp" ||
      die "assemble: internal error: changelog.d/$name was not written; nothing changed"
  done
  mv "$tmp" "$news"

  local n=0
  for name in ${FRAGMENTS[@]+"${FRAGMENTS[@]}"}; do
    rm -f "$frag_dir/$name"
    n=$((n + 1))
  done
  echo "changelog: assembled $n fragment(s) into '# ferx $version' and removed them."
  echo "Review NEWS.md - add any release-level prose (e.g. upgrade notes) under the"
  echo "new heading by hand - then commit it together with the deleted fragments."
}

# Test-only files under a user-facing prefix.
is_test_path() {
  case "$1" in
    */tests/*) return 0 ;;
  esac
  return 1
}

# Each newline-separated path of $1 that is user-facing, indented. Same
# newline split as `load_fragments`, for the same here-document reason; a
# function of its own so the IFS change cannot leak into the caller.
user_facing_of() {
  local IFS=$'\n' f
  set -f
  for f in $1; do
    if is_user_facing "$f"; then
      printf '  %s\n' "$f"
    fi
  done
  set +f
}

# $1 as `git diff --name-only` prints it under `core.quotePath=false`: plain
# for UTF-8, but a path with a quote, backslash or control character comes out
# C-quoted (`"R/a\"b.R"`). Dropping a leading quote is enough for the prefix
# test, and errs toward "user-facing", which is the safe direction.
is_user_facing() {
  local p path="${1#\"}"
  is_test_path "$path" && return 1
  for p in "${USER_FACING_PREFIXES[@]}"; do
    case "$path" in "$p"*) return 0 ;; esac
  done
  return 1
}

# 0 if fragment $2, as of commit $1, appears as ONE contiguous run - same lines,
# same order, same multiplicity - among the lines this PR ADDED to NEWS.md ($3):
# i.e. `assemble` wrote it, continuation and all. Not a first-line match
# anywhere in the file (a dropped continuation, or text that already existed,
# would pass) and not set membership (reordered or deduplicated lines would).
#
# Two checks, because each alone has a blind spot:
#  - the `+` run: blank lines are dropped from both sides, since git may align
#    a blank line inside the inserted section with a pre-existing one as
#    context, splitting the run without the text having changed. That proves
#    the lines were ADDED, in order, but is blind to a lost interior blank line;
#  - the exact body, interior blank lines and all (trailing ones trimmed, as
#    `assemble` renders it), must be one of the release section's ENTRIES
#    (directory $4, from `split_entries`) under the fragment's own category
#    heading. Each entry is consumed by the fragment it matches, so two equal
#    fragments need two entries; text inside a code block, a fix filed under
#    Breaking changes, or a bullet appended to an older release matches none.
fragment_assembled() {
  local body needle haystack stem heading h b
  body="$(git -C "$root" show "$1:$2" | awk '{ l[NR] = $0 } NF { last = NR }
    END { for (i = 1; i <= last; i++) print l[i] }')" || return 1
  [ -n "$body" ] || return 1
  needle="$(printf '%s\n' "$body" | grep -v '^[[:space:]]*$' || true)"
  haystack="$(printf '%s\n' "$3" | grep -v '^[[:space:]]*$' || true)"
  # Newline-delimited on both ends, so a match is whole lines; quoted, so the
  # needle is literal text, not a glob.
  case $'\n'"$haystack"$'\n' in
    *$'\n'"$needle"$'\n'*) ;;
    *) return 1 ;;
  esac
  stem="${2##*/}"
  stem="${stem%.md}"
  is_category "${stem##*.}" || return 1
  heading="## $(heading_of "${stem##*.}")"
  # The caller walks its deletions under `set -f`; the glob needs it lifted, and
  # put back on every way out.
  local noglob=0 found=1
  case $- in *f*) noglob=1 ;; esac
  set +f
  for h in "$4"/*.h; do
    [ -f "$h" ] || continue
    b="${h%.h}.b"
    [ "$(cat "$h")" = "$heading" ] || continue
    if [ "$(cat "$b")" = "$body" ]; then
      rm -f "$h" "$b"
      found=0
      break
    fi
  done
  [ "$noglob" -eq 0 ] || set -f
  return "$found"
}

# Split the release section on stdin into its entries, one pair of files per
# entry in directory $1: `<k>.h` the `## ` heading it sits under, `<k>.b` its
# bullet with trailing blank lines trimmed. An entry is a column-0 `- ` line
# outside any code block, plus the lines that continue it (indented, blank, or
# inside a fence it opened); any other column-0 line ends it.
split_entries() {
  awk -v dir="$1" "$FENCE_AWK"'
    function flush(   i, last, f) {
      if (!k || written[k]) return
      written[k] = 1
      last = 0
      for (i = 1; i <= n; i++) if (line[i] ~ /[^ \t]/) last = i
      f = dir "/" k ".b"
      for (i = 1; i <= last; i++) print line[i] > f
      close(f)
      f = dir "/" k ".h"
      print heading[k] > f
      close(f)
    }
    {
      was = fence_n
      fence_step($0)
      if (was) { if (cur) line[++n] = $0; next }
      if ($0 ~ /^## /) { flush(); cur = 0; h = $0; sub(/[ \t]+$/, "", h); next }
      if ($0 ~ /^- /) { flush(); k++; heading[k] = h; n = 0; line[++n] = $0; cur = 1; next }
      if ($0 ~ /^[^ \t]/) { flush(); cur = 0; next }
      if (cur) line[++n] = $0
    }
    END { flush() }
  '
}

# The text of the release section this PR adds, as `assemble` writes it: NEWS.md's
# first heading is `# ferx <DESCRIPTION's Version>`, a release heading the merge
# base ($1) did not have. Prints nothing when the PR adds no such section - then
# it is not a release, and deletes fragments for no release to carry them.
new_release_section() {
  local first heading version described base_news from to
  first="$(heading_after 0)"
  [ -n "$first" ] || return 0
  heading="$(sed -n "${first}p" "$news")"
  printf '%s' "$heading" | grep -Eq "$DEV_RE" && return 0
  version="${heading#\# ferx }"
  # The version `assemble` would accept, nothing looser: a `# ferx 0.4.0.9000`
  # heading matching an unbumped DESCRIPTION is no release.
  is_release_version "$version" || return 0
  [ -f "$description" ] || return 0
  described="$(sed -n 's/^Version:[[:space:]]*//p' "$description" | tr -d '[:space:]')"
  [ "$version" = "$described" ] || return 0
  base_news="$(mktemp)"
  git -C "$root" show "$1:NEWS.md" >"$base_news" 2>/dev/null || true
  if news_headings "$base_news" | cut -f2- | grep -Fqx -- "$heading"; then
    rm -f "$base_news"
    return 0
  fi
  rm -f "$base_news"
  from="$first"
  to="$(heading_after "$from")"
  awk -v from="$from" -v to="${to:-0}" '
    to + 0 > 0 && FNR >= to + 0 { exit }
    FNR >= from + 0             { print }
  ' "$news"
}

cmd_require() {
  local base="" pr="" opt_out="false"
  while [ $# -gt 0 ]; do
    case "$1" in
      --pr)
        [ $# -ge 2 ] || die "--pr needs a value"
        pr="$2"
        shift 2
        ;;
      --opt-out)
        [ $# -ge 2 ] || die "--opt-out needs true or false"
        opt_out="$2"
        shift 2
        ;;
      -*) die "require: unknown flag '$1'" ;;
      *)
        [ -z "$base" ] || die "require: more than one base ref given"
        base="$1"
        shift
        ;;
    esac
  done
  [ -n "$base" ] || die "require: usage: tools/changelog.sh require <base-ref> [--pr N] [--opt-out true|false]"
  case "$opt_out" in
    true | false) ;;
    *) die "require: --opt-out must be true or false, got '$opt_out'" ;;
  esac

  # Pending fragments are other PRs' entries: deleting one - or renaming it
  # away, which `--no-renames` reports as a deletion - drops it from the next
  # release. The one legitimate deletion is `assemble`'s, so a deleted fragment
  # passes only if this PR added all of it to NEWS.md. Checked before the
  # opt-out: the label waives this PR's own entry, not anyone else's.
  local mb deleted lost="" added_lines renamed
  mb="$(git -C "$root" merge-base "$base" HEAD)" ||
    die "require: no merge base between $base and HEAD (CI needs fetch-depth: 0)"

  # A rename is never a release step: `assemble` deletes, it does not move. And
  # copying the old body into NEWS.md must not excuse one - the renamed file
  # would be released a second time.
  renamed="$(git -C "$root" diff -M --name-status --diff-filter=R "$mb" HEAD -- changelog.d/ || true)"
  if [ -n "$renamed" ]; then
    {
      echo "changelog: this PR renames pending changelog fragment(s):"
      printf '%s\n' "$renamed" | sed 's/^/  /'
      echo
      echo "Leave other PRs' fragments where they are (editing one to reword it is fine)."
    } >&2
    exit 1
  fi
  deleted="$(git -C "$root" diff --no-renames --name-only --diff-filter=D "$mb" HEAD -- changelog.d/ |
    grep -Ev '^changelog\.d/(README\.md|\.gitkeep)$' || true)"
  local release=""
  if [ -n "$deleted" ]; then
    added_lines="$(git -C "$root" diff --no-renames "$mb" HEAD -- NEWS.md |
      sed -n -e '/^+++ /d' -e 's/^+//p')"
    release="$(new_release_section "$mb")"
    local entries
    entries="$(mktemp -d)"
    printf '%s\n' "$release" | split_entries "$entries"
    local IFS_SAVE="$IFS" f
    IFS=$'\n'
    set -f
    for f in $deleted; do
      if [ -z "$release" ] || ! fragment_assembled "$mb" "$f" "$added_lines" "$entries"; then
        lost="$lost  $f"$'\n'
      fi
    done
    set +f
    IFS="$IFS_SAVE"
    rm -rf "$entries"
  fi
  if [ -n "$lost" ]; then
    {
      echo "changelog: this PR deletes pending changelog fragment(s) that are not in"
      echo "a new release section of NEWS.md, so their entries would be missing from"
      echo "the next release:"
      printf '%s' "$lost"
      echo
      if [ -z "$release" ]; then
        echo "This PR adds no release section: NEWS.md's first heading must be a new"
        echo "'# ferx <version>' matching DESCRIPTION's Version, as 'assemble' writes it."
      fi
      echo "Restore them (edit one to reword it). Only 'tools/changelog.sh assemble'"
      echo "removes fragments, at release time."
    } >&2
    exit 1
  fi

  # A release PR (it assembled fragments) must leave none behind. In CI HEAD is
  # the merge with `main`, so this also catches a fragment that landed on `main`
  # after the release branch was assembled - its entry would miss the release.
  load_fragments
  if [ -n "$deleted" ]; then
    if [ ${#FRAGMENTS[@]} -gt 0 ]; then
      {
        echo "changelog: this release PR assembled fragments but leaves pending ones:"
        printf '  changelog.d/%s\n' "${FRAGMENTS[@]}"
        echo
        echo "Merge main into the release branch and run 'tools/changelog.sh assemble'"
        echo "again, so every entry is in the release section."
      } >&2
      exit 1
    fi
    echo "changelog: release PR - every pending fragment is assembled into NEWS.md."
    return 0
  fi

  # A boolean the workflow computes with `contains(labels.*.name, ...)`, not a
  # joined label list: joining loses the boundaries, so a single label named
  # `x,no-changelog` would read as the opt-out.
  if [ "$opt_out" = true ]; then
    echo "changelog: '$OPT_OUT_LABEL' label set - no fragment required."
    return 0
  fi

  local changed facing="" added=""
  # `--no-renames`: a rename reports only its destination, so moving a file
  # from R/ into tests/ would hide the production file it removed.
  changed="$(git -C "$root" -c core.quotePath=false diff --no-renames --name-only "$base"...HEAD)" ||
    die "require: git diff $base...HEAD failed (is the base fetched? CI needs fetch-depth: 0)"

  facing="$(user_facing_of "$changed")"

  # Only an ADDED fragment counts: rewording one of the pending fragments from
  # another PR does not describe this one. The opposite rename setting from the
  # diff above, and deliberately: with `-M` a renamed fragment is an `R`, not an
  # `A`, so renaming someone else's entry cannot pass for adding one.
  # Directly under changelog.d/ only: `changelog.d/.gitkeep/x.md` is no fragment.
  added="$(git -C "$root" diff -M --name-only --diff-filter=A "$base"...HEAD -- changelog.d/ |
    grep -E '^changelog\.d/[^/]+$' | grep -Ev '^changelog\.d/(README\.md|\.gitkeep)$' || true)"

  if [ -z "$facing" ]; then
    echo "changelog: no user-facing paths changed - no fragment required."
    return 0
  fi
  if [ -n "$added" ]; then
    echo "changelog: fragment(s) present:"
    printf '%s\n' "$added" | sed 's/^/  /'
    return 0
  fi

  local want="changelog.d/${pr:-<N>}.<category>.md"
  {
    echo "changelog: this PR changes user-facing code but adds no NEWS fragment."
    echo
    echo "Missing: $want"
    echo "  <category> is one of: ${CATEGORIES[*]}"
    echo "  e.g.  changelog.d/${pr:-<N>}.fix.md  containing one bullet:"
    echo "        - **What changed, in user-facing words** (#${pr:-<N>}). ..."
    echo
    echo "If the change is not user-facing (refactor, tests, CI), apply the"
    echo "'$OPT_OUT_LABEL' label to the PR instead."
    echo
    echo "User-facing paths changed:"
    printf '%s\n' "$facing"
  } >&2
  exit 1
}

[ $# -ge 1 ] || {
  usage >&2
  exit 2
}
sub="$1"
shift
case "$sub" in
  check) cmd_check "$@" ;;
  preview) cmd_preview "$@" ;;
  assemble) cmd_assemble "$@" ;;
  require) cmd_require "$@" ;;
  -h | --help | help) usage ;;
  *)
    echo "changelog: unknown subcommand '$sub' (try --help)" >&2
    exit 2
    ;;
esac
