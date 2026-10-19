# NEWS fragments

Each user-facing PR adds **one new file here** instead of editing `NEWS.md`
(#600, the same convention as ferx-core's `changelog.d/`). Two PRs adding two
different files cannot conflict, so concurrent green PRs merge back-to-back
without a rebase, a conflict resolution and a fresh CI run each.

**Name:** `<N>.<category>.md`

- `<N>` - the issue or PR number. For a second entry on the same number append
  `-2`, `-3`, ... (`541-2.fix.md`).
- `<category>` - the `NEWS.md` section it belongs in:

  | category   | section            |
  |------------|--------------------|
  | `breaking` | Breaking changes   |
  | `feature`  | New features       |
  | `fix`      | Bug fixes          |
  | `docs`     | Documentation      |
  | `internal` | Internal           |

**Body:** exactly one Markdown bullet, as it should read in `NEWS.md` -
user-facing language, with the `#NN` reference (required). The first line starts
with `- `; every further line is blank or indented by at least two spaces (a
continuation of that bullet, nested bullets and code included).

Never delete or rename another PR's fragment - CI refuses it, since the entry
would vanish from the release. Rewording one is fine.

```markdown
- **`ferx_simulate()` now honours `block_sigma`** (#672): paired rows are drawn
  from the dense residual covariance instead of independent normals.
```

Not user-facing (refactor, tests, CI)? No fragment - put the `no-changelog` label
on the PR instead.

```bash
tools/changelog.sh check            # validate the fragments
tools/changelog.sh preview          # what the next release section will look like
tools/changelog.sh assemble 0.5.0   # release time only, after setting DESCRIPTION's Version
```

`assemble` replaces the empty `# ferx X.Y.Z.9000 (development version)` heading
with `# ferx 0.5.0` and the fragments under their sections, newest first, and
deletes the fragments. The post-release version bump adds a new development
heading, which stays empty.

This README is not a fragment.
