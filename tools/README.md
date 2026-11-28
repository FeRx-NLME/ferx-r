# tools/: how the build and release scripts work

Developer notes on the scripts in this directory: the internals behind the rules
in `AGENTS.md`. `AGENTS.md` says what to do; this file says why it works and how
it is tested. Neither ships in the package (`.Rbuildignore`) nor enters an image
(`.dockerignore` is an allowlist).

## The sibling build: `sibling-cargo-build.sh`

`src/Makevars` builds through this wrapper whenever a sibling `../ferx-core`
checkout exists. It hands cargo one `--config` per crate in `FERX_CRATES`
(`ferx-core-lock-lib.sh`):

```
--config 'patch."https://github.com/FeRx-NLME/ferx-core".ferx-core.path="<abs path to ../ferx-core>"'
--config 'patch."https://github.com/FeRx-NLME/ferx-core".ferx-tools.path="<...>/crates/ferx-tools"'
--config 'patch."https://github.com/FeRx-NLME/ferx-core".ferx-search.path="<...>/crates/ferx-search"'
```

**The patch is not written to a file** (#353). `src/rust/.cargo/config.toml` is
still generated on every build (`src/rust/.cargo/` is gitignored, and Makevars
truncates and rewrites it), but it holds only `[build]`, plus, after an offline
build from a vendored tarball, the `[source]` replacement pointing at
`src/rust/vendor/` - never a patch. So nothing outside that one build sees a
patch: cargo by hand and an editor's rust-analyzer resolve the pinned revision
and leave the lock alone.

Cargo patches **per package name**: an entry for `ferx-core` alone would leave
`ferx-tools` and `ferx-search` resolving to GitHub `main` while `ferx-core` comes
from the working tree - two revisions of a workspace whose parts move together,
with no error to say so. The wrapper passes all of them, and refuses a build in
which only some land: before compiling when the manifests predict it (a crate
whose sibling version differs from the locked one), and after the resolve when
cargo produced one anyway.

### What the wrapper does around cargo

- **Snapshot and restore.** It snapshots `src/rust/Cargo.lock` before cargo runs
  and restores it afterwards, because any resolve with a patch in place rewrites
  the lock: an applied patch deletes the `source = "git+..."` lines (cargo
  records no source for a path package), and an unused one appends
  `[[patch.unused]]` tables.
- **Verdict lines.** Its last line says per crate where the build came from:
  `ferx: the sibling supplied ALL of ...`, or `ferx: WARNING the sibling was NOT
  used ...`. It also says which of three lock outcomes happened:
  `Cargo.lock untouched by this build`, `Cargo.lock restored to its pin`, or a
  `WARNING` that the lock was **already** unpinned before the build, in which
  case it is put back exactly as found and nothing is claimed about a pin.
- **Interrupts and closed pipes.** The restore runs through an interrupt and
  through a reader that goes away mid-build (`R CMD INSTALL . | head`); only
  `kill -9` gets past it.
- **One sibling build at a time per checkout.** A snapshot is only good while
  nobody else rewrites the same lock, so the wrapper holds a guard
  (`src/rust/target/.ferx-lock-guard`) across read, snapshot, build, verdict and
  restore. A second build prints `ferx: waiting for another sibling build in this
  checkout` and starts when the first is done; a guard left by a killed build is
  taken over, with a warning naming the dead process. Before the guard, two
  builds at once ended with whichever finished last putting the other's
  mid-build lock back - unpinned, both exiting 0.

## Pin checks: `check-ferx-core-pin.sh`, `update-ferx-core-lock.sh`

`check-ferx-core-pin.sh` runs no cargo. It fails if any `FERX_CRATES` entry lost
its git source, if they sit on different revisions, or if the lock carries
`[[patch.unused]]` tables, and names the three restore cases `AGENTS.md` lists.
`update-ferx-core-lock.sh` moves all of them to ferx-core `main` HEAD in one
`cargo update`, checks each `source = "git+..."` line survived on one shared
revision, and regenerates `inst/THIRD-PARTY-NOTICES.txt`.

## Test harnesses

CI never has a sibling checkout, so the scripts are driven by stub harnesses
that the `R-CMD-check` workflow runs:

| harness | drives |
|---|---|
| `test-check-ferx-core-pin.sh` | damaged locks |
| `test-sibling-cargo-build.sh` | a stub `cargo`: the applied, unused, mixed, failing, pre-stripped, closed-pipe and concurrent cases, under every shell on the box |
| `test-check-agents-md.sh` | `check-agents-md.sh` (#489) |
| `test-check-glue-*.sh` | the glue lint scripts `check-glue-*.sh` |
| `test-changelog.sh` | `changelog.sh` over scratch trees |
| `test-docker-release-tags.sh` | `docker-release-tags.sh` |
| `test-scan-image-secrets.sh` | `scan-image-secrets.sh` still catches each planted secret |

## Vendored release tarballs: `vendor-rust.sh` (#666)

A release source tarball must build offline: GxP validation rebuilds it with no
network and a read-only `CARGO_HOME`. `vendor-rust.sh` vendors every crate
`Cargo.lock` pins (the ferx-core crates at the pinned commit included) into
`src/rust/vendor.tar.xz`, together with the `[source]` replacement `cargo vendor`
printed. It fails if a vendored crate loses a required license file.

When that file is present and no sibling checkout is, `src/Makevars` unpacks it,
appends the replacement to `.cargo/config.toml` and runs
`cargo build --offline --locked`; without it the build is the usual online one.
A lock bump makes the vendored set stale: Makevars compares the lock packed into
the tarball with the current one and refuses a mismatch by name. CI's
`offline-install` job installs a vendored tarball inside `unshare --net`.

## Docker images (#723)

Two images:

- `Dockerfile.slim` - R + ferx + the `ferx` CLI. Published by
  `.github/workflows/docker.yml` to `ghcr.io/ferx-nlme/ferx-r` on every `main`
  push that can change it, and on `v*` tags (`docker-release-tags.sh` keeps the
  release tags permanent: one image per release, never replaced).
- the root `Dockerfile` - RStudio. Built and smoke-tested on PRs, never
  published.

Both build the CLI and the package from the ferx-core commit `Cargo.lock` pins,
so a pin bump needs no Dockerfile edit. What does:

- **The R version** is the rocker base tag (`rocker/r-ver` in both
  `Dockerfile.slim` stages, `rocker/tidyverse` in `Dockerfile`), written as a
  literal so Dependabot (`docker` ecosystem, monthly, one grouped PR) proposes
  new ones. Keep all three FROM lines on one R version.
- **A new system library** the package needs goes in the slim build stage's
  `apt-get` (build time) or the runtime `apt-get` (run time), and in the root
  `Dockerfile`'s. The PR smoke tests (`library(ferx)` + a fit, both images)
  catch a missing one.
- **Actions** move with the existing `github-actions` Dependabot entry.

Both Dockerfiles copy the license notices (and `ferx --licenses`) to
`/usr/local/share/doc/ferx`, which the smoke tests assert.

`docker-cleanup.yml` deletes old `:sha` images weekly (keeps `latest`, release
tags and the 30 newest). A publish run warns ("GHCR package is private") while
anonymous pulls fail; only an org owner can make the package public.

### Secret scanning

`scan-image-secrets.sh` scans every built image (config, history, filesystem)
before anything is pushed, and then the pushed bytes themselves (pulled by
digest, before any tag points at them; a failure deletes that version). It fails
on tokens, private keys, credential files, `*.lic` files and NONMEM license
banners. `known_public()` lists public upstream files it would otherwise flag.
