# AGENTS.md

This file provides guidance to AI coding agents (Claude Code, Codex, Copilot, and any other tool that reads the [AGENTS.md](https://agents.md) convention) when working with code in this repository.

## What This Is

ferx is an R package for nonlinear mixed effects (NLME) modeling. It provides a user-facing R API that delegates to a high-performance Rust backend via extendr FFI. The core computation engine lives in the sibling **ferx-core** crate at `../ferx-core` (sibling directory).

## Code style: parsimony

Prefer the simplest solution that solves the problem.

- Write the least code that does the job correctly. No speculative generality, extra options, or layers that nothing asks for yet.
- Factor shared logic into functions rather than duplicating it: two code paths computing the same thing should call the same function, or the second copy will drift. Numerics belong in ferx-core; the R layer and the glue convert and validate, they do not recompute.
- Keep functions small and single-purpose; extract a helper when a block *you are already changing* has a nameable job.
- A PR fixes its issue, not the area around it; a neighbouring defect gets its own issue - unless it is a small review finding inside the feature, which is fixed in the same PR.
- Reuse existing helpers before writing new ones - search `R/` and `src/rust/src/` first.
- Delete dead code rather than leaving it commented out or behind a flag nothing enables (code behind the `nn`/`survival`/`markov` cargo features is not dead).
- Parsimony governs production code. It is not a reason to drop test assertions, or comments recording a measured number.

## Worktree isolation

When working on a feature branch or any branch other than `main`, always use `EnterWorktree` at the start of the session. This prevents uncommitted WIP from one session contaminating another session on a different branch (a real problem when two chats share the same checkout directory).

## ferx-core dependency: do NOT edit `src/rust/Cargo.toml`

`src/rust/Cargo.toml` declares **three** git dependencies pinned to `main`, all from the ferx-core repository: `ferx-core` (the engine), `ferx-tools` (model-development tooling: the bootstrap and GAM screening) and `ferx-search` (model search: `search`, `covsearch`, `iivsearch`, `amd`, ...). The list lives once, as `FERX_CRATES` in `tools/ferx-core-lock-lib.sh`; a fourth crate from that repository goes there, in `Cargo.toml` and in `.github/dependabot.yml`'s ignore list.

**Never change any of them to a path dep**, not even locally. Reviewers and CI build the GitHub pin, so a path dep breaks their builds. The local swap is a cargo `[patch]` that `src/Makevars` hands, through `tools/sibling-cargo-build.sh`, to one cargo invocation whenever a sibling `../ferx-core` checkout exists - for all three crates, or the wrapper refuses the build. Without a sibling (CI), Makevars runs cargo directly on the revision `Cargo.lock` pins. Develop against a feature branch in `../ferx-core` freely. How the wrapper works and how it is tested: `tools/README.md`.

**A patch that does not apply is a warning, not an error**: cargo builds the pinned revision and nothing else looks wrong. Read the wrapper's last line - `ferx: the sibling supplied ALL of ...` or `ferx: WARNING the sibling was NOT used ...`. Two ways to get a "local" build that is not local:

- **The sibling's version differs from the locked one**, by any amount, in either direction (`0.4.0` locked vs `0.4.1` or `0.3.9` in the sibling). Cargo only patches a crate whose version equals the locked one. Bump the lock (below) when the sibling is at ferx-core `main`; otherwise check out a sibling revision at the locked version.
- **Working inside a `.claude/worktrees/` worktree.** Makevars looks for `../../ferx-core` relative to `src/`, i.e. `<repo>/.claude/worktrees/ferx-core`, usually absent - so the build uses the pin, as CI does, and says `ferx: no sibling checkout at ... (this is a git worktree)`. To build worktrees against the sibling, link it once from the repo root: `ln -s ../../../ferx-core .claude/worktrees/ferx-core`. The link redirects **every** worktree build, PR-review and pin-bump worktrees included: remove it before a review or a bump, or neutralise it for one command with `MAKEFLAGS="LOCAL_FERX_CORE=" R CMD INSTALL .`. Never hand-write a `[patch]` into `src/rust/.cargo/config.toml`: Makevars truncates that file on every build (#353).

### The lock, and what still rewrites it

Any cargo resolve **with a patch in place** rewrites `src/rust/Cargo.lock`, silently unpinning the crates for CI. Only the wrapper's own cargo run carries a patch, and it puts the lock back:

| command | with a sibling | without one |
|---|---|---|
| `R CMD INSTALL .`, `roxygen2::roxygenize()`, `pkgload::load_all()`, `devtools::test()` | patched build, lock restored byte-for-byte | pinned build, lock untouched |
| `cargo metadata` / `tree` / `build` in `src/rust`, rust-analyzer | pinned build, lock untouched | pinned build, lock untouched |
| `cd src && sh ../tools/sibling-cargo-build.sh <args>` | patched run, lock restored byte-for-byte | refused |
| `cargo update` in `src/rust` | moves the pin to today's `main` - use `tools/update-ferx-core-lock.sh` | same |

Run one sibling build at a time per checkout (a second one waits). Damage can still reach the lock - a stale `config.toml` from before #353, a `cargo update`, a `kill -9` - and `git status` shows it as an ordinary modified file. **Check before staging**; this runs no cargo, and it is exactly what CI runs:

```bash
tools/check-ferx-core-pin.sh
```

If it fails, read `git diff src/rust/Cargo.lock` before restoring; the script names the same three cases:

- **Only the damage:** `git checkout -- src/rust/Cargo.lock`.
- **A pin bump you meant:** re-run `tools/update-ferx-core-lock.sh`.
- **Any other change you meant**, such as a dependency newly added to `src/rust/Cargo.toml`: `git checkout -- src/rust/Cargo.lock`, then run cargo in `src/rust` (e.g. `cargo metadata --format-version 1 >/dev/null`) and check again. Not the bump script: its whole-graph `cargo update` also moves every registry crate and the pin to ferx-core `main`.

### Bumping the pinned ferx-core commit (`src/rust/Cargo.lock`)

**CI builds the commit pinned in `src/rust/Cargo.lock`**, not the latest `main`, so a new ferx-core commit is invisible to CI until the lock is bumped. Symptom: a PR using a freshly-`pub`'d ferx-core API fails CI with `error[E0603]: ... is private`.

To bump: run `tools/update-ferx-core-lock.sh` from the repo root (it moves all three crates to ferx-core `main` HEAD on one shared revision - one bump, never a second pin), then commit `src/rust/Cargo.lock` and the regenerated `inst/THIRD-PARTY-NOTICES.txt`. To pin a specific revision, run `cargo update -p ferx-core --precise <sha>` in `src/rust`, never through the wrapper.

### License notices (#773)

ferx-r's own code is MIT, but the compiled package statically links ferx-search (LGPL-3.0-or-later) and NLopt's Luksan routines (LGPL-2.1-or-later), so it ships `inst/COPYRIGHTS`, `inst/licenses/` and `inst/THIRD-PARTY-NOTICES.txt`. **Never describe the whole package as MIT** (README, DESCRIPTION, image labels). `inst/THIRD-PARTY-NOTICES.txt` is generated from `Cargo.lock` by `tools/third-party-notices.py`: the bump script regenerates it; after any other lock change (a new dependency) run the script yourself and commit the file with the lock - R-CMD-check runs `--check` and fails a stale copy.

## Build & Install

Requires a standard stable Rust toolchain.

```bash
# Install the R package (triggers Rust compilation via src/Makevars)
R CMD INSTALL .

# Build Rust library only, against the revision Cargo.lock pins
cd src/rust && cargo build --release

# Build Rust library only, against the sibling ../ferx-core checkout
cd src && sh ../tools/sibling-cargo-build.sh build --release

# Regenerate roxygen documentation
Rscript -e 'roxygen2::roxygenize()'
```

Building against the sibling is an explicit opt-in: the R-driven commands take it automatically when the checkout is there, and plain `cargo` never does. Each leaves `src/rust/Cargo.lock` as it found it; still run `tools/check-ferx-core-pin.sh` before staging.

The cargo feature set is `CARGO_FEATURES` in `src/Makevars`, assigned with `?=` so the environment wins (#288): `CARGO_FEATURES="--no-default-features --features ci,survival,markov" R CMD INSTALL .` builds without `nn` (Deep Compartment Models), the heaviest part of the dependency graph. Unset, the default `ci,nn,survival,markov` is used - that is what CI and every release build compile, so never assume a local no-`nn` build's block list matches theirs.

### Release tarballs and Docker images

Release engineering; most sessions never touch it. How it works: `tools/README.md`. The rules:

- **Never commit `src/rust/vendor.tar.xz`** (gitignored, ~25 MB per pin bump). `tools/vendor-rust.sh` builds it before `R CMD build` of an offline release tarball (#666); `tools/vendor-rust.sh --clean` removes it afterwards, or every later local build in that checkout tries the vendored set first.
- **Secrets never enter an image** (#723). `.dockerignore` is an allowlist: never widen it beyond package files. Never pass a credential as a build arg or `ENV`. Add a `known_public()` entry to `tools/scan-image-secrets.sh` only for a public upstream file you have read.
- **Outside `src/Makevars`, build ferx-core with `cargo +stable`**: its `rust-toolchain.toml` names nightly, which rustup installs unasked. Never clone it to `<package>/../ferx-core` in an image: that is the sibling path, and Makevars patches it in.
- A pin bump needs no Dockerfile edit. A new R version or system library does: see `tools/README.md`.

## Architecture

```
R API (R/*.R)
  → extendr FFI wrappers (.Call)
    → Rust glue (src/rust/src/lib.rs)
      → ferx-core crate (parser, fitting engine, PK models, ODE solver)
```

**R layer** (`R/`): the exported functions in `NAMESPACE` (e.g. `ferx_fit()`, `ferx_simulate()`, `ferx_predict()`, `ferx_example()`). They handle path normalization and input validation, and structure results into S3 objects.

**Rust glue** (`src/rust/src/lib.rs`): the `#[extendr]` FFI layer. Each function calls into ferx-core and converts its result to R lists and data frames; the fit-changing steps also read and rewrite the fit's `.fitrx` bundle (see "A fit's bytes are the fit" below). Build every ferx-core / ferx-tools `*Options` struct with a `..Default::default()` spread, never an exhaustive literal: the engine adds fields freely (ferx-core #529), and a literal then breaks this build (`error[E0063]: missing field`). R-CMD-check enforces it.

**ferx-core** (external crate): the NLME engine - model parser, data reader, estimation methods, analytical and ODE PK models, simulation and prediction. Its [README](https://github.com/FeRx-NLME/ferx-core#readme) is the overview; don't restate its method list here.

**`R/extendr-wrappers.R`** is auto-generated by extendr — do not edit by hand.

## Model DSL

Models are defined in `.ferx` text files made of sections, e.g. `[parameters]`, `[individual_parameters]`, `[structural_model]`, `[odes]`, `[error_model]`, `[fit_options]`; `ferx_rust_known_blocks()` is the list. See `inst/examples/models/` for references.

Section names are **closed-world** in the engine (ferx-core #1040): a header it does not recognise is a parse error (`E_UNKNOWN_BLOCK`), not a silently ignored block. Never hard-code the valid list in R — `ferx_rust_known_blocks()` returns it from the engine. `[initial_values]` is *not* one of them: the engine dropped it in ferx-core e5e934d (inits live inline in `[parameters]`) and silently ignored leftovers ever since; a file still carrying one parses with a `W_DEPRECATED_BLOCK` warning and the block is ignored (ferx-core #2913; it was the parse error `E_DEPRECATED_BLOCK` before). Data uses NONMEM CSV format (ID, TIME, DV, EVID, AMT, CMT, etc.); optional columns: RATE, MDV, II, SS, CENS, OCC. Use `ferx_example("warfarin")` to get paths to bundled examples.

## Bundled examples (`inst/examples/`)

`ferx_example()` discovers examples by scanning `inst/examples/models/*.ferx`, and `?ferx_example` documents each one. A bundled example is **complete only when all of these are committed together** - shipping a subset is a recurring drift source that leaves an installed package advertising examples it cannot load:

- `inst/examples/models/<name>.ferx` - the model
- its data: `inst/examples/data/<name>.csv`, **or** a `.data_aliases` entry in `R/example.R` pointing at an existing CSV
- `inst/examples/ex_<name>.R` - a runnable end-to-end script
- the `\item{<name>}{...}` entry in `ferx_example()`'s `@details` (in `R/example.R`)
- the regenerated `man/ferx_example.Rd` (run `roxygen2::roxygenize()`)
- a NEWS fragment, `changelog.d/<N>.feature.md` (see "NEWS: add a fragment" below)

Before treating example work as done:

- `git status` must show **no untracked** `inst/examples/` files - the assets have to ship, not sit as local files. (The TTE / `transit_savic` examples were authored but left untracked, so the survival feature shipped without them.)
- **Run the `ex_<name>.R` script end-to-end** against a fresh `R CMD INSTALL .`. Parse-only checks miss API-signature drift - e.g. an example calling `ferx_predict_survival(fit)` when the signature is `ferx_predict_survival(model, data, times, fit = NULL)`.

`tests/testthat/test-example.R` already asserts every `ferx_example()` data path exists on disk, but nothing runs the scripts - so the script-execution check above is manual.

## Key Return Structures

`ferx_fit()` returns an S3 "ferx_fit" list: converged, method, ofv/aic/bic, theta (named), omega (matrix), sigma, standard errors, sdtab (diagnostic data frame with PRED/IPRED/CWRES/IWRES/ETAs), and warnings.

**A fit's bytes are the fit** (ferx-core #2693). `fit$bundle` is core's `.fitrx` serialization, and the list is a view of it: `ferx_save_fit()` writes those bytes, and every step that produces a changed fit (`ferx_covariance()`, `ferx_sir()`, `ferx_calc_npde()`, `ferx_collect()`) reads them, makes new ones in the glue, and rebuilds the list with `.ferx_fit_with_bundle()`. Such a step starts from `.ferx_fit_canonical()`, the fit as its bytes say it is, and reads nothing from the caller's list. Never graft fields onto a fit: what is not in the bytes is lost on save, and a list edit does not reach the engine. In tests, put a fit in another state with `fitrx_edit_fit()`, which edits the bytes. Fits are read-only (#713): `$<-`, `[[<-` and `[<-` refuse a core field (`ferx_fit_read_only_error`), and `ferx_save_fit()` refuses a fit whose core fields no longer match the digest `.ferx_fit_seal()` took when it was built (`ferx_fit_edited_error`). Package code that builds or rebuilds a fit edits the unclassed list and seals it at the end; never assign into a sealed fit. A test probing how a reader handles a list field may edit `unsealed_fit(fit)` (helper-fitrx-bundle.R), which a save refuses. New names on a fit are the user's: saved under `r$user` in `interface_extras.json` in the allowlisted node encoding of `.ferx_user_encode()` and restored on load by `.ferx_user_decode()`. A bundle is untrusted: never decode its contents with `jsonlite::unserializeJSON()` or `unserialize()`. R-only fields go in `interface_extras.json` (`.fitrx_interface_extras()`), and the glue swaps in the wide `predictions.csv` and `sdtab_endpoint_rows.csv`, since core's table has fixed columns. `fit$model_path` and `fit$data_path` stay list fields: they say where the files are now. A subject ID's line break is escaped by the glue, since core reads its CSV entries by line: `subject_id_line_breaks.json` names the two private-use characters a bundle uses, so anything else in an ID is the ID's own.

`ferx_simulate()`, `ferx_predict()` and the other readers return data frames whose columns are documented in each function's `@return` (`?ferx_simulate`, `?ferx_predict`); read those rather than a copy here.

## Output Label Convention

All output functions must display the **bare declared variable name** — never wrap in `OMEGA()`, `SIGMA()`, or `KAPPA()`:

| Parameter type | Named | Fallback (no name) |
|---|---|---|
| IIV (omega diagonal) | `ETA_CL` | `OMEGA(1,1)` |
| Sigma | `EPS_PROP` | `SIGMA(1)` |
| IOV (kappa diagonal) | `KAPPA_CL` | `KAPPA1` |
| Theta | `TVCL` | `THETA1` |
| Shrinkage | `ETA_CL shrinkage` | `ETA1 shrinkage` |

**Off-diagonal covariances** (both IIV and IOV) use `~` and `: cov =`:
- Named: `ETA_V ~ ETA_CL : cov = 0.025  (param corr = 0.26)`
- Fallback IIV: `OMEGA(2,1) : cov = 0.025  (...)`

This applies to: `print.ferx_fit`, `ferx_estimates()`, `ferx_cor_matrix()` (via `fit$cov_matrix` dimnames), `fit$omega` row/colnames, `fit$sir_ci_omega`, `fit$sir_ci_sigma`, `summary.ferx_fit`, and any future output surface.

**Implementation rules:**
- Use `eta_names[i]` / `sigma_names[i]` / `kappa_names[i]` directly as labels.
- Always set `rownames`/`colnames` on any new matrix field holding omega/sigma values.
- `ferx_load_fit()` goes through `.ferx_fit_from_raw()`, as `ferx_fit()` does, so dimnames set there survive a `ferx_save_fit()`/`ferx_load_fit()` round trip.
- Update `roxygen2` `@return` docs when adding new output fields.

## Roxygen / Documentation

Always run `roxygen2::roxygenize()` and commit the updated `man/` files before opening or updating a PR — the `.Rd` files are checked into the repo and must stay in sync with the `#'` comments in `R/`.

**No non-ASCII characters anywhere in `R/*.R` files.** `R CMD check --as-cran` requires pure ASCII in all R source files. Violations that have burned us before:

| What to avoid | Use instead |
|---|---|
| `\uXXXX` escape sequences in `#'` comments | literal ASCII |
| em-dash `—`, en-dash `–` | `-` or `:` |
| ellipsis `…` | `...` or `etc.` |
| box-drawing chars in comment banners | `-- Section title --` |

Non-ASCII in `\preformatted{}` blocks or `@examples` code comments (`#' #`) also causes RStudio's help viewer to render `?` or silently truncate the page.

**Before every PR, verify with** (prints nothing when clean):

```r
for (f in Sys.glob("R/*.R")) if (length(tools::showNonASCIIfile(f))) message("  in ", f)
```

## Tests

testthat, under `tests/testthat/` (`test-<area>.R`, shared helpers in `helper-*.R`, fixtures in `fixtures/`). Run them (`devtools::test()` loads the package from source, so it compiles through `src/Makevars` like the other commands in the lock table):

```bash
NOT_CRAN=true Rscript -e 'devtools::test()'
NOT_CRAN=true Rscript -e 'devtools::test(filter = "ferx_simulate")'
```

Without `NOT_CRAN=true`, the `skip_on_cran()` tests are skipped. Only the coverage job sets it; the R CMD check job runs without it. CI runs `R CMD check` (`R-CMD-check.yaml`, which also runs the `tools/` guard scripts) and `covr` coverage (`test-coverage.yaml`). The Codecov targets in `codecov.yml` - 70% project, 80% patch - are informational: they report, they do not fail a PR.

**Every new feature requires a test**, and a bug fix requires a regression test that fails without the fix.

## Pull Requests

When creating a PR in this repo, always read `.github/PULL_REQUEST_TEMPLATE.md` and fill every section before calling `gh pr create`.

After reviewing a PR, always post your review on the PR, once the user has approved the draft.

## Priority labels

Issues and PRs carry GitHub labels `P0`-`P3` (correctness severity x how common the trigger): `P0` silent wrong results in common use; `P1` wrong results in realistic use, or a loud failure blocking common workflows; `P2` edge-case correctness, important features, UX/diagnostics; `P3` performance, refactors, nice-to-haves, research. The applied label is the source of truth, not a priority written in the title.

## NEWS: add a fragment, never edit `NEWS.md`

**In the same PR as any user-facing change, add a fragment under `changelog.d/`** instead of a bullet in `NEWS.md` (#600, the convention ferx-core adopted in its #1545). Every PR used to append under the development heading at the same insertion point, so merging one green PR made every other open one conflict on `NEWS.md` - each then needed a rebase, a conflict resolution and a fresh CI run for nothing but a changelog line. Separate files cannot conflict.

The file is `changelog.d/<N>.<category>.md`: `<N>` the issue or PR number (`-2`, `-3`, ... for a second entry on one number), `<category>` one of `breaking`, `feature`, `fix`, `docs`, `internal` - the `NEWS.md` sections Breaking changes, New features, Bug fixes, Documentation, Internal. Its body is exactly one bullet as it should read in `NEWS.md`: user-facing language, the `#NN` reference, continuation lines indented by two spaces. `changelog.d/README.md` has an example.

Two gates:

- `tools/changelog.sh check` (CI step in `R-CMD-check`) validates every fragment's name and shape, and fails if anything is written under `# ferx X.Y.Z.9000 (development version)` - that heading stays empty between releases.
- The `NEWS fragment` workflow (`.github/workflows/changelog.yml`) fails a PR that changes `R/`, `src/` (the glue and the ferx-core pin in `Cargo.lock`), `inst/`, `NAMESPACE` or `DESCRIPTION` without adding a fragment. A PR that is not user-facing (refactor, tests, CI) takes the **`no-changelog`** label instead. The same workflow fails a PR that deletes or renames another PR's pending fragment, label or not; rewording one is fine.

At release time only: open every issue labelled **`release-checklist`** (`gh issue list --label release-checklist`) Do their before-the-tag items before tagging and their after-the-tag items once the tag is pushed; close each issue when all its items are done, or move what is left to the next release with a comment. Then assemble `NEWS.md` as `changelog.d/README.md` describes.
