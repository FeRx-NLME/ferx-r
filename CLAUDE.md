# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ferx is an R package for nonlinear mixed effects (NLME) modeling. It provides a user-facing R API that delegates to a high-performance Rust backend via extendr FFI. The core computation engine lives in the sibling **ferx-core** crate at `../ferx-core` (sibling directory).

## Worktree isolation

When working on a feature branch or any branch other than `main`, always use `EnterWorktree` at the start of the session. This prevents uncommitted WIP from one session contaminating another session on a different branch (a real problem when two chats share the same checkout directory).

## ferx-core dependency: do NOT edit `src/rust/Cargo.toml`

`src/rust/Cargo.toml` declares **two** git dependencies pinned to `main`, both from the ferx-core repository: `ferx-core` (the engine) and `ferx-tools` (model-development tooling built on it — the bootstrap, and whatever follows). Don't change either to a path dep when working locally — the local-vs-GitHub swap is handled by `src/rust/.cargo/config.toml`:

```toml
[patch."https://github.com/FeRx-NLME/ferx-core"]
ferx-core = { path = "../../../ferx-core" }
ferx-tools = { path = "../../../ferx-core/crates/ferx-tools" }
```

**That file is generated, not committed** — `src/rust/.cargo/` is gitignored and `src/Makevars` rewrites the file from scratch on every build (`R CMD INSTALL`, `pkgload::load_all(recompile = TRUE)`, ...). So a `[patch]` line you add by hand survives only until the next build: the change belongs in `src/Makevars`.

Cargo patches **per package name**: an entry for `ferx-core` alone leaves `ferx-tools` resolving to GitHub `main` while `ferx-core` comes from your working tree — two revisions of a workspace whose halves move together, with no error to say so. Both entries, or neither.

When a sibling `../ferx-core` checkout exists, `src/Makevars` writes that `[patch]`; when it doesn't (e.g. CI without a paired checkout), nothing is patched and cargo builds the revision `Cargo.lock` pins. Writing the patch is not the same as cargo using it, so Makevars' `ferx: [patch]ing ...` line says only that the patch was written.

**A patch that does not apply is a warning, not an error.** Cargo prints `warning: patch ... was not used in the crate graph` on stderr, builds the pinned GitHub revision instead, and nothing else in the output looks wrong. Its follow-up hint says to run `cargo update` — don't: with the sibling present that strips the pin (below), and `tools/update-ferx-core-lock.sh` is the safe equivalent. Two ways to get a "local" build that is not local:

- **The sibling's version differs from the locked one — by any amount, in either direction.** Cargo uses a `[patch]` only when the patched crate's version equals the version `Cargo.lock` pins for it. A patch-level difference is enough (`0.4.0` locked, `0.4.1` in the sibling), so is an older sibling (`0.3.9` against `0.4.0`); dependents' version requirements don't enter into it, and the semver boundary #346 crossed (`0.3.1` -> `0.4.0`) was just one instance. The check is per crate, so a sibling whose `ferx-tools` version alone differs gets `ferx-core` patched and `ferx-tools` from the pin — the mixed build the per-package note above warns about. Bumping the lock (below) moves the pin to ferx-core `main`, which restores the patch when the sibling is at `main`'s version; otherwise check out a sibling revision at the locked version.
- **Working inside a `.claude/worktrees/` worktree.** Makevars looks for the sibling at `../../ferx-core` relative to `src/`, which from `<repo>/.claude/worktrees/<name>/src` is `<repo>/.claude/worktrees/ferx-core` — never there. So **no** `[patch]` is written and cargo has nothing to warn about: the build uses the pin, as CI does, and Makevars says `ferx: no sibling checkout at ... (this is a git worktree)`. Since worktrees are mandated above, that is the normal case. To build worktrees against your sibling instead, link it where the lookup lands, once, from the repo root: `ln -s ../../../ferx-core .claude/worktrees/ferx-core`. The link redirects **every** worktree build under `.claude/worktrees/` — PR-review and pin-bump worktrees included, which then stop building the pin CI builds — and each of those builds rewrites the lock (below). Remove it before a review or a bump. Don't hand-write paths into `src/rust/.cargo/config.toml` instead: Makevars truncates that file on every build.

### With the `[patch]` in place, cargo rewrites `Cargo.lock`

**Any** cargo command that resolves while the patch is in place rewrites the lock. A patch that applies deletes both `source = "git+..."` lines (cargo records no source for a path package), silently unpinning the crates for CI and everyone else; one that goes unused keeps them but appends `[[patch.unused]]` tables. Measured, from a pin of `2`:

| command | git `source` lines afterwards |
|---|---|
| `cargo metadata` | `0` |
| `cargo tree` | `0` |
| `cargo build --release` | `0`, at resolve time — long before compilation finishes |

That includes `cargo update` — both past unpinnings (`1ce7f59`, `b96c867`) were lock bumps — and every build that goes through `src/Makevars`: `R CMD INSTALL .`, `roxygen2::roxygenize()`, `pkgload::load_all()`, `devtools::test()`. It also includes an editor's rust-analyzer, which loads `src/rust` through `cargo metadata`.

So after any build or cargo run with the sibling present, check the lock before staging. This runs no cargo, so it is always safe, and it is exactly what CI runs — both crates pinned to a git source, on one revision, no `[[patch.unused]]` tables:

```bash
tools/check-ferx-core-pin.sh
```

If it fails, read `git diff src/rust/Cargo.lock` before restoring; the script names the same three cases:

- **Only the damage:** `git checkout -- src/rust/Cargo.lock`.
- **A pin bump you meant:** re-run `tools/update-ferx-core-lock.sh`.
- **Any other change you meant**, such as a dependency newly added to `src/rust/Cargo.toml`: `git checkout -- src/rust/Cargo.lock`, move `src/rust/.cargo/config.toml` aside, run cargo in `src/rust` (e.g. `cargo metadata --format-version 1 >/dev/null`), move `config.toml` back, and check again. That records the change and nothing else. Not the bump script: its whole-graph `cargo update` also moves every registry crate and the pin to ferx-core `main`.

Don't count on `git status` to flag a damaged lock: it shows up as an ordinary modified file.

This means: develop against a feature branch in `../ferx-core` freely, but never commit Cargo.toml changes that flip the dep to a path. Reviewers and CI run against the GitHub `main`, so a path dep in Cargo.toml would break their builds.

### Bumping the pinned ferx-core commit (`src/rust/Cargo.lock`)

Although `Cargo.toml` tracks `branch = "main"`, **CI builds from the commit pinned in `src/rust/Cargo.lock`** — *not* the latest `main`. The patch above only redirects local builds (which have the sibling), so a new ferx-core commit is invisible to CI until the lock is bumped. Symptom: a ferx-r PR that uses a freshly-`pub`'d ferx-core API fails CI with `error[E0603]: ... is private`, because the lock still points at a commit predating the change.

To bump: run `tools/update-ferx-core-lock.sh` from the repo root (it advances both pins to ferx-core `main` HEAD and verifies each `source = "git+..."` line survives, on one shared revision), then commit `src/rust/Cargo.lock`. `ferx-core` and `ferx-tools` come from the same repo and the same rev, so this is one bump and two lock entries — never a second pin to track. **Do not** run a bare `cargo update -p ferx-core` with the sibling present — the `[patch]` makes cargo strip the git pin and write a local path, silently unpinning it for CI and everyone else. The `R-CMD-check` workflow runs `tools/check-ferx-core-pin.sh`, which fails if either crate loses its git source, if the two land on different revisions, or if the lock carries `[[patch.unused]]` tables.

## Build & Install

Requires a standard stable Rust toolchain.

```bash
# Install the R package (triggers Rust compilation via src/Makevars)
R CMD INSTALL .

# Build Rust library only
cd src/rust && cargo build --release

# Regenerate roxygen documentation
Rscript -e 'roxygen2::roxygenize()'
```

**With a sibling `../ferx-core` patched in, all three commands rewrite `src/rust/Cargo.lock`** — `roxygenize()` too, since it compiles through `src/Makevars` — so run `tools/check-ferx-core-pin.sh` before staging (see the dependency section).

## Architecture

```
R API (R/*.R)
  → extendr FFI wrappers (.Call)
    → Rust glue (src/rust/src/lib.rs)
      → ferx-core crate (parser, fitting engine, PK models, ODE solver)
```

**R layer** (`R/`): User-facing functions — `ferx_fit()`, `ferx_simulate()`, `ferx_predict()`, `ferx_example()` — that handle path normalization, input validation, and result structuring into S3 objects.

**Rust glue** (`src/rust/src/lib.rs`): ~250 lines. Three `#[extendr]` functions that call into ferx-core and convert results to R lists/data frames via `fit_result_to_list()` and `sdtab_to_dataframe()`.

**ferx-core** (external crate): The actual NLME engine — model parser, NONMEM CSV reader, FOCE/FOCEI optimizer, analytical PK models (1/2-cpt), ODE solver (RK45), simulation, and prediction.

**`R/extendr-wrappers.R`** is auto-generated by extendr — do not edit by hand.

## Model DSL

Models are defined in `.ferx` text files with sections: `[parameters]`, `[individual_parameters]`, `[structural_model]`, `[odes]`, `[error_model]`, `[fit_options]`, and the optional `[scaling]` block (unit conversion / observation readout). See `inst/examples/models/` for references.

Section names are **closed-world** in the engine (ferx-core #1040): a header it does not recognise is a parse error (`E_UNKNOWN_BLOCK`), not a silently ignored block. Never hard-code the valid list in R — `ferx_rust_known_blocks()` returns it from the engine. `[initial_values]` is *not* one of them: the engine dropped it in ferx-core e5e934d (inits live inline in `[parameters]`) and silently ignored leftovers ever since; a file still carrying one now fails to parse with `E_DEPRECATED_BLOCK`. Data uses NONMEM CSV format (ID, TIME, DV, EVID, AMT, CMT, etc.); optional columns: RATE, MDV, II, SS, CENS, OCC. Use `ferx_example("warfarin")` to get paths to bundled examples.

## Bundled examples (`inst/examples/`)

`ferx_example()` discovers examples by scanning `inst/examples/models/*.ferx`, and `?ferx_example` documents each one. A bundled example is **complete only when all of these are committed together** - shipping a subset is a recurring drift source that leaves an installed package advertising examples it cannot load:

- `inst/examples/models/<name>.ferx` - the model
- its data: `inst/examples/data/<name>.csv`, **or** a `.data_aliases` entry in `R/example.R` pointing at an existing CSV
- `inst/examples/ex_<name>.R` - a runnable end-to-end script
- the `\item{<name>}{...}` entry in `ferx_example()`'s `@details` (in `R/example.R`)
- the regenerated `man/ferx_example.Rd` (run `roxygen2::roxygenize()`)
- a `NEWS.md` entry

Before treating example work as done:

- `git status` must show **no untracked** `inst/examples/` files - the assets have to ship, not sit as local files. (The TTE / `transit_savic` examples were authored but left untracked, so the survival feature shipped without them.)
- **Run the `ex_<name>.R` script end-to-end** against a fresh `FERX_NO_AUTODIFF=1 R CMD INSTALL .`. Parse-only checks miss API-signature drift - e.g. an example calling `ferx_predict_survival(fit)` when the signature is `ferx_predict_survival(model, data, times, fit = NULL)`.

`tests/testthat/test-example.R` already asserts every `ferx_example()` data path exists on disk, but nothing runs the scripts - so the script-execution check above is manual.

## Key Return Structures

`ferx_fit()` returns an S3 "ferx_fit" list: converged, method, ofv/aic/bic, theta (named), omega (matrix), sigma, standard errors, sdtab (diagnostic data frame with PRED/IPRED/CWRES/IWRES/ETAs), and warnings.

`ferx_simulate()` returns a data frame: SIM, ID, TIME, IPRED, DV_SIM.

`ferx_predict()` returns a data frame: ID, TIME, PRED.

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
- Use `.fitrx_named_omega()` in `persist.R` when deserialising the omega matrix so dimnames survive a `ferx_save`/`ferx_load` round-trip.
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

**Before every PR, verify with:**
```r
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  n <- sum(chartr(rawToChar(as.raw(128:255)), strrep("x", 128), readLines(f, warn = FALSE)) != readLines(f, warn = FALSE))
  if (n > 0) message(f, ": ", n, " non-ASCII lines")
}
```
Or from the shell: `python3 -c "import os; [print(f) for f in __import__('glob').glob('R/*.R') if any(b > 127 for b in open(f,'rb').read())]"`

## Pull Requests

When creating a PR in this repo, always read `.github/PULL_REQUEST_TEMPLATE.md` and fill every section before calling `gh pr create`.
