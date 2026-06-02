# Full Enzyme autodiff build — fixes & SAEM HMC verification

Tracking: PR #104

This branch carries build-system hardening plus a behavioral SAEM HMC test. The
changes are `cargo check`-clean but still need to be **verified on a machine / CI
that can complete the full fat-LTO Enzyme build** — they could not be verified on the
original dev laptop (OOMs during the fat-LTO memory peak).

## Changes in this PR

1. **`src/Makevars` — force fat LTO on the autodiff path (`override`).**
   `R CMD INSTALL` runs make with environment-overrides-makefile semantics, so a stale
   `CARGO_PROFILE_RELEASE_LTO=thin` in the environment could win on the autodiff path —
   silently building autodiff with **thin** LTO. Thin LTO keeps the differentiated
   ferx-core functions in a separate module from their ferx-r callers, breaking
   cross-crate Enzyme. Fix: `override CARGO_PROFILE_RELEASE_LTO := fat` (autodiff) /
   `:= thin` (FD) makes the makefile value authoritative regardless of the environment.

2. **`src/Makevars` — `FERX_NO_AUTODIFF=0` preflight.**
   When autodiff is forced but the `enzyme` rustup toolchain is missing, fail with a
   clear, actionable message (install Enzyme, or unset/`=1` for an FD build) instead of
   a cryptic `cargo`/`rustup` "toolchain 'enzyme' is not installed" error.

3. **`tests/testthat/test-fit.R` — behavioral SAEM HMC tests.**
   Guarded by `ferx:::ferx_rust_autodiff_enabled()`: `n_leapfrog = 3` ⇒
   `saem_n_subjects_hmc` non-NA integer `> 0` and `<= n_subjects`; `n_leapfrog = 0` ⇒
   MH only (NA / 0). Inert on FD/Tier-1 builds; CI with the Enzyme toolchain runs them.
   (Previously only a `NULL` placeholder existed for `saem_n_subjects_hmc`.)

4. **`.gitignore`** — ignore `*.Rcheck/`, `..Rcheck/`, `Rplots.pdf`.

5. **`README.md`** — split build cost into two one-time steps; CPU-dependent timing;
   fat-LTO "looks like a hang but isn't" note.

6. **`NEWS.md`** — entry for the two Makevars fixes.

(roxygen for `n_leapfrog` / `saem_n_subjects_hmc` in `R/fit.R` is already complete.)

## Build (on capable hardware / CI with the Enzyme toolchain)

```bash
FERX_NO_AUTODIFF=0 CARGO_PROFILE_RELEASE_LTO=fat R CMD INSTALL .
```

- The autodiff final link is single-threaded **fat LTO + Enzyme**: ~15–30 min on a
  fast/native CPU, up to ~1–2 h on an older laptop. Extra cores do not help.
- Fat LTO spikes memory at the end — ensure several GB free or rustc may be OOM-killed.
- Do **not** build inside a git worktree: the Makevars `[patch]` to a local
  `../ferx-core` resolves relative to the package root and breaks in a worktree.
- Ensure `ferx-r` and `ferx-core` `main` are in sync (any in-progress cross-repo
  migrations, e.g. new `FitResult` fields, must be committed in both repos before
  building against `origin/main` ferx-core).

## Verification

```r
library(ferx)                        # must NOT print "built WITHOUT autodiff"
stopifnot(isTRUE(ferx:::ferx_rust_autodiff_enabled()))

# HMC path — expect saem_n_subjects_hmc > 0:
source(system.file("examples", "ex_warfarin_saem_hmc.R", package = "ferx"))
stopifnot(fit$saem_n_subjects_hmc > 0L)

# Run the full examples suite:
for (f in list.files(system.file("examples", package = "ferx"),
                     pattern = "^ex.*\\.R$", full.names = TRUE)) {
  message("Running: ", basename(f))
  source(f)
}
```

## Before merging

- [ ] Verified `ferx_rust_autodiff_enabled()` is `TRUE` on a capable machine.
- [ ] Verified `saem_n_subjects_hmc > 0` on `ex_warfarin_saem_hmc.R`.
- [ ] Full `inst/examples/*.R` suite runs without error on the autodiff build.
- [ ] Rebase on `main` if it advanced.
