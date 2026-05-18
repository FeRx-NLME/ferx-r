# ferx-r Optimization Companion Plan (Updated 2026-05-18)

This plan is the ferx-r follow-up to `ferx-core/plans/optimization.md`.
Read that document first — the algorithmic work happens entirely in ferx-core.
ferx-r's role is to keep the R wrapper's API surface, documentation, and Rust
glue in sync with each ferx-core step as it lands.

**No new R functions or algorithms are introduced here.** Every item below is
one of: documentation, output-field surfacing in `src/rust/src/lib.rs`, or a
top-level argument addition to `ferx_fit()`.

Each section is keyed to the corresponding ferx-core step number so PRs can
be linked directly.

---

## Prerequisites

- Wait for the ferx-core step to be merged to `ferx-core/main` before opening
  the ferx-r follow-up PR. The `.cargo/config.toml` patch in
  `src/rust/.cargo/config.toml` auto-swaps in the local `../../../ferx-core`
  checkout, so you can develop and test against an unmerged branch locally —
  but the PR must target a merged ferx-core commit.
- Read `CLAUDE.md` in both repos before starting.
- Fill every section of `.github/PULL_REQUEST_TEMPLATE.md` before opening a PR.

---

## PR #22 prerequisite — fix stale default optimizer reference

**Linked to ferx-core PR #22 (`perf/cross-engine-bench-fixes`)**

PR #22 flips the ferx-core default outer optimizer from `"bobyqa"` to `"slsqp"`.
`R/fit.R` currently documents `"bobyqa"` as the default in two places:

- Line 128: `\code{"bobyqa"} (default)` in `@param settings`
- Line 417: `ferx_fit(m, d, settings = list(optimizer = "bobyqa"))  # default`

Update both to `"slsqp"`. No Rust glue changes needed — this is docs only.

**Files to touch:** `R/fit.R`

---

## Step 1 — No ferx-r work required

rayon parallelism in `likelihood.rs` and `gauss_newton.rs` is internal to
ferx-core. No API surface changes.

---

## Step 2 — No ferx-r work required

Log-Cholesky parameterization is already done in ferx-core. No API surface changes.

---

## Step 3 — No ferx-r work required

AD population gradient (`subject_nll_pop_grad`) is internal to ferx-core.
No new `FitOptions` fields, no new `FitResult` fields.

---

## Step 4 — No ferx-r work required

SAEM M-step AD gradient is internal. No API surface changes.

---

## Step 5 — Document `trust_region` optimizer (if not already done)

**Linked to ferx-core Step 5.**

`trust_region` is already accepted by ferx-core and documented in `R/fit.R`
(line 421 has the example). Verify the `@param settings` prose at line 128–134
is accurate after PR #22 merges (the default changes from bobyqa to slsqp).
No further ferx-r work expected.

**Files to touch:** `R/fit.R` (verify only, likely no change needed beyond PR #22 fix)

---

## Step 6 — Document adaptive `steihaug_max_iters`

**Linked to ferx-core Step 6.**

ferx-core changes `steihaug_max_iters: usize` (default 50) to
`steihaug_max_iters: Option<usize>` (default `None` = adaptive).

Update `R/fit.R` `@param settings` to document the new semantics:

```
\item{\code{steihaug_max_iters}}{Maximum Steihaug-CG iterations per
  trust-region outer step. Integer, or omit for the adaptive default
  (\code{sqrt(n_params)} near convergence, 5 far from it). Only relevant
  when \code{optimizer = "trust_region"}.}
```

Update the example at line 464–465 to show both forms:

```r
# Fixed budget
ferx_fit(m, d, settings = list(optimizer = "trust_region",
                               steihaug_max_iters = 50L))
# Adaptive (default — omit the key)
ferx_fit(m, d, settings = list(optimizer = "trust_region"))
```

No Rust glue changes. Settings key validation lives in Rust's `apply_fit_option`;
R passes the value as a string and Rust handles the `Option<usize>` parse.

**Files to touch:** `R/fit.R`

---

## Step 8 — Document `saem_n_leapfrog`

**Linked to ferx-core Step 8.**

ferx-core adds `saem_n_leapfrog: usize` (default 3) to `FitOptions`.

Add to `@param settings` in `R/fit.R`:

```
\item{\code{saem_n_leapfrog}}{Number of leapfrog steps per HMC proposal in
  the SAEM E-step. Integer ≥ 1, default 3. Higher values improve mixing for
  strongly correlated ETA posteriors at the cost of more gradient evaluations
  per proposal. Only used when the AD gradient path is available for the model
  (analytical PK families); falls back to Metropolis-Hastings otherwise.}
```

**Files to touch:** `R/fit.R`

---

## Step 9 — Document `sir_df`

**Linked to ferx-core Step 9.**

ferx-core adds `sir_df: f64` (default 5.0) to `FitOptions`.

Add to `@param settings` in `R/fit.R`:

```
\item{\code{sir_df}}{Degrees of freedom for the multivariate Student-t SIR
  proposal. Numeric, default 5 (following Dosne et al. 2017). Higher values
  approach the normal proposal; lower values give heavier tails and better
  coverage for skewed posteriors. Set \code{sir_df = 30} to reproduce the
  previous normal-proposal behaviour.}
```

Add a reference to Dosne 2017 in the `@references` section of `ferx_fit()` if
one does not already exist.

**Files to touch:** `R/fit.R`

---

## Step 10 — Surface `n_starts` / `start_sigma` + `start_index` output field

**Linked to ferx-core Step 10. This is the only step with Rust glue changes.**

### 10a — Document settings keys

Add to `@param settings` in `R/fit.F`:

```
\item{\code{n_starts}}{Number of independent optimization runs with perturbed
  initial estimates. Integer, default 1 (single run — no behaviour change).
  When > 1, runs are executed in parallel via rayon and the result with the
  lowest OFV is returned. Start 0 always uses the exact user-supplied initials;
  starts 1..n are log-space perturbations controlled by \code{start_sigma}.}
\item{\code{start_sigma}}{Log-space perturbation standard deviation for
  multi-start initial estimates. Numeric, default 0.3. Ignored when
  \code{n_starts = 1}.}
```

### 10b — Add `n_starts` as a top-level `ferx_fit()` argument

Multi-start is the primary robustness dial for users who hit local minima.
Expose it as a first-class argument alongside `covariance` and `optimizer_trace`,
not buried in `settings`:

```r
ferx_fit <- function(model, data, method = "focei",
                     covariance = FALSE,
                     n_starts = 1L,        # <-- new
                     optimizer_trace = FALSE,
                     settings = list(),
                     verbose = FALSE,
                     seed = NULL)
```

Merge `n_starts` into `settings` before forwarding to Rust, exactly as
`optimizer_trace` is merged at line 695–700 of `fit.R`:

```r
if (!identical(n_starts, 1L)) {
  settings <- c(list(n_starts = as.integer(n_starts)), settings)
}
```

Add validation: `n_starts` must be a positive integer.

Add a `@param n_starts` roxygen entry and examples:

```r
# Multi-start: 8 independent runs, best OFV wins
ferx_fit(m, d, n_starts = 8L)
```

### 10c — Surface `start_index` in Rust glue

ferx-core adds `start_index: usize` to `FitResult` (which start produced the
returned result; 0 = user initials). Surface it in `lib.rs`:

In `fit_result_to_list()` (around line 1271 where `uses_sde` is set):
```rust
start_index = result.start_index as i32,
```

In `default_fit_result()` (around line 884):
```rust
start_index: 0,
```

In `R/fit.F` `@returns`, add:
```
\item{start_index}{Integer (0-based). Index of the multi-start run that
  produced this result. 0 means the user-supplied initial estimates were
  used without perturbation. Always 0 when \code{n_starts = 1}.}
```

### 10d — Print `start_index` in `print.ferx_fit()`

If `fit$start_index > 0`, add a line to the printed summary:

```
Best start:        3 of 8 (starts 0–7)
```

Only print this line when `n_starts > 1` (i.e. when `fit$start_index` is
present and the fit was run with multi-start). Avoids noise for the common
single-start case.

**Files to touch:**
- `R/fit.R` (`ferx_fit()` signature, roxygen, validation, merging, `print.ferx_fit()`)
- `src/rust/src/lib.rs` (`fit_result_to_list`, `default_fit_result`)

---

## Completion Checklist

After all follow-up PRs:

- [ ] `ferx_fit()` `@param settings` documents all new keys:
      `saem_n_leapfrog`, `sir_df`, `n_starts`, `start_sigma`,
      `steihaug_max_iters` (adaptive semantics).
- [ ] Default optimizer in docs updated from `"bobyqa"` to `"slsqp"`.
- [ ] `n_starts` is a top-level `ferx_fit()` argument with validation and example.
- [ ] `start_index` is surfaced in `fit_result_to_list`, `default_fit_result`,
      `@returns`, and `print.ferx_fit()`.
- [ ] `cargo test` passes in `src/rust/`.
- [ ] `R CMD check` passes with no new warnings or notes.
- [ ] `pkgdown` site builds cleanly (`pkgdown::build_site()`).
