# ferx-r feature-parity with ferx-core (May 2026)

Branch: `feat/ferx-core-feature-parity`
ferx-core PRs covered: #167 (IOV kappa shrinkage), #166 (unused-param warnings),
#165 (flexible FIX syntax), #155 (BOBYQA default — already documented, no code
change needed).

---

## Step 1 — Extract `shrinkage_kappa_by_occ` in Rust glue

**File:** `src/rust/src/lib.rs`  
**Where:** in `fit_result_to_list()`, after `shrinkage_kappa` extraction (~line 1184)

Build a data frame from `result.shrinkage_kappa_by_occ: Vec<Vec<f64>>`:
- Columns: `occ` (1-based i32), then one `f64` column per kappa name
- When the vec is empty, return `().into()` (same pattern as `ebe_kappas_df`)
- Add `shrinkage_kappa_by_occ = shrinkage_kappa_by_occ_df` to the `list!(...)` return

**Then rebuild:**
```bash
FERX_NO_AUTODIFF=1 R CMD INSTALL --no-lock .
```

- [x] Rust glue change written
- [x] Small build passes

---

## Step 2 — R-side: postprocess, print, roxygen, tests for `shrinkage_kappa_by_occ`

**File:** `R/fit.R`

a. **Postprocess** (in `.fitrx_postprocess()`, near line 1280): after the existing
   `shrinkage_kappa` naming block, coerce the new field — set kappa column names
   from `result$kappa_names`; set to `NULL` when glue returns a zero-row frame.

b. **Print** (in `print.ferx_fit()` IOV section, near line 1872): after the
   existing per-kappa shrinkage line, add a per-occasion table when
   `!is.null(x$shrinkage_kappa_by_occ)` and `nrow(x$shrinkage_kappa_by_occ) > 1`:
   ```
     Occasion 1: ETA_CL 12.3%  ETA_V  9.1%
     Occasion 2: ETA_CL 18.7%  ETA_V 11.4%
   ```

c. **Roxygen `@return`** (near line 351): add `\item{shrinkage_kappa_by_occ}{...}`

**File:** `tests/testthat/test-fit.R`

d. Add test: fit `warfarin_iov.ferx` (use `ferx_example("warfarin_iov")`), assert:
   - `fit$shrinkage_kappa_by_occ` is a data frame
   - has column `occ` (integer)
   - has one additional column per kappa name
   - all shrinkage values are in `[0, 1]`

**File:** `tests/testthat/test-persist.R`

e. Add assertion that `shrinkage_kappa_by_occ` survives a `ferx_save`/`ferx_load`
   round-trip (check column names and nrow are identical).

**Then roxygenize:**
```r
roxygen2::roxygenize()
```

- [x] Postprocess written
- [x] Print updated
- [x] Roxygen updated
- [x] test-fit.R test written
- [x] test-persist.R assertion added
- [x] man/ regenerated

---

## Step 3 — Expose 4 diagnostic fields in Rust glue + roxygen

**File:** `src/rust/src/lib.rs`

Add four lines to the `list!(...)` in `fit_result_to_list()` (after line 1373):
```rust
saem_mu_ref_m_step_evals_saved = result.saem_mu_ref_m_step_evals_saved.map(|x| x as f64),
covariance_n_evals_estimated   = result.covariance_n_evals_estimated.map(|x| x as f64),
n_threads_used                 = result.n_threads_used as i32,
nlopt_missing_algorithms       = result.nlopt_missing_algorithms.clone(),
```

**File:** `R/fit.R`

Add four `\item{}` entries to `@return` near the other diagnostic fields:
- `saem_mu_ref_m_step_evals_saved`: M-step evaluations saved by mu-referencing; `NULL` for non-SAEM fits
- `covariance_n_evals_estimated`: estimated OFV calls for the covariance step; `NULL` when not computed or not applicable
- `n_threads_used`: integer, actual parallel threads used by the engine
- `nlopt_missing_algorithms`: character vector of NLopt algorithm names unavailable on this platform; empty on most builds

**File:** `tests/testthat/test-fit.R`

Add assertions:
- `fit$n_threads_used` is a positive integer
- `fit$nlopt_missing_algorithms` is a character vector
- `fit$saem_mu_ref_m_step_evals_saved` is NULL or numeric
- `fit$covariance_n_evals_estimated` is NULL or numeric

**Rebuild:**
```bash
FERX_NO_AUTODIFF=1 R CMD INSTALL --no-lock .
```

- [x] Rust glue lines added
- [x] Roxygen `@return` updated
- [x] test-fit.R assertions added
- [x] Small build passes (same build as Step 1)
- [x] man/ regenerated

---

## Step 4 — `unused_parameter` warning guidance + test

**File:** `R/diagnostics.R`

In `.ferx_warning_guidance()` `switch()` (~line 435), add:
```r
unused_parameter = "A declared parameter is never referenced in [individual_parameters] or [error_model]. Remove it from [parameters] or complete the expression that uses it.",
```

**File:** `tests/testthat/test-warnings.R`

Add test (no real model fit — use the fake-fit pattern from lines 47-64):
- Construct a fake `ferx_fit` whose `warnings_structured` has one row with
  `category = "unused_parameter"`, `severity = "warning"`, `message = "TVCL is declared but never used"`
- Call `ferx_warnings(fit)` and assert captured output contains guidance text
  (e.g. `"Remove it from"`)

- [x] Guidance entry added
- [x] test-warnings.R test written

---

## Step 5 — Flexible FIX syntax documentation

**File:** `R/model.R`

In the `ferx_model()` roxygen `@details` DSL reference section, add or expand
the parameter syntax examples to show:

```
  theta CL(0.134, 0.001, 10.0) (FIX)     # FIX at end (traditional)
  theta CL (FIX) (0.134, 0.001, 10.0)    # FIX anywhere (accepted since v0.x)
  theta CL(0.1, 0.001)                    # lower-bound only, no upper bound
  omega ETA_CL ~ 0.07 (FIX)              # fixed omega
  omega ETA_CL (FIX) ~ 0.07              # flexible placement
```

Also add a sentence: "A warning is emitted when a parameter is declared in
\code{[parameters]} but never referenced downstream -- this usually indicates a
commented-out expression or a typo in the parameter name."

**Roxygenize:**
```r
roxygen2::roxygenize()
```

- [x] fit.R @section added (flexible FIX + unused_parameter note)
- [x] man/ regenerated

---

## Step 6 — Final check and commit

```bash
# ASCII check
python3 -c "import glob; [print(f) for f in glob.glob('R/*.R') if any(b > 127 for b in open(f,'rb').read())]"

# R CMD check
R CMD check --as-cran --no-manual .
```

- [ ] No non-ASCII in R/*.R
- [ ] R CMD check clean (or only pre-existing NOTEs)
- [ ] All new tests pass
- [ ] Commits made, PR opened

---

## Commit message template

```
feat(fit): surface shrinkage_kappa_by_occ, diagnostic fields, unused-param
guidance, flexible FIX docs

- Extract FitResult.shrinkage_kappa_by_occ (ferx-core PR #167) in Rust
  glue as a per-occasion data frame; show in print.ferx_fit IOV section
- Extract saem_mu_ref_m_step_evals_saved, covariance_n_evals_estimated,
  n_threads_used, nlopt_missing_algorithms from FitResult
- Add guidance text for unused_parameter warning category in
  .ferx_warning_guidance()
- Document flexible FIX placement syntax and lower-bound-only theta in
  ferx_model() @details
- Tests: shrinkage_kappa_by_occ structure (IOV fit), persist round-trip,
  new diagnostic fields, unused_parameter guidance output
```
