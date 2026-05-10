# Comprehensive ferx test plan

This document covers the remaining gaps in ferx test coverage after PR #76.
It is a reference for implementing tests — not the tests themselves.

---

## Build requirements per test area

There are two build tiers. Know which one you are on before running tests.

### Tier 1 — no-autodiff build (fast, ~2 min)
```bash
FERX_NO_AUTODIFF=1 R CMD INSTALL --no-lock .
```
Covers everything that does not require gradient-based fitting. This is
**Teun's normal dev build**. Most tests below run here.

### Tier 2 — Enzyme/autodiff build (slow, ~1.5 h, requires Enzyme toolchain)
```bash
R CMD INSTALL --no-lock .
```
Required only for tests that verify **gradient correctness** or compare
autodiff vs finite-difference results. These are marked **[ENZYME ONLY]**
below and must be written with a `skip_if` guard so they are inert on
Tier 1 machines.

**Guard to use for all [ENZYME ONLY] tests:**
```r
skip_if(
  !isTRUE(ferx_rust_autodiff_enabled()),
  "Enzyme autodiff not available — skipping gradient-correctness test"
)
```

`ferx_rust_autodiff_enabled()` is already exported and returns `TRUE`/`FALSE`.

---

## Fixtures

The existing `helper-warfarin-fit.R` provides `warfarin_fit()`, a cached
FOCEI fit with `maxiter = 30L` (shape checks only, not convergence quality).
New test files that need a live fit should call `warfarin_fit()` rather than
running their own fit — this avoids redundant Rust calls in CI.

If a test needs `covariance = TRUE` (for SE/cor_matrix tests), add a second
cached helper alongside the existing one:
```r
warfarin_fit_cov <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex  <- ferx_example("warfarin")
      fit <<- ferx_fit(ex$model, ex$data,
                       method = "focei", verbose = FALSE,
                       covariance = TRUE, settings = list(maxiter = 30L))
    }
    fit
  }
})
```
Add this to `helper-warfarin-fit.R` (not a new file) so it is available
across all test files automatically.

---

## Files to create

### `tests/testthat/test-example.R` — Tier 1, no Rust call needed

`ferx_example()` just resolves bundled file paths.

| # | Test | What to check |
|---|---|---|
| 1 | known name returns list with $model and $data | `is.list(ex)`, both names present |
| 2 | $model path exists on disk | `file.exists(ex$model)` |
| 3 | $data path exists on disk | `file.exists(ex$data)` |
| 4 | $model has .ferx extension | `grepl("\\.ferx$", ex$model)` |
| 5 | $data has .csv extension | `grepl("\\.csv$", ex$data)` |
| 6 | no-arg call prints names and returns NULL invisibly | `capture.output()` non-empty, `withVisible()$visible == FALSE` |
| 7 | unknown name errors with available names in message | `expect_error(..., regexp = "warfarin")` |

---

### `tests/testthat/test-fit.R` — Tier 1 (shape/type checks) + [ENZYME ONLY] (gradient)

Uses `warfarin_fit()` and `warfarin_fit_cov()` fixtures. Only the
`ferx_check_init()` tests need a direct `ferx_fit()` call.

**Return structure (Tier 1):**

| # | Test | What to check |
|---|---|---|
| 1 | returns class ferx_fit | `expect_s3_class(warfarin_fit(), "ferx_fit")` |
| 2 | $theta is a named numeric vector | `is.numeric`, `!is.null(names)` |
| 3 | $theta names match model theta declarations | compare to `ferx_model_inspect()$theta_names` |
| 4 | $omega is a square numeric matrix | `is.matrix`, `nrow == ncol` |
| 5 | $sigma is a finite positive numeric | `is.finite`, `> 0` |
| 6 | $ofv is finite numeric | `is.finite` |
| 7 | $aic > $ofv and $bic > $ofv | arithmetic sanity |
| 8 | $method equals requested method | `expect_equal(fit$method, "focei")` |
| 9 | $sdtab is a data frame | `is.data.frame` |
| 10 | $sdtab has required columns | ID, TIME, DV, PRED, IPRED, CWRES, IWRES all present |
| 11 | $sdtab has one row per observation | nrow matches data minus dose rows |
| 12 | $se_theta absent when covariance = FALSE | `is.null(warfarin_fit()$se_theta)` |
| 13 | $se_theta present when covariance = TRUE | `!is.null(warfarin_fit_cov()$se_theta)` |
| 14 | $se_theta length matches $theta length | when covariance = TRUE |
| 15 | errors on missing model file | `expect_error(ferx_fit("no_such.ferx", ...), "not found")` |
| 16 | errors on missing data file | `expect_error(ferx_fit(ex$model, "no_such.csv"), "not found")` |

**ferx_check_init() (Tier 1):**

| # | Test | What to check |
|---|---|---|
| 17 | returns invisibly without error on valid inputs | `expect_no_error`, `withVisible()$visible == FALSE` |
| 18 | errors on missing model file | `expect_error(..., "not found")` |
| 19 | errors on missing data file | `expect_error(..., "not found")` |

**ferx_cor_matrix() (Tier 1, requires covariance = TRUE fixture):**

| # | Test | What to check |
|---|---|---|
| 20 | returns a matrix with same dimnames as $omega | `dimnames` equal |
| 21 | diagonal is all 1s | `diag(cor) == 1` |
| 22 | off-diagonal values are in [-1, 1] | `all(abs(cor[off_diag]) <= 1)` |
| 23 | returns NULL or errors gracefully when covariance was FALSE | no hard crash |

**Gradient correctness [ENZYME ONLY]:**

| # | Test | What to check |
|---|---|---|
| 24 | autodiff OFV gradient matches finite-difference gradient within tolerance | requires two fits or a gradient inspection API — skip if `!ferx_rust_autodiff_enabled()` |

---

### `tests/testthat/test-simulate.R` — Tier 1

Uses `warfarin_fit()` for simulate-from-fit mode. Fresh simulate makes one
Rust call but is fast (no optimisation loop).

| # | Test | What to check |
|---|---|---|
| 1 | returns a data frame | `is.data.frame` |
| 2 | has required columns: SIM, ID, TIME, IPRED, DV_SIM | all present |
| 3 | n_sim = 3 → SIM column has exactly 3 unique values | `length(unique(sim$SIM)) == 3` |
| 4 | n_sim = 1 → nrow matches source data obs rows | |
| 5 | same seed produces identical output | call twice, `expect_equal` |
| 6 | different seeds produce different DV_SIM | `expect_false(identical(...))` |
| 7 | DV_SIM is finite numeric (no NAs) | `all(is.finite(sim$DV_SIM))` |
| 8 | simulate-from-fit: returns data frame with same columns | using `fit = warfarin_fit()` |
| 9 | simulate-from-fit: IPRED values are finite | |
| 10 | simulate-from-fit: different from population simulation | IPRED differs when using individual vs population params |
| 11 | errors on missing model file | `expect_error` |
| 12 | errors on missing data file | `expect_error` |

---

### `tests/testthat/test-predict.R` — Tier 1

Uses `warfarin_fit()` for predict-from-fit mode.

| # | Test | What to check |
|---|---|---|
| 1 | returns a data frame | `is.data.frame` |
| 2 | has required columns: ID, TIME, PRED | all present |
| 3 | nrow matches observation rows in data (dose rows excluded) | |
| 4 | PRED is finite numeric (no NAs) | `all(is.finite(pred$PRED))` |
| 5 | PRED matches $sdtab PRED column from a full fit | same model + data → same population predictions |
| 6 | predict-from-fit: returns data frame with same columns | using `fit = warfarin_fit()` |
| 7 | predict-from-fit: PRED values are finite | |
| 8 | predict-from-fit produces different PRED than population predict | individual parameters shift predictions |
| 9 | errors on missing model file | `expect_error` |
| 10 | errors on missing data file | `expect_error` |

---

### Expand `tests/testthat/test-map-estimates.R` — Tier 1

Currently only one `ferx_eta_cov()` test (error path). Add:

| # | Test | What to check |
|---|---|---|
| + | ferx_eta_cov() returns a data frame with ID + one col per ETA | column names include "ID" + all ETA names from `warfarin_fit()$ebe_etas` |
| + | ferx_eta_cov() values are finite numerics | `all(is.finite(...))` on non-ID columns |
| + | ferx_eta_cov() nrow equals number of unique subjects | |

---

## Implementation order

1. `test-example.R` — zero Rust, done in minutes, good warm-up
2. Expand `helper-warfarin-fit.R` with `warfarin_fit_cov()`
3. `test-fit.R` — write Tier 1 tests first (1–23), add the `[ENZYME ONLY]`
   test last with the skip guard; hand off to someone with the Enzyme
   toolchain to verify #24
4. `test-simulate.R`
5. `test-predict.R`
6. Expand `test-map-estimates.R`

Tests 1–3, 5–6 can be written and verified entirely on a no-autodiff build.
**Only test #24 in test-fit.R requires the Enzyme toolchain** — all others
run cleanly with `FERX_NO_AUTODIFF=1`.
