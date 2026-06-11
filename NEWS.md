# ferx (development)

## New features

- **FREM covariate analysis** (`ferx_to_frem()`): transforms a base model and
  dataset into a Full Random Effects Model (FREM) that treats covariates as
  additional dependent variables. The extended omega block captures
  covariate-parameter relationships in a single fit, avoiding stepwise search.
  Returns the generated model/data paths plus covariate means, variances, and
  FREMTYPE mapping. Fit the result with `ferx_fit(result$model_path,
  result$data_path)`. (#194)

- `ferx_model_show()` now **syntax-highlights** `.ferx` files in colour-capable
  consoles: section headers (`[parameters]`, ...) in bold yellow, declaration
  keywords (`theta`, `omega`, `sigma`, `kappa`, ...) in cyan, and comments
  dimmed -- via the optional `cli` package. Non-colour contexts (files, pipes,
  `NO_COLOR`, or no `cli` installed) print the raw text unchanged. (#4)

- **Covariate screen** (`ferx_cov_screen()`): a quick, informal screen that
  correlates each declared covariate (from `fit$covtab`) with every
  parameter that has IIV -- against both the subject's individual parameter
  estimate and its ETA. Covariates are aggregated to one value per subject
  first (median for continuous, most-frequent level for categorical), and
  associations are reported as a signed Pearson correlation (continuous) or a
  correlation ratio (categorical), keeping only pairs above a threshold
  (default `|r| >= 0.2`). Intended to flag what is worth a formal covariate
  search, not as a covariate test itself.

- **Data-selection filtering** (`[data_selection]` block, `ferx_fit(ignore=)`,
  `ferx_selection()`): records can now be excluded from the analysis dataset at
  read time without modifying the CSV -- equivalent to NONMEM `$DATA IGNORE=` /
  `ACCEPT=`. Three entry points:

  - `[data_selection]` block in `.ferx` model files (keys `ignore`, `accept`,
    `ignore_subjects`).
  - `ferx_fit(model, data, ignore = "DV < 1.0", accept = ..., ignore_ids = ...)`
    passes conditions directly from R; conditions from both the model file and
    the R call are merged and deduplicated.
  - `ferx_selection(data, ignore = ..., accept = ..., ignore_ids = ...)` is a
    pure-R preview that returns a `ferx_data` S3 object you can inspect before
    fitting, or pass directly to `ferx_fit()` as the `data` argument.

  Exclusion counts are exposed on `fit$exclusions` (a list with
  `n_records_total`, `n_obs_excluded`, `n_dose_excluded`, `n_other_excluded`,
  `excluded_subject_ids`, `fired_ignore`, `fired_accept`).
  `print.ferx_fit()` shows a **DATA SELECTION** block when rules fired.
  `ferx_runlog()` includes an exclusion count line in the data summary.
  Exclusions survive `ferx_save_fit()` / `ferx_load_fit()` round-trips.
  The bundled `warfarin_data_selection` example demonstrates the feature.

- `ferx_selection_excluded(x)` is a new generic: called on a `ferx_data`
  object it returns the excluded rows (with a `.exclude_reason` column);
  called on a `ferx_fit` it re-reads the data file and marks records from
  excluded subjects.

- `ferx_columns(data)` prints the column headers of a NONMEM CSV dataset,
  grouped into required NONMEM columns (`ID`, `TIME`, `DV`, `EVID`, `AMT`,
  `CMT`), optional NONMEM columns (`RATE`, `MDV`, `II`, `SS`, `CENS`, `OCC`),
  and covariates / user-defined columns. Accepts a file path, a `ferx_fit`
  object (uses `fit$data_path`), or a `ferx_example()` list. Returns the
  column name vector invisibly.

- `ferx_runlog(fit)` produces a NONMEM-style `.lst` run summary: model file
  content, data summary (subject / observation counts, time range), INITIAL vs
  FINAL parameter table with SE and %RSE for every theta/omega/sigma, estimation
  settings (optimizer, max iterations, BLOQ method, NCA warm-start, random
  seeds, covariate columns present in the data), OFV / AIC / BIC with
  convergence flag, covariance-step condition number and eigenvalues, ETA/EPS
  shrinkage, Durbin-Watson autocorrelation, Shapiro-Wilk ETA normality, and
  final gradient with a convergence threshold check. Pass `verbose = FALSE` to
  capture the output as a character string.

- `ferx_runlog(fit, show_iterations = TRUE)` gains an **Iteration history**
  section when `optimizer_trace = TRUE` was used: per-iteration OFV, delta-OFV,
  and method-specific convergence metrics (GRAD_NORM / STEP_NORM for
  FOCE/FOCEI/BFGS; LM_LAMBDA + ACC for Gauss-Newton; COND_NLL + GAMMA +
  MH_ACCEPT for SAEM). Runs with more than 30 iterations are truncated to the
  first 10 and last 10. Set `show_iterations = FALSE` to suppress the section.

- `ferx_runlog_iters(fit)` is a new function that prints the complete
  untruncated per-iteration table. Accepts a `ferx_fit` object or a path to a
  trace CSV.

- `print.ferx_job(handle)` now shows a live trace snapshot (last 5 iterations)
  when the rstudio-backend handle is printed during a running fit.

- `fit$sdtab` gains a `CMT` column for multi-endpoint models (present whenever
  any observation row has `CMT != 1`). Use this column to split GOF plots by
  endpoint without rejoining the original dataset.

- `fit$sdtab` now carries the actual subject ID from the data (parsed as
  numeric where possible). Previously the ID column was a 1-based loop index
  that broke downstream joins when IDs were non-consecutive or non-numeric.
  Non-numeric IDs now trigger a `warning`-severity message.

- `ferx_fit` objects from file-based fits now carry `fit$model_text` (verbatim
  `.ferx` source), `fit$theta_init` / `fit$omega_init` / `fit$sigma_init`
  (optimizer starting values), `fit$obs_time_range`, `fit$final_gradient`,
  `fit$optimizer_label`, `fit$bloq_method_label`, `fit$n_starts`,
  `fit$inits_from_nca`, `fit$covariate_names`, and several reproducibility
  seed fields. These fields power `ferx_runlog()` and are preserved in
  `.fitrx` bundles.

## Bug fixes

- Shapiro-Wilk ETA-normality flags now fold into a single warning that lists
  every flagged ETA (with its p-value) instead of firing one warning per ETA.
  Both `fit$warnings` and the structured `eta_normality` warning shown by
  `ferx_warnings()` are affected (ferx-core#163).

- `ferx_runlog()`: theta names now resolve via `names(fit$theta)` (where
  `R/fit.R` stores them) instead of `fit$theta_names` (which is `NULL` by
  design after the R post-processing step). Fall-back chain:
  `names(fit$theta)` → `fit$theta_names` → `fit$model_structure$theta_names`
  → `THETA(i)`.

- `ferx_runlog()`: `model_text`, `inits_from_nca`, and seed fields
  (`multi_start_seed`, `saem_seed`, `sir_seed_used`, `is_seed`) now guard
  against `NA_character_` / `NA_real_` values that extendr emits for Rust
  `Option<T>::None`, preventing spurious "NA" entries in the run log.

- `ferx_runlog()`: gradient-tolerance line suppressed for derivative-free
  optimizers (BOBYQA, GN, SAEM) where a gradient tolerance is not applicable.

- `ferx_rust_fit()` (internal): `fit$model_text` was `NA` for file-based fits
  because the R binding's provenance block set `model_path` / `model_hash` but
  did not set `model_text`. Fixed by reading the model file in the same block.

## Build

- `src/Makevars`: the autodiff build now forces fat LTO with `override`. R CMD
  INSTALL runs make with `-e` semantics, so a stale
  `CARGO_PROFILE_RELEASE_LTO=thin` in the environment could previously leak past
  the makefile and silently build autodiff with thin LTO — which keeps the
  differentiated ferx-core functions in a separate module from their callers and
  breaks cross-crate Enzyme. `override` makes the makefile value authoritative
  regardless of the environment.
- `src/Makevars`: building with `FERX_NO_AUTODIFF=0` now preflights for the
  `enzyme` rustup toolchain and fails fast with an actionable message (how to
  install the toolchain, or how to fall back to finite-difference gradients)
  instead of a cryptic `toolchain 'enzyme' is not installed` error from cargo.
- `src/ad-preflight.sh` (new): the autodiff build now runs an *AD self-test*
  before the long fat-LTO link. Registering an `enzyme` toolchain is not enough
  — a standalone Enzyme plugin copied into a prebuilt nightly's sysroot passes
  the name check but then hangs forever in `llvm::Constant::getNullValue` the
  first time it differentiates anything. The self-test compiles and runs a real
  `#[autodiff_forward]` function under a portable timeout (no `timeout`/
  `gtimeout` dependency); a hang is turned into a clean failure that aborts the
  build with a pointer to the toolchain rebuild instructions, instead of leaving
  the user staring at a wedged `rustc`. The result is cached per toolchain
  fingerprint, so identical reinstalls skip it. Set `FERX_AD_PREFLIGHT_SKIP=1`
  to bypass, or `FERX_AD_PREFLIGHT_TIMEOUT=<seconds>` to adjust the budget.
  (ferx-r#111)

# ferx 0.1.5

## Documentation

- `?ferx_fit`: the `settings` parameter block is restructured into labelled
  sections, one per estimation method (Shared, FOCE/FOCEI/GN-hybrid, Trust-region,
  SAEM, Gauss-Newton, Importance Sampling, SIR, Multi-start). Each key now lists
  its default value and which methods accept it.

## New features

- `ferx_fit_async(model, data, ...)` now returns a `ferx_job` handle
  immediately so the R session stays free. Call `ferx_collect(handle)` to
  block-wait with live optimizer-trace progress; pass `verbose = FALSE` to
  suppress the display and just block until the result is ready. In RStudio
  the fit appears in the Jobs pane; elsewhere a `callr::r_bg()` background
  process is used. The returned `ferx_fit` object is identical to what
  `ferx_fit()` produces. **Breaking change from #91**: `ferx_fit_async()`
  previously blocked and returned the fit directly; it now returns a handle
  that must be passed to `ferx_collect()`.

- `print(fit)` has a new compact layout: a prominent `STATUS: CONVERGED` /
  `NOT CONVERGED` line with iteration count and wall time immediately after the
  header; OFV / AIC / BIC on one line; bold section headers with thin rules
  instead of `--- THETA Estimates ---` banners; shrinkage as a single line with
  inline `[!]` for values > 30%; a diagnostics line consolidating covariance
  status, condition number, and Durbin-Watson; and a colour-coded warning-count
  footer pointing at `ferx_warnings(fit)`. Programmatic access (`fit$theta`,
  `ferx_estimates()`, `summary()`) is unchanged.

- `ferx_warnings(fit)` pretty-prints fit warnings grouped by severity
  (critical / warning / info) with per-category remediation guidance.
  `ferx_warnings(fit, as_df = TRUE)` returns the underlying
  `fit$warnings_structured` data frame (columns: `severity`, `category`,
  `message`, `source_method`). Durbin-Watson autocorrelation guidance is
  direction-aware (positive vs negative DW) and suppresses the SDE hint
  when the model already uses a `[diffusion]` block.

- Default outer optimizer for FOCE / FOCEI changed from `slsqp` to `bobyqa`.
  BOBYQA is derivative-free (a quadratic trust-region) and re-evaluates the
  per-subject EBEs at every trial point, so it avoids the fixed-EBE FD gradient
  bias that drives SLSQP to local minima hundreds of OFV units above the true
  optimum on ODE / PD models, sparse data, and Emax-Hill identifiability
  problems. The default also flows to the FOCEI polish stage of
  `method = "gn_hybrid"` and to the polish stage of any
  `method = c(..., "focei")` chain. Pure SAEM and pure `gn` continue to ignore
  the optimizer setting. To restore the previous behaviour, pass
  `settings = list(optimizer = "slsqp")`. Requires a ferx-core build that
  includes the change.

- Log-transform-both-sides (LTBS) residual error: fit on the log scale with
  additive error, the equivalent of NONMEM's `Y = LOG(F) + EPS(1)`. Write the
  `[error_model]` block as either form:

  ```
  log(DV) ~ additive(ADD_LOG)   # DV on the natural scale; engine logs it
  DV ~ log_additive(ADD_LOG)    # DV already log-transformed in the data
  ```

  Under LTBS the fit output (`IPRED`, `PRED`, `CWRES`, `IWRES` in `sdtab`, and
  `DV_SIM` from `ferx_simulate()`) is on the log scale. `ferx_model_inspect()`
  reports the residual type as `additive (log-transformed)`. Requires a
  ferx-core build that includes the feature.

- `settings = list(reconverge_gradient_interval = N)` controls how often the
  FOCE/FOCEI population gradient re-solves each subject's inner EBE loop instead
  of holding it fixed. `0` (default) keeps the cheap fixed-EBE gradient; `1`
  reconverges every gradient evaluation; `N` reconverges every `N`-th. The
  fixed-EBE gradient can stall `slsqp` above the derivative-free (`bobyqa`)
  optimum on ill-conditioned non-IOV fits; reconverging recovers the full
  surface at ~5-6x the per-gradient cost. IOV models always reconverge and
  ignore the setting. Requires a ferx-core build that includes the option.

- Multi-endpoint (per-CMT) residual error models for simultaneous PK/PD
  fitting. A single `[error_model]` block can now assign a distinct error
  model to each observed compartment, dispatched by the dataset `CMT` column:

  ```
  [error_model]
    CMT=2: DV ~ proportional(PROP_ERR_PK)
    CMT=3: DV ~ additive(ADD_ERR_PD)
  ```

  Both endpoints contribute to one joint FOCEI/GN likelihood. ODE models only;
  supported with FOCE/FOCEI, Gauss-Newton, and SAEM.
  `ferx_model_inspect()` reports the per-CMT residual structure. New bundled
  example: `ferx_example("emax_pkpd")`.

- `[scaling]` block in `.ferx` model files for unit conversion. Form A
  (`obs_scale = <number>`) divides every model prediction by a scalar before
  residuals are computed and is safe with automatic differentiation. Form B
  (`obs_scale = <expression>`) and Form C (`y = <expr>` for ODE readout)
  support parameter and state-variable expressions but require
  `gradient = fd` in `[fit_options]`. Per-compartment variants
  (`obs_scale[CMT=N] = ...`, `y[CMT=N] = ...`) are also supported.
  New bundled example: `ferx_example("warfarin_scaled")`.

- Steady-state dosing via `SS` and `II` columns in the NONMEM CSV. Set
  `SS = 1` on a dose row and supply the dosing interval `II` (same time
  units as TIME). The engine resolves steady-state initial conditions
  analytically for 1/2/3-cpt models and via pulse expansion for ODE models.
  `SS = 2` adds the steady-state concentration to the current compartment
  state (superposition). No model-file changes are required.
  New bundled example: `ferx_example("warfarin_ss")`.

- SAEM HMC proposals: pass `settings = list(n_leapfrog = <int>)` to
  `ferx_fit()` to use Hamiltonian Monte Carlo proposals in the SAEM E-step
  (requires the Enzyme AD build). New output field `fit$saem_n_subjects_hmc`
  reports the number of subjects that used HMC proposals; `NULL` for MH-only
  or non-SAEM fits.

- SAEM fully supports inter-occasion variability (IOV / kappa) models.
  New bundled example: `ferx_example("warfarin_iov_saem")`.

- New bundled example `ferx_example("transit_2cpt")`: two-compartment ODE
  model with 3-transit-compartment absorption and allometric scaling.

- `ferx_fit()` accepts `"imp"` as a chained method
  (e.g. `method = c("focei", "imp")` or `method = c("saem", "imp")`).
  The terminal IMP stage runs an importance-sampling estimate of the
  marginal `-2 log L`, exposed on `fit$importance_sampling` (a list with
  `minus2_log_likelihood`, `mc_standard_error`, `n_samples`,
  `proposal_df`, `ess_min`/`ess_median`, `kappa_treatment`, and
  parallel `low_ess_subject_ids`/`low_ess_subject_frac` vectors).
  `print(fit)` and `summary(fit)` render the new block. New
  IMP-specific settings keys are recognized by `ferx_fit(settings = ...)`:
  `is_samples`, `is_proposal_df`, `is_seed`, `is_low_ess_threshold`.
  Requires ferx-core with importance-sampling support merged
  (FeRx-NLME/ferx-core IMP PR).

- New `stagnation_guard` key recognized by `ferx_fit(settings = ...)`.
  Pass `settings = list(stagnation_guard = FALSE)` to disable the NLopt
  outer-loop stagnation guard so SLSQP / L-BFGS run to their own xtol /
  ftol or to `maxiter` instead of short-circuiting on a numerically-flat
  OFV plateau. Useful for debugging or for problems with very slow but
  real OFV improvements below the guard's 1e-3 threshold. Consumed by
  FOCE / FOCEI / GN-hybrid only. Requires ferx-core with PR
  FeRx-NLME/ferx-core#62 merged.

## Bug fixes

- SAEM no longer collapses the random-effect variances (Omega) on sparsely
  sampled data. Previously, with few observations per subject, the
  between-subject variability could shrink toward zero during the first
  iterations while the residual error absorbed it (tiny `omega`, inflated
  additive `sigma`). A burn-in now holds Omega fixed while the MH sampler warms
  up, tunable via `settings = list(omega_burnin = <int>)` (default 20; `0`
  restores the previous behaviour). Requires the matching ferx-core update that
  adds
  the SAEM Omega burn-in.

- SIR confidence intervals now work correctly for models with `FIX`-ed
  parameters. Previously, any fixed parameter caused the proposal covariance
  to be singular, and SIR returned "All SIR samples had invalid weights".
  Sampling is now restricted to the free-parameter subspace and fixed
  parameters are held at their estimated values throughout. Requires
  ferx-core >= 0.1.0 (commit 47b48b5, ferx-core#64).

- All output functions now display the declared variable name (`ETA_CL`,
  `EPS_PROP`) rather than wrapping it in `OMEGA()`/`SIGMA()`. Affected
  surfaces: `print(fit)` OMEGA section and shrinkage, `ferx_estimates()`,
  `ferx_cor_matrix()` (via `fit$cov_matrix` dimnames), `fit$omega`
  row/column names, `fit$sir_ci_omega`, `fit$sir_ci_sigma`, and
  `summary(fit)` shrinkage. When names are absent the fallback remains
  the conventional numbered form: `OMEGA(1,1)`, `SIGMA(1)`. (#19)

## New features

- IWRES autocorrelation diagnostic: `fit$dw_statistic` (pooled Durbin-Watson)
  and `fit$iwres_lag1_r` (pooled lag-1 Pearson r) are now returned by
  `ferx_fit()`. A `--- Diagnostics ---` block is printed by `print(fit)` when
  the values are available; an actionable `message()` is emitted when DW < 1.5
  or DW > 2.5. The new `check_diagnostics(fit)` function returns a structured
  list with an `$autocorrelation` data frame and a tidy `$shrinkage` data frame
  covering both ETA and EPS components. Both fields round-trip through
  `ferx_save`/`ferx_load`; old `.fitrx` files deserialize with `NA`. Requires
  ferx-core ≥ 0.1.0 (commit 5653ddae, ferx-core#20). (#6)

- SDE support via Extended Kalman Filter: models with a `[diffusion]` block
  in the `.ferx` file now run through the EKF likelihood. `fit$uses_sde` is
  `TRUE` for these fits; diffusion variance parameters appear in `fit$theta`
  as `DIFF_<STATE>` (e.g. `DIFF_CENTRAL`). Autodiff is automatically forced
  to finite differences for SDE models; SAEM is not supported and raises a
  hard error. Requires ferx-core ≥ 0.1.0 (commit 03332951).

- Lag time parameter (`lagtime=` on the structural_model line, NONMEM-
  style `alag=` accepted as an alias) is now supported in `.ferx`
  models. Delays the effective start of every dose record by the
  parameter's value; defaults to `0.0` when omitted so existing models
  are unaffected. Random effects on lag time work the same as on any
  other PK parameter (`LAGTIME = TVLAGTIME * exp(ETA_LAGTIME)` for
  log-normal, or the additive form covered in the
  `parameter-transforms` vignette). Pairs with ferx-core#12.

- `ferx_sir(fit)` — run SIR (Sampling Importance Resampling) as a
  standalone post-fit step, without having to set `sir = TRUE` at fit
  time. Useful for expensive fits where you want to add SIR-based
  uncertainty after the fact, or when working with a fit loaded via
  `ferx_load_fit()`. `ferx_fit()` now records `model_path` /
  `data_path` and SHA-256 `model_hash` / `data_hash` on the fit; the
  hashes round-trip through `.fitrx` save/load and `ferx_sir()`
  refuses to run when either file has changed since the fit.

- `ferx_simulate_with_uncertainty()` — simulate observations while
  propagating population parameter uncertainty in addition to the usual
  between-subject (eta) and residual (eps) variability. For each
  parameter set drawn from the uncertainty distribution
  (`method = "asymptotic"` uses `fit$cov_matrix`; `method = "sir"`
  reuses SIR resamples) the standard simulator runs `n_sim_per_draw`
  replicates. Output is a long data frame with a leading `DRAW` column
  so downstream code can compute uncertainty-aware prediction bands.
  Requires `covariance = TRUE` for asymptotic mode; SIR mode requires
  `sir = TRUE` and `sir_keep_samples = TRUE` (passed via `settings`)
  at fit time. `ferx_fit()` now also exposes `sir_resamples`,
  `sir_resamples_n`, and `sir_resamples_dim` for downstream consumers.
  Pairs with ferx-core#7.

- `ferx_simulate()` output now includes a leading `DRAW` column
  (always `1` for non-uncertainty paths) for forward compatibility with
  `ferx_simulate_with_uncertainty()`. Downstream code that uses
  `subset(sim, ...)` or column selection by name is unaffected.

- `ferx_fit()` now returns `$gradient_used` — the inner-loop gradient
  method the engine actually used (`"ad"`, `"fd"`, or `"N/A"`). When
  `gradient = "auto"` it shows which branch resolved at fit time; ODE
  models silently fall back from AD to FD even when `gradient = "ad"`
  was requested, and this field surfaces that fallback without needing
  to read warnings. The raw engine labels are also exposed as
  `$gradient_method_inner` / `$gradient_method_outer`. `print()` shows
  both requested and used; `summary()` formats them as `auto
  (requested) -> ad (used)` (ferx-core#1).

## Breaking changes

- `ferx_fit()` no longer has dedicated `max_unconverged_frac` and
  `min_obs_for_convergence_check` arguments. These are estimation knobs and
  now flow through `settings`, like the other low-level fit options
  (#51):

  ```r
  # before
  ferx_fit(m, d, max_unconverged_frac = 0.1, min_obs_for_convergence_check = 2L)
  # after
  ferx_fit(m, d, settings = list(
    max_unconverged_frac = 0.1, min_obs_for_convergence_check = 2L
  ))
  ```

- `ferx_model()` argument order is now `ferx_model(data, model)` (data
  first). This enables the natural data-first pipe entry point
  `ex$data |> ferx_model(ex$model) |> ferx_fit()` (#81).

  Old-style positional calls (`ferx_model("pk.ferx")` or
  `ferx_model("pk.ferx", "data.csv")`) are detected by the `.ferx`
  extension on what is now the `data` argument and auto-corrected with a
  deprecation warning. The compatibility shim will be removed in a
  future release. Calls that pass `data` by name
  (`ferx_model("pk.ferx", data = "data.csv")`) keep working unchanged
  because R matches `data =` by name first and the remaining positional
  argument falls into the `model` slot.

## Bug fixes

- Bundled example `warfarin_additive_eta.ferx` used `tlag=TLAG` on its
  structural_model line. `tlag` was never a recognized PK parameter
  key in the engine, so the parser silently interpreted the value as
  the `cl=` parameter, overwriting the clearance value and producing
  incorrect fits. Updated to `lagtime=TLAG`. If you derived a local
  model from this example, change `tlag=` to `lagtime=` (or `alag=`).
  Pairs with ferx-core#12.

- `ferx_model_validate()` no longer flags `[initial_values]` as a missing
  required section. The block was removed from the engine in
  [ferx-core e5e934d](https://github.com/FeRx-NLME/ferx-core/commit/e5e934d)
  — theta / omega / sigma initial values are read from `[parameters]` and
  the parser silently ignores any leftover `[initial_values]` block. The
  R-side validator and `ferx_model_new()` templates hadn't been updated.
  Now: validator's `required_sections` drops `initial_values`, all five
  `ferx_model_new()` templates and every bundled `.ferx` example file
  emit the trimmed shape, and a regression test guards against the
  block creeping back into templates (#16).

- `ferx_set_section()` now applies copy-on-write when the underlying
  `ferx_model` points at a file inside the installed `ferx` package
  directory (e.g. a model returned by `ferx_example()`). The file is
  copied to `tempdir()` before editing and the returned `ferx_model`'s
  `$model` field is updated to the copy, preventing accidental
  modification of bundled examples. Plain-path callers are unaffected —
  passing a path string still edits in place (#80).

- `ferx_check_init()` now accepts a `ferx_model` as its first argument
  (in addition to a plain path), so it can be placed directly in a pipe:
  `ex$data |> ferx_model(ex$model) |> ferx_check_init()`. When a
  `ferx_model` is supplied and `data` is not, the data path on the
  object is used (#79).

## Documentation

- New vignette *"Editing ferx model files programmatically"* covering
  `ferx_model_new()`, `ferx_model_section()`, and `ferx_model_set_section()`:
  skeleton creation, section inspection, read-modify-write patterns, a
  console-only fit workflow, and overwrite-guard behaviour.

# ferx 0.1.2

## New features

- New `ferx_model_validate(path)` checks a `.ferx` file for syntax errors and
  missing required sections without running the optimizer. Prints a section
  presence report and returns `TRUE`/`FALSE` invisibly. Required sections are
  `[parameters]`, `[individual_parameters]`, `[structural_model]`,
  `[error_model]`, and `[initial_values]`; `[odes]` and `[fit_options]` are
  optional.

- `ferx_fit()` now returns `$eigenvalues` (sorted descending) and
  `$condition_number` (ratio of largest to smallest eigenvalue) for the
  covariance correlation matrix. Both are `NULL` when the covariance step was
  not run or failed. `condition_number = Inf` signals a non-positive eigenvalue.
  A warning is appended when `condition_number > 1000`. The condition number is
  shown on the `Covariance:` line in `print()` and `summary()` output.

- `ferx_model_inspect(path)` parses a `.ferx` file without fitting and prints
  a compact structural summary (model type, IIV, IOV, residual error). Pass a
  `ferx_fit` object instead of a path to inspect the structure post-fit without
  re-supplying the file.

- `ferx_fit()` now attaches `$model_structure` to every result: a named list
  with fields `theta_names`, `model_type`, `iiv`, `iov`, and `residual`.
  The same summary is shown in `print()` and `summary()` output.

- `ferx_model_section(path, section)` — extract and print the body of a named
  section from a `.ferx` model file; returns lines invisibly for scripted use.

- `ferx_model_set_section(path, section, lines)` — replace the body of a named
  section in-place; the write complement to `ferx_model_section()`.

## Documentation

- New vignette **"Model workflow: inspect, edit, fit"** (`vignette("model-workflow", package = "ferx")`)
  demonstrates the pre-fit inspection workflow: `ferx_model_inspect()` before
  `ferx_fit()` to catch structural mistakes cheaply, and re-inspecting the
  fitted result via `ferx_model_inspect(fit)` (closes #57).

## Changes to existing functions

- `ferx_fit()` now returns `$sigma_names` and `$sigma_types` (parallel to
  `$sigma`), and `print()` displays each sigma with its declared name, the
  derived variance (`sigma^2`), and — for proportional components — the
  CV% (`sigma * 100`). Sigma is on the standard-deviation scale for both
  proportional and additive components, matching the new YAML output added
  in ferx-core#57. Closes #59.

- `result$model_structure` is now sourced from the Rust engine (built from the
  parsed `CompiledModel`) instead of an R-side regex re-parse of the `.ferx`
  file (closes ferx-core#49). The shape is unchanged — `theta_names`,
  `model_type`, `iiv`, `iov`, `residual` — so `ferx_model_inspect()` callers
  see the same fields. `model_type` now distinguishes IV bolus from infusion
  (e.g. `"1-cpt IV infusion"`) and adds 3-cpt variants; the pre-fit
  `ferx_model_inspect(path)` parser was updated to the same label set so both
  the file-based and fit-based inspection paths report identical strings.

- `ferx_model_edit()` gains `overwrite = FALSE`. Previously the function
  silently skipped the file copy when the destination already existed; it now
  errors with a clear message. Callers that relied on the silent-skip must add
  `overwrite = TRUE`.

- `ferx_model_new()` gains `print = FALSE`. When `print = TRUE` the skeleton
  is printed to the console without writing any file or opening an editor;
  `path` becomes optional in that mode. Five templates are available:
  `"1cpt_oral"` (default), `"1cpt_iv"`, `"2cpt_oral"`, `"2cpt_iv"`, `"ode"`.

## Bug fixes

- `print.ferx_fit()` now uses the exact coefficient of variation formula for
  `EXP(OMEGA)` log-normal parameters: `sqrt(exp(omega) - 1) * 100` instead of
  the approximation `sqrt(omega) * 100` (doi:10.1002/psp4.12507). Applied only
  when `eta_param_types == "log_normal"` (defaults to log-normal when the field
  is absent). Display for logit, additive, and custom ETA types is deferred to
  #53.

- `ferx_model_section()`: fixed a descending-index bug where an empty section
  body (header immediately followed by another header) returned lines in reverse
  instead of `character(0)`.

- `ferx_model_set_section()`: fixed a last-section splice bug where
  `seq.int(end+1, length)` produced a descending sequence and appended `NA`
  plus a duplicate line when replacing the last section in a file.

- `ferx_estimates()`: `estimate_natural` is now `NA` when SE is unavailable,
  matching the documented contract that all natural-scale columns are `NA`
  when the covariance step was not run. Previously the back-transform was
  always populated for `log` and `logit` thetas regardless of SE.

- `.ferx_model_type()` (used by `ferx_model_inspect()`): now returns
  `"X-cpt IV infusion"` for `*_infusion` PK models, matching the label the
  Rust engine attaches post-fit. Previously the pre-fit label dropped the
  `IV` token.

- `print.ferx_fit()`: the logit-ETA `+/-1SD` summary line is now ASCII; the
  previous `±` rendered as `<U+00B1>` under non-UTF-8 locales.
