# ferx (development)

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
