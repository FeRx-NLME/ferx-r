# ferx (development)

## New features

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
