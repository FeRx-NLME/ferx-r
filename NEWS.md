# ferx (development)

## Documentation

- New vignette *"Editing ferx model files programmatically"* covering
  `ferx_model_new()`, `ferx_model_section()`, and `ferx_model_set_section()`:
  skeleton creation, section inspection, read-modify-write patterns, a
  console-only fit workflow, and overwrite-guard behaviour.

# ferx 0.1.2

## New features

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

## Changes to existing functions

- `ferx_fit()` now returns `$sigma_names` and `$sigma_types` (parallel to
  `$sigma`), and `print()` displays each sigma with its declared name, the
  derived variance (`sigma^2`), and — for proportional components — the
  CV% (`sigma * 100`). Sigma is on the standard-deviation scale for both
  proportional and additive components, matching the new YAML output added
  in ferx-nlme#57. Closes #59.

- `result$model_structure` is now sourced from the Rust engine (built from the
  parsed `CompiledModel`) instead of an R-side regex re-parse of the `.ferx`
  file (closes ferx-nlme#49). The shape is unchanged — `theta_names`,
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
