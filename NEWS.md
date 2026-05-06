# ferx 0.1.1

## New features

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
