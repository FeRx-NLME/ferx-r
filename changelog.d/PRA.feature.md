- **A multi-start fit reports what every start did** (#PRA,
  [ferx-core #2867](https://github.com/FeRx-NLME/ferx-core/issues/2867)). With
  `n_starts > 1`, the new `fit$starts` is a data frame with one row per start:
  `start` (0 is the model file's own initial values), `ofv`, `penalized`,
  `objective`, `converged`, `status`, `status_reasons`, `n_iterations`,
  `wall_time_secs`, `theta` (that start's final theta, named), `selected` (the
  start the fit returned), and `error` / `error_code` for a start that failed.
  `print()` shows the engine's one-line verdict under the OFV, for example
  `Multi-start: 8 starts, 8 converged; best OFV -280.1234 (start 3); 7 of 8
  converged starts within 1.0 of best; 1 converged above it: start 4 (+5.27)`.
  The rows survive `ferx_save_fit()` / `ferx_load_fit()`. A single-start fit has
  zero rows.
