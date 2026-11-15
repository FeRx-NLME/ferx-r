- **Results and output that change with ferx-core `a959f4b6..4f3ce9b4`** (#PRA):

  - **Adaptive dosing: every decision is a record**
    ([ferx-core #2895](https://github.com/FeRx-NLME/ferx-core/issues/2895)).
    The interval that ends at a decision now runs under that decision's
    parameters, whether the controller doses or holds: the kappa of its own
    occasion under IOV, `TIME` at the decision, and a covariate row chosen by
    the new `[adaptive_dosing]` key `decision_covariates = locf | nocb`
    (default `locf`, the row the controller reads). Before, such an interval
    ran under the next *data* record. Trajectories, decisions and metrics of
    `ferx_simulate_adaptive()` move on layouts with no data record at some
    decision (by up to 2e-2 in core's fixtures); a layout with a data record on
    every decision, including every bundled example, does not change.
  - **Durbin-Watson values move under a `correlation ... by` line whose blocks
    follow one another in time**, such as `by OCC` over visits
    ([ferx-core #2858](https://github.com/FeRx-NLME/ferx-core/issues/2858)):
    `fit$dw_statistic` and `fit$iwres_lag1_r` now treat those blocks as one
    series, as occasions already were. No estimate or OFV changes.
  - **`method = "gn"` no longer reports `converged = TRUE` when its trust
    radius has collapsed**
    ([ferx-core #2738](https://github.com/FeRx-NLME/ferx-core/issues/2738)); the
    collapse warning now reads "... trust radius collapsed without
    convergence" and is critical.
  - **An `[error_model]` magnitude reads `[individual_parameters]` names**
    ([ferx-core #2915](https://github.com/FeRx-NLME/ferx-core/issues/2915)): a
    name assigned once without a random effect (`ISPK = 1`, then
    `DV ~ additive(S * ISPK)`) is substituted. A theta that
    `[individual_parameters]` reassigns under its own name
    (`CL = CL * exp(ETA_CL)`) now means the assigned value in the error model
    too, where it used to read the raw theta silently, so such a model's fit
    changes.
  - **`[output] MIXNUM` on a `[mixture]` model tables the class**
    ([ferx-core #2784](https://github.com/FeRx-NLME/ferx-core/issues/2784)):
    the simulated class, else the fitted `MIXEST`, as NONMEM's `$TABLE MIXNUM`
    does, where it used to echo a data column of that name. A data column
    named `MIXNUM` now draws `W_OUTPUT_MIXNUM_IS_CLASS`, and `MIXNUM` /
    `MIXEST` are refused as `[derived]` names.
  - **A leftover `[initial_values]` block is a warning, not an error**
    ([ferx-core #2913](https://github.com/FeRx-NLME/ferx-core/issues/2913)): the
    model parses with `W_DEPRECATED_BLOCK` and the block is ignored, as its
    contents always were. `E_DEPRECATED_BLOCK` no longer exists, and
    `ferx_model_validate()` reports such a model as valid.
  - **The multi-start info line is new**
    ([ferx-core #2867](https://github.com/FeRx-NLME/ferx-core/issues/2867)):
    every multi-start fit carries one `Multi-start: ...` line saying whether
    the starts agree. It replaces `Multi-start: best result from start k/n`,
    which was printed only when start 0 lost.
  - **`W_TTE_DEGENERATE_HAZARD` is reported once per subject and CMT**
    ([ferx-core #1648](https://github.com/FeRx-NLME/ferx-core/issues/1648)),
    with the number of replicates (or uncertainty draws) affected and which
    ones, instead of once per replicate. No simulated row changes.
  - **Finite-difference steps change** on a theta declared at
    `0 < |theta| < 1` with an O(1) effect (such as a covariate exponent near
    0), on `[covariate_nn]` weights and on a small sigma under M3 censoring
    ([ferx-core #2406](https://github.com/FeRx-NLME/ferx-core/issues/2406)).
    Gradients and FD standard errors on those coordinates are more accurate,
    so fits that reach them can stop at slightly different values.
