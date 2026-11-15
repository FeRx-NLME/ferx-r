- **Engine fixes and new warnings with ferx-core `a959f4b6..4f3ce9b4`** (#PRA):

  - An evaluation-only fit (`maxiter = 0`) now warns
    `W_THETA_INIT_OUTSIDE_BOUNDS` when a theta starts outside its own declared
    bounds, naming the bound and the value the OFV was scored at; it used to
    clamp silently. `gn_hybrid` and `vi` at `maxiter = 0` get the start checks
    every searching method gets
    ([ferx-core #2865](https://github.com/FeRx-NLME/ferx-core/issues/2865)).
  - The positive Durbin-Watson warning on a model whose occasions follow one
    another in time now suggests `correlation ar1(RHO) by ID`, a line whose
    block is the whole subject, instead of a line `by` the occasion column,
    which cannot correlate rows across occasions; a second warning reports
    autocorrelation between consecutive blocks
    ([ferx-core #2858](https://github.com/FeRx-NLME/ferx-core/issues/2858)).
  - A `[mixture]` fit that ends with two classes merged or a class empty warns
    `W_MIXTURE_MERGED_CLASSES`
    ([ferx-core #2907](https://github.com/FeRx-NLME/ferx-core/issues/2907)), and
    a theta used only in a mixing expression no longer draws `W_UNUSED_PARAM`
    ([ferx-core #2766](https://github.com/FeRx-NLME/ferx-core/issues/2766)).
  - Adaptive dosing: a run with `round_to` warns `W_ADAPTIVE_STEP_ABSORBED`
    when the dose grid absorbs the smallest `increase N%` at a dose a subject
    can be at, so that rule cannot raise the dose
    ([ferx-core #2810](https://github.com/FeRx-NLME/ferx-core/issues/2810)). A
    run under a time-varying covariate or IOV no longer fails with
    `frozen-schedule replay verification failed` on a compartment off row
    inside an `ADDL` train
    ([ferx-core #2811](https://github.com/FeRx-NLME/ferx-core/issues/2811)).
  - An `ADDL`-implied dose takes a route's own lag from the row that carried
    `ADDL`, as NONMEM does; with a covariate in the lag that changed along the
    train it used to arrive at the wrong time
    ([ferx-core #2367](https://github.com/FeRx-NLME/ferx-core/issues/2367)).
  - `ode_method = "auto"` no longer switches a stiff segment back to `rk45`
    mid-way
    ([ferx-core #2845](https://github.com/FeRx-NLME/ferx-core/issues/2845)), and
    `sqrt()` at a zero argument has a zero derivative on every gradient path
    instead of sending the subject to finite differences
    ([ferx-core #2892](https://github.com/FeRx-NLME/ferx-core/issues/2892)).
  - `[data] path` is optional, so a model whose data is always passed as
    `data =` parses and can be copied away from its dataset
    ([ferx-core #2914](https://github.com/FeRx-NLME/ferx-core/issues/2914)).
  - Faster fits: subjects with an exactly-zero EBE no longer rerun the cold
    inner multi-start
    ([ferx-core #2688](https://github.com/FeRx-NLME/ferx-core/issues/2688)), a
    trial the EBE guard is certain to reject stops solving subjects
    ([ferx-core #2709](https://github.com/FeRx-NLME/ferx-core/issues/2709)),
    linear ODE models with idle states skip them in the matrix exponential
    ([ferx-core #2760](https://github.com/FeRx-NLME/ferx-core/issues/2760)),
    closed-form event walks no longer rescan the dose list
    ([ferx-core #2754](https://github.com/FeRx-NLME/ferx-core/issues/2754)),
    `sqrt` in `[odes]` compiles to native code
    ([ferx-core #2844](https://github.com/FeRx-NLME/ferx-core/issues/2844)), and
    one mixture subject outside the analytic gradient no longer sends the whole
    gradient to finite differences
    ([ferx-core #2906](https://github.com/FeRx-NLME/ferx-core/issues/2906)).
    Estimates and OFVs are unchanged by these.
