# ferx 0.2.0

## Breaking changes

Public functions renamed for verb-clarity and naming consistency (part of the
API cleanup in #223; naming rule + hard-break policy decided in #224). Old
names are removed - no deprecation shims. Update calls as follows:

- `ferx_npde()` -> `ferx_calc_npde()`
- `ferx_selection()` -> `ferx_apply_selection()`
- `ferx_to_frem()` -> `ferx_model_to_frem()` (moves into the `ferx_model_*` family)
- `ferx_warnings()` -> `ferx_get_warnings()`
- `ferx_columns()` -> `ferx_get_columns()`
- `ferx_plot_trace(fit)` -> `plot(fit)` (#229). New S3 methods `plot.ferx_fit()`
  and `plot.ferx_job()` replace it; `plot.ferx_job()` plots the trace
  accumulated so far by an in-progress `ferx_fit_async()` job, not just a
  completed fit. FOCE/FOCEI traces now show the running-minimum OFV by
  default (`monotonic = TRUE`), since the raw per-evaluation trace includes
  rejected line-search trial steps that can transiently increase OFV; pass
  `monotonic = FALSE` for the raw trace.

`ferx_selection_excluded()` is removed. To retrieve excluded records, pass
`excluded = TRUE` to `ferx_apply_selection()`, which now also accepts a
`ferx_data` or `ferx_fit` object as its `data` argument:

```r
# before
sel  <- ferx_selection(data, ignore = "DV < 1")
excl <- ferx_selection_excluded(sel)
excl <- ferx_selection_excluded(fit)

# after
excl <- ferx_apply_selection(data, ignore = "DV < 1", excluded = TRUE)
excl <- ferx_apply_selection(fit, excluded = TRUE)
```

`ferx_cor_matrix()`, `ferx_estimates()`, and `ferx_eta_cov()` are removed and
replaced with fields computed automatically at the end of `ferx_fit()` (and
recomputed by `ferx_load_fit()`), part of the fit-accessor cleanup in #226:

- `ferx_cor_matrix(fit)` -> `fit$cor_matrix`
- `ferx_estimates(fit)` -> `fit$estimates`
- `ferx_eta_cov(fit, data)` -> `fit$eta_cov` (no longer takes a `data` argument;
  it is computed from the dataset used to fit the model)

```r
# before
ferx_cor_matrix(fit)
ferx_estimates(fit)
ferx_eta_cov(fit, read.csv(ex$data))

# after
fit$cor_matrix
fit$estimates
fit$eta_cov
```

## Added

- `ferx_stop()` terminates a background fit started by `ferx_fit_async()`
  without waiting for it to finish (#235). Previously the only way to stop a
  running job was to send a kill signal manually.
- `ferx_model_to_frem()` gains a `fit` argument (#239). Pass a `ferx_fit`
  result from fitting the base model and its theta/omega estimates seed the
  generated FREM model's PK theta inits and PK-PK omega block, so a subsequent
  fit of the FREM model warm-starts from converged parameters instead of the
  base model's declared inits. Optional; `NULL` (default) is unchanged
  behaviour.
- `ferx_model_inspect()` now reports covariate-selected residual error models as
  `covariate-selected (...)` in `model_structure$residual`, matching the new
  `[error_model]` `if/else` selector in ferx-core
  ([ferx-core #658](https://github.com/FeRx-NLME/ferx-core/issues/658)).
`ferx_model_new()` is removed (#231). Scaffolding a new model from a template
is now a mode of the `ferx_model()` constructor, selected by passing
`template =` (or `print = TRUE` to preview a skeleton without writing a file).
Unlike the old function, which returned the file path, scaffold mode returns a
`ferx_model` object, so it pipes straight into `ferx_fit()`. The output path
moves from the first positional argument to the named `path =` argument:

```r
# before
ferx_model_new("m.ferx", template = "1cpt_oral", edit = FALSE)
ferx_model_new(print = TRUE)

# after
ferx_model(template = "1cpt_oral", path = "m.ferx", edit = FALSE)
ferx_model(print = TRUE)

# scaffold + fit in one pipe
ferx_model(template = "1cpt_oral", path = "m.ferx", edit = FALSE) |>
  ferx_fit(data)
```

# ferx 0.1.6

## Added

- **Analytic Savic transit absorption is now available in model files** — a
  `pk one_cpt_transit(cl, v, n, mtt)` structural model: Savic transit-compartment
  absorption fed straight into a one-compartment disposition as a fast analytical
  closed form (exponential tilting; no ODE solve), with exact FOCE/FOCEI
  sensitivities and a continuous, estimable number of transit compartments `N`
  (via ferx-core #611 / #386). It is the closed-form counterpart to the ODE
  `transit(n, mtt)` input rate (`ferx_example("transit_savic")`) and is much
  faster; with `N = 0` it reduces to first-order oral absorption. New bundled
  example `ferx_example("one_cpt_transit")` with a runnable
  `inst/examples/ex_one_cpt_transit.R`.

- **State-reactive (adaptive / feedback) dosing simulation** — `ferx_simulate_adaptive()`
  runs a forward simulation whose dosing regimen is decided at run time by the
  model file's `[adaptive_dosing]` block: a declarative first-matching-rule
  controller that titrates the next dose from the simulated (optionally
  assay-noised) trough at each decision time (via ferx-core #585, epic #391). The
  base subjects are dose-free — the controller supplies every dose. Returns the
  concentration trajectories, the realized dose ledger, the per-decision log
  (including holds), and per-subject outcome metrics — cumulative dose,
  realized dose-change counts, holds, discontinuation, the observed-signal
  summary, and the fraction of monitored values inside the model's
  `target_window` when one is declared (metrics via ferx-core #605); the
  frozen-schedule replay verifier runs on every replicate. New bundled example
  `ferx_example("adaptive_tdm")` (a vancomycin-style TDM trough titration) with a
  runnable `inst/examples/ex_adaptive_tdm.R`.
- **Joint PK-TTE (drug-driven hazard) is now available in model files** — an
  `[event_model] hazard = <expr>` that references the ODE PK state (e.g.
  `H0 * exp(BETA * (central / V))`) is accumulated as a cumulative-hazard ODE
  compartment and estimated jointly with the PK by FOCEI/SAEM, with shared random
  effects (via ferx-core #564). Mutually exclusive with the analytic `family`
  hazard; requires an ODE model. Validated three-way (ferx vs NONMEM vs nlmixr2).
  New bundled example `ferx_example("pktte_joint")` with a runnable
  `inst/examples/ex_pktte_joint.R` (`ferx_fit()` + `ferx_predict_survival()`).
- **Joint PK-TTE event-time simulation** — `ferx_simulate()` gains a `horizon`
  argument and now samples drug-driven (ODE-accumulated) time-to-event endpoints
  (via ferx-core #564, Slice 2.2). With a finite `horizon`, a joint PK-TTE model
  yields, per subject, its continuous PK rows plus a TTE row on the event CMT
  carrying the sampled event/censor `TIME` and an `OBSERVED` flag (1 = event
  before the horizon, 0 = right-censored at it; `NA` for continuous rows). The
  simulation output gains `CMT` and `OBSERVED` columns.
- **Zero-order absorption is now available in model files** — the built-in
  `zero_order(dur)` `[odes]` input rate (a constant-rate / modeled-duration input,
  NONMEM `RATE=-2`/`D1`), via the ferx-core update (ferx-core #504). Two new
  bundled examples: `ferx_example("zero_order_absorption")` (constant-rate input
  into central) and `ferx_example("sequential_absorption")` (zero-order fill of a
  depot, then first-order `ka` to central).
- **Biphasic / parallel absorption in model files** — an `[odes]` input-rate term
  can now be scaled by a declared pathway fraction (`FR*igd(...)`) and more than one
  term can feed a compartment, so the Freijer & Post biphasic inverse-Gaussian model
  is `d/dt(central) = FR1*igd(...) + FR2*igd(...)` (via ferx-core #388). New bundled
  example `ferx_example("biphasic_igd_absorption")`. Validated against a NONMEM
  `$DES` biphasic run (ferx FOCEI objective vs `#OBJV` to ~1e-5).
- **Parallel / mixed dual-pathway absorption in model files** — the new built-in
  `first_order(ka)` `[odes]` input rate (classic first-order / Bateman absorption,
  exposed as a composable input rate) plus a pathway fraction on `zero_order(...)`
  let two absorption pathways be split by a dose fraction: `parallel`
  (`FR1*first_order(ka=KA1) + FR2*first_order(ka=KA2)`) and `mixed`
  (`FZO1*first_order(ka=KA) + FZO*zero_order(dur=DUR)`) (via ferx-core #505). Two
  new bundled examples: `ferx_example("parallel_absorption")` and
  `ferx_example("mixed_absorption")`. Validated against NONMEM `$DES` runs (ferx
  FOCEI objective vs `#OBJV` to ~1e-5 parallel / ~1e-4 mixed).
- **`ferx_predict_survival()`** — survival-function predictions (`S(t)`, `H(t)`,
  `h(t)`, plus median and mean survival) on a user-supplied time grid for
  `[event_model]` (time-to-event) endpoints, for every subject and TTE CMT.
  Mirrors `ferx_predict()`; optionally uses a `fit`'s estimated `theta`. For
  competing risks (multiple TTE CMTs) it also returns the cause-specific
  cumulative incidence `cif` and all-cause survival `survival_all`, with
  `sum(cif) + survival_all = 1` (ferx-core #501).
- **Time-to-event support is now compiled into the package.** The Rust backend
  enables ferx-core's `survival` feature by default, so `[event_model]` blocks
  and the TTE datareader routing are active in shipped builds (previously the
  feature was off, making that routing a no-op).
- **Bundled examples for time-to-event and Savic transit absorption.** New
  `ferx_example()` models, each with a runnable `inst/examples/ex_*.R` script:
  `tte_exponential`, `tte_weibull`, `tte_gompertz`, and `tte_competing_risks`
  (standalone `[event_model]` TTE, paired with `ferx_predict_survival()`), plus
  `transit_savic` (Savic transit-compartment absorption via the built-in
  `transit(n, mtt)` input rate).
- **New `outer_xtol` / `outer_ftol` fit settings** — expose the derivative-free
  `bobyqa` outer optimizer's step / objective stop tolerances (NLopt
  `xtol_rel` / `ftol_rel`), settable via `ferx_fit(settings = list(...))` or the
  model file's `[fit_options]`. `outer_ftol` defaults to an automatic per-model
  value (tighter for time-to-event, where the objective is exact). See
  `?ferx_fit` (ferx-core #469).

## Breaking changes

- **Importance-sampling result/settings names use the `imp_*` prefix.**
  `fit$is_seed` is now `fit$imp_seed`, and IMP settings now use
  `imp_samples`, `imp_proposal_df`, `imp_seed`, and
  `imp_low_ess_threshold`. The short-lived `is_*` names are no longer accepted,
  matching ferx-core (FeRx-NLME/ferx-core#422).

## Fixed

- **Time-to-event frailty variance now matches NONMEM / nlmixr2.** A
  weakly-identified `omega^2` on a *nonlinear* hazard parameter (e.g. a Weibull
  shape frailty) previously read high because the derivative-free outer optimizer
  stopped short on the near-flat objective ridge. It now converges onto the
  NONMEM LAPLACIAN / nlmixr2 FOCEI consensus (the reference Weibull dataset moves
  `omega^2` 0.204 → 0.176). Automatic for `[event_model]` fits — no model change
  needed (ferx-core #469).
- **`ferx_fit()` no longer overrides model-file `[fit_options]` with accepted
  defaults.** `covariance`, `verbose`, `mu_referencing`, `sir`, and `gradient`
  now default to `NULL`, meaning "use the model file's value" (falling back to
  the engine default when the model file is silent). Previously their non-`NULL`
  R defaults (e.g. `covariance = TRUE`) silently overrode a model file that set
  the opposite — so a model with `covariance = false` still ran the covariance
  step. Pass the argument explicitly to override the model file
  (FeRx-NLME/ferx-core#558).
- **`ferx_fit()` no longer overrides the model file's estimation method.**
  `method` now defaults to `NULL`, meaning "use the `[fit_options] method` from
  the model file" (falling back to FOCEI only when the model file sets none).
  Previously the R-side default `method = "focei"` silently overrode a model
  file that specified e.g. `method = saem`. Pass `method` explicitly to override
  the model file as before (FeRx-NLME/ferx-core#558).
- **`ferx_selection()` preview recognizes the bare `ignore = C` shorthand**
  (NONMEM `IGNORE=C`). The pure-R preview parser previously returned no match for
  an operator-less clause, so the preview reported zero exclusions while the Rust
  fit dropped the flagged comment rows. A lone column name now expands to
  `C == C` (case preserved on the value to match the raw cell), `Inf`/`-Inf`
  now join `NaN` in being treated as label strings rather than numeric values,
  and an ordered comparison against a non-numeric value (e.g. `BW < abc`) yields
  no exclusion instead of a lexical string compare - all matching ferx-core
  (FeRx-NLME/ferx-core#536).

- **`ferx_to_frem()` now warns about estimated parameters with no random effect**
  and carries the base model's scaling over to the FREM model. A non-fixed
  parameter without an `ETA` is estimated poorly by IMP/IMPMAP (the
  importance-weighted M-step is biased for weakly-identified fixed effects), so
  `ferx_to_frem()` now emits a warning at conversion time recommending an `ETA` be
  added (ferx mu-references automatically), the parameter be held fixed, or FOCEI
  be used. The base model's `[scaling]`/`[odes]` blocks (e.g. `obs_scale`) are also
  now transferred to the generated FREM model instead of being dropped — dropping
  `obs_scale` rescaled every prediction and collapsed a PK typical value during
  FREM fits. Requires ferx-core with these fixes
  (FeRx-NLME/ferx-core#406, #407).
  
## Changed

- **`optimizer` now defaults to `"auto"`** (FeRx-NLME/ferx-core#490). The new
  `"auto"` choice picks the population optimizer per model: `"nlopt_lbfgs"` when
  the exact analytic FOCE/FOCEI gradient is available, and `"bobyqa"` when only
  finite differences are. Pass `settings = list(optimizer = "auto")` (or omit it)
  to get the automatic choice; the fit reports the resolved optimizer as
  `"auto (<resolved>)"`. Set `optimizer = "bobyqa"` for the previous fixed
  default.
- **No special Rust toolchain is needed to build.** ferx-core now uses
  hand-rolled analytic sensitivities (FeRx-NLME/ferx-core#381), so the package
  builds with the stable Rust toolchain. The former custom-toolchain build
  switches and preflight check are gone. The legacy build-mode probe is retained
  but now always returns `FALSE`. Gradients are unchanged (exact analytic
  sensitivities). Also bumps the bundled `nalgebra` to 0.35 to match ferx-core.

- **`method = "imp"` is now an estimator by default** (NONMEM `METHOD=IMP`): it
  updates the population parameters by importance-sampling Monte-Carlo EM instead
  of only evaluating the marginal `-2 log L` at fixed parameters. **Breaking:**
  calls that used `method = "imp"` (or `c("focei", "imp")`) purely to *score* a
  fit now re-estimate — pass `settings = list(imp_eval_only = TRUE)` (NONMEM
  `EONLY=1`) to recover the old evaluation-only behaviour. New `settings`:
  `imp_iterations`, `imp_averaging`, `imp_eval_only`; `imp_proposal_df` now also
  accepts `"normal"`/`"mvn"`. The estimating `"imp"` may lead or sit mid-chain;
  the evaluation-only `"imp"` must still be terminal. Plain `"imp"` is fragile on
  rich data (warm-start with `c("focei", "imp")`, or use `"impmap"`). Requires
  ferx-core with the `METHOD=IMP` estimator (FeRx-NLME/ferx-core#402). (#181)

## Performance

- **`ferx_fit()` no longer pays a ~100 ms per-call latency floor.** The R-interrupt
  poll loop in the Rust binding slept a fixed 100 ms between checks, so any fit
  that finished in between (single-subject MAP/posthoc, small datasets, quick
  refits) still took ~0.1 s of wall time regardless of the engine's actual
  runtime. The worker now signals completion on a channel, so the call returns
  the instant the fit finishes; `POLL_MS` bounds only Ctrl-C latency. A
  single-subject fixed-parameter MAP fit drops from ~0.118 s to ~0.005–0.009 s
  (~13–24×); estimates and interrupt behaviour are unchanged. (#178)

## New features

- **Weibull absorption — `weibull(td, beta)`**: a new built-in absorption input
  rate for `[odes]` models, alongside `transit(...)` and `igd(...)`. It adds a
  Weibull absorption-time distribution — scale `Td`, shape `beta` — fed straight
  into the central compartment, modelling the entire absorption delay in one term
  (no first-order `ka`). The shape selects the profile: `beta > 1` a delayed
  interior peak, `beta = 1` first-order absorption (`ka = 1/Td`), `beta < 1` fast
  early uptake. The dose feeds the density over time (∫ R_in dt = F·Dose), not as a
  bolus, exactly like `transit(...)`/`igd(...)`, and drives exact analytic
  FOCE/FOCEI/Bayes gradients. New example `weibull_absorption`. Anchored against a
  NONMEM `$DES` Weibull run (ferx FOCEI matches NONMEM `#OBJV` to ~1e-6). Requires
  ferx-core with FeRx-NLME/ferx-core#497.
- **IIV on residual error** (`iiv_on_ruv`): a `.ferx` model can now place a random
  effect on the residual error, matching NONMEM `Y = IPRED + EPS*EXP(ETA)`. Declare
  an `omega` and reference it from `[error_model]` with `iiv_on_ruv = NAME`; each
  subject then gets a log-normally scaled residual SD. Supported under FOCEI, IMP,
  IMPMAP, and SAEM. Validated against NONMEM 7.5.1 (ΔOFV 0.017). Requires ferx-core
  with this feature (FeRx-NLME/ferx-core#409).
- **M3 LOQ censoring supports upper limits**: datasets may now use
  `CENS = -1` to mark observations censored above an upper limit of
  quantification, with `DV` carrying the ULOQ value. Existing `CENS = 1`
  lower-limit handling is unchanged. Requires ferx-core with
  FeRx-NLME/ferx-core#416.
- **Inverse-Gaussian (Freijer & Post) absorption — `igd(mat, cv2)`**: a new
  built-in absorption input rate for `[odes]` models, alongside `transit(...)`.
  It adds an inverse-Gaussian absorption-time distribution — mean absorption time
  `MAT`, relative dispersion `CV2` (= Var/mean²) — fed straight into the central
  compartment, modelling the entire absorption delay in one term (no first-order
  `ka`). The dose feeds the density over time (∫ R_in dt = F·Dose), not as a
  bolus, exactly like `transit(...)`. New example `igd_inverse_gaussian`. Anchored
  against a NONMEM `$DES` inverse-Gaussian run. (Requires ferx-core with
  FeRx-NLME/ferx-core#347; the biphasic Freijer sum-of-two is a planned follow-up,
  FeRx-NLME/ferx-core#388.)

- **FREM covariate analysis** (`ferx_to_frem()`): transforms a base model and
  dataset into a Full Random Effects Model (FREM) that treats covariates as
  additional dependent variables. The extended omega block captures
  covariate-parameter relationships in a single fit, avoiding stepwise search.
  Covariates (and their continuous/categorical kind) are taken from the model's
  `[covariates]` block; the `covariates` argument is an optional subset filter
  to FREM only some of them. Returns a `ferx_model` referencing the generated
  model and data files, so it composes directly: `ferx_fit(ferx_to_frem(...))`.
  (#194)

- **IMPMAP estimator**: `ferx_fit(..., method = "impmap")` (alias
  `"importance_sampling_map"`) runs the NONMEM `METHOD=IMPMAP` Monte-Carlo EM
  estimator — importance sampling assisted by mode-a-posteriori re-centering.
  Runs standalone or as a chain stage (`c("focei", "impmap")`). Tuned via
  `settings` keys `impmap_iterations`, `impmap_samples`, `impmap_proposal_df`
  (`"normal"` for the MVN proposal, or a Student-t DoF), `impmap_averaging`,
  `impmap_seed`, `impmap_low_ess_threshold`. Requires a mu-referenced
  parameterization; IOV is not yet supported. Needs a ferx-core that provides
  the `impmap` method (separate `Cargo.lock` bump). (ferx-core #270)
- **Modeled infusion duration (`RATE = -2`)**: a NONMEM
  `RATE = -2` dose now infuses `AMT` over a *modeled* duration — declare an
  individual parameter `D{n}` for the dose compartment `n` and ferx infuses at
  rate `AMT / D{n}`, resolved per iteration and occasion (so it carries
  covariate and IOV effects), on **both** the analytical `pk(...)` engine and
  `ode(...)` models. Composes with `F{n}` and `ALAG{n}`, steady state,
  multi-dose, and system resets. A `RATE = -2` dose with no matching `D{n}`
  parameter is a clear error rather than a silent
  bolus, and a `D{n}` that is non-positive at the initial estimate is flagged.
  Handled entirely in the data reader and model parser, so no R-side change is
  needed. (Requires ferx-core with FeRx-NLME/ferx-core#384.)

- **Modeled infusion rate (`RATE = -1`)**: a NONMEM `RATE = -1` dose now infuses
  `AMT` at a *modeled* rate — declare an individual parameter `R{n}` for the dose
  compartment `n` and ferx infuses at rate `R{n}` (duration `AMT / R{n}`),
  resolved per iteration and occasion (so it carries covariate and IOV effects).
  The mirror of the modeled-duration `RATE = -2`, supported on **both** the
  analytical `pk(...)` engine and `ode(...)` models. Composes with `F{n}` and
  `ALAG{n}`, steady state, multi-dose, and system resets. A `RATE = -1` dose with
  no matching `R{n}` parameter is a clear error rather than a silent bolus, and an
  `R{n}` that is non-positive at the initial estimate is flagged. Handled entirely
  in the data reader and model parser, so no R-side change is needed. (Requires
  ferx-core with FeRx-NLME/ferx-core#418.)

- **`ferx_npde(fit, nsim, seed)`**: compute simulation-based NPDE (Normalized
  Prediction Distribution Errors, decorrelated within subject) and NPD
  (Normalized Prediction Discrepancies) post-hoc from an existing fit, without
  re-running `ferx_fit()`. Useful when a model was fitted without
  `[fit_options] npde_nsim`. Returns the `fit` with `NPDE`/`NPD` columns added
  to `fit$sdtab`, so `ferx_xpose()` and goodness-of-fit plots pick them up
  automatically; model/data default to the paths recorded on the fit.
  (ferx-r #172, requires ferx-core #377)
  
- **Bayesian estimation (`method = "bayes"`)**: full MCMC posterior sampling
  (Gibbs-within-HMC, NONMEM `METHOD=BAYES` parity). Returns posterior means
  with 95% credible intervals and convergence diagnostics (split-R-hat, ESS) on
  `fit$bayes` instead of a point estimate; `print()` shows a posterior-summary
  table. Tuning via `settings = list(bayes_warmup=, bayes_iters=, bayes_chains=,
  bayes_thin=, bayes_seed=)`. Supports BSV and zero-mean inter-occasion
  variability (per-occasion `kappa`; the IOV variance posterior appears as
  `OMEGA_IOV(...)`). Validated against FOCEI and NONMEM `METHOD=BAYES` on
  warfarin (ferx-core #380).

- **`ode_template` — generate the disposition ODE**: `ode_template NAME(...)`
  in `[structural_model]` writes the standard disposition ODE for a named model
  (`one`/`two`/`three_cpt` `iv`/`oral`) for you — the same states, micro-constant
  RHS, and `obs_scale` the analytical `pk NAME(...)` uses, but as an explicit ODE
  you can extend. It takes the same parameters as `pk NAME(...)` (including `ka`
  for oral routes). Re-declaring a `d/dt(X)` in `[odes]` **overrides** the
  generated equation for compartment `X` (undeclared compartments keep theirs) —
  the standard way to attach a built-in absorption input such as `transit(...)`.
  Combining an ODE-only absorption function with an analytical `pk NAME(...)` is
  now a clear error pointing at `ode_template`, never a silent conversion. New
  example `two_cpt_oral_cov_ode_template` (verified identical to its analytical
  and hand-ODE siblings in `test-ode-analytical-equivalence.R`). (Requires
  ferx-core with FeRx-NLME/ferx-core#363.)

- **Xpose interoperability**: `ferx_xpose(fit)` turns a fit into a ready-to-use
  Xpose object in memory (no NONMEM table files written to disk), so all
  downstream Xpose goodness-of-fit, covariate, and parameter diagnostics work
  out-of-the-box. Supports both the modern tidyverse `xpose` package
  (`backend = "xpose"`, default) and the classic S4 `xpose4`
  (`backend = "xpose4"`). Continuous vs categorical covariates are split using
  the model's `[covariates]` types, overridable via the `continuous` /
  `categorical` arguments. `RES`/`IRES` are derived and `WRES` is `NA` (ferx
  does not compute the FO-weighted residual). The estimation-iteration trace is
  not populated, so `xpose::prm_vs_iteration()` / `grd_vs_iteration()` are not
  supported (pending an engine change); use `ferx_plot_trace()` for OFV over
  iterations. When the fit carries simulation-based `NPDE`/`NPD` columns (from
  `[fit_options] npde_nsim > 0`), they are mapped to the Xpose residual role, so
  residual diagnostics (e.g. `xpose::res_vs_idv(xpdb, res = "NPDE")`) work on
  them out-of-the-box. (ferx-r #165)

- **Configurable ODE solver tolerance**: ODE models accept `ode_reltol`
  (default `1e-4`), `ode_abstol` (default `1e-6`), and `ode_max_steps`
  (default `10000`) in the model file's `[fit_options]` block or via
  `ferx_fit(settings = list(ode_reltol = ...))`. Defaults are unchanged, so
  existing fits are unaffected. PRED reproduces the analytical closed form to
  about `1e-4`, but the FOCE objective amplifies solver error, so the OFV of an
  ODE-form model could differ from its analytical equivalent by several units;
  a tighter `ode_reltol` lets the two agree. The shipped `*_ode` examples now
  set `ode_reltol = 1e-10`, and `test-ode-analytical-equivalence.R` checks the
  OFV agrees within a tolerance band in addition to PRED. (Requires ferx-core
  with FeRx-NLME/ferx-core#334.)

- **Standard PK models in ODE form**: every standard analytical model
  (`one_cpt_iv`, one-compartment oral = `warfarin`, `two_cpt_iv`,
  `two_cpt_oral_cov`, `three_cpt_iv`, `three_cpt_oral`) now ships an ODE-form
  example alongside its analytical counterpart (`*_ode`, plus new analytical
  `one_cpt_iv` / `three_cpt_oral` examples and datasets). The ODE forms use an
  amount-based convention (states are amounts; observed concentration via
  `[scaling] obs_scale = V`/`V1`), with bioavailability `F` and lag time
  applied by the engine at the dose rather than baked into the `[odes]` RHS.
  A new test (`test-ode-analytical-equivalence.R`) asserts each shipped pair
  gives identical predictions; the exhaustive cross-check across all dosing
  modes (bolus, infusion, multi-dose, steady state, lag, F) lives in
  ferx-core (`tests/analytical_ode_equivalence.rs`). Also fixes
  `bioavailability_ode`, which double-counted `F` (it was both declared as an
  individual parameter -- applied at the dose by the engine -- and baked into
  the absorption flux). (#127)

- **Propensity-score-matched simulation**: `ferx_simulate(..., match = ...)`
  reassigns each replicate's drawn etas to subjects by Mahalanobis matching
  (under the model omega) against the subjects' fitted (posthoc) etas, so a
  subject's observed dosing/sampling design is paired with a similar drawn eta.
  This corrects VPC bias from treatment adaptation in real-world data (e.g.
  longer dosing intervals for high-clearance patients). `match` accepts
  `FALSE`/`"none"` (off), `"optimal"` (or `TRUE`; global linear-assignment
  minimum, best on average and recommended), `"nearest"` (greedy
  nearest-neighbour), or `"rank"` (pair by Mahalanobis-norm rank). Requires
  observed data; the posthoc etas use the fitted parameters when a `fit` is
  supplied. Needs a ferx-core that provides `simulate_with_options` with the
  `match_method` option (separate `Cargo.lock` bump). (ferx-core #288, #396)

- **Standalone importance sampling**: `ferx_fit(..., method = "imp")` now runs
  without a preceding estimator, scoring the model's initial parameters
  (`fit$importance_sampling` is populated, `fit$method_chain` is `"IMP"`). The
  R-side guard that rejected a lone `"imp"` has been removed; the at-most-once
  and must-be-terminal checks remain. Needs a ferx-core that allows standalone
  IMP (separate `Cargo.lock` bump). (ferx-core #269)

- **Covariance estimator & non-PD fallback options**, forwarded via `settings`:
  `covariance_method` (`"r"` inverse-Hessian / `"s"` score cross-product /
  `"rsr"` Huber-White sandwich standard errors) and `covariance_fallback`
  (`"sir"` runs SIR with an absolute-eigenvalue-rectified proposal when the
  finite-difference Hessian is not positive definite). `fit$covariance_status`
  can now be `"sir_fallback"`, which is documented and labelled. Needs a
  ferx-core that provides these options (separate `Cargo.lock` bump).
  (ferx-core #245, #248)

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

- Oral models with a depot-bypassing infusion (`RATE > 0` into the central
  compartment) now return correct concentrations for subjects fit through the
  event-driven analytical path (those with time-varying covariates, reset
  records, or IOV); the infusion input was previously dropped, giving ~0
  predictions for those subjects while no-covariate subjects were unaffected.
  Delivered by bumping the ferx-core pin; no wrapper change (ferx-core#351).

- A `[structural_model]` PK parameter that maps to a name not defined in
  `[individual_parameters]` (e.g. `pk one_cpt_oral(cl=CL, ...)` with no `CL`)
  is now a clear parse error instead of silently fitting a structurally broken
  model (every prediction floored, 100% shrinkage). An unrecognized PK-parameter
  key (a typo such as `clx=`) is likewise rejected, and a numeric-literal value
  (e.g. `ka=1.0`) is honored as a constant rather than silently zeroed. Delivered
  by bumping the ferx-core pin; no wrapper change (ferx-core#261).

- Datasets without an `EVID` column no longer silently fit a dose-free model.
  ferx now infers a dose from a nonzero `AMT` when `EVID` is absent (matching
  NONMEM), so a NONMEM dataset that marks doses only by `AMT`/`MDV=1` administers
  correctly instead of dropping every dose. The reader also warns when `AMT != 0`
  rows are not treated as doses, or when a population parses zero doses despite
  having observations. Delivered by bumping the ferx-core pin; no wrapper change
  (ferx-core#262).

- IOV models: the `sdtab` diagnostic table (`fit$sdtab`) now reports each
  observation's **occasion** individual parameters -- `CL`, `V`, `KA`, any
  `[derived]`/`[output]` column, and `TAD` -- instead of silently using
  `kappa = 0` for every row, so a parameter with inter-occasion variability no
  longer looks identical across occasions. `TAD` additionally shifts each dose
  by its own occasion (and covariate) absorption lag. Delivered by bumping the
  ferx-core pin; no wrapper change (ferx-core#238).

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
  (`multi_start_seed`, `saem_seed`, `sir_seed_used`, `imp_seed`) now guard
  against `NA_character_` / `NA_real_` values that extendr emits for Rust
  `Option<T>::None`, preventing spurious "NA" entries in the run log.

- `ferx_runlog()`: gradient-tolerance line suppressed for derivative-free
  optimizers (BOBYQA, GN, SAEM) where a gradient tolerance is not applicable.

- `ferx_rust_fit()` (internal): `fit$model_text` was `NA` for file-based fits
  because the R binding's provenance block set `model_path` / `model_hash` but
  did not set `model_text`. Fixed by reading the model file in the same block.

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
  residuals are computed. Form B
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
  `ferx_fit()` to use Hamiltonian Monte Carlo proposals in the SAEM E-step.
  New output field `fit$saem_n_subjects_hmc` reports the number of subjects
  that used HMC proposals; `NULL` for MH-only or non-SAEM fits.

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
  `imp_samples`, `imp_proposal_df`, `imp_seed`, `imp_low_ess_threshold`.
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
  `gradient = "auto"` it shows which branch resolved at fit time. The raw
  engine labels are also exposed as
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
