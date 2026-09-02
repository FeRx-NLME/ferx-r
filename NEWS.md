# ferx 0.3.0.9000 (development version)

## Breaking changes

- **An unrecognised `[section]` in a `.ferx` file is now an error** (ferx-core #1040).
  Sections were read by name lookup, so one the engine did not know was never read
  and never reported: a misspelled `[fit_option]` left the model validating clean
  while the fit ran with the default method and no covariance step, returning
  without standard errors. Section names are now closed-world — an unknown header
  fails the parse with `E_UNKNOWN_BLOCK`, naming the offender, its line, the valid
  set and a did-you-mean.

  `ferx_model_validate()` follows suit on two counts. An unknown section now
  counts against the returned `$ok`; it was printed as `[unknown section]` and
  then left out of the status, so `res$ok` was `TRUE` for a model carrying a
  section the engine would ignore. And the list of valid sections comes from the
  engine (`ferx_rust_known_blocks()`) instead of a copy maintained in R — the copy
  had drifted, so `ferx_model_validate(ferx_example("two_cpt_oral_cov")$model)`
  reported `covariates [unknown section]` for a perfectly valid model, and it
  still advertised `[initial_values]`, which the engine dropped years ago
  (ferx-core e5e934d).

  **If your model file carries `[initial_values]`, delete it.** It has been dead
  weight since initial estimates moved inline into `[parameters]`, and it now
  fails the parse with `E_DEPRECATED_BLOCK`, naming the replacement. The three
  files under `examples/models/` that still had one have been fixed.

- **A dose attribute your model also reads is now an error** (ferx-core #993).
- **A dose attribute your model also reads is now an error**
  ([ferx-core #993](https://github.com/FeRx-NLME/ferx-core/issues/993),
  [ferx-core #1004](https://github.com/FeRx-NLME/ferx-core/issues/1004)).
  `F`, `LAGTIME`/`ALAG` and the compartment-indexed `F{n}`/`ALAG{n}`/`LAGTIME{n}`
  are applied by the engine **at the dose event**. A model that declared one and
  *also* referenced it in `[odes]` (the right-hand side or an `init(...)` seed),
  the `[scaling]` readout, or the `[adaptive_dosing] observe` signal was applying
  it twice — silently, and by exactly that factor. Such a model now fails to parse
  with `E_DOSE_ATTR_DOUBLE_USE`, naming both readings and the fix. `D{n}`/`R{n}`
  carry the same reservation but are consulted only for a coded `RATE=-2`/`-1`
  dose, so that collision is reported against the **dataset** instead and a model
  whose data never codes `RATE` is untouched. Reads from `[derived]` / `[output]`
  are post-solve reporting and stay silent.

  **Analytical (`pk ...`) models are covered too** (#1004). The first pass
  rejected ODE models only, on the argument that an explicit `pk(..., f=F)`
  mapping made a second use "stated rather than silent". It did not: nothing in
  the model file says the value is applied twice, and a `[scaling]` or
  `[adaptive_dosing] observe` expression that reads a mapped
  `f=`/`lagtime=`/`alag=` parameter applied it once at the dose and once where it
  was read — on the **default** engine, with no diagnostic, measured at exactly
  `F` on the prediction. Both engines now reject it with the same
  `E_DOSE_ATTR_DOUBLE_USE` code.

  The **remedy differs by engine**, and the message says which applies. On an
  ODE model the *name* routes the parameter, so renaming it fixes the model. On
  an analytical model the name is inert and the **mapping** binds it, so the fix
  is to drop the `f=`/`lagtime=`/`alag=` argument from the `pk(...)` call —
  renaming changes nothing, because the mapping follows the parameter. The
  message quotes the mapping as you wrote it, `alag=` spelling included.

  **If you map a dose attribute and also read it** — e.g.
  `pk one_cpt_oral(cl=CL, v=V, ka=KA, f=F)` together with
  `[scaling] obs_scale = V / F` — that model now fails to parse. Drop whichever
  half was not meant.

  **If you have an ODE model that folds `F` into the absorption flux** — the
  pre-dose-entry convention, e.g. `d/dt(central) = F * KA * depot / V - ...` — it
  will now fail to parse instead of quietly computing `F²`. Drop the `F` from the
  right-hand side; if the parameter was never bioavailability, rename it (the name
  is what routes it).

  Two shapes stay accepted. A parameter merely *named* `F` that no `pk(...)`
  argument maps is an ordinary parameter, so the usual `CL/F`, `V/F`
  apparent-parameter convention is unaffected. And an **analytical**
  `[initial_conditions]` read is fine: an initial condition is not an absorbed
  dose, so the engine seeds the amount with `F = 1` and no lag, and
  `init(depot) = F * 500` — the bioavailable residue of a pre-study dose —
  applies `F` exactly once. Note the scope of that second one: the same reasoning
  has not yet been carried over to the ODE engine, where an `init(...)` seed
  reading a dose attribute is still rejected
  ([ferx-core #1046](https://github.com/FeRx-NLME/ferx-core/issues/1046)).

  **One analytical shape still escapes the check.** A parameter assigned only
  inside an `if` with no `else` is bound by `pk(..., f=F)` but resolves to nothing
  when `[scaling]` reads it, giving `NaN` predictions on a clean parse — tracked
  as [ferx-core #1026](https://github.com/FeRx-NLME/ferx-core/issues/1026).
  Do not read a mapped dose attribute back in
  `[scaling]` on the assumption that the parser will stop you.

  Both readings make ferx **stricter than NONMEM**, measured on NONMEM 7.6.0:
  `$PK` defining `F1` *and* `S2 = V/F1` runs clean under `ADVAN2` and returns
  predictions scaled by exactly `F1`, and a `$PK` `F1` referenced in `$DES` under
  `ADVAN13` quietly computes `F²` — neither draws a diagnostic (ferx-core's
  `nonmem_anchor/analytical_dose_attr_double_use_*` and
  `nonmem_anchor/dose_attr_double_use_*`). A control stream translated literally
  can therefore fail to parse in ferx even though it ran in NONMEM.

- **An unrecognised `[block]` name is now an error**
  ([ferx-core #1040](https://github.com/FeRx-NLME/ferx-core/issues/1040)).
  Blocks were read by name lookup, so a header the parser did not know was never
  read and never reported: a misspelled `[fit_option]` left the model validating
  clean while the fit ran with the default method, the default iteration cap and
  **no covariance step** — returning without standard errors and no indication
  why. The same went for `[scalings]`, `[outputs]` and friends. Block names are
  now closed-world, like the keys inside a block already were: an unknown header
  is `E_UNKNOWN_BLOCK`, listing every offender with its line, the full valid set,
  and a did-you-mean for a near match. Three neighbouring silent drops go with
  it — an instance name on a block that takes none (`[fit_options DOSE]`) or
  missing where one is required (`[covariate_nn]`) is `E_BLOCK_INSTANCE_NAME`; a
  block whose cargo feature this build lacks (`[event_model]`, `[markov_model]`)
  is `E_BLOCK_FEATURE_DISABLED`; and `[initial_values]` — ferx's own former
  spelling for initial estimates, unread since they moved inline into
  `[parameters]` — is `E_DEPRECATED_BLOCK`. **If a model of yours still carries
  an `[initial_values]` block, delete it**: it has done nothing for several
  releases and now stops the parse. No bundled model under `inst/` carries one.
  One rough edge to know about: `ferx_model_validate()` still prints its own
  "Sections present" table from a list that predates this change, so on an
  unusual block that table can disagree with the verdict below it — the engine
  diagnostics are the authoritative half.

- **`obs_scale = TIME` is now a parse error, and an undefined name in `[scaling]`
  no longer reads as zero**
  ([ferx-core #1028](https://github.com/FeRx-NLME/ferx-core/issues/1028)). `obs_scale` is
  subject-static — evaluated once at `t = 0` — so `obs_scale = TIME` (or `= T`)
  could only ever have read `0`; the error names the Form C `y = <expr>` readout
  as the place for a time-dependent term. Separately, `[scaling]` expressions
  never registered their covariate references as required data columns, and
  `ferx_predict()` ran no covariate check at all, so an unresolvable identifier
  reached the predictor as the covariate map's `0.0` default. Both halves now
  register their references, and **`ferx_predict()` reports
  `E_MISSING_COVARIATE`** for a missing column just as `ferx_fit()` and
  `ferx_simulate()` already did — so a `ferx_predict()` call that silently used
  `0` for a missing column now errors. One more narrow break: a `[scaling]`
  expression referencing an *undeclared* data column named `T` now reads the
  model-time built-in instead, matching `[odes]` where that name has always been
  reserved. Declare `T` in `[covariates]` to keep it a data column; whenever the
  fold does happen ferx warns and names both escapes, so it is never silent.
  `TAFD` / `TAD` are unaffected and stay ordinary covariate references.

- **`method = "imp"`, `"impmap"` and `"bayes"` are rejected on a model with no
  random effects, and pure `"gn"` warns**
  ([ferx-core #1006](https://github.com/FeRx-NLME/ferx-core/issues/1006),
  [ferx-core #1007](https://github.com/FeRx-NLME/ferx-core/issues/1007)).
  All three already refused `n_eta = 0` at run
  time, so a `method = c("focei", "imp")` chain ran its whole FOCEI stage before
  failing. `E_METHOD_NO_RANDOM_EFFECTS` now fires up front, anywhere in a chain.
  One consequence when upgrading: a chain with `imp_eval_only = TRUE` on a
  fixed-effects-only model previously returned a fit result with the IMP failure
  downgraded to a warning, and now returns this error instead. Pure Gauss-Newton
  is start-sensitive at `n_eta = 0` — with no inner EBE loop to absorb a poor
  `sigma` start the BHHH step can collapse far from the optimum, with only
  `Converged: NO` as the signal — so it warns (`W_GN_NO_RANDOM_EFFECTS`) and
  points at `"gn_hybrid"` / `"focei"`. Both are suppressed when a later stage
  re-optimises the GN result. The bundled fixed-effects-only examples
  (`one_cpt_iv_pooled`, `binary_logistic`) request `focei` and are unaffected.

## New features

- **`ferx_bootstrap(resume = TRUE)`: continue an interrupted run**
  ([#317](https://github.com/FeRx-NLME/ferx-r/issues/317)). A 200-replicate
  bootstrap that died partway - a killed session, a full disk, a laptop closing
  - had to be started over from R: the engine has had `resume` since ferx-core
  #1143 and the CLI has exposed it as `--resume`, but the R entry point did not.
  It now does. `resume = TRUE` (which needs `directory`) refits only the sample
  indices that directory's `raw_results.csv` does not already carry, and reuses
  the base fit rather than refitting it.

  A resumed run is not an approximation of the uninterrupted one, it is the same
  run: a replicate's draw is a pure function of `(seed, index)`, so a reused
  replicate is bit-for-bit the one a fresh run would have produced. The engine
  refuses to resume from a directory belonging to a different run - the model
  and data hashes, the parameter names and the run settings are all recorded
  next to the replicates and checked first.

  `retry_failed = TRUE` (needs `resume = TRUE`) refits the replicates a previous
  run recorded as *failed* instead of carrying the failure forward. Off by
  default, matching PsN: a fit that failed usually fails again, so this is for a
  transient failure - an out-of-memory kill, a full disk - not a model one.

  ```r
  # Interrupted at replicate 137 of 200
  bs <- ferx_bootstrap(ex$model, ex$data, samples = 200, seed = 1,
                       directory = "warfarin-bootstrap", resume = TRUE)
  ```

- **`ferx_bootstrap(progress = TRUE)`: a progress bar while the replicates fit.**
  A 200-sample bootstrap is minutes to hours behind one call, and until it
  returned there was nothing to say whether it was working or wedged. The engine
  now reports each fit as it completes and R draws it - a cli progress bar when
  the cli package is installed, `utils::txtProgressBar()` otherwise. Default
  `interactive()`, so a script or a knitted document is unaffected.

  The bar reports what the run will actually do: the base model gets a spinner
  of its own (a bar sitting at 0/200 through a single fit reads as stuck),
  `dofv = TRUE` gets a second bar for its second pass over the replicates, and a
  resumed run counts only the replicates still to fit. Fitting happens on a
  worker thread while this R session draws, so nothing in the run waits on the
  bar, and a run watched and the same run unwatched produce the same numbers
  from the same seed. `ferx bootstrap` on the command line grows the same bar
  (`--no-progress` to suppress it).

  ```r
  bs <- ferx_bootstrap(ex$model, ex$data, samples = 200, threads = 8)
  #> Bootstrap replicates [=========>          ]  97/200 | ETA 1m 28s
  ```

- `ferx_gam_screen()`: GAM-based covariate pre-screening. For each ETA x covariate pair
  fits `eta ~ f(cov)` and ranks covariates by AIC improvement over the null model
  (delta_aic = AIC_null - AIC_best). Functional forms tried: linear, natural cubic
  spline (df = 2 and 3 by default), and one-hot categorical. Warns when ETA shrinkage
  exceeds 30%. All numerical work happens in the engine, via
  `ferx_tools::gam::gam_screen()`.

- **`ferx_bootstrap()`: the non-parametric case bootstrap, at PsN's
  `bootstrap` feature parity** ([ferx-core #1144](https://github.com/FeRx-NLME/ferx-core/issues/1144),
  engine side ferx-core #1140). Resample whole subjects with replacement, refit
  the model to each replicate, and get bias, bootstrap standard errors and both
  intervals - the percentile one (`ci_lower` / `ci_upper`) and the
  normal-approximation one built from the bootstrap SE (`ci_lower_normal` /
  `ci_upper_normal`) - side by side, so the disagreement between them is
  visible. Unlike the covariance step it survives a failed or non-positive-
  definite `R^-1`, and it does not assume a symmetric interval.

  ```r
  bs <- ferx_bootstrap(ex$model, ex$data, samples = 200, seed = 12345,
                       threads = 8)
  bs$parameters   # one row per parameter
  bs$raw          # one row per fit, the original dataset first
  bs$diagnostics  # run counts, exclusion tallies, diagnostic means
  plot(bs)        # histogram per parameter (+ the dofv panel when dofv = TRUE)
  ```

  Supports stratified resampling (`stratify_on`), PsN's per-stratum
  `-sample_size` spelled as an R named vector
  (`sample_size = c("1001" = 12, "1002" = 24)`), `dofv`, `keep_covariance`, and
  the four exclusion filters. **`directory` defaults to `NULL`**, where the CLI
  writes `{model}-bootstrap/`: an R user has the data frames in hand and rarely
  wants eight CSVs appearing in the working directory. Set it to get the
  artefacts - and to be able to call `ferx_bootstrap_summarize()` later.

- **`ferx_bootstrap_summarize()`** re-computes a finished run's statistics from
  its `raw_results.csv` under different exclusion criteria - PsN's
  `-summarize`. Refits nothing: it is the recovery path for a run where too
  many replicates were filtered out to resolve a percentile interval.

- **Sample-size-weighted inter-occasion variability is now declared on the
  `kappa` itself** ([ferx-core #1031](https://github.com/FeRx-NLME/ferx-core/issues/1031),
  [ferx-core #1062](https://github.com/FeRx-NLME/ferx-core/pull/1062)).
  Between-treatment-arm variability (BTAV) - the arm-level random effect in
  every published longitudinal MBMA - is distributed
  `kappa_ik ~ N(0, gamma^2 / N_ik)`: a 400-subject arm's mean wanders a quarter
  as far as a 25-subject arm's. Until now the only way to express that in ferx
  was to write the divisor into a structural expression by hand
  (`... + KAPPA_EMAX / sqrt(NARM)`), where `/ NARM` instead of `/ sqrt(NARM)`
  produces a plausible wrong answer rather than an error. It is now declared
  where the rest of the variance structure is:

  ```
  [parameters]
    kappa KAPPA_EMAX ~ 2.0 (sd) weight = NARM

  [individual_parameters]
    LEMAX = LEMAX0 + log(OR_ABATA) * ABATA + ETA_EMAX + KAPPA_EMAX
  ```

  The engine applies the weight by reparameterisation, so the estimate in
  `fit$omega_iov` stays the **unweighted** `gamma^2` - the quantity a published
  analysis reports - and so do the kappa EBEs and their shrinkage. On the R
  side:

  - `print()` on the fit annotates the weighted kappa with the effective
    between-arm SD at a typical arm, `gamma / sqrt(W)`:

    ```
    KAPPA_CL = 1.840240  (CV% = 230.2)  SE = N/A  Shrinkage = 37.4%
                         weight = NARM  ->  SD = 0.0678 at NARM = 400.0000
    ```

  - `ferx_model_inspect()` reports the weight both pre-fit (from the model
    file) and post-fit (from what the engine parsed), as
    `IOV:  KAPPA_CL (weight = NARM)`, and returns it in `$iov_weights`.
  - `ferx_fit()` results carry `kappa_weights` (the weight expression per
    kappa, `NA` where unweighted) and `kappa_weight_typical` (the median
    weight over the dataset's subject-occasions), and both survive a
    `ferx_save_fit()` / `ferx_load_fit()` round-trip.

  Models with no weighted kappa are unaffected in every one of those places:
  both fields are `NULL` and the `.fitrx` bundle is byte-identical to before.

  Three follow-ups closed the gaps where the weight was still invisible:

  - `ferx_model_inspect()` on a **model file** now recognises the modifier by
    the same rules the engine's parser uses (case-insensitive, whole-word, at
    bracket depth 0, on a single `=`). A model written `WEIGHT = NARM` fits
    with the weight applied but was reported pre-fit as unweighted, and a
    weight expression containing a comparison (`weight = NARM * (FLAG == 1)`)
    was truncated to `1)`.
  - `summary()` annotates the weighted kappa on its `IOV:` line, as `print()`
    and `ferx_model_inspect()` already did. All three now share one helper, so
    they cannot disagree about whether a model is weighted.
  - `fit$estimates` gains a `weight` column carrying the weight expression on a
    weighted kappa row and `NA` everywhere else. The `estimate` on such a row
    is the unweighted `gamma^2`, so a reader of the table alone had no way to
    tell it from an ordinary kappa and would take `sqrt(estimate)` for the
    between-occasion SD instead of `sqrt(estimate / W)`.

- **`ferx_model_inspect()` no longer depends on the case of a declaration
  keyword.** The engine's `theta` / `omega` / `sigma` / `kappa` declaration
  regexes are all case-insensitive, so a model written `THETA TVCL(1.0, 0.001,
  100.0)` fits exactly like the lowercase spelling - but the R-side reader
  matched case-sensitively and reported *no* population parameters, no IIV and
  no IOV for it: an entirely blank structure for a model the engine parses
  fine.

  As a guard against the next such drift, a `[parameters]` block that has
  content but not one recognised declaration in it now warns instead of
  silently returning an empty structure. It is keyword-agnostic, so it fires on
  whatever the next divergence turns out to be. No bundled example triggers it.

- **`ferx_simulate()` and `ferx_predict()` now accept a design template whose
  `DV` column is empty** (#286, ferx-core #957). Dosing plus sampling times with
  `DV = "."` / `NA` - the natural way to write a design, and what NONMEM's
  `$SIMULATION` accepts - previously returned **zero rows**: every observation
  record with a missing `DV` was dropped as a forgotten `MDV = 1`, which is the
  right reading only when the `DV` is an input. Simulation and prediction both
  *produce* that column, so such a record is now read as a sampling time whose
  value is about to be generated, and no placeholder number is needed.
  `MDV = 1` still excludes a record, and fitting is unchanged. Applies to
  `ferx_simulate()` (both the default-parameter and the `fit = ` path),
  `ferx_simulate_adaptive()`, `ferx_simulate_with_uncertainty()`, and
  `ferx_predict()`.

  Because the simulate and fit readers now disagree on the same file, a kept
  empty-DV record is **counted and reported**: the returned object carries the
  count in its `simulation_warnings` attribute, and it is re-emitted as an R
  warning. A dataset with an *accidental* missing `DV` (rather than a deliberate
  design) would otherwise silently gain simulated rows at times `ferx_fit()`'s
  `sdtab` has no observation for - biasing a VPC built by overlaying the two.

- **Models with no random effects - fixed-effects-only (naive-pooled) fits**
  (ferx-core #989). A continuous model may now omit every `omega` declaration.
  With `n_eta = 0` there is no inner empirical-Bayes problem and no `log|Omega|`
  term, so FOCE/FOCEI collapse to plain maximum likelihood: every subject shares
  one set of parameters and `sigma` alone carries the spread. `PRED` equals
  `IPRED` and `CWRES` equals `IWRES`. An `[error_model]` and its `sigma` are
  still required - only Omega is optional. New bundled example
  `ferx_example("one_cpt_iv_pooled")`, anchored against a NONMEM 7.6.0 run of
  the same model written as `$OMEGA 0 FIX` (OFV -269.6370, TVCL 4.8408,
  TVV 52.834).

- **`TIME` now works in a `[scaling]` Form C readout**
  ([ferx-core #1028](https://github.com/FeRx-NLME/ferx-core/issues/1028)).
  A `y = <expr>` / `y[CMT=N] = <expr>` readout referencing `TIME` parsed fine but
  was never bound to the observation, so it read `0` at every row and the whole
  time-dependent term vanished — a response-versus-time readout such as
  `y[CMT=1] = EMAX * TIME / (TIME + T50)` therefore fit, converged, and reported
  plausible parameters for a structural model nobody wrote. `TIME` (and the `T`
  alias) now resolves to each observation's own time on both the ODE and
  analytical paths, and to each decision's time in `[adaptive_dosing] observe`.
  The dummy `d/dt(clock) = 1` workaround is no longer needed, and is better
  dropped: `clock` starts at the subject's first record, not at `t = 0`. The
  readout's `TIME` is the raw data-file clock — the one sdtab, `ferx_predict()` /
  `ferx_simulate()` and `[derived]` windows report.

- **New `mstep_damping` setting for SAEM**
  ([ferx-core #1011](https://github.com/FeRx-NLME/ferx-core/issues/1011)).
  SAEM's numerical theta/sigma M-step assigned the maximiser outright,
  re-maximising against a *single* MCMC eta draw each iteration — `argmax` of one
  draw rather than the stochastic-approximation average — a Monte-Carlo bias that
  does not decay with iteration count. It bit only a theta left to the numerical
  M-step without a mu-reference; a log-mu-referenced theta was already exact. The
  M-step result is now blended in as `theta <- theta + gamma * (theta* - theta)`,
  with `gamma` capped at `0.03` during exploration. Pass
  `settings = list(mstep_damping = ...)` to change it; **`1.0` disables the
  damping entirely** and restores the previous behaviour exactly. A fit whose
  every estimated theta is mu-referenced, `FIX`ed or pinned out is bit-identical
  to before — both bundled SAEM examples (`warfarin_saem`, `warfarin_iov_saem`)
  are log-mu-referenced throughout and are unchanged. SAEM also now warns when an
  estimated theta carries no ETA at all, which is the shape the bias was measured
  on; attaching an ETA or holding the parameter `FIX` remains the better fix.

## Bug fixes

- **An infusion into a built-in absorption compartment is no longer delivered
  twice** ([ferx-core #1187](https://github.com/FeRx-NLME/ferx-core/issues/1187)).
  A `RATE > 0` or `RATE = -2` dose into the absorption compartment of an ODE model
  had its rate applied both through the absorption kernel and, a second time, as a
  plain input rate straight into the target compartment - so the dose was double
  counted *and* the second copy bypassed absorption entirely. On a mass-balance
  readout the excess is exactly `2 * F * amt`; on a concentration readout it is
  worst early (214x at t = 0.5 on a transit model) and settles near 1.86x once
  both copies have distributed.

  Six surfaces were affected: ODE `IPRED` and the state columns in `sdtab`, the
  cumulative hazard and `ferx_simulate()` event times of a joint PK-TTE model,
  the CTMM likelihood, and `[derived]` grid integrals. The event-driven path
  returned a correct `IPRED` beside wrong states, so an sdtab row could disagree
  with itself. Bolus dosing, and any model whose infusion targets a central
  compartment, were never affected.

- **`ferx_get_warnings()` now answers a failed covariance step with the advice
  it was written to give.** Every targeted branch of the covariance guidance -
  a non-positive-definite Hessian with its eigenvalue list, ill-conditioned
  Hessian entries naming the parameter, a near-singular omega, a non-finite
  base OFV, and the minor/moderate/severe regularisation tiers - sat behind
  `category == "covariance_step"`. ferx-core does not code those messages that
  way: it codes them `covariance_failed` and `covariance_regularized`, and
  reserves `covariance_step` for an informational note about the step's cost.
  So the whole block was unreachable and a failed covariance step printed no
  guidance at all, while the one message that did reach it - the benign cost
  note - was answered with "Standard errors unavailable. Check identifiability",
  reporting a failure that had not happened.

  Guidance is now selected by the message, which is what those branches were
  always keyed on. Four categories reach it: ferx-core's three, plus the
  `covariance` category `ferx_covariance()` assigns to the same engine messages,
  so the post-hoc covariance step gets guidance too. A covariance message
  arriving under `general` is answered as well, which covers every fit read back
  with `ferx_load_fit()` - it does not restore the structured table, so all its
  warnings arrive under that category. Admission is anchored to messages that
  *begin* with "Covariance step", so SIR's diagnostics, which mention the step
  in passing, keep their own guidance instead of being told to re-run with
  `covariance = FALSE` - the one setting that removes the matrix SIR needs.

  Three messages that were being answered wrongly now have their own advice: an
  off-diagonal FD stencil that could not be evaluated, which ferx-core reports
  on its *success* path (the standard errors exist and are merely
  over-optimistic, not missing); a covariance step cancelled part-way, which
  produced no standard errors but diagnosed nothing about the model; and the
  informational cost note. A near-singular omega now reaches the omega branch
  too - ferx-core picks that descriptor from the sign of the smallest
  eigenvalue, and the guidance had matched only the other one.

- **The unused-parameter warning gets its guidance back.** ferx-core has no
  `unused_parameter` code - its unused-declaration messages classify to
  `general` - so the guidance written for them was unreachable in the same way
  the covariance block was. They are now routed by message. ferx-core's
  flat-theta warning contains the same "computed but never used" phrase and is
  excluded explicitly, so it is not answered with the unused-parameter text.

- `ferx_sir()` now reports the engine's SIR-step warnings instead of dropping
  them (ferx-core #1021). The binding returned only the CIs and the ESS, so the
  proposal diagnostics ferx-core emits - a covariance that is rank-deficient
  beyond its `FIX`ed parameters, or a proposal direction shrunk to keep draws
  inside the parameter bounds - never reached the fit. They now land in
  `fit$warnings` and in `fit$warnings_structured` under the `sir` category, and
  `ferx_get_warnings()` answers them with guidance about the *model* (the named
  parameters are not identified by the data, and their SIR intervals understate
  the uncertainty) rather than the old "raise `sir_samples`" advice, which does
  not help in that case. The same ferx-core release stops those fits from
  failing outright with `All SIR samples had invalid weights` - which is what
  every model-based meta-analysis hit, since fixing the residual variance is the
  inverse-variance weighting scheme and cannot be dropped.

- `ferx_covariance()` and `ferx_sir()` no longer reject a fit that has no random
  effects (#290). Both required a non-empty `fit$ebe_etas`, which a
  fixed-effects-only fit does not have - there are no EBEs to warm-start from -
  so both stopped with "fit$ebe_etas is empty; cannot warm-start the inner
  loop". They now pass an empty warm start through. This matters for
  naive-pooled models in particular: NONMEM's `$COVARIANCE` default is the RSR
  sandwich, and `ferx_covariance(fit, covariance_method = "rsr")` is how you
  reproduce it from R. On such a model the sandwich runs about twice as wide as
  ferx's default `"r"`, because the model ignores within-subject correlation by
  construction.

- `print()` on a fit with no random effects no longer prints an empty `OMEGA`
  section header (#290). An empty header over an empty body reads as an
  estimation that failed rather than one that was never requested.

- **`ferx_simulate(..., fit = f)` works again for models with inter-occasion
  variability** (ferx-core #1019). The fitted IOV covariance (`fit$omega_iov`) was
  never passed to the engine, so any model declaring a `kappa` crashed the R session
  with `omega_iov is present whenever the model declares kappa (n_kappa > 0)` instead
  of simulating. `fit$omega_iov` is now threaded through every from-fit entry point
  (`ferx_simulate()`, `ferx_simulate_with_uncertainty()`, `ferx_calc_npde()`,
  `ferx_predict()`, `ferx_predict_survival()`). `ferx_calc_npde()` simulates
  internally and hit the same crash; `ferx_simulate_with_uncertainty()` silently drew
  around the model file's *initial* IOV variance instead of the fitted one, and now
  uses the fitted value.

- **An adaptive-dosing `dv` monitor no longer floors a negative Form C
  `[scaling]` readout at zero**
  ([ferx-core #1039](https://github.com/FeRx-NLME/ferx-core/issues/1039),
  [ferx-core #1020](https://github.com/FeRx-NLME/ferx-core/issues/1020)). The
  assay floor ("an assay cannot read below zero") was applied unconditionally
  after the residual draw. A Form C
  `y = <expr>` readout is an arbitrary expression — a change from baseline, a
  difference from a comparator, a z-score — so the *same model* read correctly
  under `mode = "ipred"` and came back as exactly `0` under `mode = "dv"` for
  every negative sample, silently: a controller thresholding a
  change-from-baseline signal saw `0` over precisely the region it was written to
  react to, and dosed accordingly. The floor is now gated on the same predicate
  as the prediction path, so it applies only to the bare-state readout and to
  Forms A/B, which keep it.

- **Engine warnings now carry an `absorption_twin_declined` category**
  ([ferx-core #1008](https://github.com/FeRx-NLME/ferx-core/issues/1008)).
  A transit or inverse-Gaussian
  model that cannot build the ODE twin it falls back on now says so at parse time, with its
  own reason, instead of declining silently. It reaches `fit$warnings` and
  `fit$warnings_structured` like any other engine warning. Note
  `ferx_get_warnings()` prints the category but has no remediation text for it
  yet, so it shows without the guidance block other categories get.

## Internal

- **`ferx-tools` is now a second git dependency**, from the same ferx-core
  repository and the same revision as `ferx-core` - one extra `Cargo.lock`
  entry, no second pin to track. `src/Makevars` writes the matching `[patch]`
  line for it, so a local `../ferx-core` checkout patches *both* crates: with
  only the `ferx-core` entry a local build silently mixed a working-tree
  `ferx-core` with a GitHub-`main` `ferx-tools`. The `R-CMD-check` pin guard
  and `tools/update-ferx-core-lock.sh` now check both crates and that they sit
  on one revision.

- Bumped the pinned ferx-core commit and updated the extendr glue for the
  `mixture` / `pmix` / `mixest` fields added by ferx-core #977/#985 (#291).
- Bumped the pinned ferx-core commit to `25b5f473` for ferx-core #993 (the
  dose-attribute double-use rejection above). No glue change: the new items on
  the Rust side are additive and unused here, so `src/rust/src/lib.rs` is
  untouched.
- Bumped the pinned ferx-core commit to `0f571d83` (#304). Beyond ferx-core
  #1004 above, that range carries #1040, #1028/#1042/#1045, #1039/#1020,
  #1011/#1012, #1006/#1007, #1008, #1019 and #1021 — every user-visible one
  is written up in the sections above. No glue change.

# ferx 0.3.0

## Breaking changes

- **`CMT=0` on an ODE dataset now predicts differently** (ferx-core #899). `CMT=0` is
  NONMEM's *default dose compartment* and resolves to compartment 1. The ODE engine
  previously did four different things with it depending on which internal driver a
  subject took — including dropping the dose silently — so the same dataset could
  produce different answers, and a fit could differentiate a different dosing history
  than it predicted. Every site now resolves it to compartment 1, and compartment-indexed
  dose attributes (`F1`, `ALAG1`) are read correctly. **If you have ODE datasets written
  with `CMT=0`, earlier results were wrong and should be regenerated.**
- `method = "agq"` has been removed (ferx-core #251). Adaptive Gauss-Hermite
  quadrature is not a separate estimator — it is the single-point method with more
  nodes, so the node count is now an argument and the method name selects the
  Hessian anchor: `ferx_fit(..., method = "laplace", settings = list(n_agq = N))`
  is the exact-anchor quadrature (`n_agq = 1` is Laplace), and
  `method = "focei"` with `n_agq > 1` is the new Gauss-Newton-anchored quadrature.
  The old `"agq"` / `"gauss_hermite"` tokens now error with a pointer to the
  replacement.

Public functions renamed for verb-clarity and naming consistency (part of the
API cleanup in #223; naming rule + hard-break policy decided in #224). Old
names are removed - no deprecation shims. Update calls as follows:

## Added

- **`?ferx_fit` now documents every fit option the engine accepts** (99 of 101; the
  two exceptions are the FREM structural maps that `ferx_model_to_frem()` writes
  for you, and the help says so). Previously 32 keys were reachable through
  `settings` but documented nowhere, so the only way to find them was to read the
  engine source. Newly documented: `inner_restarts`, `inner_optimizer`,
  `cov_inner_tol`, `parameter_scaling`, `ebe_warm_start`, `checkpoint` /
  `checkpoint_interval_secs`, `iov_column` / `iov_occasion`, `npde_nsim` /
  `npde_seed`, `sir_df` / `sir_keep_samples`, `conddist` and its three companions,
  `imp_auto` / `impmap_auto`, `imp_defensive_alpha`, `iscale_min` / `iscale_max`,
  `frem_rao_blackwell`, `impmap_mceta` and `impmap_sobol`. Two of these change how
  an existing documented option behaves and are worth knowing: `imp_auto` /
  `impmap_auto` default to `TRUE`, which makes `imp_samples` / `impmap_samples` a
  *starting* count that ramps up rather than a fixed one; and the applicability
  headings were wrong — `method = "laplace"` accepts the whole outer-optimizer and
  iteration-cap block, and the inner-loop keys apply to `"imp"`, `"impmap"` and
  `"bayes"` too.
- **Stiff and high-order ODE solvers via `settings = list(ode_method = ...)`** —
  `"rk45"` (default), `"vern7"`, `"rosenbrock23"`, `"rodas4"` and `"rodas5p"`.
  These cover two independent problems that want opposite fixes: a
  *stability*-limited (stiff) model — fast reversible binding / TMDD,
  Michaelis-Menten with `KM` far below observed concentrations, long transit
  chains — takes tiny steps whatever `ode_reltol` asks for, and wants one of the
  linearly implicit Rosenbrock methods; an *accuracy*-limited model accepts
  nearly every step and only slows down as `ode_reltol` tightens, where a stiff
  method buys nothing and `vern7`'s higher order is the lever (~2.3× at `1e-9`
  on ferx-core's transit benchmark, but ~1.4× *slower* at default tolerances).
  Every method is a full peer — analytic sensitivities, time-to-event and
  categorical endpoints, simulation and adaptive dosing work with all of them.
  Also settable in the model file's `[fit_options]` block. Delivered via the
  ferx-core pin bump (ferx-core #952 / #387).

- **Exact analytic covariance R-matrix, on by default** — the covariance step now
  assembles the observed information from third-order sensitivities of the
  closed-form prediction rather than differencing the objective, for models in
  scope (plain analytical Gaussian; no IOV, LTBS, `[scaling]`, M3, FREM or
  non-Gaussian endpoint). Both `method = "focei"` and `method = "foce"` are
  served, from two separate assemblies — the non-interaction one is built on the
  Sheiner–Beal gradient and carries no `log|H~|` term — so neither falls back to
  finite differences (ferx-core #954 pins both end-to-end). This removes the
  `eps/h^2` differencing noise and the
  `fd_hessian_step` tuning knob, and costs `2 * (n_theta + n_eta) + 1`
  sensitivity evaluations per subject instead of roughly `2 * n_free^2` objective
  evaluations that each re-solve every inner loop. Out-of-scope models keep the
  finite-difference stencil unchanged. **Standard errors on in-scope models may
  shift slightly** — they are now exact rather than finite-difference
  approximations; set `settings = list(analytic_cov_hessian = FALSE)` to
  reproduce pre-bump values. Note `fd_hessian_step` is inert on in-scope models
  for the same reason. Delivered via the ferx-core pin bump (ferx-core #953 /
  #436).

- **Adaptive dosing now accepts a pre-scheduled base regimen (loading /
  maintenance dose)** — `ferx_simulate_adaptive()` no longer requires dose-free
  base subjects. Ordinary dose rows in the data (`EVID = 1/4`, with `AMT`, and
  optionally `RATE`, `SS`, `II`) are integrated as a standing prescription and the
  `[adaptive_dosing]` controller augments them, the real therapeutic-drug-monitoring
  / model-informed-precision-dosing workflow. System resets (`EVID = 3/4`) are also
  honored. Note the returned dose ledger and `metrics$CUM_DOSE` count
  controller-issued doses only, so a pre-scheduled base dose is excluded from them
  (it is still reflected in the trajectories, `PCT_TIME_IN_WINDOW`, and the
  `auc_target` metric). New bundled example
  `ferx_example("adaptive_vanco_loading")` — a vancomycin loading-dose +
  maintenance titration — with a runnable `inst/examples/ex_adaptive_vanco_loading.R`.
  Delivered via the ferx-core pin bump (#276; ferx-core #702 / #716 / #929 and
  follow-ups).

- **New bundled examples `ferx_example("ss_absorption")` and
  `ferx_example("infusion_absorption")`** — steady-state dosing (`SS=1`, `II`) and
  infusion (`RATE>0`) into a built-in absorption compartment (`first_order(ka)`
  forcing central), the two dosing routes that ferx-core #719 (gaps 1 and 2)
  added for the pointwise density absorption kernels. Both were previously
  rejected at parse time. Each ships a runnable `inst/examples/ex_*.R` and is
  anchored to NONMEM 7.6.0 (ADVAN2 TRANS2): the dataset `DV` is the NONMEM
  population prediction and `ferx_predict()` reproduces it to < 1e-4. Steady-state
  equilibrates the periodic dosing through the absorption kernel; the infusion
  becomes the zero-order source feeding the kernel.

- **New bundled example `ferx_example("binary_logistic")`** — a fixed-effects
  `[binary_model]` (logistic) endpoint: the 0/1 outcome on CMT 3 is Bernoulli with
  `logit P(DV = 1) = TH0 + THX * X + THT * TIME`, the exact analogue of base-R
  `glm(DV ~ X + TIME, family = binomial)`. Ships with a runnable
  `inst/examples/ex_binary_logistic.R` that fits it and shows `ferx_simulate()`
  returning 0/1 `DV_SIM` on the binary CMT (#271). Delivered via the ferx-core pin
  bump (#900).

- **Per-route absorption lag in model files** — every built-in input-rate function
  now takes an optional `lag=` argument (`first_order(ka=KA, lag=L)`,
  `zero_order(dur=DUR, lag=L)`, `transit`, `igd`, `weibull`), giving each parallel /
  mixed pathway its own onset delay on top of any compartment lagtime — the classic
  immediate-release + delayed-release picture that a single per-dose lagtime cannot
  express (ferx-core #856). New bundled example
  `ferx_example("per_route_lag_absorption")` with a runnable
  `inst/examples/ex_per_route_lag_absorption.R`. Delivered via the ferx-core pin bump.

- **`ferx_xpose()` now populates the estimation-iteration trace**, so
  `xpose::prm_vs_iteration()` (parameter value vs iteration) and
  `xpose::grd_vs_iteration()` (gradient vs iteration) work on the returned
  object (#168). When the fit was run with `optimizer_trace = TRUE`, the
  per-parameter value and gradient trajectories are written into the xpose
  `$files` slot as synthetic NONMEM `.ext` / `.grd` tables. A new
  `iterations` argument (default `TRUE`) gates this; the `.grd` table is only
  emitted for gradient-based methods, and when no trace is present the slot is
  left empty (the iteration plots then raise xpose's usual "no files" message
  while the goodness-of-fit / covariate plots are unaffected). Only the
  `"xpose"` backend is affected. Builds on the ferx-core optimizer trace now
  carrying per-parameter estimates and gradients per iteration (ferx-core #640).

- **Analytic inverse-Gaussian (IG) absorption is now available in model files** -
  `pk one_cpt_ig(cl, v, mat, cv2)` and `pk two_cpt_ig(cl, v1, q, v2, mat, cv2)`
  structural models: Freijer & Post inverse-Gaussian absorption fed straight into
  a one- or two-compartment disposition as an analytic closed form (ferx-core
  #790), with exact FOCE/FOCEI sensitivities that do not depend on any ODE-solver
  tolerance and a uniform pk-line interface matching the analytic transit models.
  It is the closed-form counterpart to the ODE `igd()` input rate
  (`ferx_example("igd_inverse_gaussian")`); `f` and `lagtime` are supported. Outside
  the closed form's convergence domain a plain model transparently reroutes to its
  ODE `igd()` twin (a model that also maps `f`/`lagtime` has no twin and is
  rejected, rather than rerouted). New bundled examples `ferx_example("one_cpt_ig")` and
  `ferx_example("two_cpt_ig")` with runnable `inst/examples/ex_one_cpt_ig.R` /
  `ex_two_cpt_ig.R`.

- **Bundled analytic two-compartment transit example** - a new
  `ferx_example("two_cpt_transit")` for the
  `pk two_cpt_transit(cl, v1, q, v2, n, mtt)` closed form (ferx-core #634): Savic
  transit-compartment absorption superposed bi-exponentially onto a
  two-compartment disposition, the 2-cpt analytic counterpart to
  `one_cpt_transit` and the closed-form counterpart to the ODE `transit_2cpt`.
  Paired with a model-simulated 2-cpt transit-truth dataset (a genuine
  parameter-recovery example, unlike `one_cpt_transit`'s shared anchor) and a
  runnable `inst/examples/ex_two_cpt_transit.R` (closes #251).

- **Inter-occasion variability (IOV) now composes with the analytic absorption
  closed forms** - `pk one_cpt_transit` / `two_cpt_transit` / `one_cpt_ig` /
  `two_cpt_ig` previously rejected a `kappa` random effect at parse time. A
  subject carrying IOV is now transparently rerouted, per subject, to the model's
  exact `transit()` / `igd()` ODE twin, which integrates the cross-occasion dose
  carryover the closed-form superposition cannot express - no switch to a
  hand-written ODE model is needed, just the same `pk ...` line plus a `kappa`
  random effect and `iov_column` (via ferx-core #719). New bundled example
  `ferx_example("one_cpt_transit_iov")` - analytic `one_cpt_transit` with IOV on
  CL, paired with an 8-subject subset of the ferx-core transit+IOV NONMEM anchor
  dataset (simulated from the model, so the fit recovers the data-generating
  parameters) and a runnable `inst/examples/ex_one_cpt_transit_iov.R`. Steady-state
  dosing and infusions
  under IOV on the analytic path remain unsupported (ferx-core #719). Requires the
  bumped ferx-core.

- `ferx_simulate()` now surfaces per-subject simulation diagnostics from
  ferx-core (#762 / #763): a degenerate or pathological hazard that would
  otherwise censor a subject with no event is raised as an R warning and
  attached to the returned data frame as a `simulation_warnings` attribute
  (a character vector, empty for a clean run). Requires the bumped ferx-core
  (`simulate_with_options_diag`).

- New `ferx_covariance(fit)` runs the finite-difference-Hessian covariance step
  against an existing fit without re-estimating (#738), the covariance-step
  analogue of `ferx_sir()`. Add standard errors to a fit produced with
  `covariance = FALSE`, or re-run the step with a different `covariance_method`
  (e.g. the `"rsr"` sandwich), including on a fit loaded from a `.fitrx` bundle.
  It re-reads the model/data from the fit's recorded paths with SHA-256
  integrity checks (refusing stale inputs) and refreshes `cov_matrix`,
  `cor_matrix`, `se_theta`/`se_omega`/`se_sigma`/`se_kappa`,
  `covariance_status`, `eigenvalues`, and `condition_number`. The numerics
  closely match `ferx_fit()`'s inline covariance step (the same engine step; the
  standalone re-reads the data and cold-starts the inner EBE loop, so agreement
  is close but not bit-exact); a step that runs but fails
  (non-PD / unusable Hessian) is non-fatal, reporting
  `covariance_status = "failed"` with a diagnostic warning.

- Model files may now declare a `[data]` block (`path = ...`, resolved relative
  to the model file's directory). When `data` is omitted, `ferx_fit()`,
  `ferx_model()`, `ferx_simulate()`, `ferx_predict()`,
  `ferx_predict_survival()`, `ferx_simulate_adaptive()`, `ferx_check_init()`,
  and `ferx_inits_from_nca()` fall back to the declared dataset; an explicit
  `data` argument still overrides it (#254).

## Fixed

- **`ferx_example("warfarin_scaled")` could not be fitted at all.** Its model file
  carried `gradient = ad` in `[fit_options]`; the Enzyme automatic-differentiation
  path was retired in ferx-core in favour of the analytic `Dual2` sensitivities, so
  the engine now rejects that token outright (`E_AD_RETIRED`) and the bundled
  `ex_warfarin_scaled.R` failed immediately. Now `gradient = auto`.
- **`fit$gradient_used` never reported the analytic gradient.** The internal label
  map still translated the retired `"Enzyme AD"` string and had no case for the
  engine's actual `"analytic (Dual2)"`, so the long string fell through unmapped —
  meaning `fit$gradient_used == "ad"` was permanently `FALSE` and `print()` /
  `summary()` rendered the raw engine string. It now reports `"analytic"`.
- **Documentation corrections found by an audit against the engine.** The
  `[scaling]` help told users that expression (`obs_scale = V`) and Form C
  (`y = <expr>`) readouts force finite-difference gradients and to set
  `gradient = fd`. Both are differentiated exactly under the default
  `gradient = auto` (ferx-core #486), so that advice lost the analytic gradient
  *and* switched the outer optimizer from L-BFGS to BOBYQA; only the per-CMT
  variants fall back, and they do so silently. `bloq_method = "drop"` was described
  as discarding BLOQ rows when it keeps them, fitting each at its limit value.
  Corrected defaults: `inner_tol` (`1e-4` → `1e-5`), `n_mh_steps` (`10` → `20`),
  `impmap_proposal_df` (documented as `"normal"`, actually Student-t `4`), and
  `threads` (documented as one worker per logical CPU, actually cores − 1 capped at
  8). Also: `max_unconverged_frac` rejects an outer step rather than relaxing the
  converged flag; `covariance_method = "s"`/`"rsr"` work under FOCE, not FOCEI only;
  `optimizer = "bfgs"` is a deprecated alias for `nlopt_lbfgs`, not a distinct
  algorithm; `gradient = "ad"` now errors rather than being tolerated; and
  `method = "laplace"` corresponds to NONMEM `LAPLACIAN INTER` (plain `LAPLACIAN`
  differs by ~9 OFV units).
- **ODE-form models fit with `method = "foce"` now match their analytical
  closed-form equivalent's marginal objective** (via ferx-core #378). When a
  subject's per-subject (inner EBE) objective was multimodal, the analytical and
  ODE forms could condition on different modes, so their FOCE marginal OFV
  diverged - by up to ~18 units on some models/platforms (most visibly the
  3-compartment IV `three_cpt_iv` example on Linux). ferx-core now keeps the
  better inner estimate on the analytical path, matching the ODE path, so the two
  forms agree to solver round-off. Delivered via the ferx-core pin bump.

- **Simulated binary / categorical outcomes are no longer `NA`** (#271). With a
  `[binary_model]` endpoint (ferx-core #900), `ferx_simulate()` mapped every
  simulated row through `continuous_value()`, whose categorical arm returns `NaN`,
  so each binary draw came back as `DV_SIM = NA` and was indistinguishable from a
  PK row that failed to predict. The simulate frame now folds a categorical draw
  into `DV_SIM` as its numeric 0/1 outcome (matching how the input CSV codes DV);
  combined with the existing `CMT` column, simulated binary outcomes are now
  usable from R. (TTE `Event` rows are unchanged: `DV_SIM = NA`, event time in
  `TIME`, `OBSERVED` flag set.)

- **`ferx_fit()` on a model with no random effects (`n_eta = 0`)** - e.g. a
  fixed-effects `[binary_model]` logistic regression - no longer errors with
  "missing value where TRUE/FALSE needed" (#271). An eta-less fit returns an empty
  omega; the R post-processing left it as a bare `numeric(0)` instead of a 0x0
  matrix, so `nrow(omega)` was `NULL` and poisoned the eta-metadata guards with
  `NA`. `omega` is now always a matrix, so `n_eta = 0` models fit cleanly.

- **`ferx_save_fit()` / `ferx_load_fit()` round-trip on extension-less paths**
  (#268). Saving to a path with no file extension (e.g.
  `ferx_save_fit(fit, "results/run1_base_diag")`) previously landed at
  `results/run1_base_diag.zip` - because Info-ZIP appends `.zip` to a
  suffix-less archive name - so the mirroring
  `ferx_load_fit("results/run1_base_diag")` failed with "File does not exist".
  `ferx_save_fit()` now appends the conventional `.fitrx` extension when the
  output path has none, and `ferx_load_fit()` falls back to `<path>.fitrx`, so
  the bare-path round-trip works. Explicit extensions are still honoured as
  given.

- **Closed-form transit / inverse-Gaussian absorption under IOV, time-varying
  covariates, or a `TIME` switch now honors a call-time ODE tolerance and
  converges its per-subject estimates correctly** (via ferx-core #814, a #719
  follow-up). These models serve such subjects on an internally generated ODE
  "twin"; before this a `ferx_fit(settings = list(ode_reltol = ...))` (or
  `ode_abstol` / `ode_max_steps`) override was silently dropped on the twin path,
  and an estimate that fell back to the finite-difference inner gradient could
  stop short of convergence. The `one_cpt_transit_iov` example above is the main
  beneficiary. Requires the bumped ferx-core.

# ferx 0.2.0

## Breaking changes

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

The four section-editing functions collapse into one get/set pair, part of
the API cleanup in #223 (#233; decision recorded on #227). Both accept
either a `ferx_model` object or a plain path:

- `ferx_model_section()`, `ferx_get_section()` -> `ferx_model_get_section()`
  (returns the section's lines; always a data-return, not a pipe passthrough)

- `ferx_set_section()` -> `ferx_model_set_section()` (unchanged behaviour:
  returns `x` for piping, with copy-on-write for bundled package models)

The `ferx_get_section()` mid-pipe peek (printing a section and passing the
`ferx_model` object through unchanged) is gone - there is no replacement that
both prints and continues the pipe. Call `ferx_model_get_section()` on its
own line before the pipe, or use `ferx_model_show()` to peek at the whole
file:

```r
# before
fit <- ferx_model(ex$data, ex$model) |>
  ferx_get_section("parameters") |>
  ferx_fit()

# after
ferx_model_get_section(ex$model, "parameters")
fit <- ferx_model(ex$data, ex$model) |>
  ferx_fit()
```

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

- New `ferx_conddist(fit)` exposes the SAEM conditional-distribution results
  (`settings = list(conddist = TRUE)`) to R (#244): per-subject/per-eta
  conditional mean, SD, and mode (`fit$cond_dist`), with distribution-based
  eta-shrinkage as an attribute. Previously `cond_dist` was computed by
  ferx-core but never reached R, for either in-process fits or `.fitrx`
  bundles; it now survives `ferx_save_fit()` / `ferx_load_fit()` too.

- `ferx_fit(..., optimizer_trace = TRUE)` now stores the per-iteration trace
  on the fit object itself as `fit$trace` (a data frame), not just its temp
  file path (`fit$trace_path`) (#228). `fit$trace` survives
  `ferx_save_fit()` / `ferx_load_fit()`, and `ferx_trace()`,
  `ferx_runlog()`, `ferx_runlog_iters()` all read it directly when present
  instead of re-reading a temp file that may have since been deleted.
  `fit$impmap_trace` is now only ever populated when `impmap_trace = TRUE`
  was actually requested (via `settings =` or `[fit_options]`), guarding
  against it leaking from an intermediate stage of a method chain.
  `ferx_job` handles (from `ferx_fit_async()`) gain a computed
  `trace_path` field alongside the existing `sidecar_path`.

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

## Fixed

- Datasets that reuse a subject ID in a non-contiguous block (e.g. a second
  cohort reusing IDs 12/13/14) are now handled correctly by the per-subject
  joins in `ferx_xpose()`, `ferx_save_fit()`, `fit$eta_cov`, and
  `ferx_cov_screen()`. ferx-core (like NONMEM) processes records sequentially,
  so each block is a distinct subject even when two share a textual ID; these
  functions previously keyed their joins on the raw ID, which silently gave the
  second block the first subject's ETAs/parameters (xpose), dropped it and wrote
  an `ebes.csv` that disagreed with `predictions.csv` so the loader rejected the
  bundle (save_fit), or double-weighted the shared covariate row (eta_cov /
  cov_screen). Joins now key on subject order instead (#252, alongside
  ferx-core #743).

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

- **Laplace estimator (`method = "laplace"`)**: alias `"laplacian"`. The Laplace
  approximation with the **exact** Hessian — NONMEM `$EST METHOD=1 LAPLACIAN`,
  which it reproduces to six significant figures. This is *not* the same estimator
  as `"focei"`, which builds its Gaussian from the Gauss-Newton Hessian and reports
  a different OFV. Internally it is `"agq"` with the node count pinned to 1 (a
  bit-identical OFV), making it the cheapest member of that family — on warfarin it
  converges *faster than FOCEI*. It supports IOV at any occasion count.
  (ferx-core #251)
- **AGQ and `laplace` now support inter-occasion variability (`[iov]`)**: the
  integral runs over the stacked `(eta, kappa_1..kappa_K)` vector. AGQ's grid grows
  with the occasion count (`n_agq^(n_eta + K*n_kappa)`) and is capped; `laplace` is a
  single node regardless, so it is always tractable under IOV. (ferx-core #251)
- **Adaptive Gaussian quadrature (`method = "agq"`)**: `ferx_fit(..., method =
  "agq")` (aliases `"aghq"`, `"gauss_hermite"`) selects the new AGQ estimator,
  with `settings = list(n_agq = 3)` setting the Gauss-Hermite nodes per random
  effect. AGQ generalises Laplace — instead of a single Gaussian at each
  subject's empirical-Bayes mode it evaluates the *exact* conditional likelihood
  on a Gauss-Hermite grid around that mode, so `n_agq = 1` reproduces Laplace
  identically and more nodes refine the marginal. Because it makes no
  Gaussian-residual assumption it covers non-Gaussian endpoints (time-to-event,
  categorical) that FOCE/FOCEI structurally cannot, and unlike SAEM/IMP its
  objective is deterministic (the OFV is bit-identical run to run). It carries an
  exact analytic outer gradient, so a converged warfarin fit is *faster* than
  FOCEI. Validated against NONMEM `$EST METHOD=1 LAPLACIAN`, which `n_agq = 1`
  reproduces to six significant figures. Cost is `n_agq^n_eta` per subject per
  iteration, so it suits models with few random effects. IOV is supported (see
  above).
  (ferx-core #251)
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
