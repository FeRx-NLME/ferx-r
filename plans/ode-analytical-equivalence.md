# Standard models: analytical + ODE form, with equivalence verification

Tracking: [FeRx-NLME/ferx-r#127](https://github.com/FeRx-NLME/ferx-r/issues/127)
"make the standard analytical models also available in their ode form" -
all standard analytical models as an ODE version.

Branch (worktree): `worktree-ode-analytical-equiv`.

## Goal

1. Every standard analytical PK model type has an example model **and** an
   equivalent **ODE-form** example, all backed by a dataset.
2. A verification that each ODE form produces the **same predictions** as its
   analytical counterpart - in **both** ferx-core (Rust integration test, the
   source of truth) and ferx-r (testthat over the shipped example pairs).

Scope decided with user: all 6 standard model types need data; each gets an
analytical example and an ODE example; all combinations (model x dosing mode)
must run.

---

## Phase 1 - Inventory of analytical models

The analytical PK engine lives in ferx-core (`src/pk/`, dispatched in
`src/pk/mod.rs`; closed forms in `one_compartment.rs`, `two_compartment.rs`,
`three_compartment.rs`). The `PkModel` enum (`src/types.rs:700`) is the
complete set - **6 model types**:

| `PkModel`        | DSL name (`pk ...`)              | Required params                  | Compartments |
|------------------|----------------------------------|----------------------------------|--------------|
| `OneCptIv`       | `one_cpt_iv(cl, v)`              | cl, v                            | central |
| `OneCptOral`     | `one_cpt_oral(cl, v, ka)`       | cl, v, ka                        | depot, central |
| `TwoCptIv`       | `two_cpt_iv(cl, v1, q, v2)`     | cl, v1, q, v2                    | central, periph |
| `TwoCptOral`     | `two_cpt_oral(cl, v1, q, v2, ka)` | cl, v1, q, v2, ka              | depot, central, periph |
| `ThreeCptIv`     | `three_cpt_iv(cl, v1, q2, v2, q3, v3)` | cl, v1, q2, v2, q3, v3    | central, periph1, periph2 |
| `ThreeCptOral`   | `three_cpt_oral(cl, v1, q2, v2, q3, v3, ka)` | + ka                | depot, central, periph1, periph2 |

Common, model-independent dosing features (all flow through the same engine for
analytical and ODE paths):
- **Bolus vs infusion**: driven by `RATE` per dose (IV models accept both;
  `DoseEvent::is_infusion`). No separate model variant.
- **Bioavailability `f`** (oral): optional, default 1.0. Engine applies
  `F * AMT` at dose entry (and `F * RATE` for infusions) - issue #122.
- **Lag time** (`lagtime`/`alag`): optional, default 0.0. Engine shifts the
  effective dose time.

Reserved-slot routing (`ode_param_slots`, `model_parser.rs:4248`): for ODE
models, an individual parameter named `f`/`lagtime`/`alag` (case-insensitive,
via `name_to_index(name.to_lowercase())`) is routed to the engine-reserved
`PK_IDX_F` / `PK_IDX_LAGTIME` slots and applied **by the engine at the dose** -
*not* via the `[odes]` RHS.

### Existing example coverage (ferx-r `inst/examples/models/`)

| Model type      | Analytical example                  | Dataset                | ODE example today |
|-----------------|-------------------------------------|------------------------|-------------------|
| one_cpt_iv      | **none**                            | **none**               | none |
| one_cpt_oral    | `warfarin.ferx`                     | `warfarin.csv`         | partial (`warfarin_ode_time` adds extra states; not a clean 1:1) |
| two_cpt_iv      | `two_cpt_iv.ferx`                   | `two_cpt_iv.csv`       | none |
| two_cpt_oral    | `two_cpt_oral_cov.ferx`             | `two_cpt_oral_cov.csv` | none |
| three_cpt_iv    | `three_cpt_iv.ferx`                 | `three_cpt_iv.csv`     | none |
| three_cpt_oral  | **none**                            | **none**               | none |
| (oral + F)      | `bioavailability.ferx`              | `bioavailability.csv`  | `bioavailability_ode.ferx` - **stale, see findings** |

Gaps to fill: analytical examples + datasets for **one_cpt_iv** and
**three_cpt_oral**; ODE siblings for all 6.

---

## Phase 2 - Transformation rules (analytical -> ODE)

### Unit convention (decision)

Use a **uniform amount-based** convention for every ODE example:
- All compartment states are **amounts**; doses add amount directly
  (`u[cmt-1] += F * AMT`).
- Output concentration via the `[scaling]` block: **`obs_scale = V`** (or `V1`
  for multi-compartment) - the observed central *amount* is divided by the
  central volume. (Equivalently `[scaling] y = central / V`, "Form C".)

Rationale: this is the only convention that works **identically for IV and
oral** (an IV bolus adds amount to `central`, so `central` must be in amount
units), it matches NONMEM's `S2 = V` idiom, and it keeps F / lag / infusion
handling entirely in the engine. The existing concentration-baked oral examples
(`bioavailability_ode`, `warfarin_ode_time`) use a different per-model
convention; we standardize.

### Per-model ODE templates

Micro-constants: `k10 = CL/V1`, `k12 = Q2/V1`, `k21 = Q2/V2`,
`k13 = Q3/V1`, `k31 = Q3/V3` (single `Q` -> `q2/v2` for 2-cpt).

**one_cpt_iv**
```
ode(obs_cmt=central, states=[central])
d/dt(central) = -(CL/V) * central
[scaling] obs_scale = V
```

**one_cpt_oral** (F via `F` indiv param; *not* in the flux)
```
ode(obs_cmt=central, states=[depot, central])
d/dt(depot)   = -KA * depot
d/dt(central) =  KA * depot - (CL/V) * central
[scaling] obs_scale = V
```

**two_cpt_iv**
```
ode(obs_cmt=central, states=[central, periph])
d/dt(central) = -(CL/V1 + Q/V1) * central + (Q/V2) * periph
d/dt(periph)  =  (Q/V1) * central - (Q/V2) * periph
[scaling] obs_scale = V1
```

**two_cpt_oral** = two_cpt_iv + `depot` with `-KA*depot` / `+KA*depot` into central.

**three_cpt_iv**
```
ode(obs_cmt=central, states=[central, periph1, periph2])
d/dt(central)  = -(CL/V1 + Q2/V1 + Q3/V1)*central + (Q2/V2)*periph1 + (Q3/V3)*periph2
d/dt(periph1)  =  (Q2/V1)*central - (Q2/V2)*periph1
d/dt(periph2)  =  (Q3/V1)*central - (Q3/V3)*periph2
[scaling] obs_scale = V1
```

**three_cpt_oral** = three_cpt_iv + `depot` absorption into central.

### Transformation rules summary

| Analytical feature        | ODE handling |
|---------------------------|--------------|
| `cl, v[1]`                | central elimination `CL/V[1]` |
| `q[2], v2; q3, v3`        | inter-compartmental micro-constants (table above) |
| `ka`                      | `depot` compartment, `-KA*depot` out / `+KA*depot` into central |
| `f=F` (bioavailability)   | declare `F` as individual parameter; engine applies at dose. **Do NOT bake F into the flux.** |
| `lagtime`                 | declare `lagtime`/`alag` as individual parameter; engine applies. |
| bolus / infusion / SS / ADDL | unchanged - same `[odes]`, driven by data columns. |
| concentration readout     | `[scaling] obs_scale = V` (or `V1`). |

---

## Phase 3 - Verification (ODE == analytical)

The engine's public `predict()` (population PRED, eta=0) is the comparison
point; the established pattern already exists in
`ferx-core/tests/bioavailability_ode_nonmem.rs::ode_bioavailability_matches_analytical`:
parse both models, `predict()` the same population, assert pointwise relative
error.

**Tolerance**: RK45 defaults are `abstol = 1e-6`, `reltol = 1e-4`
(`src/ode/solver.rs:59`). The existing analytical-vs-ODE test passes at
**`rel < 1e-4`**; use that as the headline threshold (tightening solver
tolerances in the test can reach ~1e-6 if needed).

### 3a. ferx-core: `tests/analytical_ode_equivalence.rs` (source of truth)

Parametric integration test. For each of the **6 model types**, embed the
analytical and amount-based ODE model strings (etas fixed to 0 -> population
PRED; plus one case with non-zero eta exercising IPRED), build synthetic
populations covering **all dosing combinations**:

- single **bolus**
- **infusion** (RATE > 0)
- **multiple doses**
- **steady state** (SS + II)
- **lag time** (lagtime > 0)
- **bioavailability** F < 1 (oral models)

Assert `predict()` agrees pointwise (`rel < 1e-4`) between analytical and ODE.
Runs in the default/CI feature set (no autodiff, FD-safe) and fast (no `fit()`),
so it is not gated behind `slow-tests`. This is the rigorous, exhaustive check.

### 3b. ferx-r: `tests/testthat/test-ode-analytical-equivalence.R`

For each shipped analytical/ODE example **pair**, load both via
`ferx_example()`, run `ferx_predict()` on the shared dataset, and assert the
`PRED` columns agree (`expect_equal(..., tolerance = 1e-4)`). This directly
guards the **shipped example files** (which the core test, using embedded
strings, does not), and catches regressions like the stale `bioavailability_ode`
(below).

---

## Deliverables / file changes

ferx-r (`inst/examples/`):
- New analytical examples + datasets: `one_cpt_iv.ferx` + `one_cpt_iv.csv`,
  `three_cpt_oral.ferx` + `three_cpt_oral.csv` (datasets simulated from the
  model so the examples also fit sensibly).
- New ODE siblings for all 6: `one_cpt_iv_ode.ferx`, `one_cpt_oral_ode.ferx`,
  `two_cpt_iv_ode.ferx`, `two_cpt_oral_ode.ferx`, `three_cpt_iv_ode.ferx`,
  `three_cpt_oral_ode.ferx` - each reusing its analytical sibling's dataset via
  the `.data_aliases` map in `R/example.R:153`.
- Runnable `ex_*.R` scripts for the new examples (matching existing style), and
  update the "Available bundled examples" roxygen list in `R/example.R`.
- `tests/testthat/test-ode-analytical-equivalence.R` (Phase 3b).
- **Fix** `inst/examples/models/bioavailability_ode.ferx` (finding below) and
  reconcile / re-point `warfarin_ode_time` if we want it as the clean 1cpt-oral
  ODE pair (or leave it as the TIME/TAFD demo and add a separate clean pair).
- `NEWS.md` entry.

ferx-core:
- `tests/analytical_ode_equivalence.rs` (Phase 3a). No engine code change
  expected - this is verification only.

Docs: note the ODE-form availability in the model-file docs if the analytical
models are documented there.

---

## Findings / risks discovered during evaluation

1. **`bioavailability_ode.ferx` double-counts F (latent bug).** It declares `F`
   as an individual parameter *and* bakes `F*KA*depot/V` into the central flux.
   Post-issue #122 the engine already applies `F*AMT` at dose entry for ODE
   models (reserved `PK_IDX_F` slot, `event_driven.rs:313`,
   `ode/predictions.rs:117`), so F is applied twice -> central is ~F^2 of
   correct. It is not currently covered by any test (the NONMEM test embeds its
   own corrected string). Fix to plain flux as part of this work; the new R test
   would have caught it.

2. **Two unit conventions exist in shipped ODE examples** - concentration-baked
   (`bioavailability_ode`, `warfarin_ode_time`) vs amount + `obs_scale`. We
   standardize on amount + `obs_scale=V/V1` for the new equivalence pairs;
   decide whether to also migrate the two legacy examples for consistency.

3. **Tolerance is solver-bound** (~1e-4). If a stricter bound is wanted, tighten
   `abstol/reltol` in the test harness rather than loosening the assertion.

4. **`warfarin_ode_time` is not a clean 1cpt-oral pair** (it adds `AUC_D1` /
   `TAM_INT` states for a TIME/TAFD/TAD demo). Plan adds a dedicated
   `one_cpt_oral_ode.ferx` rather than overloading it.

---

## Suggested implementation order

1. ferx-core `tests/analytical_ode_equivalence.rs` (all 6 x dosing modes) -
   proves the transforms are correct before shipping any example file.
2. Create the 2 missing analytical examples + simulated datasets.
3. Create the 6 ODE sibling examples (+ `.data_aliases`, `ex_*.R`, roxygen).
4. Fix `bioavailability_ode.ferx`.
5. ferx-r `test-ode-analytical-equivalence.R`.
6. Roxygenize, ASCII check, `NEWS.md`, PR (fill template).
