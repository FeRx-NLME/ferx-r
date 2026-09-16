# Review: ferx-core #520 — loose default ODE tolerances and the covariance step

Tracking: [FeRx-NLME/ferx-core#520](https://github.com/FeRx-NLME/ferx-core/issues/520)
"Loose default ODE tolerances can yield ill-conditioned covariance
(regularized/unreliable SEs) on ODE models".

Review posted on the issue:
[#520 (comment)](https://github.com/FeRx-NLME/ferx-core/issues/520#issuecomment-5694083528).

Reviewed 2026-09-16 against ferx-core `main` @ `d66046e3` (2026-09-15), which
is also the revision `src/rust/Cargo.lock` pins for this package, so every
finding below about the engine applies to the ferx-r build as shipped. This is
an evaluation only; no engine or package change is proposed here beyond the
follow-up list at the end.

---

## 1. Where the issue stands

- Filed 2026-06-25 by hiddevandebeek out of the NONMEM cross-check for the
  analytic ODE lagtime gradient (ferx-core PR #472). Open, assigned to
  roninsightrx.
- One comment, 2026-09-15, by roninsightrx: a re-measurement on `main` @
  `d66046e3` with a remediation proposal (section 4 below).
- Triaged 2026-09-11 into the 0.4.0 rollout epic ferx-core #1350 as Tier 0
  row 3, glyph "decision needed" between a tighter default, a
  covariance-step auto-tighten, or a targeted warning. The row notes that
  #1291 "helps SEs in scope, not the OFV".
- Cited by ferx-core #960 (2026-08-04, knife-edge FD Hessian on warfarin) as
  one of three fix routes; #960 closed when the analytic R-matrix went
  default-on.
- Not mentioned by the #486 analytic-sensitivity umbrella or the #1217 bug
  order. #1217 does record one trap for reproducing #520 (see #1221 below).

Engine changes since filing that bear on it:

| Date | Change | Effect on #520 |
|---|---|---|
| 2026-07-29 | #953 (closes #436): exact analytic covariance R-matrix for closed-form FOCE/FOCEI, `analytic_cov_hessian = true` by default | None for `[odes]` models: #436's CHANGELOG entry lists ODE first among the out-of-scope models that keep the FD covariance |
| 2026-09-03 | #1221 (closes #1212): ODE settings passed on a Rust `FitOptions` now reach the solver, `run_covariance` included | Any Rust-side reproduction before this date silently ran at the default; model-file `[fit_options]` was always honoured |
| 2026-09-08 | #1291: analytic covariance for in-scope `[odes]` models (FOCE, FOCEI, FOCEI-anchored AGQ; IOV and M3 included) | The FD-of-OFV Hessian is no longer used for the issue's model class |
| 2026-09-11 | #1350 triage | Decision pending |

The default is unchanged: `ode_reltol = 1e-4`, `ode_abstol = 1e-6`
(`src/types.rs`, `FitOptions::default`). It was inherited from a hard-coded
constant when #127 made it configurable on 2026-06-14, and that commit kept it
"so existing fits are unchanged". It was never chosen on accuracy grounds.

## 2. The issue's claims, scrutinised

**The mechanism holds up.** The FD stencil perturbs by
`fd_hessian_step * (1 + |x|)` with `fd_hessian_step = 1e-2` and divides by
`h²`, so objective noise is amplified by roughly `4e4`. A resulting
non-positive-definite Hessian goes through `invert_psd_with_floor`, which
floors eigenvalues at `max_eig * 1e-10`; the SE that comes out of a floored
direction is of order `1/sqrt(floor)`. That is exactly the "%RSE in the
thousands" signature, and it is the same artifact ferx-core #1021 hit in SIR
proposals.

**Independent corroboration exists.** The #127 commit body says an ODE-form
OFV can differ from its analytical twin by several units at the default, up to
about 15 on `three_cpt_iv`. `docs/model-file/absorption.qmd` reports transit
omega-squared inflated at the default and tells users to tighten toward
`1e-9`. #960 found the FD Hessian fragile even without integrator noise.

**Project practice already contradicts the default.**

| Corpus | Count |
|---|---|
| ferx-core NONMEM anchor control streams, all `TOL=9` | 83 of 83 |
| ferx-core `[odes]` model files pinning `1e-9` or `1e-10` | 34 of 75 |
| ferx-r bundled `[odes]` examples pinning `1e-10` | 9 of 28, and every standard PK one |

**The evidence is thin.** One simulated model. The tight-tolerance run is
recorded in `docs/model-file/lagtime.qmd`; the default-tolerance run that
motivates the issue has no numbers, model, or data anywhere in the repo. The
reproduction section is generic. The tight settings used (`1e-12` / `1e-14`)
demonstrate convergence, not what a sensible default is.

**The framing understates the problem.** The issue says point estimates are
fine and only SEs suffer. But the reported 1.9-unit OFV gap is a perturbed
optimum, the absorption docs show variance components off by about 15 %, and
the `outer_ftol` auto default was pinned at `1e-6` rather than `1e-8`
precisely because ODE noise makes the tighter stop unreachable. The default
affects estimation, not just covariance. "Optimization averages over noise"
is not how BOBYQA or the gradient optimizers behave; noise sets a convergence
floor.

**The proposed absolute tolerance ignores units.** `ode_abstol` applies to
state amounts, so a fixed `1e-11` means different things in milligrams and in
nanograms.

## 3. What the analytic ODE route does at the default

For in-scope models the third-order blocks come from central differences of
the second-order sensitivity jet. `src/sens/provider.rs::third_order_fd_step`
sizes the step as `cbrt(reltol).min(1e-2) * (1 + |x|)`:

| `ode_reltol` | cube root | step used |
|---|---|---|
| `1e-4` (default) | 0.046 | **0.01, cap binds** |
| `1e-6` | 0.010 | 0.010 |
| `1e-9` | 0.001 | 0.001 |

So at the default the step is the same 1 % the FD stencil uses, the cube-root
balance the design rests on is not in effect, and `ode_abstol` is not
consulted. The CHANGELOG line "the ODE step accounts for `ode_reltol`" is true
below `1e-6` only. PR #1291's self-review listed the cap, the missing `abstol`
term, and the lag-arrival jump (an observation within `±h` of a lagged arrival
is differenced across the jump) as items to resolve before merge; there is no
follow-up commit or issue for any of the three. The only default-tolerance
test on this path is a one-state IV fixture asserting agreement within 2 % of
a `1e-10` reference (`ode_cov_hessian_is_stable_at_default_solver_tolerance`).
No test runs the analytic ODE covariance on a lag model.

Scope, for reference: a bare estimated lagtime on an ODE **is** admitted on the
non-IOV path (the covariance gate only excludes lag under IOV, and
`ode_analytical_supported` only restricts lag for `Weibull` forcings). Models
that still take the FD-of-OFV Hessian: `method = laplace` (exact anchor needs
fourth-order sensitivities), non-Gaussian endpoints (TTE, categorical, CTMM),
mixtures, `gradient = fd`, LTBS, any `[scaling]` other than the Form C `y =`
readout (see section 4), analytic readouts, `iiv_on_ruv`, custom or
correlated residuals, FREM, covariate-selected error models.

## 4. The 2026-09-15 re-measurement, scrutinised

What roninsightrx measured and proposed:

- The lag anchor `nonmem_anchor/first_order_alag_fit.ferx` (30 subjects) is
  clean and agrees to 4-5 figures across `1e-4` to `1e-12` on the default
  analytic path. Attributed to #436.
- With `analytic_cov_hessian = false`, `three_cpt_iv_ode` at 180 subjects
  regularises at the default (`min eig = -3.2e1`, SE(TVV2) 5829, %RSE 24982)
  and is clean one order of magnitude tighter.
- Cost on that model: 1.7x at `1e-6` / `1e-8`, 2.5x at `1e-8`, 5.5x at
  `1e-10`; past ~`1e-10` the solve becomes stability-limited (minimum-step
  clamps).
- Shipped examples on the FD route: `warfarin_ode_lagtime.ferx`,
  `emax_pkpd.ferx`; `mm_oral.ferx` is analytic.
- Defect: regularisation severity is graded on the *count fraction* of
  clipped eigenvalues (`covariance.rs` ~1308), so 1 of 13 clipped printed
  "minor. Standard errors are likely reliable." above a 4400x-inflated SE.
- Proposal: **C1** grade severity by magnitude; **C2** name
  `ode_reltol`/`ode_abstol` when the FD route regularised on an ODE model at
  loose tolerance; then **A** global default `1e-6` / `1e-8`, or **B** a
  `cov_ode_reltol` / `cov_ode_abstol` pair defaulting to
  `min(ode_reltol, 1e-8)` applied only when the analytic R-matrix is not in
  use. Recommends C + B.

Where the code or history reads differently:

1. **Dating.** For `[odes]` models the analytic route opened with #1291 on
   2026-09-08, not #436 on 2026-07-29; #436's changelog lists ODE as out of
   scope and #1291's body opens with "previously rejected every `[odes]`
   model". #1350 row 3 has it right.
2. **Different model.** The issue's model is `ode(states=[depot, central])`,
   60 subjects, IIV on CL/V/KA, a state jump at the lagged arrival. The
   re-measurement used the one-state `first_order(ka)` forcing anchor, 30
   subjects, IIV on CL/KA, a rate-on saltation at arrival. Both integrate: the
   indexed `ALAG1` populates `dose_attr_map`, which keeps `mr_scope` from
   serving the anchor in closed form. But the depot form is the one the
   #1291 self-review flagged for differencing across the jump, and nothing
   in the table speaks to it. Whichever fixture becomes the regression test
   should be the two-state depot form and should keep the indexed `ALAG1` on
   purpose: a bare `ALAG` / `LAGTIME` on a `first_order` forcing is
   `mr_scope`-eligible, is served in closed form, and makes a tolerance sweep
   trivially flat.
3. **Why the analytic route is flat at the default** is section 3 above: by
   construction of the capped step, and validated on one IV fixture plus this
   one lag measurement. The measurement should be committed as a test.
4. **Two of the three FD examples are FD for other reasons.**
   `emax_pkpd.ferx` sets `gradient = fd` itself (its header says the per-CMT
   Form C readout requires it), which is the gate's explicit opt-out.
   `warfarin_ode_lagtime.ferx` is FD because of `[scaling] obs_scale = V`,
   which parses to `ScalingSpec::ExpressionScale` and trips the covariance
   gate's `!matches!(model.scaling, ScalingSpec::None)` clause
   (`provider.rs` ~5282). It is not the lagtime and not the IIV on it. The
   same model written with `[scaling] y = central / V` should be admitted.
   Two consequences: admitting `ExpressionScale` to the analytic ODE
   covariance (the gradient path already differentiates it via its `deriv`
   program, #367) would move the `obs_scale` population off the FD route
   without touching tolerances; and C2 should name *which gate clause* put
   the fit on the FD route, because an `obs_scale` user can rewrite the
   readout and be done.
5. **The `1e-6` plateau is shown on one non-stiff, closed-form-equivalent
   model.** The absorption docs' transit case and #127's `three_cpt_iv` OFV
   gap both point tighter. Confirm on a transit fixture, watching OFV and
   omega-squared rather than SE, before `1e-6` becomes a default or a
   covariance floor.
6. **"Bias, not ill-conditioning" understates it.** The issue's 1.9-unit OFV
   gap is a perturbed optimum; the transit omega bias is an error in a
   reported parameter. Option B leaves both exactly where they are.
7. **C1 has an R-side twin** (section 5 below).
8. **A vs B.** B is mechanism-targeted, has precedent in `cov_inner_tol`, and
   is defensible against #1221's "same surface" rationale because a central
   second difference does not need stationarity; but its `min(ode_reltol,
   1e-8)` and A's `1e-6` are two different plateau claims from one
   measurement, and B applies only off the analytic route, so the capped step
   in section 3 stays. A at `1e-6` costs the measured 1.7x (rk45 cost;
   `ode_method = auto` does not select `vern7` by tolerance), moves every
   default-tolerance fixture, and is the only option that touches the OFV /
   omega problem. Honest framing: C1 + C2 fix the *diagnostic* failure; B
   fixes the *SE* failure on the FD route cheaply; only A (or A + B)
   addresses what the issue's evidence actually measured, which was an OFV
   gap.

## 5. ferx-r follow-ups

Nothing here is blocked on the engine decision except the last item.

- **Paired change for C1.** `R/ferx_get_warnings.R::.ferx_warning_guidance()`
  mirrors ferx-core's count-graded tiers and, for "minor", adds "this is
  common on smooth OFV surfaces and is usually benign". When ferx-core grades
  severity by magnitude, the R tiers must follow, and the "minor" wording
  must stop asserting benignity it cannot know. Until then the same inflated
  SE is called benign twice.
- **Docs.** The `ode_reltol` entry in `?ferx_fit` (`R/ferx_fit.R`, mirrored
  in `man/ferx_fit.Rd`) says the default affects the OFV against an
  analytical reference and nothing about standard errors or the FD route.
  Same gap in ferx-core's `docs/model-file/fit-options.qmd` and
  `docs/warnings.qmd` `covariance_regularized` row.
- **Cosmetic, both sides.** The regularisation message says "eigenvalue floor
  applied to FD Hessian" on the analytic route too, and the
  `ill-conditioned entries` branch recommends tuning `fd_hessian_step` when
  no stencil was used. R guidance repeats the `fd_hessian_step` advice.
- **`obs_scale` in bundled examples.** 9 of the 28 `[odes]` examples under
  `inst/examples/models/` use `obs_scale`, so a third of the shipped ODE
  examples take the FD covariance route for the reason in section 4 item 4.
  If ferx-core admits `ExpressionScale`, nothing changes here; if it does
  not, consider rewriting those readouts as `y = <expr>` where the two are
  equivalent, and say so in `?ferx_example`.
- **Pin.** The current pin already includes #1291 and #1221, so R users get
  the analytic ODE covariance and honoured call-time tolerances now. A
  default change (option A) or a new `cov_ode_reltol` key (option B) would
  need the usual `[fit_options]` key plumbing check in `ferx_fit()` and a
  `NEWS.md` entry.

## 6. Suggested order (as posted on the issue)

1. C1 + C2 (with the declining gate clause named) + the three doc rows, plus
   the paired ferx-r guidance change. Independent of the tolerance decision.
2. Commit the lag measurement as a regression test on the two-state depot
   model, analytic on, default tolerance, asserting SE parity with the
   `1e-9` run.
3. Re-run a transit fixture at `1e-6` / `1e-8` for OFV and omega-squared,
   not SE, to confirm the plateau.
4. Decide A vs B on that evidence. Separately, consider admitting
   `ExpressionScale` to the analytic ODE covariance to shrink the FD
   population.
