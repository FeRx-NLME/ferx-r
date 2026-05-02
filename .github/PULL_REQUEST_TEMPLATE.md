<!-- Title format: type(scope): short description  [closes #N] -->
<!--
  type  : feat | fix | perf | refactor | docs | test | chore
  scope : fit | simulate | predict | output | install | api
  e.g.  : feat(output): expose CWRES column from sdtab in fit() result [closes #7]
-->

## Why
<!-- What problem does this solve for R users?
     Link to a ferx-nlme PR if this exposes a new engine capability. -->

## What changed
<!-- Function signatures, new arguments, output columns, parsing logic. -->

## Alternatives considered
<!-- What else was tried or evaluated, and why this approach won. Omit if obvious. -->

## Cross-repo dependency
| Repo | PR | Status | Must merge first |
|------|----|--------|-----------------|
| ferx-nlme (engine) | InsightRX/ferx-nlme#___ | open / merged / not needed | yes / no / — |

<!-- If the engine PR is not yet merged, keep this PR as Draft. -->

## Breaking changes (R users)
- [ ] Function signature changed
- [ ] Output column names / types changed (add migration note below)
- [ ] New minimum ferx-nlme version required — `DESCRIPTION` git dep updated to commit `______`
- [ ] None

<details>
<summary>Migration notes (if breaking)</summary>

<!-- What do users need to change in their R scripts? Before/after example. -->

</details>

## Output changes
<!-- For changes that affect what fit(), simulate(), or predict() return:
     paste a before/after of the relevant list elements or data frame columns. -->

<details>
<summary>Before / after (if applicable)</summary>

```r
# before

# after
```

</details>

## Tests
- [ ] `testthat` tests added / updated
- [ ] `R CMD check` passes (no ERRORs or WARNINGs)
- [ ] Test against real ferx-nlme binary (not mocked)

## Example (user-facing features only)
<!-- If this adds or changes any user-visible behaviour (new argument, new output column,
     new function), add a minimal reproducible R example. This becomes the basis for
     the roxygen @examples block and can be copy-pasted by users. -->
- [ ] `@examples` block added / updated in the relevant roxygen docs
- [ ] Standalone example added to `inst/examples/` or inline below
- [ ] Not applicable (internal / refactor / fix with no new API surface)

<details>
<summary>Example (if not a standalone file)</summary>

```r
# minimal reproducible example demonstrating the new feature

```

</details>

## Docs
- [ ] `roxygen2` docs updated (`devtools::document()` run)
- [ ] Vignette updated if relevant
- [ ] `NEWS.md` entry added

## Checklist
- [ ] `DESCRIPTION` version bumped if warranted
- [ ] `DESCRIPTION` ferx-nlme git dep points to merged engine commit (not a branch)

## Reviewer hints
<!-- Where to focus. What's subtle. What can be skimmed. -->

## Open questions
<!-- Things you're uncertain about and want input on. -->
