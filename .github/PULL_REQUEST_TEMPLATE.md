<!-- Title format: type(scope): short description  [closes #N] -->
<!--
  type  : feat | fix | perf | refactor | docs | test | chore
  scope : fit | simulate | predict | output | install | api
  e.g.  : feat(output): expose CWRES column from sdtab in fit() result [closes #7]
-->

## Why
<!-- What problem does this solve for R users?
     Link to a ferx-core PR if this exposes a new engine capability. -->

## What changed
<!-- Function signatures, new arguments, output columns, parsing logic. -->

## Alternatives considered
<!-- What else was tried or evaluated, and why this approach won. Omit if obvious. -->

## Cross-repo dependency
| Repo | PR | Status | Must merge first |
|------|----|--------|-----------------|
| ferx-core (engine) | FeRx-NLME/ferx-core#___ | open / merged / not needed | yes / no / — |

<!-- If the engine PR is not yet merged, keep this PR as Draft. -->

## Breaking changes (R users)
- [ ] Function signature changed
- [ ] Output column names / types changed (add migration note below)
- [ ] New minimum ferx-core version required — `DESCRIPTION` git dep updated to commit `______`
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
- [ ] Test against real ferx-core binary (not mocked)

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
- [ ] `DESCRIPTION` ferx-core git dep points to merged engine commit (not a branch)
- [ ] `src/rust/Cargo.lock` updated to pin ferx-core to the correct commit
  <!-- Run from src/rust/ with the local [patch] override removed:
       cp .cargo/config.toml .cargo/config.toml.bak
       printf '[build]\n' > .cargo/config.toml
       cargo fetch
       cp .cargo/config.toml.bak .cargo/config.toml
       Then commit Cargo.lock. Verify with:
       grep -A2 'name = "ferx-core"' Cargo.lock | grep source -->
- [ ] Not applicable (no ferx-core dependency change in this PR)

## Docs & examples

### ferx-site
- [ ] New function / argument / output field → example page exists in `ferx-site/examples/` or a ferx-site PR is opened
- [ ] New `inst/examples/*.R` script → matching `.qmd` page in `ferx-site/examples/` covers it
- [ ] New bundled model in `inst/examples/models/` → referenced in at least one ferx-site example page
- [ ] New bundled data in `inst/examples/data/` → listed in ferx-site reference page and ferx-book
- [ ] `roxygen2` `@examples` blocks match the standalone scripts in `inst/examples/`

### ferx-book
- [ ] New user-visible feature → relevant book chapter updated or a ferx-book PR is opened
- [ ] New bundled example → `ferx_example("name")` call added/updated in the relevant chapter

### Example execution (run locally before marking ready for review)
- [ ] Rebuilt ferx-r against local ferx-core: `FERX_NO_AUTODIFF=1 R CMD INSTALL .` (triggers Rust compile without autodiff)
- [ ] All affected `inst/examples/*.R` scripts run without error on the local build
- [ ] All affected ferx-site example `.qmd` pages render cleanly: `quarto render examples/<page>.qmd`
- [ ] All affected ferx-book chapters render cleanly: `quarto render chapters/<chapter>.qmd`
- [ ] No example execution step needed (internal refactor / docs-only / no R API surface changed)

## Reviewer hints
<!-- Where to focus. What's subtle. What can be skimmed. -->

## Open questions
<!-- Things you're uncertain about and want input on. -->
