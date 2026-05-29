# ferx

[![codecov](https://codecov.io/gh/FeRx-NLME/ferx-r/branch/main/graph/badge.svg)](https://codecov.io/gh/FeRx-NLME/ferx-r)

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend with [Enzyme](https://enzyme.mit.edu/) automatic differentiation for exact gradients.

## Features

- **FOCE/FOCEI estimation** with automatic differentiation
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **BLOQ handling**: Beal's M3 likelihood for observations below the LLOQ
- **Model DSL**: define models in `.ferx` text files

## Installation

A Rust installation with the Enzyme AutoDifferentiation engine is required for FeRx to compute gradients. Most likely you will need to build Rust from source. Expect a few hours end-to-end the first time: ~45-60 min to build the Enzyme toolchain, plus ~1-2 hours for the first ferx-r compile.

See [documentation](https://ferx-nlme.github.io/ferx-core/installation.html) for installation instructions.

> **Heads-up - two gotchas that bite people building the toolchain from source:**
> - If the autodiff verification step makes `rustc` *hang* (one core busy, memory flat, never returns), you have an LLVM/Enzyme mismatch. Rebuild with `--set llvm.download-ci-llvm=false`.
> - The toolchain is initially linked into `/tmp`, which macOS/Linux purge - relocate it to a permanent path (e.g. `~/.local/share/enzyme-toolchain`) and re-link, or it will silently break later. See the install docs for the exact steps.

### Install the package

After installing Rust, in R run:

```r
pak::pak("FeRx-NLME/ferx-r")
```

Or from a local clone:
```bash
R CMD INSTALL .
```

### Mac/Linux without Enzyme (finite-difference gradients)

No Enzyme toolchain? No problem. The build auto-detects whether Enzyme is present; if it is not, it falls back to finite-difference gradients automatically. Just install normally:

```r
pak::pak("FeRx-NLME/ferx-r")
```

To skip the probe and force the FD path explicitly (e.g. in scripts or CI):

```r
Sys.setenv(FERX_NO_AUTODIFF = "1")
pak::pak("FeRx-NLME/ferx-r")
```

First install takes ~1-2 hours (Rust compilation from scratch). Subsequent installs are fast.

All estimation methods (FOCE/FOCEI/SAEM/IMP) work. When ferx loads it prints a startup message listing the specific AD-only features that are limited.

### Windows (supported with finite-difference gradients)

Native Windows installs are supported. The build automatically uses finite-difference gradients on Windows (no `FERX_NO_AUTODIFF=1` needed) - all estimation methods work, but AD-only features are limited. Native Enzyme autodiff is **not** available on Windows; for that, use the Docker image below (or WSL2).

Prerequisites:

- [R](https://cran.r-project.org/bin/windows/base/) and [Rtools44](https://cran.r-project.org/bin/windows/Rtools/) (Rtools44 ships the MinGW gcc that R uses to link the package)
- [rustup](https://www.rust-lang.org/tools/install) with the GNU-ABI toolchain:

```powershell
rustup toolchain install stable-x86_64-pc-windows-gnu
```

You do **not** need to `rustup default` it — the package's build pins this toolchain automatically on Windows. (The rustup default on Windows is the MSVC ABI, which is not link-compatible with Rtools' MinGW linker.)

Then install in R:

```r
pak::pak("FeRx-NLME/ferx-r")
```

When ferx loads it prints a startup message listing the specific AD-only features that are limited on Windows.

### Docker

A Docker image is available that bundles the Enzyme toolchain (built from source), ferx CLI, the ferx R package, and RStudio Server — no local Rust/Enzyme setup required. On Windows, this is the recommended path if you need Enzyme autodiff.

```bash
# Build (first build takes ~45-60 min; cached after that)
docker build -t ferx:latest .

# Run RStudio Server
docker run --rm -p 8787:8787 -e PASSWORD=ferx ferx:latest
# -> http://localhost:8787   user: rstudio   password: ferx
```

## Quick Start

```r
library(ferx)

# Get bundled example paths
ex <- ferx_example("warfarin")

# Fit a one-compartment oral PK model
result <- ferx_fit(ex$model, ex$data, method = "focei")
result

# Simulate at the fitted estimates (typical VPC flow)
sim <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42, fit = result)

# Population predictions at the fitted estimates
preds <- ferx_predict(ex$model, ex$data, fit = result)
```

Pass `fit = <ferx_fit result>` to `ferx_simulate()` / `ferx_predict()` to use
the fitted theta / omega / sigma. Omit it to use the model file's initial
values.

## BLOQ handling (M3 method)

For observations below the lower limit of quantification, flag them with a
`CENS` column in the data (1 = censored, with `DV` carrying the LLOQ value)
and pass `bloq_method = "m3"` to `ferx_fit()`. Each censored observation then
contributes `P(y < LLOQ | θ, η) = Φ((LLOQ − f)/√V)` to the likelihood instead
of a Gaussian residual, avoiding the terminal-phase bias that comes from
simply dropping BLOQ rows.

```r
bloq <- ferx_example("warfarin_bloq")
result <- ferx_fit(bloq$model, bloq$data, method = "focei", bloq_method = "m3")
sim <- ferx_simulate(bloq$model, bloq$data, n_sim = 100, seed = 42, fit = result)
```

See `inst/examples/ex1a_warfarin_bloq.R` for a full fit + VPC walkthrough.

## Model Specification

Models are defined in `.ferx` files:

```
[parameters]
  theta TVCL(0.2, 0.001, 10.0)   # name(initial, lower, upper)
  theta TVV(10.0, 0.1, 500.0)
  theta TVKA(1.5, 0.01, 50.0)

  omega ETA_CL ~ 0.09            # between-subject variability (variance)
  omega ETA_V  ~ 0.04
  omega ETA_KA ~ 0.30

  sigma PROP_ERR ~ 0.02

[individual_parameters]
  CL = TVCL * exp(ETA_CL)
  V  = TVV  * exp(ETA_V)
  KA = TVKA * exp(ETA_KA)

[structural_model]
  pk one_cpt_oral(cl=CL, v=V, ka=KA)

[error_model]
  DV ~ proportional(PROP_ERR)
```

See `ferx_example()` for available bundled examples.

## API Reference

| Function | Description |
|---|---|
| `ferx_fit()` | Fit a NLME model (FOCE/FOCEI). `bloq_method = "m3"` enables M3. |
| `ferx_simulate()` | Simulate replicates with BSV and residual error. Pass `fit =` to use fitted estimates. |
| `ferx_predict()` | Population predictions (ETA = 0). Pass `fit =` to use fitted theta. |
| `ferx_example()` | Get paths to bundled example models and data |

## License

MIT — see [LICENSE.md](LICENSE.md).
