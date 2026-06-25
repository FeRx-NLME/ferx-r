# ferx R package

[![R-CMD-check](https://github.com/FeRx-NLME/ferx-r/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FeRx-NLME/ferx-r/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/FeRx-NLME/ferx-r/branch/main/graph/badge.svg)](https://codecov.io/gh/FeRx-NLME/ferx-r)
[![CodeFactor](https://www.codefactor.io/repository/github/ferx-nlme/ferx-r/badge)](https://www.codefactor.io/repository/github/ferx-nlme/ferx-r)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Dependencies](https://img.shields.io/badge/dependencies-2-brightgreen.svg)](DESCRIPTION)

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend.

## Features

- **FOCE/FOCEI estimation** 
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **BLOQ handling**: Beal's M3 likelihood for observations below the LLOQ
- **Model DSL**: define models in `.ferx` text files

## Installation

Install a standard stable Rust toolchain with [rustup](https://www.rust-lang.org/tools/install), then install the R package. No custom Rust toolchain is required.

### Install the package

In R, run:

```r
pak::pak("FeRx-NLME/ferx-r")
```

Or from a local clone:
```bash
R CMD INSTALL .
```

First install takes ~1-2 hours (Rust compilation from scratch). Subsequent installs are fast.

All estimation methods (FOCE/FOCEI/SAEM/IMP) work with the standard build.

### Windows

Native Windows installs are supported.

Prerequisites:

- [R](https://cran.r-project.org/bin/windows/base/) and the matching [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for your R version — Rtools45 for R 4.5, Rtools44 for R 4.4. (Rtools ships the MinGW gcc that R uses to link the package.)
- [rustup](https://www.rust-lang.org/tools/install) with the GNU-ABI toolchain:

```powershell
rustup toolchain install stable-x86_64-pc-windows-gnu
```

You do **not** need to `rustup default` it — the package's build pins this toolchain automatically on Windows. (The rustup default on Windows is the MSVC ABI, which is not link-compatible with Rtools' MinGW linker.)

Then install in R:

```r
pak::pak("FeRx-NLME/ferx-r")
```

### Docker

A Docker image is available that bundles ferx CLI, the ferx R package, and RStudio Server.

```bash
# Build
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
