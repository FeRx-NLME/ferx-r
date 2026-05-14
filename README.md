# ferx

Fast nonlinear mixed effects (NLME) modeling in R, powered by a Rust backend with [Enzyme](https://enzyme.mit.edu/) automatic differentiation for exact gradients.

## Features

- **FOCE/FOCEI estimation** with automatic differentiation
- **Analytical PK models**: 1- and 2-compartment (oral/IV)
- **ODE-based models**: Dormand-Prince RK45 solver for general ODEs
- **NONMEM-compatible**: reads standard NONMEM CSV datasets
- **BLOQ handling**: Beal's M3 likelihood for observations below the LLOQ
- **Model DSL**: define models in `.ferx` text files

## Installation

A Rust installation with the Enzyme AutoDifferentiation engine is required for FeRx to compute gradients. Most likely
you will need to build Rust from source, which may take an hour or so.

See [documentation](https://ferx-nlme.github.io/ferx-core/installation.html) for installation instructions. 

### Install the package

After installing Rust, in R run:

```r
devtools::install_github("FeRx-NLME/ferx-r")
```

Or from a local clone:
```bash
R CMD INSTALL .
```

### Docker

A Docker image is available that bundles the Enzyme toolchain (built from source), ferx CLI, the ferx R package, and RStudio Server — no local Rust/Enzyme setup required. If you're on Windows, this is currently the only supported way of running FeRx.

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
