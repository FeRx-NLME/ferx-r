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
- **Compartment-free models**: direct response equations (`$PRED`-equivalent)
- **ODE-based models**: general `[odes]` models, solved with `ode_method = "auto"` (default), which solves a linear system exactly by the matrix exponential and otherwise picks the stepper per integration segment, or with a named stepper -- Dormand-Prince RK45, high-order Verner 7(6), or one of three stiff Rosenbrock methods -- selected with `settings = list(ode_method = ...)`
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

#### Building without Deep Compartment Models

The default build enables `nn` (Deep Compartment Models / low-dim neural ODEs,
the `[covariate_nn NAME]` block). It is the heaviest part of the dependency
graph, and the final link of the Rust static library is where a machine short
on memory tends to give up. Drop it by setting `CARGO_FEATURES` in the
environment:

```bash
CARGO_FEATURES="--no-default-features --features ci,survival" R CMD INSTALL .
```

Keep `ci` and `survival`: they carry the covariance step and time-to-event
support. A build without `nn` rejects a model declaring `[covariate_nn]` with
`E_NN_FEATURE_DISABLED`; everything else is unchanged.

#### Installing without network access

An install from GitHub or a local clone downloads the Rust dependencies
during the build: ferx-core at the commit `src/rust/Cargo.lock` pins, and
the crates.io crates. To build where there is no network (a validation
image, an air-gapped server), make a source tarball that carries them all:

```bash
tools/vendor-rust.sh                 # writes src/rust/vendor.tar.xz (~25 MB)
R CMD build --no-build-vignettes .   # the tarball now includes it
```

`R CMD INSTALL ferx_<version>.tar.gz` then builds offline, with no network
and a read-only `CARGO_HOME`. The vendored crates are the same pinned
revisions an online install would fetch. `vendor.tar.xz` is not committed
(it is gitignored), so a GitHub install is unaffected by any of this.

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

#### Slim image (pipelines, CI, clusters)

A slim image with R, the ferx package and the `ferx` CLI - no RStudio, no Rust
toolchain - is published to the GitHub Container Registry for `linux/amd64`
and `linux/arm64`, so you skip the Rust build:

```bash
docker pull ghcr.io/ferx-nlme/ferx-r:latest
```

Pick the tag that fits:

| Tag | What it is | Moves? |
|---|---|---|
| `latest` | the newest `main` build (development) | with every merge that changes the image |
| `release` | the newest release | when a newer version is released |
| `0.5.0` / `v0.5.0` | exactly that release (from the first release after 0.4.0 on) | never |
| `<sha>` | one `main` build, by ferx-r commit (full or 7-character short) | never |

Release tags are published once and never replaced: the workflow refuses to
overwrite an existing one, and the weekly cleanup never deletes them, `latest`
or `release`. SHA tags are kept for the 30 most recent published commits, so
pipelines that must stay reproducible for longer (validated work) should pin a
release, and record the image digest (`docker buildx imagetools inspect
ghcr.io/ferx-nlme/ferx-r:0.5.0`; pull by `ferx-r@sha256:<digest>`), which
names the exact bytes. Not every `main` commit has an image: merges that only
touch docs, tests or CI publish nothing, and when several PRs merge in quick
succession the commits in between may be skipped.

The image carries ferx's hard dependencies only (no Suggests such as testthat
or xpose), and the default user cannot write to the site library. Add packages
your pipeline needs (targets, crew, ...) in a derived image. Use pak, which
also installs the system libraries a binary package needs (targets pulls in
igraph, which needs `libglpk`):

```dockerfile
FROM ghcr.io/ferx-nlme/ferx-r:<sha>
USER root
RUN R -q -e 'install.packages("pak"); pak::pkg_install(c("targets", "crew"))' \
    && rm -rf /var/lib/apt/lists/*
USER ferx
```

Mount your project at `/work` (the image's working directory) and run R or the
CLI there:

```bash
docker run --rm -v "$PWD:/work" ghcr.io/ferx-nlme/ferx-r:latest Rscript fit.R
docker run --rm -v "$PWD:/work" ghcr.io/ferx-nlme/ferx-r:latest \
  ferx model.ferx --data data.csv
docker run --rm -it -v "$PWD:/work" ghcr.io/ferx-nlme/ferx-r:latest   # R console
```

The container runs as an unprivileged user `ferx` (uid 1000). If your host uid
differs and the outputs must be yours, add `--user "$(id -u):$(id -g)"`.

**Threads.** By default the engine uses the available cores minus one, capped
at 8. Inside a container "available" honours a CPU quota (`docker run
--cpus 4`, a Kubernetes CPU *limit*) and CPU pinning (`--cpuset-cpus`), but not
a Kubernetes CPU *request* without a limit, or a Slurm allocation that is not
enforced through cgroups: there the engine sees every core of the node. When
several workers share a node, set the count explicitly - `ferx_fit(...,
threads = N)` in R, `--threads N` for the CLI - to match what each job was
given.

To build it yourself: `docker build -f Dockerfile.slim -t ferx-r:slim .`

#### RStudio image

A Docker image is available that bundles ferx CLI, the ferx R package, and RStudio Server.

```bash
# Build
docker build -t ferx:latest .

# Run RStudio Server, reachable from this machine only
docker run --rm -p 127.0.0.1:8787:8787 -e PASSWORD='<your password>' ferx:latest
# -> http://localhost:8787   user: rstudio
```

Pick your own password: the image contains none, and RStudio Server sets the
`rstudio` user's password from `PASSWORD` when the container starts. Keep the
`127.0.0.1:` prefix unless you mean to serve RStudio to your network; without
it, Docker publishes the port on every interface, so anyone who can reach the
machine can try to log in.

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

### Compartment-free structural models

Dose-response curves, response time courses, and meta-regression models can
declare their prediction directly in `[structural_model]`, without a `pk`,
`ode`, or `ode_template` line:

```
[structural_model]
  EFF = EMAX * TIME / (ET50 + TIME)
  y   = E0 - EFF
```

Here `EFF` is a named intermediate and `y` is the required prediction. The
equation is evaluated once per observation, so there are no compartments to
dose and nothing to integrate; the dataset does not need an `AMT` column. See
`ferx_example("emax_timecourse")` for a complete synthetic example equivalent
to a NONMEM `$PRED` model.

## API Reference

| Function | Description |
|---|---|
| `ferx_fit()` | Fit a NLME model (FOCE/FOCEI). `bloq_method = "m3"` enables M3. |
| `ferx_simulate()` | Simulate replicates with BSV and residual error. Pass `fit =` to use fitted estimates. |
| `ferx_predict()` | Population predictions (ETA = 0). Pass `fit =` to use fitted theta. |
| `ferx_example()` | Get paths to bundled example models and data |

## License

MIT — see [LICENSE.md](LICENSE.md).
