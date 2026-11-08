# syntax=docker/dockerfile:1.7
#
# ferx: R package + Rust (ferx-core) engine + RStudio Server
#
# Build (from this directory):
#   docker build -t ferx:latest .
#
# Run RStudio Server:
#   docker run --rm -p 8787:8787 -e PASSWORD=ferx ferx:latest
#   -> http://localhost:8787   user: rstudio   password: ferx
#
# Run the ferx CLI directly:
#   docker run --rm -v "$PWD:/work" -w /work ferx:latest ferx model.ferx --data data.csv
#
# ferx-core's gradients are pure stable Rust, so the image builds against the
# stock stable toolchain - no custom rustc/LLVM build is needed.

# ===========================================================================
# Final image based on rocker/tidyverse.
# ===========================================================================
# Pinned to a specific R release tag for reproducible builds (avoid :latest drift).
FROM rocker/tidyverse:4.6.1

# ---------------------------------------------------------------------------
# 1. System build + runtime deps.
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake ninja-build lld libssl-dev pkg-config \
        python3 build-essential curl git ca-certificates \
        libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. rustup + stable toolchain. CARGO_HOME/RUSTUP_HOME under /opt so the
#    `rstudio` user can use the toolchain too.
# ---------------------------------------------------------------------------
ENV CARGO_HOME=/opt/cargo \
    RUSTUP_HOME=/opt/rustup \
    PATH=/opt/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain stable --profile minimal \
    && chmod -R a+rX /opt/cargo /opt/rustup

# ---------------------------------------------------------------------------
# 3. Build the `ferx` CLI at the ferx-core commit src/rust/Cargo.lock pins -
#    the same engine the R package builds below - with the package's feature
#    set (CARGO_FEATURES in src/Makevars). `+stable` because ferx-core's own
#    rust-toolchain.toml names nightly, which rustup would auto-install.
#    The source goes to /tmp/ferx-core and is deleted: a checkout at
#    /opt/ferx-core would sit exactly where src/Makevars looks for a sibling
#    (../../ferx-core from /opt/ferx/src) and be patched into the package build.
#    Same recipe as Dockerfile.slim.
# ---------------------------------------------------------------------------
COPY src/rust/Cargo.lock /tmp/ferx-Cargo.lock
RUN set -eux; \
    rev="$(sed -n 's|^source = "git+https://github.com/FeRx-NLME/ferx-core?branch=main#\([0-9a-f]\{40\}\)"$|\1|p' \
             /tmp/ferx-Cargo.lock | sort -u)"; \
    case "$rev" in \
        ''|*[!0-9a-f]*) \
            echo "Dockerfile: expected exactly one ferx-core revision in src/rust/Cargo.lock, got: '$rev'" >&2; \
            exit 1;; \
    esac; \
    git init -q /tmp/ferx-core; \
    cd /tmp/ferx-core; \
    git fetch -q --depth 1 https://github.com/FeRx-NLME/ferx-core "$rev"; \
    git checkout -q FETCH_HEAD; \
    cargo +stable build --release --locked -p ferx-cli \
        --features ferx-core/ci,ferx-core/nn,ferx-core/survival; \
    install -m 0755 target/release/ferx /usr/local/bin/ferx; \
    cd /; \
    rm -rf /tmp/ferx-core /tmp/ferx-Cargo.lock \
           /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git

# ---------------------------------------------------------------------------
# 4. Copy the R package source and install it. With no sibling checkout, the
#    Rust staticlib builds from the ferx-core commit src/rust/Cargo.lock pins
#    (on stable: src/Makevars writes the toolchain file). Then clean all
#    build/caching state.
# ---------------------------------------------------------------------------
COPY . /opt/ferx
RUN set -eux; \
    R -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes')"; \
    R -e "remotes::install_deps('/opt/ferx', dependencies=TRUE, upgrade='never')"; \
    R CMD INSTALL --no-multiarch /opt/ferx; \
    rm -rf /opt/ferx/src/rust/target \
           /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git \
           /tmp/Rtmp* /tmp/downloaded_packages

# ---------------------------------------------------------------------------
# 5. Make the toolchain visible to interactive R sessions inside RStudio.
# ---------------------------------------------------------------------------
RUN printf 'PATH=/opt/cargo/bin:${PATH}\n' \
        > /home/rstudio/.Renviron \
    && chown rstudio:rstudio /home/rstudio/.Renviron

# RStudio Server entrypoint, port 8787, default user `rstudio` are all
# inherited from rocker/tidyverse - nothing more to do here.
