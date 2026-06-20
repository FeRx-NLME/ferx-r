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
# ferx-core's gradients are pure stable Rust (the Enzyme autodiff path was
# retired in FeRx-NLME/ferx-core#381), so the image builds against the stock
# stable toolchain - no custom rustc/LLVM build is needed.

# ===========================================================================
# Final image based on rocker/tidyverse.
# ===========================================================================
# Pinned to a specific R release tag for reproducible builds (avoid :latest drift).
FROM rocker/tidyverse:4.6.0

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
# 3. Clone ferx-core from GitHub, build the `ferx` CLI binary, keep the
#    source tree at /opt/ferx-core (the R package's Cargo.toml has a
#    relative path dep). Drop build artifacts and cargo caches.
# ---------------------------------------------------------------------------
RUN set -eux; \
    git clone --depth 1 https://github.com/FeRx-NLME/ferx-core /opt/ferx-core; \
    cd /opt/ferx-core; \
    cargo build --release; \
    install -m 0755 target/release/ferx /usr/local/bin/ferx; \
    rm -rf target .git; \
    rm -rf /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git

# ---------------------------------------------------------------------------
# 4. Copy the R package source, install it (which rebuilds the Rust staticlib
#    against the retained /opt/ferx-core via the relative path dep in
#    src/rust/Cargo.toml), then clean all build/caching state.
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
