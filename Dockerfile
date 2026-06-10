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
# NOTE: The first build takes ~45-60 min (builds Rust+LLVM+Enzyme from source).
#       Docker caches this layer, so subsequent rebuilds are fast.

# ===========================================================================
# Stage 1: Build the Enzyme-enabled Rust toolchain from source.
# This ensures rustc, LLVM, and Enzyme are all compiled together with
# matching ABIs. The ~20 GB build tree stays in this throwaway stage.
# ===========================================================================
FROM ubuntu:24.04 AS enzyme-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake ninja-build build-essential curl git ca-certificates \
        python3 pkg-config libssl-dev libzstd-dev \
    && rm -rf /var/lib/apt/lists/*

# Bootstrap rustup (needed for `rustup toolchain link` later)
ENV CARGO_HOME=/opt/cargo \
    RUSTUP_HOME=/opt/rustup \
    PATH=/opt/cargo/bin:/usr/local/bin:/usr/bin:/bin

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain nightly --profile minimal

# Clone the Rust source and build a stage-1 compiler with Enzyme.
# download-ci-llvm=false is required: Enzyme cannot be built against CI LLVM.
RUN set -eux; \
    git clone --depth 1 https://github.com/rust-lang/rust /tmp/rust-src; \
    cd /tmp/rust-src; \
    ./configure \
        --release-channel=nightly \
        --enable-llvm-enzyme \
        --enable-llvm-link-shared \
        --enable-ninja \
        --disable-docs \
        --set llvm.download-ci-llvm=false; \
    ./x build --stage 1 library; \
    rustup toolchain link enzyme build/host/stage1; \
    # Verify the toolchain works
    rustc +enzyme --version --verbose

# ===========================================================================
# Stage 2: Final image based on rocker/tidyverse.
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
# 2. rustup + nightly (for cargo). The enzyme toolchain is copied from the
#    builder stage.  CARGO_HOME/RUSTUP_HOME under /opt so the `rstudio`
#    user can use the toolchain too.
# ---------------------------------------------------------------------------
ENV CARGO_HOME=/opt/cargo \
    RUSTUP_HOME=/opt/rustup \
    PATH=/opt/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --default-toolchain nightly --profile minimal \
    && rustup component add rust-src --toolchain nightly \
    && chmod -R a+rX /opt/cargo /opt/rustup

# ---------------------------------------------------------------------------
# 3. Copy the enzyme toolchain from the builder stage and register it.
# ---------------------------------------------------------------------------
COPY --from=enzyme-builder /tmp/rust-src/build/host/stage1 /opt/enzyme-toolchain
RUN rustup toolchain link enzyme /opt/enzyme-toolchain \
    && rustup default enzyme \
    && chmod -R a+rX /opt/enzyme-toolchain \
    && rustc +enzyme --version --verbose

ENV RUSTUP_TOOLCHAIN=enzyme

# ---------------------------------------------------------------------------
# 4. Clone ferx-core from GitHub, build the `ferx` CLI binary, keep the
#    source tree at /opt/ferx-core (the R package's Cargo.toml has a
#    relative path dep). Drop build artifacts and cargo caches.
# ---------------------------------------------------------------------------
RUN set -eux; \
    git clone --depth 1 https://github.com/FeRx-NLME/ferx-core /opt/ferx-core; \
    cd /opt/ferx-core; \
    RUSTFLAGS="-Z autodiff=Enable" cargo build --release; \
    install -m 0755 target/release/ferx /usr/local/bin/ferx; \
    rm -rf target .git; \
    rm -rf /opt/cargo/registry/cache /opt/cargo/registry/src /opt/cargo/git

# ---------------------------------------------------------------------------
# 5. Copy the R package source, install it (which rebuilds the Rust staticlib
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
# 6. Make the toolchain visible to interactive R sessions inside RStudio.
# ---------------------------------------------------------------------------
RUN printf 'PATH=/opt/cargo/bin:${PATH}\nRUSTUP_TOOLCHAIN=enzyme\n' \
        > /home/rstudio/.Renviron \
    && chown rstudio:rstudio /home/rstudio/.Renviron

# RStudio Server entrypoint, port 8787, default user `rstudio` are all
# inherited from rocker/tidyverse — nothing more to do here.
