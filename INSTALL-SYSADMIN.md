# Sysadmin Install Guide: Shared Rust + Enzyme toolchain for ferx

This guide is for system administrators setting up a shared Rust + Enzyme autodiff toolchain at `/opt/rust-nightly` so multiple users can install the ferx R package.

If you just want ferx on your own dev machine, see the main [README](README.md) — you can link your personal nightly as the `enzyme` toolchain without this.

## Why this is needed

ferx depends on Rust's experimental autodiff, which uses LLVM's Enzyme plugin. As of 2026:
- The `std::autodiff` macros and `-Z autodiff` flag are in upstream rustc nightly
- **But** the Enzyme shared library (`libEnzyme-<N>.so`) is not yet shipped via rustup — it must be built separately

Full rustc rebuild from the EnzymeAD fork is ~2 hours. This guide does the ~30-minute shortcut: upstream nightly from rustup + build only Enzyme against matching LLVM.

## Prerequisites

- Ubuntu 22.04+ (adapt for other distros)
- At least 5 GB free disk
- Sudo access

## Steps

### 1. Install apt dependencies

```bash
sudo apt install -y cmake ninja-build clang libssl-dev pkg-config \
                    python3 build-essential curl git libzstd-dev
```

### 2. Install rustup + upstream nightly

Don't use snap's rustup — its filesystem confinement breaks on non-standard home directories.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
rustup toolchain install nightly
rustup component add rust-src --toolchain nightly
```

### 3. Find your nightly's LLVM major version

```bash
rustc +nightly --version --verbose | grep LLVM
# e.g. "LLVM version: 22.1.2" — major is 22
```

Use this `<MAJOR>` in the rest of this guide.

### 4. Install matching LLVM

Ubuntu's default repos usually lag behind current nightly LLVM. Use apt.llvm.org:

```bash
wget https://apt.llvm.org/llvm.sh
chmod +x llvm.sh
sudo ./llvm.sh <MAJOR>

# Fix GPG keyring permissions if apt warns:
sudo chmod 644 /etc/apt/trusted.gpg.d/apt.llvm.org.asc
sudo apt update
sudo apt install -y llvm-<MAJOR>-dev clang-<MAJOR>
```

### 5. Build the Enzyme plugin against LLVM

```bash
git clone https://github.com/EnzymeAD/Enzyme /tmp/enzyme-build
cd /tmp/enzyme-build/enzyme
mkdir build && cd build

cmake -G Ninja .. \
  -DLLVM_DIR=/usr/lib/llvm-<MAJOR>/lib/cmake/llvm \
  -DENZYME_CLANG=OFF \
  -DENZYME_FLANG=OFF

ninja
# 15–30 min. Produces Enzyme/LLVMEnzyme-<MAJOR>.so
```

### 6. Stage the nightly toolchain in /opt

```bash
sudo mkdir -p /opt/rust-nightly
sudo cp -a ~/.rustup/toolchains/nightly-x86_64-unknown-linux-gnu/. /opt/rust-nightly/
sudo chown -R root:root /opt/rust-nightly
sudo chmod -R a+rX /opt/rust-nightly
sudo chmod a+rx /opt/rust-nightly/bin/*
```

Use `cp -a` (not `-r`) to preserve symlinks. `chmod a+rX` (capital X) adds execute only on directories and files that are already executable — safe for the whole tree.

### 7. Drop Enzyme into the target-specific lib directory

**Important:** rustc searches for Enzyme in `<sysroot>/lib/rustlib/<target>/lib/`, not in `<sysroot>/lib/`. Despite the error message saying "folder", it's actually looking for a file named `libEnzyme-<MAJOR>.so` in that target-specific location.

```bash
TARGET=x86_64-unknown-linux-gnu
sudo cp /tmp/enzyme-build/enzyme/build/Enzyme/LLVMEnzyme-<MAJOR>.so \
  /opt/rust-nightly/lib/rustlib/$TARGET/lib/libEnzyme-<MAJOR>.so
sudo chmod a+r /opt/rust-nightly/lib/rustlib/$TARGET/lib/libEnzyme-<MAJOR>.so
```

Note the filename rewrite: `LLVMEnzyme-<N>.so` → `libEnzyme-<N>.so` (with `lib` prefix).

### 8. Verify

```bash
/opt/rust-nightly/bin/rustc -Z autodiff=Enable - </dev/null 2>&1 | head
# Should print: "error[E0601]: `main` function not found"
# That's the expected success indicator — rustc + Enzyme loaded fine.
# If you see "autodiff backend not found in the sysroot" instead, step 7 didn't take.
```

### 9. Cleanup (optional)

```bash
rm -rf /tmp/enzyme-build
```

## Per-user instructions (distribute to users)

Each user who wants to install ferx runs these steps once:

```bash
# Install their own rustup (no default toolchain)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain none
source "$HOME/.cargo/env"

# Link the shared toolchain under the name `enzyme`
rustup toolchain link enzyme /opt/rust-nightly

# Verify
rustc +enzyme --version
```

Then in `~/.Renviron`:
```
PATH=/opt/rust-nightly/bin:${HOME}/.cargo/bin:${PATH}
RUSTUP_TOOLCHAIN=enzyme
```

Restart R, then:
```r
devtools::install_github("FeRx-NLME/ferx-r")
```

## Troubleshooting

**"autodiff backend not found in the sysroot: failed to find a `libEnzyme-<N>` folder"**
- The `.so` is in `<sysroot>/lib/` instead of `<sysroot>/lib/rustlib/<target>/lib/` — move it to the right place (step 7)
- LLVM version mismatch — check `rustc --version --verbose | grep LLVM`, rebuild Enzyme against that LLVM
- Filename must be `libEnzyme-<MAJOR>.so` with the `lib` prefix

**"not a directory: '/<path>/lib'"** from `rustup toolchain link`
- Permission issue — non-root user can't read into the target path. Run `sudo chmod -R a+rX /opt/rust-nightly`.

**"error: the option `Z` is only accepted on the nightly compiler"** during R install
- R is resolving a stable rustc. Check `Sys.which("rustc")` in R, ensure `~/.Renviron` `PATH` includes `/opt/rust-nightly/bin` first.

**`"Enzyme: cannot handle (forward) unknown intrinsic llvm.maximumnum"`** during ferx build
- Not a sysadmin issue — this is a ferx-core code limitation. Recent rustc lowers `.max()`/`.min()` to intrinsics Enzyme doesn't differentiate. Report upstream or use an older commit of ferx-core if urgent.

## Refresh after rustc nightly updates

Upstream nightly rolls forward, but `/opt/rust-nightly` is frozen at whatever you installed. To refresh:

1. Install a fresh nightly to your own `~/.rustup`
2. Check if the LLVM major version changed
3. If it did — rebuild Enzyme against the new LLVM (step 5), update the `<MAJOR>` in step 7
4. Re-stage `/opt/rust-nightly` (step 6) and drop the new Enzyme `.so` (step 7)

Plan a refresh roughly quarterly, or whenever a user reports `autodiff_forward` or similar std items aren't found — a sign ferx-core moved ahead of your cached toolchain.
