# Security Policy

## Supported versions

ferx (the R package) is pre-1.0 and under active development. Security fixes are
applied to the **latest release** (and `main`); older versions are not
maintained. Please update to the latest version before reporting an issue.

## Reporting a vulnerability

**Please do not open a public GitHub Issue for security vulnerabilities.**

Report privately via GitHub's
[**Report a vulnerability**](https://github.com/FeRx-NLME/ferx-r/security/advisories/new)
flow (Security → Advisories → Report a vulnerability). If you can't use that,
contact a maintainer directly.

Please include:

- a description of the issue and its impact,
- steps to reproduce (a minimal model / data subset if relevant),
- the affected version or commit, and
- any suggested fix.

We aim to acknowledge a report within a few business days, agree on a disclosure
timeline, and credit reporters who wish to be named once a fix ships.

## Scope

ferx is a thin R wrapper over the **ferx-core** Rust engine (via the extendr
FFI). Its main security surface is:

- **Native / FFI code** — the Rust glue (`src/rust/`) and the extendr boundary;
  a memory-safety bug there can crash the host R session.
- **Dependencies** — Rust crates monitored by Dependabot.
- **Untrusted input** — parsing `.ferx` model files and NONMEM-format CSV data.

The core engine has its own policy in
[ferx-core's `SECURITY.md`](https://github.com/FeRx-NLME/ferx-core/security/policy);
report engine-level issues there.
