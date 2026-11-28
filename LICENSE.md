# MIT License

Copyright (c) 2026 the FeRx core developer team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Code linked into the compiled package

The MIT License above covers the ferx R package's own code: the R code, the
Rust glue in `src/rust/src/`, the build scripts, documentation and examples.

The compiled package also statically links code under other licenses, which
the MIT grant does not replace:

- **ferx-core** and **ferx-tools** (MIT), from
  [FeRx-NLME/ferx-core](https://github.com/FeRx-NLME/ferx-core).
- **ferx-search** (LGPL-3.0-or-later), the model-search engine behind
  `ferx_modelsearch()`, `ferx_covsearch()`, `ferx_amd()` and the other search
  tools. Parts of it follow the behaviour of
  [Pharmpy](https://github.com/pharmpy/pharmpy) (LGPL-3.0-or-later).
- **NLopt**, whose Luksan routines (L-BFGS among them) are
  LGPL-2.1-or-later.
- The other Rust crates the build links, each under its own license.

[`inst/COPYRIGHTS`](inst/COPYRIGHTS) summarises the copyright holders and
licenses and says how to rebuild the package with a modified ferx-search;
[`inst/THIRD-PARTY-NOTICES.txt`](inst/THIRD-PARTY-NOTICES.txt) reproduces the
license files of every linked crate, and [`inst/licenses/`](inst/licenses) holds
the LGPL-3.0, GPL-3.0 and LGPL-2.1 texts. All of them are installed with the
package (`system.file("COPYRIGHTS", package = "ferx")`) and shipped in the
Docker images.
