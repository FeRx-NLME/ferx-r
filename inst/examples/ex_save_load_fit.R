library(ferx)

# Portable .fitrx fit bundle: a cross-language zip of JSON (scalars,
# vectors, matrices) and CSV (sdtab, ebes). Roundtrips through ferx-r
# and ferx-core; can also be inspected with any zip tool.
ex <- ferx_example("warfarin")

# Two ways to produce a .fitrx:
#
# (1) Pass `output` to ferx_fit() and the fit is auto-saved before return.
#     `include_data = TRUE` embeds the input CSV in the bundle so it is
#     fully self-contained (otherwise the bundle references the original
#     CSV by path).
tmp <- tempfile(fileext = ".fitrx")
fit <- ferx_fit(
  ex$model, ex$data,
  output       = tmp,
  include_data = TRUE
)

# (2) Or save explicitly after the fact:
fit2 <- ferx_fit(ex$model, ex$data)
ferx_save_fit(fit2, tempfile(fileext = ".fitrx"))

# Roundtrip: load the bundle back into a ferx_fit. The reloaded object
# carries the same theta / omega / sigma / sdtab as the original, including
# named dimnames on omega.
fit_loaded <- ferx_load_fit(tmp)
isTRUE(all.equal(fit$theta, fit_loaded$theta))
identical(dimnames(fit$omega), dimnames(fit_loaded$omega))

# Use the loaded fit anywhere the original would have worked:
print(fit_loaded)
fit_loaded$estimates
