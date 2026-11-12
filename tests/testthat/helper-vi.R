# A small VI fit of a bundled example through a model-file `method = vi`
# (#507 / #508): the example's own `[fit_options]` minus its method and
# covariance step, 30 Adam iterations, a fixed seed. `final_ofv` sets
# `vi_final_ofv`; NULL leaves the engine default (`none`, so `fit$ofv` is NaN).
vi_fit <- function(final_ofv = NULL, example = "warfarin") {
  ex <- ferx_example(example)
  src <- readLines(ex$model)
  at <- grep("^\\[fit_options\\]", src)
  # Keep the example's own options (an IOV model's `iov_column`), not its
  # method or covariance step.
  opts <- src[-seq_len(at)]
  opts <- opts[!grepl("^\\s*(method|covariance)\\s*=", opts)]
  src <- c(
    src[seq_len(at - 1L)],
    "[fit_options]",
    opts,
    "  method     = vi",
    "  vi_iters   = 30",
    "  vi_seed    = 42",
    if (!is.null(final_ofv)) paste0("  vi_final_ofv = ", final_ofv),
    "  covariance = false"
  )
  model <- tempfile(fileext = ".ferx")
  writeLines(src, model)
  suppressWarnings(ferx_fit(model, ex$data, verbose = FALSE))
}
