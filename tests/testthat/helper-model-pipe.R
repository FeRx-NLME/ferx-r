# Shared fixtures for the ferx_model() pipe-family tests, split out of the
# old test-model-pipe.R when R/model.R was split one-function-per-file.
# testthat sources helper-*.R once before any test-*.R file, so a single
# copy here is visible to all of them.

write_pipe_test_model <- function() {
  path <- tempfile(fileext = ".ferx")
  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  path
}

modifying_editor <- function(p, ...) {
  lines <- readLines(p)
  writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
}
