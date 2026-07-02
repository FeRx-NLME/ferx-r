# Shared fixtures for the ferx_model()-family tests (test-ferx_model*.R),
# split out of the old test-model.R when R/model.R was split one-function-
# per-file. testthat sources helper-*.R once before any test-*.R file, so a
# single copy here is visible to all of them.

# Write a minimal .ferx file from a named list of section -> lines.
write_test_model <- function(sections) {
  path  <- tempfile(fileext = ".ferx")
  lines <- character(0)
  for (nm in names(sections)) {
    lines <- c(lines, paste0("[", nm, "]"), sections[[nm]])
  }
  writeLines(lines, path)
  path
}

# Build a minimal ferx_fit stub for ferx_model_inspect() dispatch tests.
make_ferx_fit_stub <- function(model_structure = NULL, model_name = "test_model") {
  obj <- list(model_structure = model_structure, model_name = model_name)
  class(obj) <- "ferx_fit"
  obj
}

null_editor <- function(...) invisible(NULL)
