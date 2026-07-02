
# ---- header from test-model-pipe.R ----
# Tests for the ferx_model S3 pipe object and pipe-friendly wrappers added in #47:
#   ferx_model()       — constructor
#   print.ferx_model() — console summary
#   ferx_set_section() — pipe-friendly section replacement
#   ferx_get_section() — pipe-friendly section display
#   ferx_fit()         — ferx_model dispatch (inline, see #52)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

write_pipe_test_model <- function() {
  path <- tempfile(fileext = ".ferx")
  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  path
}

modifying_editor <- function(p, ...) {
  lines <- readLines(p)
  writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
}

# ---------------------------------------------------------------------------
# Block 1 — ferx_model() constructor
# ---------------------------------------------------------------------------



test_that("ferx_get_section() on ferx_model prints the section to console", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(ferx_get_section(ferx_model(model = path), "parameters"))
  expect_true(any(grepl("parameters", out)))
})
test_that("ferx_get_section() on ferx_model returns the ferx_model invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m      <- ferx_model(model = path)
  result <- withVisible(ferx_get_section(m, "parameters"))
  expect_s3_class(result$value, "ferx_model")
  expect_false(result$visible)
})
test_that("ferx_get_section() on a plain path returns the path invisibly", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  result <- withVisible(ferx_get_section(path, "parameters"))
  expect_equal(result$value, path)
  expect_false(result$visible)
})
