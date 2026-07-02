
# write_pipe_test_model()/modifying_editor() come from helper-model-pipe.R

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
