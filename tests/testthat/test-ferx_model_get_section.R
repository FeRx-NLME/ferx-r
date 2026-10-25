# ferx_model_get_section
# ---------------------------------------------------------------------------



test_that("ferx_model_get_section() returns body lines of a named section", {
  path <- write_test_model(list(
    parameters  = c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- ferx_model_get_section(path, "parameters")
  expect_equal(result, c("  theta TVCL(1.0, 0.001, 100.0)", "  omega ETA_CL ~ 0.09"))
})
test_that("ferx_model_get_section() returns the last section correctly", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = c("  method = focei", "  maxiter = 300")
  ))
  on.exit(unlink(path))

  result <- ferx_model_get_section(path, "fit_options")
  expect_equal(result, c("  method = focei", "  maxiter = 300"))
})
test_that("ferx_model_get_section() returns character(0) for an empty section", {
  path <- write_test_model(list(
    parameters  = character(0),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- suppressMessages(ferx_model_get_section(path, "parameters"))
  expect_equal(result, character(0))
})
test_that("ferx_model_get_section() errors with available names when section is missing", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  expect_error(
    ferx_model_get_section(path, "odes"),
    regexp = "parameters, fit_options"
  )
})
test_that("ferx_model_get_section() errors on non-existent file", {
  expect_error(
    ferx_model_get_section(file.path(tempdir(), "no_such.ferx"), "parameters"),
    regexp = "File not found"
  )
})
test_that("ferx_model_get_section() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_get_section(path, "parameters"), regexp = "\\.ferx")
})
test_that("ferx_model_get_section(strip = FALSE) preserves leading whitespace (default)", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)", "    omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  expect_equal(
    ferx_model_get_section(path, "parameters"),
    c("  theta TVCL(1.0, 0.001, 100.0)", "    omega ETA_CL ~ 0.09")
  )
})
test_that("ferx_model_get_section(strip = TRUE) trims only leading whitespace", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)  ", "    omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  result <- ferx_model_get_section(path, "parameters", strip = TRUE)
  # Leading whitespace removed
  expect_equal(
    result,
    c("theta TVCL(1.0, 0.001, 100.0)  ", "omega ETA_CL ~ 0.09")
  )
  # Internal spacing untouched (between tokens)
  expect_true(grepl("theta TVCL", result[1]))
  expect_true(grepl("omega ETA_CL ~ 0.09", result[2]))
  # Trailing whitespace not trimmed
  expect_true(grepl("  $", result[1]))
})
test_that("ferx_model_get_section() prints '# [section]' header to console", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_get_section(path, "parameters"))
  expect_true(any(grepl("# [parameters]", out, fixed = TRUE)))
})
test_that("ferx_model_get_section() returns body invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(suppressMessages(ferx_model_get_section(path, "parameters")))
  expect_false(result$visible)
})
test_that("ferx_model_get_section() preserves blank lines within section body", {
  path <- write_test_model(list(
    parameters = c("  theta TVCL(1.0, 0.001, 100.0)", "", "  omega ETA_CL ~ 0.09")
  ))
  on.exit(unlink(path))

  result <- ferx_model_get_section(path, "parameters")
  expect_equal(result, c("  theta TVCL(1.0, 0.001, 100.0)", "", "  omega ETA_CL ~ 0.09"))
})
test_that("ferx_model_get_section() parses header with surrounding whitespace", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("  [parameters]  ", "  theta TVCL(1.0, 0.001, 100.0)"), path)
  on.exit(unlink(path))

  result <- ferx_model_get_section(path, "parameters")
  expect_equal(result, "  theta TVCL(1.0, 0.001, 100.0)")
})
test_that("ferx_model_get_section(strip = TRUE) on empty section returns character(0)", {
  path <- write_test_model(list(
    parameters  = character(0),
    fit_options = "  method = foce"
  ))
  on.exit(unlink(path))

  result <- suppressMessages(ferx_model_get_section(path, "parameters", strip = TRUE))
  expect_equal(result, character(0))
})
test_that("ferx_model(template=) created file is readable by ferx_model_get_section()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  ferx_model(template = "1cpt_oral", path = path, edit = FALSE)
  expect_no_error(ferx_model_get_section(path, "parameters"))
})

# ---------------------------------------------------------------------------
# ferx_model object input
# ---------------------------------------------------------------------------

test_that("ferx_model_get_section() on a ferx_model returns the section lines", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  m <- ferx_model(model = path)

  result <- ferx_model_get_section(m, "parameters")
  expect_equal(result, ferx_model_get_section(path, "parameters"))
})
test_that("ferx_model_get_section() on a ferx_model prints the section to console", {
  path <- write_pipe_test_model()
  on.exit(unlink(path))
  out <- capture.output(ferx_model_get_section(ferx_model(model = path), "parameters"))
  expect_true(any(grepl("parameters", out)))
})

# ---------------------------------------------------------------------------
# Instance-named sections (#445)
# ---------------------------------------------------------------------------

test_that("ferx_model_get_section() addresses each instance of [event_model] by name", {
  p <- ferx_example("tte_competing_risks")$model

  a <- ferx_model_get_section(p, "event_model cause_a", strip = TRUE)
  b <- ferx_model_get_section(p, "event_model cause_b", strip = TRUE)
  expect_true("cmt    = 2" %in% a)
  expect_true("cmt    = 3" %in% b)
  # Each instance stops at the next header, not at the next instance's body.
  expect_false(any(grepl("cmt    = 3", a, fixed = TRUE)))
})
test_that("ferx_model_get_section() refuses a bare type with several instances, naming them", {
  p <- ferx_example("tte_competing_risks")$model
  expect_error(
    ferx_model_get_section(p, "event_model"),
    regexp = "ambiguous.*cause_a, cause_b.*event_model cause_a"
  )
})
test_that("ferx_model_get_section() resolves a bare type with a single instance", {
  p <- ferx_example("warfarin_dcm")$model

  by_type <- ferx_model_get_section(p, "covariate_nn")
  by_name <- ferx_model_get_section(p, "covariate_nn TYPICAL_PK")
  expect_identical(by_type, by_name)
  expect_true(length(by_type) > 0L)
  expect_null(attributes(by_type))
  # The printed header is the one resolved to, instance included.
  out <- capture.output(ferx_model_get_section(p, "covariate_nn"))
  expect_true(any(out == "# [covariate_nn TYPICAL_PK]"))
})
test_that("ferx_model_get_section() matches the type case-insensitively, the instance exactly", {
  p <- ferx_example("tte_competing_risks")$model
  expect_identical(
    ferx_model_get_section(p, "EVENT_MODEL cause_a"),
    ferx_model_get_section(p, "event_model cause_a")
  )
  expect_error(ferx_model_get_section(p, "event_model CAUSE_A"), regexp = "not found")
  expect_error(ferx_model_get_section(p, "event_model cause_c"),
               regexp = "event_model cause_a, event_model cause_b")
})
test_that("ferx_section_headers() reads headers with the engine's grammar", {
  hdr <- ferx_section_headers(c(
    "[parameters]",
    "  block_omega (ETA_CL, ETA_V) = [",
    "    0.1,",
    "    [0.01, 0.1]",
    "  ]",
    "[event_model cause_a]   # trailing comment",
    "  states = [central]",
    "[Fit_Options]"
  ))
  expect_equal(hdr$positions, c(1L, 6L, 8L))
  expect_equal(hdr$types, c("parameters", "event_model", "fit_options"))
  expect_equal(hdr$instances, c(NA_character_, "cause_a", NA_character_))
  expect_equal(hdr$names, c("parameters", "event_model cause_a", "fit_options"))
})
test_that("ferx_model_get_section() rejects a non-string section", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))
  expect_error(ferx_model_get_section(path, c("a", "b")), regexp = "single string")
})
test_that("ferx_model_get_section() reaches an unnamed block beside named instances of its type", {
  # The engine accepts `[event_model]` and `[event_model NAME]` together.
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "[event_model]",
    "  cmt = 2",
    "[event_model cause_a]",
    "  cmt = 3",
    "[event_model cause_b]",
    "  cmt = 4"
  ), path)

  expect_equal(ferx_model_get_section(path, "event_model"), "  cmt = 2")
  expect_equal(ferx_model_get_section(path, "event_model cause_a"), "  cmt = 3")
  expect_equal(ferx_model_get_section(path, "event_model cause_b"), "  cmt = 4")

  ferx_model_set_section(path, "event_model", "  cmt = 5")
  expect_equal(ferx_model_get_section(path, "event_model"), "  cmt = 5")
  expect_equal(ferx_model_get_section(path, "event_model cause_a"), "  cmt = 3")
})
test_that("ferx_model_get_section() refuses a duplicated header rather than picking one", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))
  writeLines(c("[event_model]", "  cmt = 2", "[event_model]", "  cmt = 3",
               "[event_model cause_a]", "  cmt = 4", "[event_model cause_a]", "  cmt = 5"),
             path)
  expect_error(ferx_model_get_section(path, "event_model"), regexp = "2 headers")
  expect_error(ferx_model_get_section(path, "event_model cause_a"), regexp = "2 headers")
})
