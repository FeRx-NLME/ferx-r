# ferx_model_section
# ---------------------------------------------------------------------------



test_that("ferx_model_show() prints '# model: <filename>' header", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl(paste0("# model: ", basename(path)), out, fixed = TRUE)))
})
test_that("ferx_model_show() prints all file lines", {
  path <- write_test_model(list(
    parameters  = "  theta TVCL(1.0, 0.001, 100.0)",
    fit_options = c("  method = foce", "  maxiter = 300")
  ))
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl("[parameters]",             out, fixed = TRUE)))
  expect_true(any(grepl("theta TVCL",               out, fixed = TRUE)))
  expect_true(any(grepl("[fit_options]",             out, fixed = TRUE)))
  expect_true(any(grepl("method = foce",             out, fixed = TRUE)))
  expect_true(any(grepl("maxiter = 300",             out, fixed = TRUE)))
})
test_that("ferx_model_show() returns path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_show(path))
  expect_equal(result$value, path)
  expect_false(result$visible)
})
test_that("ferx_model_show() errors on missing file", {
  expect_error(
    ferx_model_show(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})
test_that("ferx_model_show() errors on wrong extension", {
  path <- tempfile(fileext = ".txt")
  writeLines("hello", path)
  on.exit(unlink(path))

  expect_error(ferx_model_show(path), regexp = "\\.ferx")
})
test_that("ferx_model_show() handles an empty .ferx file without error", {
  path <- tempfile(fileext = ".ferx")
  writeLines(character(0), path)
  on.exit(unlink(path))

  out <- capture.output(ferx_model_show(path))
  expect_true(any(grepl("# model:", out, fixed = TRUE)))
})
test_that("ferx_model_show() handles a file containing only comment lines", {
  path <- tempfile(fileext = ".ferx")
  writeLines(c("# this is a comment", "# another comment"), path)
  on.exit(unlink(path))

  expect_no_error(out <- capture.output(ferx_model_show(path)))
  expect_true(any(grepl("# this is a comment", out, fixed = TRUE)))
})
test_that("ferx_model_show() pipes: ferx_model(template=)$model |> ferx_model_show()", {
  path <- tempfile(fileext = ".ferx")
  on.exit(unlink(path))

  result <- ferx_model(template = "1cpt_oral", path = path, edit = FALSE)$model |>
    ferx_model_show()
  expect_equal(result, path)
})

# ---- header from test-model-show-color.R ----
# Console syntax highlighting for ferx_model_show() (issue #4, option 3).



test_that(".ferx_use_color() follows cli's colour detection", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 1)
  expect_false(ferx:::.ferx_use_color())
  withr::local_options(cli.num_colors = 256)
  expect_true(ferx:::.ferx_use_color())
})
test_that(".ferx_highlight_line styles headers, keywords, and comments", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 256)

  hdr <- ferx:::.ferx_highlight_line("[parameters]")
  kw  <- ferx:::.ferx_highlight_line("  theta TVCL(1.0)  # clearance")
  cmt <- ferx:::.ferx_highlight_line("# a comment")

  # Each is actually colourised...
  expect_true(cli::ansi_has_any(hdr))
  expect_true(cli::ansi_has_any(kw))
  expect_true(cli::ansi_has_any(cmt))

  # ...and stripping the ANSI recovers the exact input (no text mangling).
  expect_identical(cli::ansi_strip(hdr), "[parameters]")
  expect_identical(cli::ansi_strip(kw), "  theta TVCL(1.0)  # clearance")
  expect_identical(cli::ansi_strip(cmt), "# a comment")
})
test_that(".ferx_highlight_line leaves non-keyword / blank lines untouched", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 256)
  # No leading keyword and no comment -> returned unchanged (no ANSI).
  expect_false(cli::ansi_has_any(ferx:::.ferx_highlight_line("  CL = TVCL * exp(ETA_CL)")))
  expect_identical(ferx:::.ferx_highlight_line(""), "")
})
test_that("ferx_model_show prints plain, unchanged text when colour is off", {
  skip_if_not_installed("cli")
  withr::local_options(cli.num_colors = 1)
  ex  <- ferx_example("warfarin")
  out <- capture.output(ferx_model_show(ex$model))

  expect_false(any(cli::ansi_has_any(out)))         # no escape codes anywhere
  raw <- readLines(ex$model, warn = FALSE)
  expect_true(all(raw %in% out))                    # every source line present verbatim
})
test_that("ferx_model_show validates its input", {
  expect_error(ferx_model_show(tempfile(fileext = ".ferx")), "File not found")
  f <- tempfile(fileext = ".txt")
  file.create(f)
  on.exit(unlink(f), add = TRUE)
  expect_error(ferx_model_show(f), "must be a .ferx file")
})
test_that("ferx_model_show emits ANSI when colour is on, same content when stripped", {
  skip_if_not_installed("cli")
  ex  <- ferx_example("warfarin")

  plain <- withr::with_options(
    list(cli.num_colors = 1),  capture.output(ferx_model_show(ex$model)))
  col   <- withr::with_options(
    list(cli.num_colors = 256), capture.output(ferx_model_show(ex$model)))

  expect_true(any(cli::ansi_has_any(col)))          # colour actually applied
  expect_identical(cli::ansi_strip(col), plain)     # identical once ANSI removed
})
