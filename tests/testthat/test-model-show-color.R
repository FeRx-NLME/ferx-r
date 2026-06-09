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
