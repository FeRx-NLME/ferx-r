# ferx_model_get_section
# ---------------------------------------------------------------------------



test_that("ferx_model_edit() errors when dest file exists and overwrite = FALSE", {
  ex   <- ferx_example("warfarin")
  dest <- tempfile()
  dir.create(dest)
  on.exit(unlink(dest, recursive = TRUE))

  file.copy(ex$model, file.path(dest, basename(ex$model)))

  expect_error(
    ferx_model_edit(ex$model, dest = dest),
    regexp = "already exists"
  )
})
test_that("ferx_model_edit() rejects malformed save_as before opening the editor", {
  # Use a user-owned (non-package) file so we exercise the save_as validation
  # path without triggering the in-package copy step. The error must fire
  # before utils::file.edit() is reached.
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  expect_error(
    ferx_model_edit(path, save_as = c("a.ferx", "b.ferx")),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
  expect_error(
    ferx_model_edit(path, save_as = 42),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
  expect_error(
    ferx_model_edit(path, save_as = NA_character_),
    regexp = "must be NULL, TRUE, FALSE, or a single character string"
  )
})
test_that("ferx_model_edit() save_as = FALSE collapses to NULL (no copy, returns path)", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_edit(path, save_as = FALSE, .editor = null_editor))
  expect_equal(result$value, path)
  expect_false(result$visible)
})
test_that("ferx_model_edit() save_as = NULL returns path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  on.exit(unlink(path))

  result <- withVisible(ferx_model_edit(path, save_as = NULL, .editor = null_editor))
  expect_equal(result$value, path)
  expect_false(result$visible)
})
test_that("ferx_model_edit() save_as = 'dest' copies file to that path", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  ferx_model_edit(path, save_as = dest, .editor = null_editor)

  expect_true(file.exists(dest))
  expect_equal(readLines(dest), readLines(path))
})
test_that("ferx_model_edit() save_as returns the dest path invisibly", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  result <- withVisible(ferx_model_edit(path, save_as = dest, .editor = null_editor))
  expect_equal(result$value, dest)
  expect_false(result$visible)
})
test_that("ferx_model_edit() save_as + overwrite = FALSE errors when dest exists", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  writeLines("existing", dest)
  on.exit({ unlink(path); unlink(dest) })

  expect_error(
    ferx_model_edit(path, save_as = dest, overwrite = FALSE, .editor = null_editor),
    regexp = "already exists"
  )
})
test_that("ferx_model_edit() save_as + overwrite = TRUE silently overwrites existing dest", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  writeLines("old content", dest)
  on.exit({ unlink(path); unlink(dest) })

  expect_no_error(
    ferx_model_edit(path, save_as = dest, overwrite = TRUE, .editor = null_editor)
  )
  expect_true(file.exists(dest))
})
test_that("ferx_model_edit() copies pkg file to dest before editing", {
  ex   <- ferx_example("warfarin")
  dest <- tempfile()
  dir.create(dest)
  on.exit(unlink(dest, recursive = TRUE))

  expect_message(
    ferx_model_edit(ex$model, dest = dest, .editor = null_editor),
    regexp = "Copied to"
  )
  expect_true(file.exists(file.path(dest, basename(ex$model))))
})
test_that("ferx_model_edit() post-edit copy: .editor can modify the file before save_as", {
  path <- write_test_model(list(parameters = "  theta TVCL(1.0, 0.001, 100.0)"))
  dest <- tempfile(fileext = ".ferx")
  on.exit({ unlink(path); unlink(dest) })

  modifying_editor <- function(p, ...) {
    lines <- readLines(p)
    writeLines(c(lines, "  theta TVV(10.0, 0.1, 1000.0)"), p)
  }

  ferx_model_edit(path, save_as = dest, .editor = modifying_editor)

  dest_lines <- readLines(dest)
  expect_true(any(grepl("TVV", dest_lines)))
})
test_that("ferx_model_edit() errors on missing file", {
  expect_error(
    ferx_model_edit(file.path(tempdir(), "no_such.ferx")),
    regexp = "File not found"
  )
})
