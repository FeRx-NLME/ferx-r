# .fitrx bundles in the tests. A fit's bytes (`fit$bundle`) are the fit (#712): the
# post-hoc steps and ferx_save_fit() read them and not the list, so a test that
# needs a fit in another state edits the bytes, with `fitrx_edit_fit()`. An edit to
# the list reaches neither.

# `fit` with its bytes edited: `edit` takes and returns `fit.json` as lists, and the
# archive entries `drop` go. The list is rebuilt from the new bytes, the way a step
# that returns a new fit builds it, so the paths stay the fit's.
fitrx_edit_fit <- function(fit, edit = identity, drop = character()) {
  path <- tempfile(fileext = ".fitrx")
  writeBin(fit$bundle, path)
  out <- fitrx_edit_bundle(path, edit, drop)
  ferx:::.ferx_fit_with_bundle(fit, readBin(out, "raw", n = file.size(out)))
}

# A bundle's fit.json, as lists, from its path or its bytes.
fitrx_fit_json <- function(path) {
  if (is.raw(path)) {
    bytes <- path
    path <- tempfile(fileext = ".fitrx")
    writeBin(bytes, path)
  }
  staging <- tempfile("fitrx_json_")
  dir.create(staging)
  utils::unzip(path, exdir = staging, junkpaths = TRUE)
  jsonlite::read_json(file.path(staging, "fit.json"), simplifyVector = FALSE)
}

# A copy of the bundle at `path` with fit.json rewritten through `edit`, in
# the writer's own JSON settings, and the archive entries `drop` left out.
fitrx_edit_bundle <- function(path, edit, drop = character()) {
  staging <- tempfile("fitrx_edit_")
  dir.create(staging)
  utils::unzip(path, exdir = staging, junkpaths = TRUE)
  unlink(file.path(staging, drop))
  json <- file.path(staging, "fit.json")
  wire <- edit(jsonlite::read_json(json, simplifyVector = FALSE))
  jsonlite::write_json(wire, json, auto_unbox = TRUE, pretty = TRUE,
                       digits = I(17), null = "null", na = "null")
  out <- tempfile(fileext = ".fitrx")
  utils::zip(out, list.files(staging, full.names = TRUE), flags = "-j -q")
  out
}

# A copy of the bundle at `path` with the entry `name` replaced by `file`.
fitrx_replace_entry <- function(path, name, file) {
  staging <- tempfile("fitrx_edit_")
  dir.create(staging)
  utils::unzip(path, exdir = staging, junkpaths = TRUE)
  file.copy(file, file.path(staging, name), overwrite = TRUE)
  out <- tempfile(fileext = ".fitrx")
  utils::zip(out, list.files(staging, full.names = TRUE), flags = "-j -q")
  out
}

# The columns of `fit$theta_levels` the wire carries.
fitrx_layout_cols <- c("block", "label", "group", "contrast")
