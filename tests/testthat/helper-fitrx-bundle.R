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

# `fit` with its bytes' `interface_extras.json` edited: `edit` takes and returns this
# interface's fields (`call_settings`, `model_file_settings`, ...) as lists.
fitrx_edit_extras <- function(fit, edit) {
  text <- ferx:::ferx_rust_load_fit(fit$bundle)$interface_extras
  extras <- jsonlite::parse_json(text, simplifyVector = FALSE)
  extras$r <- edit(extras$r)
  json <- jsonlite::toJSON(extras, auto_unbox = TRUE, digits = I(17), null = "null",
                           na = "null")
  ferx:::.ferx_fit_with_bundle(
    fit, ferx:::ferx_rust_bundle_with_extras(fit$bundle, as.character(json))
  )
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

# A copy of the bundle at `path` as ferx-core alone writes it: `predictions.csv`
# narrowed to core's own columns (`write_predictions_csv` in core's
# src/io/fitrx.rs, in that order; no CMT, no TTE or [output] columns), and
# without the R package's `sdtab_endpoint_rows.csv`. The discrete and
# time-to-event rows are then only in core's `predictions_discrete.csv` and
# `predictions_tte.csv` (ferx-core #911).
fitrx_core_only_bundle <- function(path) {
  staging <- tempfile("fitrx_core_")
  dir.create(staging)
  utils::unzip(path, exdir = staging, junkpaths = TRUE)
  unlink(file.path(staging, "sdtab_endpoint_rows.csv"))
  preds <- file.path(staging, "predictions.csv")
  wide <- utils::read.csv(preds, colClasses = "character", check.names = FALSE,
                          na.strings = character())
  core <- c("ID", "TIME", "DV", "PRED", "IPRED", "CWRES", "IWRES", "NPDE", "NPD",
            "EBE_OFV", "N_OBS", "CENS", "OCC")
  utils::write.table(wide[intersect(core, names(wide))], preds, sep = ",",
                     quote = FALSE, row.names = FALSE)
  out <- tempfile(fileext = ".fitrx")
  utils::zip(out, list.files(staging, full.names = TRUE), flags = "-j -q")
  out
}

# The columns of `fit$theta_levels` the wire carries.
fitrx_layout_cols <- c("block", "label", "group", "contrast")

# `fit` without its read-only seal (#713): a plain list a test may edit, to probe
# how a reader handles a field the engine never writes that way. A save refuses
# it (`ferx_fit_edited_error`); a step reads the bytes, not the edit.
unsealed_fit <- function(fit) {
  attr(fit, "ferx_core_digest") <- NULL
  fit
}
