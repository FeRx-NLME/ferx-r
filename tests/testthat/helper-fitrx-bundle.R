# .fitrx bundles through the engine's own reader and writer (#466), for the
# tests of the native `data_bindings` slot (test-theta-levels.R T17,
# test-covariate-stats.R C16). The oracle is ferx-core's `load_fit` /
# `save_fit` (`ferx_rust_fitrx_engine_resave()`), never ferx_load_fit()
# reading back what ferx_save_fit() wrote: R's two halves agreeing with each
# other is what let #379 survive.

# Save `fit` from R (data bundled, which the engine's writer needs), read it
# with the engine's loader and write it back with the engine's writer. Returns
# the R bundle, the engine-written one, and the level layout the engine read
# (`block`, `label`, `group`, `contrast`; blocks in sorted order).
fitrx_engine_trip <- function(fit) {
  r_path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, r_path, include_data = TRUE)
  core_path <- tempfile(fileext = ".fitrx")
  read <- ferx:::ferx_rust_fitrx_engine_resave(r_path, core_path)
  list(r_path = r_path, core_path = core_path, read = read)
}

# A bundle's fit.json, as lists.
fitrx_fit_json <- function(path) {
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

# The columns of `fit$theta_levels` the wire carries.
fitrx_layout_cols <- c("block", "label", "group", "contrast")
