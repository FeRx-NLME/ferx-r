# ----- save -----------------------------------------------------------------

#' Save a ferx fit to a portable \code{.fitrx} bundle
#'
#' Writes the fit's bundle, the bytes ferx-core serialized it to, as a single,
#' self-describing zip archive that any language with stdlib JSON + CSV + zip
#' readers can open (R, Python, Julia, Rust). The archive includes parameter
#' estimates, per-subject EBEs, per-observation predictions, and the verbatim
#' \code{.ferx} model source. The file is the same one the ferx-core command
#' line reads and writes: \code{ferx summary run1.fitrx} works on it.
#'
#' A fit's bundle (\code{fit$bundle}) is what is saved, not its list: the
#' engine wrote it when the fit was made, and \code{\link{ferx_covariance}},
#' \code{\link{ferx_sir}} and \code{\link{ferx_calc_npde}} return a new fit with
#' new bytes. A fit with no
#' bundle, such as one from an \code{.rds} written by an older ferx, is refused
#' with an error of class \code{ferx_no_bundle_error}: refit it.
#'
#' Every number is written so that it reads back exactly: every
#' double value in the estimates and in the per-subject and per-observation
#' tables reads back as the same double, and the predictions and simulations
#' a reloaded fit drives are \code{identical()} to the original's.
#'
#' SIR draws kept with \code{sir_keep_samples = TRUE} (\code{fit$sir_resamples})
#' are written to \code{sir.resamples_packed} in \code{fit.json}, so
#' \code{ferx_simulate_with_uncertainty(method = "sir")} works on the reloaded
#' fit. They are the bulk of such a bundle: \code{n_resamples x n_packed x 8}
#' bytes as doubles, and about the same again in the compressed archive (1000
#' draws of warfarin's 7 packed parameters, 56 kB of doubles, add 60 kB to
#' the \code{.fitrx}). Each row is one draw in the engine's
#' packed parameter space, not on the natural scale: how each column is
#' packed depends on the model (a theta with a negative lower bound is not
#' log-packed, for one), and the bundle does not yet name the columns or
#' their packing.
#'
#' The bundle's \code{model.ferx} holds the bytes that were fitted, the ones
#' \code{fit$model_hash} names: editing the model file after
#' \code{\link{ferx_fit}} does not change what is saved.
#'
#' A fit is read-only (see \dQuote{Read-only fits} in \code{\link{ferx_fit}}).
#' Its core fields are compared with the digest taken when it was built, and a
#' fit that does not match, such as one edited through \code{unclass()}, or that
#' carries no digest, is refused with an error of class
#' \code{ferx_fit_modified_error}. What was added under a new name is written to
#' the bundle's \code{interface_extras.json} and comes back with
#' \code{\link{ferx_load_fit}}; an addition that is not plain data (a function,
#' an environment, a fitted \code{lm}) cannot be written as JSON and is left out
#' with a warning of class \code{ferx_fit_addition_dropped_warning}, whose
#' \code{fields} names it.
#'
#' The schema is ferx-core's; see its \code{docs/src/file-formats/fitrx.md} for
#' the full field reference.
#'
#' @param fit A \code{ferx_fit} object returned by \code{\link{ferx_fit}}.
#' @param output Path to write. When the path has no file extension, the
#'   conventional \code{.fitrx} extension is appended automatically (so
#'   \code{"results/run1"} is written as \code{"results/run1.fitrx"});
#'   \code{\link{ferx_load_fit}} accepts either the bare or the extended
#'   path. An explicit extension is honoured as given. The bundle replaces the
#'   file in one step, after anything else the call writes, so a failed save
#'   leaves what was at \code{output}.
#' @param include_data Logical. When \code{TRUE}, the input NONMEM CSV used
#'   to fit the model is embedded verbatim inside the archive (as
#'   \code{data.csv}). Requires that the file is still accessible at the
#'   path captured at fit time (\code{fit$data_path}), unchanged since the
#'   fit: a file whose SHA-256 no longer matches \code{fit$data_hash} is not
#'   the data that was fitted, and \code{ferx_save_fit()} refuses it with an
#'   error of class \code{ferx_data_edited_error} rather than bundle it.
#'   Default \code{FALSE}.
#'   A fit whose \code{data} was a data.frame (\code{fit$data_source} set)
#'   always embeds it: its file lives only in a temporary directory, so a
#'   bundle without it could not be re-read.
#' @param write_model Logical. When \code{TRUE}, also write the fitted model
#'   with its final estimates as initial values (\code{\link{ferx_update_model}})
#'   next to the bundle, as \code{<name>_final.ferx} for an \code{output} of
#'   \code{<name>.fitrx}. An existing file of that name is replaced, unless it
#'   is the fit's own model file (\code{fit$model_path}), which is an error.
#'   Default \code{FALSE}.
#' @return Invisibly returns \code{fit}, so the function can be used inside a
#'   pipe without breaking the chain:
#'   \preformatted{
#'   fit |> ferx_save_fit("run1.fitrx") |> (\(f) f$estimates)()
#'   }
#' @family utilities
#' @export
#' @seealso \code{\link{ferx_load_fit}}, \code{\link{ferx_fit}} (and its
#'   \code{output} argument for save-during-fit).
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' tmp <- tempfile(fileext = ".fitrx")
#' ferx_save_fit(fit, tmp)
#' fit2 <- ferx_load_fit(tmp)
#' identical(fit$theta, fit2$theta)
ferx_save_fit <- function(fit, output, include_data = FALSE,
                          write_model = FALSE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (returned by ferx_fit()).")
  }
  if (!is.character(output) || length(output) != 1L || is.na(output)) {
    stop("`output` must be a single character path.")
  }
  if (!is.logical(include_data) || length(include_data) != 1L || is.na(include_data)) {
    stop("`include_data` must be TRUE or FALSE.")
  }
  if (!is.logical(write_model) || length(write_model) != 1L || is.na(write_model)) {
    stop("`write_model` must be TRUE or FALSE.")
  }
  bundle <- .ferx_fit_bundle(fit, "ferx_save_fit")
  .ferx_fit_check_saveable(fit, "ferx_save_fit")
  bundles_data <- isTRUE(include_data) || !is.null(fit$data_source)
  # What goes beside the bytes is read from them too, so a hand edit of the list
  # cannot make a `_final.ferx` or a data check disagree with the bundle.
  view <- if (write_model || bundles_data) .ferx_fit_canonical(fit, "ferx_save_fit")

  # Default the archive extension to `.fitrx` when the caller passes a bare
  # path with no extension, so `ferx_save_fit(fit, "results/run1")` lands at
  # "results/run1.fitrx" - the documented convention - rather than at a
  # bare "results/run1" or (via Info-ZIP's own habit) "results/run1.zip".
  # `ferx_load_fit()` mirrors this by falling back to `<path>.fitrx`.
  if (!nzchar(tools::file_ext(output))) {
    output <- paste0(output, ".fitrx")
  }

  # The updated model is rendered before anything is written, so a fit the
  # engine cannot seed fails without leaving a bundle behind. Its name never
  # takes the fit's own model file: `run1.fitrx` saved beside `run1.ferx`
  # writes `run1_final.ferx`, and a path that still lands on the model file
  # is refused rather than overwriting the source.
  if (write_model) {
    model_out_path <- paste0(tools::file_path_sans_ext(output), "_final.ferx")
    src <- fit$model_path
    if (is.character(src) && length(src) == 1L && !is.na(src) &&
        file.exists(src) && file.exists(model_out_path) &&
        identical(normalizePath(src), normalizePath(model_out_path))) {
      stop("write_model = TRUE would overwrite the fit's own model file ",
           src, "; save the bundle under another name.")
    }
    updated_model <- .ferx_updated_model_text(view)
  }

  # data.csv: when asked, and always for a fit of a data.frame (#419). A bundled
  # copy of the user's own CSV goes when this save did not ask for it.
  if (bundles_data) {
    data_path <- .ferx_fit_data_path(view)
    if (is.null(data_path) || is.na(data_path) || !file.exists(data_path)) {
      stop(if (is.null(view$data_source)) {
        "include_data = TRUE but fit$data_path is missing or does not exist."
      } else {
        "The fit's data (fit$data_source) could not be written to a file to bundle."
      })
    }
    bundle <- ferx_rust_bundle_with_data(
      bundle, .fitrx_fitted_data_bytes(view, data_path), data_path
    )
  } else {
    bundle <- ferx_rust_bundle_without_data(bundle)
  }
  bundle <- .fitrx_bundle_with_additions(bundle, fit)

  out_abs <- normalizePath(dirname(output), mustWork = FALSE)
  if (!dir.exists(out_abs)) dir.create(out_abs, recursive = TRUE)
  out_abs <- file.path(out_abs, basename(output))

  # Both files are written beside their targets before either moves, so a failed
  # write changes nothing. A target that cannot be replaced is refused first, and
  # the bundle moves last: whatever goes wrong, a failed save leaves what was at
  # `output`.
  targets <- c(if (write_model) model_out_path, out_abs)
  if (any(dir.exists(targets))) {
    stop("ferx_save_fit: ", targets[dir.exists(targets)][[1L]],
         " is a directory, so the save cannot replace it.", call. = FALSE)
  }
  stage <- .ferx_private_dir(dirname(out_abs))
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(targets))
  if (write_model) writeBin(charToRaw(updated_model), staged[[1L]])
  writeBin(ferx_rust_bundle_for_file(bundle), staged[[length(staged)]])
  for (i in seq_along(targets)) {
    # An existing file keeps its permissions.
    if (file.exists(targets[[i]])) {
      Sys.chmod(staged[[i]], file.mode(targets[[i]]), use_umask = FALSE)
    }
    ferx_rust_replace_file(staged[[i]], targets[[i]])
  }

  invisible(fit)
}

# ----- internals ------------------------------------------------------------

# A new directory in `parent` only this user can enter. `tempfile()` picks a name
# without creating anything, so a file staged at such a name in a shared directory
# could meet a symlink another user put there first, and the write would follow it.
# mkdir does not follow one and fails on any name already taken, which is why its
# result decides here.
.ferx_private_dir <- function(parent, path = tempfile("fitrx_", tmpdir = parent)) {
  if (!dir.create(path, showWarnings = FALSE, mode = "0700")) {
    stop("ferx_save_fit: cannot create a staging directory in ", parent, ".",
         call. = FALSE)
  }
  path
}

# A fit's bundle, the bytes every step that reads or replaces them starts from.
# A fit has none when it comes from an `.rds` written by an older ferx, or when
# the engine could not serialize it; nothing falls back to the list.
.ferx_fit_bundle <- function(fit, who) {
  bundle <- fit$bundle
  if (is.raw(bundle) && length(bundle) > 0L) return(bundle)
  stop(errorCondition(paste0(
    who, ": this fit carries no .fitrx bundle: it was made by an older ferx, or ",
    "the engine could not serialize it. Old fits are not supported. Refit with ",
    "ferx_fit()."), class = "ferx_no_bundle_error"))
}

# The bytes to bundle as data.csv. A fit of a data.frame carries the rows it
# fitted in `fit$data_source`; when that text hashes to `fit$data_hash` (or
# the fit recorded no hash) it is bundled without touching the cache file,
# which may since have become unreadable. Otherwise one read of `data_path`,
# when it hashes to `fit$data_hash` or the fit recorded no hash. Anything
# else is refused with a `ferx_data_edited_error` naming the file, before the
# bundle is written: the fit's data is gone, and a data.csv its own
# `data_hash` does not name would make from-fit calls on the reloaded fit
# refuse (or, without a hash check, use the wrong rows).
.fitrx_fitted_data_bytes <- function(fit, data_path) {
  hash <- fit$data_hash
  has_hash <- is.character(hash) && length(hash) == 1L && !is.na(hash) &&
    nzchar(hash)
  if (has_hash) hash <- tolower(hash)
  src <- fit$data_source
  if (is.character(src) && length(src) == 1L && !is.na(src)) {
    src_bytes <- charToRaw(src)
    if (!has_hash || identical(ferx_rust_sha256_raw(src_bytes), hash)) {
      return(src_bytes)
    }
  }
  bytes <- .fitrx_read_bytes(data_path)
  if (!has_hash || identical(ferx_rust_sha256_raw(bytes), hash)) return(bytes)
  stop(errorCondition(sprintf(paste0(
    "ferx_save_fit(): the data file %s has changed since the fit (its ",
    "SHA-256 no longer matches fit$data_hash), so it is not the data that ",
    "was fitted and is not bundled. Restore the file and save again, save ",
    "with include_data = FALSE, or refit."), data_path),
    class = "ferx_data_edited_error", data_path = data_path))
}

# A file's bytes, read to EOF in one connection (no `file.size()` race: a file
# that grows between a size query and the read would be cut short).
.fitrx_read_bytes <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  chunks <- list()
  repeat {
    chunk <- readBin(con, "raw", n = 1048576L)
    if (length(chunk) == 0L) break
    chunks[[length(chunks) + 1L]] <- chunk
  }
  if (length(chunks) == 0L) raw() else do.call(c, chunks)
}

# What an R fit holds and core's schema has no slot for, as the text of the
# bundle's interface_extras.json, which core stores and never reads. Under one
# key per interface, so a bundle another interface wrote reads as having none.
.fitrx_interface_extras <- function(fit) {
  out <- list()
  for (k in c(.FITRX_INTERFACE_KEYS, "model_structure")) {
    if (!is.null(fit[[k]])) out[[k]] <- fit[[k]]
  }
  # Core's conddist.csv keeps neither count.
  cd <- fit$cond_dist
  if (!is.null(cd)) out$cond_dist <- list(nsamp = cd$nsamp, burnin = cd$burnin)
  as.character(jsonlite::toJSON(
    list(r = out), auto_unbox = TRUE, digits = I(17), null = "null", na = "null"
  ))
}
