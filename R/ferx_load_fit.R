# ----- load -----------------------------------------------------------------

#' Load a fit from a \code{.fitrx} bundle
#'
#' Reads a portable \code{.fitrx} archive written by \code{\link{ferx_save_fit}}
#' (or by the ferx-core Rust CLI's \code{--output} flag) with ferx-core's own
#' loader, and builds the list with the converter \code{\link{ferx_fit}} uses,
#' so a loaded fit has the fields a fresh one has, plus the R-specific fields
#' the writer kept in the bundle's \code{interface_extras.json}.
#' \code{fit$individual_estimates} needs the data, so it is \code{NULL} for a
#' bundle without \code{data.csv}. Only bundles of format version 2 load: an
#' older one fails with a version error, and the model must be refitted.
#'
#' A bundle carrying SIR draws (\code{sir.resamples_packed}, written when the
#' fit kept them with \code{sir_keep_samples = TRUE}) restores
#' \code{sir_resamples} (flat, one packed-space draw after another),
#' \code{sir_resamples_n} and \code{sir_resamples_dim}, so
#' \code{ferx_simulate_with_uncertainty(method = "sir")} runs on the loaded
#' fit. A bundle without them leaves all three \code{NULL}. The
#' \code{sir_ci_*} matrices keep their row names (the declared parameter
#' names).
#'
#' To call \code{\link{ferx_predict}} on a loaded fit, the embedded model
#' source can be re-parsed through the existing pipeline.
#'
#' A level-block fit's \code{theta_levels} is rebuilt by laying the bundle's
#' level bindings (\code{data_bindings.levels}) out on its model file. When
#' this ferx cannot do that - the model file is missing or was edited since
#' the fit, its engine no longer parses the model, or no longer accepts the
#' bindings - the bundle still loads: \code{theta_levels} is \code{NULL},
#' \code{unavailable_theta_levels} records why (\code{code}, the engine's
#' diagnostic code or \code{NA}; \code{message}; and the bundle's
#' \code{levels}, which \code{\link{ferx_save_fit}} writes back unchanged),
#' and one warning of class \code{ferx_levels_unavailable_warning} says so.
#' Estimates, \code{fit$estimates}, \code{print()} and the runlog work
#' as usual; \code{\link{ferx_predict}}, \code{\link{ferx_simulate}},
#' \code{\link{ferx_calc_npde}}, \code{\link{ferx_sir}} and
#' \code{\link{ferx_covariance}} refuse the fit with that reason (class
#' \code{ferx_levels_unavailable_error}). Level bindings that are themselves
#' malformed, or two copies of them that disagree, still refuse the load.
#'
#' A \code{scoring_settings} or \code{sir.settings} record that names a
#' setting this ferx does not know - a bundle written by a newer ferx-core -
#' is not loaded: \code{fit$scoring_settings} / \code{fit$sir_settings} is
#' \code{NULL}, one warning of class
#' \code{ferx_settings_record_dropped_warning} names the keys (\code{keys})
#' and the record (\code{record}), and the rest of the bundle loads.
#' \code{\link{ferx_covariance}} and \code{\link{ferx_sir}} then run as on a
#' fit that records none. A setting the record lacks (one written before the
#' engine recorded it) is filled with the value that engine ran under.
#'
#' An \code{iov_occasion} this ferx does not recognise (not
#' \code{"column"}, \code{"per_dose"} or \code{{"time_windows": [...]}}) -
#' a newer ferx-core's variant, or a hand-edited key - is not loaded either:
#' \code{fit$iov_occasion} is \code{NULL}, one warning of class
#' \code{ferx_iov_occasion_dropped_warning} carries the wire value
#' (\code{value}), and the rest of the bundle loads.
#' \code{\link{ferx_covariance}} and \code{\link{ferx_sir}} then derive
#' occasions from the model file's rule, as for a bundle saved before
#' ferx-core #1783.
#'
#' @param path Path to a \code{.fitrx} file. If the given path does not exist
#'   and has no extension, a \code{.fitrx} extension is tried, mirroring
#'   \code{\link{ferx_save_fit}}'s default (so a fit saved to
#'   \code{"results/run1"} loads back from the same bare path).
#' @return A \code{ferx_fit}-classed list.
#' @examples
#' ex  <- ferx_example("warfarin")
#' fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
#' tmp <- tempfile(fileext = ".rds")
#' ferx_save_fit(fit, tmp)
#' fit2 <- ferx_load_fit(tmp)
#' identical(fit$theta, fit2$theta)
#' @family utilities
#' @export
#' @seealso \code{\link{ferx_save_fit}}.
ferx_load_fit <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("`path` must be a single character path to a .fitrx file.")
  }
  if (!file.exists(path)) {
    # Mirror ferx_save_fit()'s `.fitrx` default: a caller who saved to a bare
    # "results/run1" (stored as "results/run1.fitrx") can load it back with
    # the same bare path.
    alt <- paste0(path, ".fitrx")
    if (!nzchar(tools::file_ext(path)) && file.exists(alt)) {
      path <- alt
    } else {
      stop("File does not exist: ", path)
    }
  }

  bundle <- readBin(path, "raw", n = file.size(path))
  loaded <- tryCatch(
    ferx_rust_load_fit(bundle),
    error = function(e) {
      stop("Failed to read .fitrx archive '", path, "': ", conditionMessage(e),
           call. = FALSE)
    }
  )
  .fitrx_warn_dropped_records(loaded)
  raw <- loaded$fit
  extras <- .fitrx_read_interface_extras(loaded$interface_extras)

  # The model and the data are staged to files outside any cleanup, and the
  # wire's `model_path` / `data_path` overridden with them, so `ferx_sir(fit)`
  # and `ferx_predict(fit)` work on this machine even if the original paths do
  # not exist. The stored hashes still match: the bundle holds the source bytes.
  model_path <- tempfile("fitrx_model_", fileext = ".ferx")
  writeBin(charToRaw(loaded$model_text), model_path)
  wire_data_path <- .fitrx_unwrap_opt_chr(raw$data_path)
  raw$model_path <- model_path
  raw$data_path <- .fitrx_stage_data(loaded$data_csv, wire_data_path)

  # What needs the model file, which the converter leaves to this caller: a
  # failure here costs these fields, not the fit (#531).
  levels <- .fitrx_theta_levels_from_bindings(raw, loaded$levels, bundle)
  raw$theta_levels <- levels$theta_levels
  raw$cov_labels <- .fitrx_cov_labels(raw, levels$theta_levels)
  raw$model_structure <- raw$model_structure %||% extras$model_structure
  if (!nzchar(raw$model_name %||% "")) {
    raw$model_name <- tools::file_path_sans_ext(basename(path))
  }

  # A fit saved with a trace carries it in the bundle, where the path the
  # engine wrote it to has usually gone.
  trace <- .fitrx_read_trace_text(loaded$trace_csv)

  result <- .ferx_fit_from_raw(
    raw,
    model           = model_path,
    data            = wire_data_path %||% raw$data_path %||% path,
    settings_used   = extras$call_settings %||% list(),
    model_file_opts = extras$model_file_settings %||% list(),
    loaded          = TRUE
  )
  # Never a name that starts with `theta_levels`: `$` matches partially, so
  # with `theta_levels` absent, `fit$theta_levels` would read that record.
  if (!is.null(levels$unavailable)) result["theta_levels"] <- list(NULL)
  result$unavailable_theta_levels <- levels$unavailable
  result$trace <- trace
  result$bundle <- bundle
  for (key in intersect(names(extras), .FITRX_INTERFACE_KEYS)) {
    result[[key]] <- extras[[key]]
  }
  # "Unknown" (NA), not 0, for a bundle that does not record the counts.
  if (!is.null(result$cond_dist)) {
    result$cond_dist$nsamp <- as.integer(extras$cond_dist$nsamp %||% NA_integer_)
    result$cond_dist$burnin <- as.integer(extras$cond_dist$burnin %||% NA_integer_)
  }
  result
}

# The keys of `interface_extras.json` a loaded fit takes as they are.
.FITRX_INTERFACE_KEYS <- c("call_settings", "model_file_settings", "data_name",
                           "gradient", "gradient_used")

# Warn about what `ferx_rust_load_fit()` dropped from a bundle before core read it.
#
# A key this ferx does not know in a settings record - a bundle written by a newer
# ferx-core, whose record gained a setting - drops that record with a
# `ferx_settings_record_dropped_warning` naming it (#535): the fit is then one
# that records none. ferx-core's own reader ignores such a key and keeps the
# record; R does not, since re-scoring under the record without a setting the
# writer used would not repeat the writer's run. An `iov_occasion` of another
# shape drops the field with a `ferx_iov_occasion_dropped_warning` naming the
# value (#524), and the post-hoc steps fall back to the model file's rule.
.fitrx_warn_dropped_records <- function(loaded) {
  records <- list(
    list(kind = "scoring_settings", block = "scoring_settings",
         unknown = loaded$dropped_scoring),
    list(kind = "sir_settings", block = "sir.settings", unknown = loaded$dropped_sir)
  )
  for (r in records) {
    if (length(r$unknown) == 0L) next
    warning(warningCondition(paste0(
      "ferx_load_fit: the bundle's `", r$block, "` record names ",
      if (length(r$unknown) == 1L) "a setting" else "settings",
      " this ferx does not know (", paste0("\"", r$unknown, "\"", collapse = ", "),
      "), so fit$", r$kind, " is not loaded; it was likely written by a newer ",
      "ferx-core. The rest of the fit loads, and ferx_covariance() / ferx_sir() ",
      "run as on a fit that records no ", r$kind, ". Update ferx to keep it."),
      class = "ferx_settings_record_dropped_warning", keys = r$unknown,
      record = r$kind))
  }
  value <- loaded$dropped_iov_occasion
  if (!is.null(value)) {
    warning(warningCondition(paste0(
      "ferx_load_fit: the bundle's `iov_occasion` is ", value, ", not \"column\", ",
      "\"per_dose\" or {\"time_windows\": [...]}, so fit$iov_occasion is not ",
      "loaded; it was likely written by a newer ferx-core. The rest of the fit ",
      "loads, and ferx_covariance() / ferx_sir() derive occasions from the model ",
      "file's rule, as on a fit that records none. Update ferx to keep it."),
      class = "ferx_iov_occasion_dropped_warning",
      value = jsonlite::parse_json(value, simplifyVector = FALSE)))
  }
  invisible()
}

# This interface's part of a bundle's `interface_extras.json`, which core stores
# and never reads: an empty list for a bundle another writer made.
.fitrx_read_interface_extras <- function(text) {
  if (is.null(text)) return(list())
  extras <- jsonlite::parse_json(text, simplifyVector = FALSE)$r %||% list()
  # JSON has one number type, and a whole number parses as an integer; settings
  # are doubles on a fresh fit.
  for (key in c("call_settings", "model_file_settings")) {
    if (length(extras[[key]])) {
      extras[[key]] <- rapply(extras[[key]], as.numeric, classes = "integer",
                              how = "replace")
    }
  }
  extras
}

# `fit$trace` from the text of a bundle's `trace.csv`, NULL when it has none.
.fitrx_read_trace_text <- function(text) {
  if (is.null(text)) return(NULL)
  path <- tempfile("fitrx_trace_", fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  writeBin(charToRaw(text), path)
  .ferx_read_trace_csv(path)
}

# Stage a bundle's `data.csv` (raw bytes, or NULL) and return the path the fit
# reads its data from. Data fitted from a data.frame (#419) was read from a file
# named for its md5 and exists nowhere but in the bundle, so it goes back to a
# file of that kind, which is what restores `fit$data_source` and makes the next
# ferx_save_fit() bundle it again. A bundled copy of the user's own CSV stays an
# ordinary `include_data` copy, so re-saving without `include_data` leaves it
# out. With nothing bundled the recorded path stands: if it resolves here,
# `ferx_sir()` uses it, and the hash check still applies.
.fitrx_stage_data <- function(data_csv, wire_path) {
  if (is.null(data_csv)) {
    if (is.null(wire_path) || !nzchar(wire_path)) return(NA_character_)
    return(wire_path)
  }
  staged <- tempfile("fitrx_data_", fileext = ".csv")
  writeBin(data_csv, staged)
  frame_md5 <- .ferx_data_cache_md5(wire_path)
  if (!is.null(frame_md5) && identical(unname(tools::md5sum(staged)), frame_md5)) {
    text <- .ferx_read_data_text(staged)
    unlink(staged)
    return(.ferx_write_data_text(text))
  }
  staged
}

# `theta_levels` for a loaded fit, as the glue's list of columns, laid out on the
# bundle's model file staged at `raw$model_path` (#466), and checked against the
# fit's `model_hash` first: the layout places theta by position, so a file that
# is not the fitted one would give every level a neighbour's value. `levels` is
# the bindings' model-free columns.
#
# Returns `list(theta_levels =, unavailable =)`; `unavailable` is NULL when the
# layout was made, else a list of the engine's `code` (NA when it gave none), its
# `message`, and the bundle's `levels` wire, which `ferx_save_fit()` writes back
# as it found it, so a re-saved bundle loses nothing.
.fitrx_theta_levels_from_bindings <- function(raw, levels, bundle) {
  if (is.null(levels) || nrow(levels) == 0L) {
    return(list(theta_levels = NULL, unavailable = NULL))
  }
  frame <- .ferx_theta_levels_frame(
    levels$block, levels$index, levels$label, levels$group, levels$contrast,
    rep(NA_character_, nrow(levels))
  )
  cs <- raw$covariate_stats
  args <- .ferx_fit_binding_args(list(
    theta_levels = frame,
    covariate_stats = .ferx_covariate_stats_frame(
      cs$covariate, cs$median, cs$mean, cs$min, cs$max, cs$mode, cs$levels
    ),
    profiled_blocks = raw$profiled_blocks,
    profiled_levels_at_bound = raw$profiled_levels_at_bound
  ))
  has_cov <- (raw$cov_matrix_dim %||% 0L) > 0L
  # The layout call, returning the frame or the condition it failed with. The
  # engine's record of a coded refusal is read at once, before anything else
  # calls the engine: an edited model file is coded E_FIT_INPUT_CHANGED
  # (ferx-core #2027), and the retry below would replace that record.
  lay_out <- function(hash) {
    tryCatch(
      ferx_rust_theta_levels_from_fit(
        raw$model_path, hash, as.numeric(raw$theta),
        as.character(raw$theta_names), args$fit_bindings,
        cov_matrix_flat = if (has_cov) as.numeric(raw$cov_matrix) else numeric(0),
        cov_matrix_dim = if (has_cov) raw$cov_matrix_dim else 0L,
        sir_resamples_flat = as.numeric(raw$sir_resamples %||% numeric(0)),
        sir_resamples_n = as.integer(raw$sir_resamples_n %||% 0L),
        sir_resamples_dim = as.integer(raw$sir_resamples_dim %||% 0L)
      ),
      error = function(e) .ferx_engine_coded_error(e, conditionMessage(e)) %||% e
    )
  }
  hash <- .ferx_hash_arg(.fitrx_unwrap_opt_chr(raw$model_hash) %||% NA_character_)
  tl <- lay_out(hash)
  if (!inherits(tl, "error")) return(list(theta_levels = tl, unavailable = NULL))
  # Whether the hash is what refused it: the same call without the check (an
  # empty hash skips it) fails the same way only when it is not. Then - and only
  # then - is the model's own validation the cause, and its code the failure's.
  # Read from the engine's behaviour, not from its wording.
  hash_refused <- nzchar(hash) && {
    again <- lay_out("")
    !inherits(again, "error") ||
      !identical(conditionMessage(again), conditionMessage(tl))
  }
  wire <- ferx_rust_fitrx_wire_keys(bundle, "data_bindings")
  .fitrx_unavailable_theta_levels(
    jsonlite::parse_json(wire, simplifyVector = FALSE)$levels, tl,
    if (hash_refused) NULL else raw$model_path
  )
}

# The warn tier of `.fitrx_theta_levels_from_bindings()` (#531): warn once,
# naming the engine's text and code, and return the record the loaded fit
# carries as `unavailable_theta_levels`. The code is the one validation ties
# to the failure (`.ferx_engine_error()`), so a model this engine no longer
# parses gets its parse code (`E_UNKNOWN_BLOCK`, `E_DEPRECATED_BLOCK`, ...).
# `model_path` is NULL when there is no model to validate, or when the hash is
# what refused it (the caller tells): an edited file's own errors are not why
# the layout was refused. An edited file arrives already coded
# E_FIT_INPUT_CHANGED by the engine (ferx-core #2027), and a coded refusal is
# not re-validated. A missing model file is not coded, nor is a binding the
# engine refuses, which the glue raises without one; those carry NA, and the
# text says "no diagnostic code".
.fitrx_unavailable_theta_levels <- function(levels, e, model_path) {
  msg <- conditionMessage(e)
  if (!is.null(model_path) && !inherits(e, "ferx_engine_error")) {
    e <- .ferx_engine_error(e, model_path, NULL)
    msg <- conditionMessage(e)
  }
  code <- if (inherits(e, "ferx_engine_error")) e$code else NA_character_
  why <- list(code = code, message = msg, levels = levels)
  warning(warningCondition(
    paste0(
      "ferx_load_fit: the bundle's theta level bindings (data_bindings$levels) ",
      "could not be laid out on its model file, so fit$theta_levels is ",
      "unavailable. The estimates load; ferx_predict(), ferx_simulate(), ",
      "ferx_calc_npde(), ferx_sir() and ferx_covariance() refuse this fit. ",
      .ferx_levels_unavailable_reason(why)
    ),
    class = "ferx_levels_unavailable_warning", code = code
  ))
  list(theta_levels = NULL, unavailable = why)
}

# The engine's method table (#508): each method's `.fitrx` token and the label
# a live fit reports (ferx-core's `EstimationMethod::label()`), read once from
# the glue, whose exhaustive match fails to compile when the engine gains a
# method. Cached in `.ferx_state`.
.fitrx_method_table <- function() {
  tab <- .ferx_state$method_table
  if (is.null(tab)) {
    tab <- ferx_rust_method_labels()
    .ferx_state$method_table <- tab
  }
  tab
}

# The engine's labels for `raw$cov_matrix`, read from the bundle's model file
# laid out with the fit's bindings and checked against its `model_hash`, as
# ferx_fit() gets them. None - an unnamed matrix - when the engine cannot read
# the model or lay it out (an edited file, a level layout the load could not
# make), or the layout does not have this matrix's size.
.fitrx_cov_labels <- function(raw, theta_levels) {
  d <- raw$cov_matrix_dim %||% 0L
  if (d == 0L) return(character())
  tl <- theta_levels
  cs <- raw$covariate_stats
  tryCatch(
    ferx_rust_cov_labels(
      raw$model_path,
      .ferx_hash_arg(.fitrx_unwrap_opt_chr(raw$model_hash) %||% NA_character_),
      .ferx_fit_binding_args(list(
        theta_levels = .ferx_theta_levels_frame(
          tl$block, tl$index, tl$label, tl$group, tl$contrast, tl$theta_name,
          tl$value, tl$se, tl$sir_lo, tl$sir_hi, tl$profiled
        ),
        covariate_stats = .ferx_covariate_stats_frame(
          cs$covariate, cs$median, cs$mean, cs$min, cs$max, cs$mode, cs$levels
        ),
        profiled_blocks = raw$profiled_blocks,
        profiled_levels_at_bound = raw$profiled_levels_at_bound
      ))$fit_bindings,
      d
    ),
    error = function(e) character()
  )
}
