# ----- read-only fits (#713) ---------------------------------------------------
#
# A fit's core fields are the engine's: the list is a view of the fit's bytes
# (`fit$bundle`), and an edit of the list reaches neither the bytes nor the
# engine. So they may not be changed, while new objects may be added beside them.
# R cannot make that airtight, so there are two layers:
#
# 1. The assignment methods below refuse `fit$x <-`, `fit[["x"]] <-` and
#    `fit["x"] <-` on a core field, with a `ferx_fit_read_only_error`. A new name
#    is a user addition and goes through.
# 2. Every fit is sealed when it is built (`.ferx_fit_seal()`): it carries a
#    digest of each core field. `ferx_save_fit()` recomputes them and refuses a
#    fit whose core fields changed anyway (`unclass()`, edit, reclass) with a
#    `ferx_fit_edited_error`.
#
# User additions are saved in the bundle's `interface_extras.json`, under
# `r$user`, one jsonlite::serializeJSON() text per name, and come back on load.
#
# Internal code that builds or rebuilds a fit works on the unclassed list and
# seals it at the end; nothing in the package assigns into a sealed fit.

# The attribute that holds the digest: a named character vector, one SHA-256 per
# core field.
.FERX_FIT_DIGEST_ATTR <- "ferx_core_digest"

# Fields a user may set as their own: they say where the fit's files are now, not
# what the fit is (the bytes hold the model, and data is hash-checked).
.FERX_FIT_FREE_FIELDS <- c("model_path", "data_path")

# One field's digest: SHA-256 of its serialization, version 2 so an ALTREP vector
# serializes as its values whether or not it has been expanded since, without the
# 14-byte header that names the R version which wrote it.
.ferx_fit_field_digest <- function(value) {
  bytes <- serialize(value, NULL, xdr = TRUE, version = 2L)
  ferx_rust_sha256_raw(bytes[-seq_len(14L)])
}

# `x` (a fit's list) as a sealed `ferx_fit`: every field but `user` (the user's
# additions) and the free fields is core, and digested.
.ferx_fit_seal <- function(x, user = character()) {
  x <- unclass(x)
  attr(x, .FERX_FIT_DIGEST_ATTR) <- NULL
  core <- setdiff(names(x), c(user, .FERX_FIT_FREE_FIELDS, ""))
  digest <- vapply(core, function(k) .ferx_fit_field_digest(x[[k]]), "")
  attr(x, .FERX_FIT_DIGEST_ATTR) <- digest
  class(x) <- "ferx_fit"
  x
}

# A fit's core field names, NULL for an unsealed fit (one built by hand, or by
# an older ferx).
.ferx_fit_core_fields <- function(fit) {
  digest <- attr(fit, .FERX_FIT_DIGEST_ATTR, exact = TRUE)
  if (is.null(digest)) NULL else names(digest)
}

# The names a user added to a sealed fit; none for an unsealed one.
.ferx_fit_user_fields <- function(fit) {
  core <- .ferx_fit_core_fields(fit)
  if (is.null(core)) return(character())
  setdiff(names(unclass(fit)), c(core, .FERX_FIT_FREE_FIELDS, ""))
}

# The names an assignment `x[i] <-` or `x[[i]] <-` writes to.
.ferx_fit_target_names <- function(x, i, missing_i) {
  nms <- names(unclass(x)) %||% character()
  if (missing_i) return(nms)
  if (is.character(i)) return(i)
  if (is.logical(i)) {
    return(nms[rep_len(i, length(nms)) %in% TRUE])
  }
  if (is.numeric(i)) {
    i <- i[!is.na(i)]
    if (any(i < 0)) return(nms[i])
    return(nms[i[i >= 1 & i <= length(nms)]])
  }
  character()
}

# Refuse an assignment that writes to a core field of a sealed fit.
.ferx_fit_guard <- function(x, targets) {
  core <- .ferx_fit_core_fields(x)
  if (is.null(core)) return(invisible())
  hit <- unique(intersect(targets, core))
  if (length(hit) == 0L) return(invisible())
  one <- length(hit) == 1L
  stop(errorCondition(paste0(
    "ferx_fit objects are read-only: ", paste0("fit$", hit, collapse = ", "),
    if (one) " is a field" else " are fields", " of the fit, a view of its ",
    ".fitrx bytes (fit$bundle), and cannot be assigned to. A changed fit comes ",
    "from the step that makes one (ferx_covariance(), ferx_sir(), ",
    "ferx_calc_npde()) or from a new ferx_fit(). Objects of your own can be ",
    "added under a new name (fit$my_notes <- ...); they are saved with the fit."),
    class = "ferx_fit_read_only_error", fields = hit, call = NULL))
}

#' Fits are read-only
#'
#' A \code{ferx_fit} is a view of the bytes the engine serialized it to
#' (\code{fit$bundle}), so its fields cannot be changed: assigning to one with
#' \code{$<-}, \code{[[<-} or \code{[<-} is an error of class
#' \code{ferx_fit_read_only_error}. Objects of your own can be added under a new
#' name; \code{\link{ferx_save_fit}} keeps them, and \code{\link{ferx_load_fit}}
#' restores them. \code{fit$model_path} and \code{fit$data_path}, which say
#' where the fit's files are now, may be set.
#'
#' @param x A \code{ferx_fit}.
#' @param name,i The field to set.
#' @param j Unused.
#' @param ... Unused.
#' @param value The value to assign.
#' @return The fit with the new object added.
#' @name ferx_fit-read-only
#' @keywords internal
#' @export
`$<-.ferx_fit` <- function(x, name, value) {
  .ferx_fit_guard(x, name)
  NextMethod()
}

#' @rdname ferx_fit-read-only
#' @export
`[[<-.ferx_fit` <- function(x, i, ..., value) {
  target <- .ferx_fit_target_names(x, if (length(i)) i[[1L]] else i,
                                   missing_i = FALSE)
  .ferx_fit_guard(x, target)
  NextMethod()
}

#' @rdname ferx_fit-read-only
#' @export
`[<-.ferx_fit` <- function(x, i, j, ..., value) {
  target <- .ferx_fit_target_names(x, if (!missing(i)) i, missing_i = missing(i))
  .ferx_fit_guard(x, target)
  NextMethod()
}

# Refuse to save a fit whose core fields are not the ones it was built with: the
# list was edited around the guards (`unclass()`, edit, reclass), or it carries
# no digest at all (one built by hand, by an older ferx, or with the digest
# stripped).
.ferx_fit_check_sealed <- function(fit, who) {
  digest <- attr(fit, .FERX_FIT_DIGEST_ATTR, exact = TRUE)
  if (is.null(digest)) {
    stop(errorCondition(paste0(
      who, ": this fit carries no record of the fields it was built with, so ",
      "it cannot be checked against its bytes: it was built by hand or by an ",
      "older ferx, or the record was removed. Save a fit as ferx_fit() (or ",
      "ferx_load_fit()) returned it, or refit."),
      class = "ferx_fit_edited_error", fields = character(), call = NULL))
  }
  x <- unclass(fit)
  changed <- names(digest)[vapply(names(digest), function(k) {
    !(k %in% names(x)) || !identical(.ferx_fit_field_digest(x[[k]]), digest[[k]])
  }, NA)]
  if (length(changed)) {
    stop(errorCondition(paste0(
      who, ": the fit's ", paste0("`", changed, "`", collapse = ", "),
      if (length(changed) == 1L) " was" else " were",
      " changed after it was built. Fits are read-only: the list is a view of ",
      "the fit's bytes, and a save writes the bytes, so the edit would be lost. ",
      "Save the fit as it was returned, or make a changed fit with the step ",
      "that makes one (ferx_covariance(), ferx_sir(), ferx_calc_npde())."),
      class = "ferx_fit_edited_error", fields = changed, call = NULL))
  }
  invisible()
}

# ----- user additions in the bundle --------------------------------------------

# The user's additions to `fit` as `name -> serializeJSON text`. One that does not
# come back from JSON as it went in (an environment, an external pointer, a
# closure over local state, ...) is left out with a
# `ferx_fit_user_field_dropped_warning`.
.ferx_fit_user_json <- function(fit, who) {
  user <- .ferx_fit_user_fields(fit)
  out <- list()
  dropped <- character()
  x <- unclass(fit)
  for (k in user) {
    value <- x[[k]]
    text <- tryCatch({
      s <- jsonlite::serializeJSON(value, digits = I(17))
      if (identical(jsonlite::unserializeJSON(s), value)) as.character(s) else NULL
    }, error = function(e) NULL)
    if (is.null(text)) dropped <- c(dropped, k) else out[[k]] <- text
  }
  if (length(dropped)) {
    warning(warningCondition(paste0(
      who, ": ", paste0("fit$", dropped, collapse = ", "),
      if (length(dropped) == 1L) " does" else " do",
      " not survive a round trip through JSON, so ",
      if (length(dropped) == 1L) "it is" else "they are",
      " not saved with the fit. Other objects you added are saved."),
      class = "ferx_fit_user_field_dropped_warning", fields = dropped,
      call = NULL))
  }
  out
}

# `bundle` with its `interface_extras.json` holding `user` (from
# `.ferx_fit_user_json()`) under `r$user`, in place of whatever the bundle held
# there. The rest of the file is kept, other interfaces' keys included; a bundle
# with no user additions either way is returned as it is.
.ferx_bundle_with_user_fields <- function(bundle, user) {
  text <- ferx_rust_load_fit(bundle)$interface_extras
  extras <- if (is.null(text)) list() else
    jsonlite::parse_json(text, simplifyVector = FALSE)
  had <- length(extras$r$user) > 0L
  if (!had && length(user) == 0L) return(bundle)
  if (is.null(extras$r)) extras$r <- structure(list(), names = character())
  extras$r["user"] <- list(if (length(user)) user else NULL)
  json <- jsonlite::toJSON(extras, auto_unbox = TRUE, digits = I(17),
                           null = "null", na = "null")
  ferx_rust_bundle_with_extras(bundle, as.character(json))
}

# The user additions a bundle's interface extras (`.fitrx_read_interface_extras()`)
# hold, as a named list of values. One that no longer parses is left out with a
# `ferx_fit_user_field_dropped_warning`.
.ferx_fit_user_values <- function(extras) {
  user <- extras$user
  if (!length(user)) return(list())
  out <- list()
  dropped <- character()
  for (k in names(user)) {
    value <- tryCatch(list(jsonlite::unserializeJSON(user[[k]])),
                      error = function(e) NULL)
    if (is.null(value)) dropped <- c(dropped, k) else out[k] <- value
  }
  if (length(dropped)) {
    warning(warningCondition(paste0(
      "ferx_load_fit: the bundle's saved ", paste0("fit$", dropped, collapse = ", "),
      " could not be read back, so ", if (length(dropped) == 1L) "it is" else
        "they are", " left out."),
      class = "ferx_fit_user_field_dropped_warning", fields = dropped,
      call = NULL))
  }
  out
}

# `x` (a fit's unclassed list) with `user` (a named list) added as the user's
# additions, then sealed. A name the fit already has is the fit's: such an
# addition, saved by a ferx that did not have that field, is left out with a
# warning.
.ferx_fit_seal_with_user <- function(x, user) {
  x <- unclass(x)
  taken <- intersect(names(user), c(names(x), .FERX_FIT_FREE_FIELDS))
  if (length(taken)) {
    warning(warningCondition(paste0(
      "ferx_load_fit: the saved ", paste0("fit$", taken, collapse = ", "),
      " added by the user ", if (length(taken) == 1L) "is" else "are",
      " now a field of the fit, so the fit's own value is kept."),
      class = "ferx_fit_user_field_dropped_warning", fields = taken,
      call = NULL))
    user <- user[setdiff(names(user), taken)]
  }
  for (k in names(user)) x[k] <- list(user[[k]])
  .ferx_fit_seal(x, user = names(user))
}
