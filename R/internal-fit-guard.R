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
# `r$user`, one node of this file's own JSON encoding per name (see
# `.ferx_user_encode()`), and come back on load.
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

# The assignment's result, `new` (the unclassed list after base R's own
# replacement), as a fit; refused when it changed a core field of the sealed fit
# `x`. Judged by what the replacement did, not by reading the subscript: R
# accepts names, positions, negatives, logicals, factors (by code), Dates and
# other classed numbers (as positions), and a resolver of its own would miss
# one. `identical()` returns at once for a field the replacement left as it was
# (the same object), so the check costs nothing for an untouched fit. An
# assignment that puts back the value a field already has changes nothing, and
# goes through.
.ferx_fit_assigned <- function(x, new) {
  core <- .ferx_fit_core_fields(x)
  if (!is.null(core)) {
    old <- unclass(x)
    have <- names(new) %||% character()
    changed <- core[vapply(core, function(k) {
      !(k %in% have) || !identical(old[[k]], new[[k]])
    }, NA)]
    if (length(changed)) .ferx_fit_read_only_error(changed)
  }
  class(new) <- class(x)
  new
}

.ferx_fit_read_only_error <- function(hit) {
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
  new <- unclass(x)
  new[[name]] <- value
  .ferx_fit_assigned(x, new)
}

#' @rdname ferx_fit-read-only
#' @export
`[[<-.ferx_fit` <- function(x, i, ..., value) {
  new <- unclass(x)
  new[[i, ...]] <- value
  .ferx_fit_assigned(x, new)
}

#' @rdname ferx_fit-read-only
#' @export
`[<-.ferx_fit` <- function(x, i, j, ..., value) {
  new <- unclass(x)
  if (missing(i)) {
    new[] <- value
  } else if (missing(j)) {
    new[i, ...] <- value
  } else {
    new[i, j, ...] <- value
  }
  .ferx_fit_assigned(x, new)
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

# ----- the user-addition encoding ----------------------------------------------
#
# A bundle is untrusted input, so what `ferx_load_fit()` decodes from it must be
# data and nothing else: no `unserialize()`, no namespace loaded, no S4 object
# made, no method dispatched. jsonlite::unserializeJSON() does all of those, so
# user additions use this encoding instead, which only base vector types,
# allowlisted attributes and allowlisted classes get into, and which is decoded
# with primitives alone (`vapply()`, `lapply()`, `attr<-`).
#
# A node is `{"type": t, "values": [...], "attributes": {name: node, ...}}`:
# `t` is one of "NULL", "logical", "integer", "double", "character", "list";
# `values` holds the elements, NA as null, and for "double" the strings "NaN",
# "Inf" and "-Inf"; "list" values are nodes. `attributes` is left out when there
# are none.

.FERX_USER_TYPES <- c("NULL", "logical", "integer", "double", "character", "list")
.FERX_USER_ATTRS <- c("names", "dim", "dimnames", "levels", "row.names", "tzone",
                      "units", "class")
.FERX_USER_CLASSES <- c("data.frame", "factor", "ordered", "Date", "POSIXct",
                        "POSIXt", "difftime", "table", "noquote")
.FERX_USER_MAX_DEPTH <- 64L

.ferx_user_bad <- function(what) stop(what, call. = FALSE)

# `x` as a node (an R list for jsonlite::toJSON(auto_unbox = TRUE)); an error for
# anything outside the encoding.
.ferx_user_encode <- function(x, depth = 0L) {
  if (depth > .FERX_USER_MAX_DEPTH) .ferx_user_bad("nested too deeply")
  if (isS4(x)) .ferx_user_bad("an S4 object")
  type <- typeof(x)
  if (!type %in% .FERX_USER_TYPES) .ferx_user_bad(paste("a value of type", type))
  if (type == "NULL") return(list(type = "NULL"))
  a <- attributes(x)
  bare <- x
  attributes(bare) <- NULL
  values <- switch(type,
    double = lapply(bare, function(v) {
      if (is.nan(v)) "NaN" else if (is.na(v)) NULL else if (v == Inf) "Inf" else
        if (v == -Inf) "-Inf" else v
    }),
    list = lapply(bare, .ferx_user_encode, depth = depth + 1L),
    as.list(bare)
  )
  node <- list(type = type, values = values)
  if (length(a)) {
    extra <- setdiff(names(a), .FERX_USER_ATTRS)
    if (length(extra)) .ferx_user_bad(paste("attribute", extra[[1L]]))
    cls <- setdiff(a$class, .FERX_USER_CLASSES)
    if (length(cls)) .ferx_user_bad(paste("class", cls[[1L]]))
    node$attributes <- lapply(a, .ferx_user_encode, depth = depth + 1L)
  }
  node
}

# The value a node (as jsonlite::parse_json(simplifyVector = FALSE) reads it)
# stands for; an error for anything that is not a node of the encoding.
.ferx_user_decode <- function(node, depth = 0L) {
  if (depth > .FERX_USER_MAX_DEPTH) .ferx_user_bad("nested too deeply")
  if (!is.list(node) || is.null(names(node)) ||
      !all(names(node) %in% c("type", "values", "attributes"))) {
    .ferx_user_bad("not a node")
  }
  type <- node$type
  if (!is.character(type) || length(type) != 1L || !type %in% .FERX_USER_TYPES) {
    .ferx_user_bad("an unknown type")
  }
  if (type == "NULL") {
    if (length(node$values) || length(node$attributes)) .ferx_user_bad("NULL")
    return(NULL)
  }
  vals <- node$values %||% list()
  if (!is.list(vals) || !is.null(names(vals))) .ferx_user_bad("values")
  scalar <- function(v, ok) is.null(v) || (length(v) == 1L && ok(v))
  check <- function(ok) {
    if (!all(vapply(vals, scalar, NA, ok = ok))) .ferx_user_bad(type)
  }
  x <- switch(type,
    logical = {
      check(is.logical)
      vapply(vals, function(v) if (is.null(v)) NA else v, NA)
    },
    integer = {
      check(function(v) is.numeric(v) && is.finite(v) && v == trunc(v) &&
              abs(v) <= .Machine$integer.max)
      vapply(vals, function(v) if (is.null(v)) NA_integer_ else as.integer(v),
             NA_integer_)
    },
    double = {
      check(function(v) is.numeric(v) ||
              (is.character(v) && v %in% c("NaN", "Inf", "-Inf")))
      vapply(vals, function(v) {
        if (is.null(v)) NA_real_ else if (is.character(v))
          c("NaN" = NaN, "Inf" = Inf, "-Inf" = -Inf)[[v]] else as.numeric(v)
      }, NA_real_)
    },
    character = {
      check(is.character)
      vapply(vals, function(v) if (is.null(v)) NA_character_ else v, "")
    },
    list = lapply(vals, .ferx_user_decode, depth = depth + 1L)
  )
  attrs <- node$attributes
  if (length(attrs)) {
    if (!is.list(attrs) || is.null(names(attrs)) ||
        !all(names(attrs) %in% .FERX_USER_ATTRS) || anyDuplicated(names(attrs))) {
      .ferx_user_bad("attributes")
    }
    attrs <- lapply(attrs, .ferx_user_decode, depth = depth + 1L)
    cls <- attrs$class
    if (!is.null(cls) && (!is.character(cls) || !all(cls %in% .FERX_USER_CLASSES))) {
      .ferx_user_bad("class")
    }
    # dim before dimnames, class last: `attr<-` is a primitive and dispatches
    # nothing, so setting a class here runs no code.
    order <- c(intersect(c("dim", "dimnames"), names(attrs)),
               setdiff(names(attrs), c("dim", "dimnames", "class")),
               intersect("class", names(attrs)))
    for (nm in order) attr(x, nm) <- attrs[[nm]]
  }
  x
}

# A value's node, when it comes back through JSON as it went in; NULL otherwise.
.ferx_user_node <- function(value) {
  tryCatch({
    node <- .ferx_user_encode(value)
    json <- jsonlite::toJSON(node, auto_unbox = TRUE, digits = I(17),
                             null = "null", na = "null")
    back <- .ferx_user_decode(jsonlite::parse_json(json, simplifyVector = FALSE))
    if (identical(back, value)) node else NULL
  }, error = function(e) NULL)
}

# The user's additions to `fit` as `name -> node`. One the encoding cannot carry
# exactly (an environment, a function, an S4 object, a class or attribute
# outside the allowlist, ...) is left out with a
# `ferx_fit_user_field_dropped_warning`.
.ferx_fit_user_json <- function(fit, who) {
  user <- .ferx_fit_user_fields(fit)
  out <- list()
  dropped <- character()
  x <- unclass(fit)
  for (k in user) {
    node <- .ferx_user_node(x[[k]])
    if (is.null(node)) dropped <- c(dropped, k) else out[[k]] <- node
  }
  if (length(dropped)) {
    warning(warningCondition(paste0(
      who, ": ", paste0("fit$", dropped, collapse = ", "),
      if (length(dropped) == 1L) " is" else " are",
      " not plain data that survives a round trip through JSON (base vectors, ",
      "lists, data frames, factors, dates and times), so ",
      if (length(dropped) == 1L) "it is" else "they are",
      " not saved with the fit. Other objects you added are saved."),
      class = "ferx_fit_user_field_dropped_warning", fields = dropped,
      call = NULL))
  }
  out
}

# `bundle` with its `interface_extras.json` holding `user` (the nodes from
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
# hold, as a named list of values, decoded by `.ferx_user_decode()`, never by a
# deserializer that could run code. One that is not a valid node is left out
# with a `ferx_fit_user_field_dropped_warning`.
.ferx_fit_user_values <- function(extras) {
  user <- extras$user
  if (!length(user)) return(list())
  out <- list()
  dropped <- character()
  keys <- names(user)
  if (!is.list(user) || is.null(keys)) {
    warning(warningCondition(
      "ferx_load_fit: the bundle's saved user objects are malformed, so none is loaded.",
      class = "ferx_fit_user_field_dropped_warning", fields = character(),
      call = NULL))
    return(list())
  }
  for (i in seq_along(user)) {
    k <- keys[[i]]
    value <- if (nzchar(k) && !k %in% names(out)) {
      tryCatch(list(.ferx_user_decode(user[[i]])), error = function(e) NULL)
    }
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
