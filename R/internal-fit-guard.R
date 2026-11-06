# Read-only fits (#713): what a user may not change on a `ferx_fit`, and how
# that is enforced. R cannot make a list read-only, so there are two layers:
#
#   * `$<-`, `[[<-` and `[<-` refuse a core field, which stops the honest mistake
#     at the line that makes it.
#   * `ferx_save_fit()` compares the core fields with a digest taken when the fit
#     was built, which catches what those methods never see (`unclass()`, `attr<-`,
#     `names<-`), so an edited fit cannot leave the session as a result.
#
# The record is the `ferx_seal` attribute: `fields`, the names that were core
# when the fit was built, and `digest`. A list with the class and no record (one
# built by hand, or read from an `.rds` older than this) has no list of what was
# fitted, so every name it holds counts as core, and `ferx_save_fit()` refuses it.
#
# Code in this package that builds or refreshes a fit works on the plain list:
# `.ferx_fit_unseal()`, then assign, then `.ferx_fit_seal()`.

# The fields a fit can lack. A fit records the names it has when it is built,
# which covers everything else; these are the ones a build can leave out (no
# covariance step, no IOV, no SIR, ...), and they are core all the same, or
# `fit$se_theta <- 0.1` on a fit without a covariance step would pass as a new
# name. A new field the builders can leave out belongs here.
.FERX_FIT_OPTIONAL_FIELDS <- c(
  "bundle", "call_settings", "cond_dist", "condition_number", "cor_matrix",
  "cov_fixed", "cov_matrix", "covariance_method", "covariance_status",
  "covariance_wall_time_secs", "covariate_stats", "data_hash", "data_name",
  "data_path", "data_source", "dw_statistic", "ebe_kappas",
  "ebe_kappas_by_level", "eigenvalues", "estimates", "eta_cov",
  "eta_linked_theta", "eta_log_transformed", "eta_normality",
  "eta_param_types", "gradient", "gradient_used", "impmap_trace",
  "individual_estimates", "iov_occasion", "kappa_is_diagonal", "kappa_levels",
  "kappa_names", "kappa_param_types", "kappa_weight_typical", "kappa_weights",
  "left_init", "max_abs_correlation", "model_file_settings", "model_hash",
  "model_name", "model_path", "model_source", "model_structure", "model_text",
  "n_threads_used", "omega", "omega_iov", "omega_iov_param_corr",
  "omega_is_diagonal", "omega_param_corr", "population_fingerprint",
  "profiled_blocks", "profiled_levels_at_bound", "reader_settings",
  "residual_correlations", "scoring_settings", "se_kappa", "se_omega",
  "se_sigma", "se_theta", "shrinkage_eps", "shrinkage_kappa",
  "shrinkage_kappa_by_occ", "sigma_types", "sir_ci_kappa", "sir_ci_omega",
  "sir_ci_sigma", "sir_ci_theta", "sir_ess", "sir_resamples",
  "sir_resamples_dim", "sir_resamples_n", "sir_seed_used", "sir_settings",
  "stalled_at_init", "theta_levels", "theta_transforms", "trace",
  "trace_path", "unavailable_theta_levels", "vi", "warning_details",
  "warnings", "warnings_structured"
)

# Names a finished fit never holds and the package reads all the same: the
# engine's wire names, which `.ferx_fit_from_raw()` turns into the fields above
# and drops. The writer takes `fit$cov_condition_number` before
# `fit$condition_number`, for one, so an addition under such a name would reach
# the bundle in place of the fitted value. Anything the package reads off a fit
# belongs here or in the list above.
.FERX_FIT_WIRE_NAMES <- c(
  "bundle_error", "cov_condition_number", "cov_eigenvalues", "cov_labels",
  "cov_matrix_dim", "omega_dim", "omega_iov_dim", "residual_correlation_fixed",
  "residual_correlation_i", "residual_correlation_j",
  "residual_correlation_names", "residual_correlation_rho",
  "se_residual_correlation", "theta_names"
)

# `typeof()`s of what a user addition may hold, at any depth and in its
# attributes: data, which JSON carries exactly. This is the writer's test and the
# reader's: `jsonlite::unserializeJSON()` also rebuilds functions, namespaces, S4
# objects and `unserialize()`s builtins, none of which a `.fitrx` from someone
# else may reach.
.FERX_FIT_ADDITION_TYPES <- c("NULL", "logical", "integer", "double", "complex",
                              "character", "raw", "list")

# ----- the record -------------------------------------------------------------

# The record, or NULL for a fit with none. `attr<-` can leave anything there, and
# what is not a record counts as none: refused at save, like any fit without one.
.ferx_fit_seal_of <- function(x) {
  seal <- attr(x, "ferx_seal", exact = TRUE)
  if (is.list(seal)) seal
}

# SHA-256 of the named fields of `x`. Serialization version 2 writes ALTREP
# vectors out in full, so the digest follows the values and not how R holds them,
# and the 14 bytes dropped are its header, which names the running R version: a
# fit read from an `.rds` under a newer R is the same fit.
.ferx_fit_digest <- function(x, fields) {
  bytes <- serialize(unclass(x)[fields], NULL, version = 2L)
  ferx_rust_sha256_raw(bytes[-seq_len(14L)])
}

# What no addition may be called, whether or not this fit has it.
.ferx_fit_reserved_names <- function() {
  c(.FERX_FIT_OPTIONAL_FIELDS, .FERX_FIT_WIRE_NAMES)
}

# The names added to a sealed fit since it was built.
.ferx_fit_user_fields <- function(x) {
  seal <- .ferx_fit_seal_of(x)
  if (is.null(seal)) return(character())
  setdiff(names(x), c(seal$fields, .ferx_fit_reserved_names()))
}

# The list `x` (or a fit) as a sealed `ferx_fit`. On a fit being refreshed, what
# the user added stays theirs and outside the digest; everything else is core.
.ferx_fit_seal <- function(x) {
  user <- .ferx_fit_user_fields(x)
  x <- unclass(x)
  fields <- setdiff(names(x), user)
  attr(x, "ferx_seal") <- list(fields = fields,
                               digest = .ferx_fit_digest(x, fields))
  class(x) <- "ferx_fit"
  x
}

# The plain list of a fit, for the code that refreshes one. Hand it back with
# `.ferx_fit_reseal()`.
.ferx_fit_unseal <- function(fit) unclass(fit)

# `.ferx_fit_seal()` for a fit that came through `.ferx_fit_unseal()`. A fit with
# no record stays without one: refreshing it must not make it saveable.
.ferx_fit_reseal <- function(x) {
  if (!is.null(.ferx_fit_seal_of(x))) return(.ferx_fit_seal(x))
  class(x) <- "ferx_fit"
  x
}

# Do the core fields still match the record? FALSE for a fit with none. The
# digest covers the names the fit was built with, so a name that arrived since
# must be one the assignment methods would have let in: `unclass()` can otherwise
# add `se_theta` to a fit with no covariance step, and the writer would read it.
.ferx_fit_intact <- function(fit) {
  seal <- .ferx_fit_seal_of(fit)
  is.character(seal$digest) && length(seal$digest) == 1L &&
    is.character(seal$fields) &&
    identical(.ferx_fit_digest(fit, seal$fields), seal$digest) &&
    all(vapply(setdiff(names(fit), seal$fields), function(name) {
      is.null(.ferx_fit_addition_refusal(fit, name))
    }, logical(1L)))
}

.ferx_fit_modified_stop <- function(who) {
  stop(errorCondition(paste0(
    who, ": this fit's core fields no longer match the ones it was built with, ",
    "so it is not the result of a fit and is refused. A ferx fit is read-only: ",
    "refit, or keep your own values under a new name (fit$my_value <- ...)."),
    class = "ferx_fit_modified_error", call = NULL))
}

# For a step that builds on a fit: an edited one is refused, since the step
# would otherwise seal the edit in. A fit with no record has nothing to compare,
# and stays unsaveable (`.ferx_fit_reseal()`).
.ferx_fit_check_unmodified <- function(fit, who) {
  if (!is.null(.ferx_fit_seal_of(fit)) && !.ferx_fit_intact(fit)) {
    .ferx_fit_modified_stop(who)
  }
  invisible(fit)
}

# For `ferx_save_fit()`, where no record is a refusal too.
.ferx_fit_check_saveable <- function(fit, who) {
  if (is.null(.ferx_fit_seal_of(fit))) {
    stop(errorCondition(paste0(
      who, ": this fit carries no digest of its core fields: it was made by an ",
      "older ferx, built by hand, or had the digest removed. Refit with ",
      "ferx_fit()."), class = "ferx_fit_modified_error", call = NULL))
  }
  if (!.ferx_fit_intact(fit)) .ferx_fit_modified_stop(who)
  invisible(fit)
}

# ----- assignment -------------------------------------------------------------

.ferx_fit_readonly_stop <- function(message, field) {
  stop(errorCondition(message, class = "ferx_fit_readonly_error", field = field,
                      call = NULL))
}

# Why `name` cannot be added to the fit `x`, or NULL if it can. `$` matches
# partially, so with `se_theta` absent, `se_theta_note` would be what
# `fit$se_theta` reads: in `print()` and everywhere else in this package.
.ferx_fit_addition_refusal <- function(x, name) {
  seal <- .ferx_fit_seal_of(x)
  core <- union(seal$fields %||% names(x), .ferx_fit_reserved_names())
  if (name %in% core) {
    return(paste0("fit$", name, " is read-only: the core fields of a ferx fit ",
                  "cannot be changed. Refit, or keep your own value under a new ",
                  "name (fit$my_", name, " <- ...)."))
  }
  shadowed <- setdiff(core, names(x))
  shadowed <- shadowed[startsWith(name, shadowed)]
  if (length(shadowed) > 0L) {
    return(paste0("fit$", name, " cannot be added: `$` matches partially, so ",
                  "fit$", shadowed[[1L]], ", which this fit does not have, would ",
                  "read it. Choose a name that does not start with `",
                  shadowed[[1L]], "`."))
  }
  NULL
}

# Refuse the assignment to `names` (NA for a position past the end, or a name
# that is not one).
.ferx_fit_guard <- function(x, names) {
  if (anyNA(names) || !all(nzchar(names))) {
    .ferx_fit_readonly_stop(paste0(
      "A new field on a ferx fit needs a name: use fit$my_value <- ... or ",
      "fit[[\"my_value\"]] <- ..."), NA_character_)
  }
  for (name in names) {
    why <- .ferx_fit_addition_refusal(x, name)
    if (!is.null(why)) .ferx_fit_readonly_stop(why, name)
  }
  invisible()
}

# The field names an index selects, NA for a position past the end. Positions
# go through `seq_along()`, which reads logical and negative indices as `[` does.
.ferx_fit_index_names <- function(x, i) {
  if (is.character(i)) return(i)
  names(x)[seq_along(x)[i]]
}

#' Read-only fit fields
#'
#' Assigning to a core field of a \code{ferx_fit} is an error of class
#' \code{ferx_fit_readonly_error}; a new name is allowed. See the
#' \dQuote{Read-only fits} section of \code{\link{ferx_fit}}.
#'
#' @param x A \code{ferx_fit}.
#' @param name,i The field.
#' @param value The value.
#' @param ... Passed on to the list method.
#' @return \code{x} with the new field.
#' @name ferx_fit-assign
#' @keywords internal
#' @export
`$<-.ferx_fit` <- function(x, name, value) {
  .ferx_fit_guard(x, name)
  NextMethod()
}

#' @rdname ferx_fit-assign
#' @export
`[[<-.ferx_fit` <- function(x, i, ..., value) {
  # `fit[[c("theta", "TVCL")]] <- 1` reaches into the first of them.
  .ferx_fit_guard(x, .ferx_fit_index_names(x, i[[1L]]))
  NextMethod()
}

#' @rdname ferx_fit-assign
#' @export
`[<-.ferx_fit` <- function(x, i, ..., value) {
  .ferx_fit_guard(x, if (missing(i)) names(x) else .ferx_fit_index_names(x, i))
  NextMethod()
}

# ----- user additions in a bundle ---------------------------------------------

# Is `x` data, all the way down, attributes included?
.ferx_fit_is_data <- function(x) {
  if (!typeof(x) %in% .FERX_FIT_ADDITION_TYPES || isS4(x)) return(FALSE)
  if (is.list(x) && !all(vapply(unclass(x), .ferx_fit_is_data, logical(1L)))) {
    return(FALSE)
  }
  all(vapply(attributes(x), .ferx_fit_is_data, logical(1L)))
}

.ferx_fit_addition_warning <- function(message, fields) {
  warning(warningCondition(message, class = "ferx_fit_addition_dropped_warning",
                           fields = fields, call = NULL))
}

# A fit's user additions for `interface_extras.json`: name -> the text of
# `jsonlite::serializeJSON()`, which keeps the type and attributes plain JSON
# would lose and, at 17 digits, every double. Text and not a nested object,
# because the reader parses the file once with a parser that does not keep them.
# What is not data is left out, with one warning naming it.
.fitrx_user_additions_wire <- function(fit) {
  fields <- .ferx_fit_user_fields(fit)
  ok <- vapply(fields, function(f) .ferx_fit_is_data(fit[[f]]), logical(1L))
  if (!all(ok)) {
    dropped <- fields[!ok]
    .ferx_fit_addition_warning(paste0(
      "ferx_save_fit: ", paste0("fit$", dropped, collapse = ", "),
      if (length(dropped) == 1L) " is" else " are",
      " not plain data (numbers, text, lists, data frames, ...), so cannot be ",
      "written as JSON and is left out of the bundle."), dropped)
  }
  lapply(stats::setNames(fields[ok], fields[ok]), function(f) {
    as.character(jsonlite::serializeJSON(fit[[f]], digits = I(17)))
  })
}

# Is this parsed `serializeJSON()` tree data? Read before it is unpacked, which
# would otherwise run on whatever a bundle holds.
.fitrx_addition_wire_is_data <- function(node) {
  if (!is.list(node) || !is.character(node$type) || length(node$type) != 1L ||
        !node$type %in% .FERX_FIT_ADDITION_TYPES) {
    return(FALSE)
  }
  children <- c(if (identical(node$type, "list")) node$value, node$attributes)
  all(vapply(children, .fitrx_addition_wire_is_data, logical(1L)))
}

# One addition from its wire text. An error for anything that is not data.
.fitrx_addition_from_wire <- function(text) {
  if (!is.character(text) || length(text) != 1L ||
        !.fitrx_addition_wire_is_data(jsonlite::parse_json(text))) {
    stop("not plain data")
  }
  jsonlite::unserializeJSON(text)
}

# The sealed `fit` with the additions a bundle holds (`interface_extras.json`'s
# `user_additions`). One that cannot come back - it is not data, or a newer ferx
# made its name a core field - costs that addition, with one warning.
.fitrx_restore_user_additions <- function(fit, wire) {
  dropped <- character()
  for (name in names(wire)) {
    value <- tryCatch(
      if (is.null(.ferx_fit_addition_refusal(fit, name))) {
        list(.fitrx_addition_from_wire(wire[[name]]))
      },
      error = function(e) NULL
    )
    if (is.null(value)) {
      dropped <- c(dropped, name)
    } else {
      fit[name] <- value
    }
  }
  if (length(dropped) > 0L) {
    .ferx_fit_addition_warning(paste0(
      "ferx_load_fit: the bundle's user addition",
      if (length(dropped) > 1L) "s", " ", paste0("`", dropped, "`", collapse = ", "),
      " could not be restored: not plain data, or the name of a core field ",
      "of this fit. The rest of the fit loads."), dropped)
  }
  fit
}
