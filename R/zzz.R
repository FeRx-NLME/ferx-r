# Package-internal state: holds the trace path of the most recent ferx_fit()
# call (when optimizer_trace = TRUE), so ferx_trace() can be called with no
# argument to inspect the last run.
`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

# Coerce a fit's recorded file hash into the string the Rust bindings expect,
# where "" means "skip the integrity check". Three states are handled:
#   1. Non-empty hex string -> forwarded as-is; Rust enforces equality.
#   2. NULL (older binary that never populated the field) -> "".
#   3. NA_character_ (sha256_file failed at fit time; api.rs `.ok()` turned the
#      Err into None, stored R-side as NA) -> "". This coercion is essential:
#      `%||%` only catches NULL, and NA at the FFI boundary stringifies to "NA",
#      which compares unequal to every real digest and would trigger a spurious
#      "hash mismatch" error. Shared by ferx_sir() and ferx_covariance().
.ferx_hash_arg <- function(x) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x[[1L]]))) {
    ""
  } else {
    as.character(x)
  }
}

# The interaction flag the fit's objective ran under, handed to the engine by
# ferx_covariance() and ferx_sir().
#
# `fit$method` is the last *estimating* stage (the engine skips an
# evaluation-only IMP when it sets it), so it decides first. FOCE and FOCEI fix
# the flag (FALSE / TRUE), whatever a fit records: a bundle ferx-r wrote before
# #591 stores FALSE for every fit, because fit lists did not carry the flag
# then. Every other method (SAEM, IMP, IMPMAP, VI, Gauss-Newton, Laplace)
# inherits it from the options it ran under, so a recorded TRUE is taken as
# it stands (#591): an SAEM fit, or the estimating IMP of a chain
# `method = [foce, imp]`, chosen in a model file runs with interaction. A
# recorded FALSE cannot be told from that legacy FALSE, so it falls through
# to the chain reading below, as does a fit that records no method.
#
# The chain reading (legacy): the last stage that is not IMP, so
# c("focei", "imp") is interaction = TRUE - right for an evaluation-only IMP,
# and for an estimating one that inherited FOCEI's flag through ferx_fit().
.ferx_fit_interaction <- function(fit) {
  method <- toupper(as.character(fit$method %||% character()))
  if (length(method) == 1L && !is.na(method) && nzchar(method)) {
    if (method %in% c("FOCE", "FOCEI")) {
      return(identical(method, "FOCEI"))
    }
    recorded <- fit$interaction
    if (is.logical(recorded) && length(recorded) == 1L && isTRUE(recorded)) {
      return(TRUE)
    }
  }
  mchain <- toupper(as.character(fit$method_chain %||% fit$method))
  if (length(mchain) == 0L) return(FALSE)
  estimating <- mchain[mchain != "IMP"]
  last_method <- if (length(estimating) > 0L) {
    estimating[length(estimating)]
  } else {
    mchain[length(mchain)]
  }
  identical(last_method, "FOCEI")
}

# The thetas the fit held, as ferx_covariance() and ferx_sir() hand them to the
# engine (#613): its FIX flags plus any theta the engine froze during the fit
# for having no effect on the objective, which the model file does not record.
# The glue holds them on the step's skeleton as the fit's own step did; NA
# reads as not held, and logical() (a fit without the field) leaves the model
# file's FIX alone.
.ferx_fit_theta_fixed <- function(fit) {
  as.logical(fit$theta_fixed %||% logical())
}

# The fit's method and method chain in the engine's labels ("SAEM", "FOCEI",
# ...), as ferx_covariance() and ferx_sir() hand them to the engine (#591).
# The engine scores the objective of the fit's method and routes the
# covariance step on it; "" / character() leave the old FOCE / FOCEI reading.
.ferx_fit_method_args <- function(fit) {
  method <- as.character(fit$method %||% "")
  list(
    method = if (length(method) == 1L && !is.na(method)) method else "",
    method_chain = as.character(fit$method_chain %||% character())
  )
}

# `fit$ebe_kappas` flattened for the engine (#591): the kappa columns row-major,
# one row per (subject, occasion) in fit order, and the ID of each row. The
# engine starts the covariance step from the EBEs the fit reports only when
# their kappas come too; without them it re-solves cold. Empty for a fit
# without IOV.
.ferx_ebe_kappa_args <- function(fit) {
  k <- fit$ebe_kappas
  if (is.null(k) || !is.data.frame(k) || nrow(k) == 0L) {
    return(list(ebe_kappas_flat = numeric(0), ebe_kappa_ids = character(0)))
  }
  cols <- setdiff(names(k), c("ID", "OCC"))
  mat <- as.matrix(k[, cols, drop = FALSE])
  storage.mode(mat) <- "double"
  list(
    ebe_kappas_flat = as.numeric(t(mat)),
    ebe_kappa_ids = as.character(k$ID)
  )
}

.ferx_state <- new.env(parent = emptyenv())
.ferx_state$last_trace_path  <- NULL
.ferx_state$last_trace_time  <- NULL
.ferx_state$last_trace_model <- NULL

# The worker-pool width a standalone ferx_covariance() / ferx_sir() step
# runs on (#527), as the `threads` the glue takes: 0 is the engine default,
# n > 0 pins. NULL takes the width the fit ran on, `fit$n_threads_used`, so
# a fit pinned to 2 threads is not refined on every core; a fit that records
# no width (NULL / NA, e.g. a hand-built list) gets the engine default.
.ferx_step_threads <- function(threads, fit) {
  if (is.null(threads)) {
    used <- fit$n_threads_used
    if (is.numeric(used) && length(used) == 1L && is.finite(used) && used >= 1) {
      return(as.integer(used))
    }
    return(0L)
  }
  # The range before the whole-number test, as in ferx_sir() (#536).
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
        threads < 0 || threads > .Machine$integer.max ||
        threads != floor(threads)) {
    stop("`threads` must be NULL or a non-negative integer scalar")
  }
  as.integer(threads)
}
