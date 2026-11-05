# The variational-inference result, `fit$vi` (method = "vi", #507 / #508).
#
# One constructor, `.ferx_vi_from_raw()`, builds `fit$vi` from the engine's
# layout - per-subject vectors, matrices flattened row-major - for both the
# live fit and a reloaded one, both the glue's `vi` list, so the two cannot drift. `.fitrx_build_vi_wire()` writes it back in the
# layout ferx-core's own `.fitrx` uses (`ViWire` in src/io/fitrx.rs, the YAML
# writer's `vi:` section, ferx-core #2058), so either side loads the other's
# bundle.

# `raw` holds:
#   scalars        neg_two_elbo, data_term, kl_term, n_iterations, converged,
#                  family, n_mc_samples, kl, n_kl_fallback_subjects,
#                  n_fd_subjects, elbo_tightness_ratio, superseded_by
#                  ("" or NA when unset)
#   elbo_trace, objective_trace   numeric; an empty objective_trace is the
#                  engine's "identical to elbo_trace"
#   ids            the subject of each eta_means entry
#   eta_means      list of per-subject numeric vectors
#   eta_covs, kappa_means   lists of row-major numeric vectors, with
#                  eta_covs_nrow / kappa_means_nrow giving each one's rows
.ferx_vi_from_raw <- function(raw, eta_names, kappa_names) {
  if (is.null(raw) || !length(raw)) return(NULL)
  num <- function(x) .fitrx_num_or_nan(x)
  count <- function(x) as.numeric(unlist(x, use.names = FALSE) %||% NA_real_)[[1L]]
  ids <- as.character(unlist(raw$ids, use.names = FALSE) %||% character())
  means <- lapply(raw$eta_means %||% list(), as.numeric)

  n_eta <- if (length(means)) length(means[[1L]]) else length(eta_names)
  eta_labels <- if (length(eta_names) == n_eta) {
    as.character(eta_names)
  } else {
    paste0("ETA", seq_len(n_eta))
  }

  eta_means <- NULL
  if (length(means)) {
    if (any(lengths(means) != n_eta)) {
      stop("fit$vi: the per-subject posterior means differ in length.",
           call. = FALSE)
    }
    m <- matrix(unlist(means), nrow = length(means), byrow = TRUE)
    colnames(m) <- eta_labels
    eta_means <- data.frame(ID = ids, m, check.names = FALSE,
                            stringsAsFactors = FALSE)
  }

  # Row-major vectors back to matrices, named by subject. A list shorter than
  # `ids` is one the engine stopped short of; an entry with no rows is a
  # subject it has nothing for.
  as_mats <- function(vals, nrows, dimnames_for) {
    vals <- vals %||% list()
    nrows <- as.integer(unlist(nrows, use.names = FALSE) %||% integer())
    out <- lapply(seq_along(vals), function(i) {
      v <- as.numeric(vals[[i]])
      nr <- nrows[[i]]
      nc <- if (nr > 0L) length(v) %/% nr else 0L
      m <- matrix(v, nrow = nr, ncol = nc, byrow = TRUE)
      dimnames(m) <- dimnames_for(nr, nc)
      m
    })
    names(out) <- ids[seq_along(out)]
    out
  }
  eta_covs <- as_mats(raw$eta_covs, raw$eta_covs_nrow, function(nr, nc) {
    if (nr == n_eta && nc == n_eta) list(eta_labels, eta_labels) else NULL
  })
  kappa_means <- as_mats(raw$kappa_means, raw$kappa_means_nrow, function(nr, nc) {
    if (nr > 0L && nc == length(kappa_names)) {
      list(NULL, as.character(kappa_names))
    } else {
      NULL
    }
  })

  elbo_trace <- as.numeric(unlist(raw$elbo_trace, use.names = FALSE) %||% numeric())
  objective_trace <- as.numeric(unlist(raw$objective_trace, use.names = FALSE) %||% numeric())
  if (!length(objective_trace)) objective_trace <- elbo_trace
  superseded_by <- as.character(unlist(raw$superseded_by, use.names = FALSE))
  if (!length(superseded_by) || is.na(superseded_by[[1L]]) ||
        !nzchar(superseded_by[[1L]])) {
    superseded_by <- NA_character_
  }

  list(
    neg_two_elbo = num(raw$neg_two_elbo),
    data_term = num(raw$data_term),
    kl_term = num(raw$kl_term),
    n_iterations = count(raw$n_iterations),
    converged = isTRUE(as.logical(unlist(raw$converged))),
    family = as.character(unlist(raw$family) %||% NA_character_)[[1L]],
    n_mc_samples = count(raw$n_mc_samples),
    kl = as.character(unlist(raw$kl) %||% NA_character_)[[1L]],
    n_kl_fallback_subjects = count(raw$n_kl_fallback_subjects),
    elbo_trace = elbo_trace,
    objective_trace = objective_trace,
    eta_means = eta_means,
    eta_covs = eta_covs,
    kappa_means = kappa_means,
    n_fd_subjects = count(raw$n_fd_subjects),
    elbo_tightness_ratio = num(raw$elbo_tightness_ratio),
    superseded_by = superseded_by[[1L]]
  )
}

# A matrix as the list of its rows, the wire's `Vec<Vec<f64>>`.
.fitrx_vi_rows <- function(m) {
  m <- as.matrix(m)
  lapply(seq_len(nrow(m)), function(i) unname(as.numeric(m[i, ])))
}

# `fit$vi` as ferx-core's `ViWire`. NULL for a fit without one, so the key is
# left out, as the engine leaves it out for every other method.
.fitrx_build_vi_wire <- function(fit) {
  v <- fit$vi
  if (is.null(v)) return(NULL)
  em <- v$eta_means
  ids <- if (is.null(em)) character() else as.character(em$ID)
  eta_cols <- if (is.null(em)) character() else setdiff(names(em), "ID")
  covs <- v$eta_covs %||% list()
  kappas <- v$kappa_means %||% list()
  posterior <- lapply(seq_along(ids), function(i) {
    p <- list(
      id = ids[[i]],
      mean = unname(as.numeric(unlist(em[i, eta_cols], use.names = FALSE)))
    )
    if (i <= length(covs)) p$cov <- .fitrx_vi_rows(covs[[i]])
    if (i <= length(kappas)) p$kappa_means <- .fitrx_vi_rows(kappas[[i]])
    p
  })
  elbo <- as.numeric(v$elbo_trace %||% numeric())
  obj <- as.numeric(v$objective_trace %||% numeric())
  out <- list(
    neg_two_elbo = as.numeric(v$neg_two_elbo),
    data_term = as.numeric(v$data_term),
    kl_term = as.numeric(v$kl_term),
    family = as.character(v$family),
    n_iterations = as.integer(v$n_iterations),
    n_mc_samples = as.integer(v$n_mc_samples),
    kl = as.character(v$kl),
    n_kl_fallback_subjects = as.integer(v$n_kl_fallback_subjects),
    converged = isTRUE(v$converged),
    n_fd_subjects = as.integer(v$n_fd_subjects),
    superseded_by = .fitrx_opt_chr(v$superseded_by),
    elbo_tightness_ratio = as.numeric(v$elbo_tightness_ratio),
    elbo_trace = elbo,
    # The engine stores an objective trace identical to the bound as empty
    # and omits the key; `fit$vi` fills it back in, so leave it out here.
    objective_trace = if (identical(obj, elbo)) NULL else obj
  )
  if (length(posterior)) {
    out$eta_names <- as.character(fit$eta_names %||% eta_cols)
    out$eta_posterior <- posterior
  }
  out[!vapply(out, is.null, logical(1L))]
}
