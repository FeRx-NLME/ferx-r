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

# -- Reporting a VI fit (#507) ------------------------------------------------
#
# VI estimates a bound, not an objective: under `vi_final_ofv = none` the fit
# has no OFV, so `fit$ofv`, `fit$aic` and `fit$bic` are NaN by design. Every
# surface that would show those NaNs - print(), summary(), ferx_bic(),
# ferx_runlog() - says why, in the same words, through these helpers.

# ferx-core's `ElboTightness::IMPLAUSIBLE_RATIO` (src/estimation/vi/elbo.rs):
# above it the engine warns `W_VI_BAD_BASIN` and reports `converged: false`.
# Not reachable through ferx-core's public API, so mirrored here.
.ferx_vi_tightness_limit <- 25

# TRUE when `fit` is a VI fit whose OFV is missing because VI reports a bound.
# A VI result with a finite OFV (vi_final_ofv = laplace, an eval-only Laplace
# stage, or a later estimating stage that superseded VI) has an OFV to show.
.ferx_vi_ofv_missing <- function(fit) {
  if (is.null(fit$vi)) return(FALSE)
  ofv <- suppressWarnings(as.numeric(fit$ofv %||% NA_real_))
  length(ofv) != 1L || !is.finite(ofv)
}

# The explanation, one sentence on why and one on what to do instead.
.ferx_vi_ofv_note <- function() {
  paste0(
    "OFV is NaN because VI estimates a bound (-2*ELBO) rather than an ",
    "objective value, so the fit has no OFV, AIC or BIC. For a comparable ",
    "OFV, set vi_final_ofv = \"laplace\", or chain method = c(\"vi\", ",
    "\"laplace\") with agq_eval_only = TRUE."
  )
}

# The VI block's lines, unindented and unstyled, shared by print.ferx_fit(),
# print.ferx_summary() and ferx_runlog(). NULL for a fit without `fit$vi`.
.ferx_vi_report_lines <- function(fit, width = 72L) {
  v <- fit$vi
  if (is.null(v)) return(NULL)
  fmt <- function(x, spec = "%.4f") {
    x <- suppressWarnings(as.numeric(x %||% NA_real_))
    if (length(x) != 1L || is.na(x)) "NA" else sprintf(spec, x)
  }
  int <- function(x) fmt(x, "%.0f")
  lines <- sprintf("-2*ELBO: %s  (a bound, not an OFV; data term %s, KL term %s)",
                   fmt(v$neg_two_elbo), fmt(v$data_term), fmt(v$kl_term))
  lines <- c(lines, sprintf(
    "Converged: %s    Iterations: %s    Family: %s    MC samples: %s    KL: %s",
    if (isTRUE(v$converged)) "YES" else "NO", int(v$n_iterations),
    v$family %||% "NA", int(v$n_mc_samples), v$kl %||% "NA"
  ))
  ratio <- suppressWarnings(as.numeric(v$elbo_tightness_ratio %||% NA_real_))
  if (length(ratio) == 1L && !is.na(ratio)) {
    flag <- if (is.finite(ratio) && ratio > .ferx_vi_tightness_limit) {
      sprintf("  [!] above %g: the bound is not usable here (bad basin)",
              .ferx_vi_tightness_limit)
    } else {
      ""
    }
    lines <- c(lines, sprintf("ELBO tightness ratio: %s%s", fmt(ratio, "%.2f"), flag))
  }
  n_kl <- suppressWarnings(as.numeric(v$n_kl_fallback_subjects %||% 0))
  if (length(n_kl) == 1L && !is.na(n_kl) && n_kl > 0) {
    lines <- c(lines, sprintf("KL sampled (no closed form) for %s subject(s)", int(n_kl)))
  }
  n_fd <- suppressWarnings(as.numeric(v$n_fd_subjects %||% 0))
  if (length(n_fd) == 1L && !is.na(n_fd) && n_fd > 0) {
    lines <- c(lines, sprintf(
      "Finite-difference eta gradients for %s subject(s) (slower, not wrong)", int(n_fd)))
  }
  # A later estimating stage re-estimated after VI: its EBEs are the fit's,
  # and everything in `fit$vi` describes the point VI ended on. Otherwise -
  # an eval-only readout included, which keeps the preceding stage's EBEs -
  # the reported EBEs are VI's posterior means, not modes.
  sup <- v$superseded_by
  if (length(sup) == 1L && !is.na(sup) && nzchar(sup)) {
    lines <- c(lines, strwrap(sprintf(
      "Superseded by %s: these VI numbers describe the estimates VI ended on, not the reported ones.",
      sup), width = width))
  } else {
    lines <- c(lines, "EBEs, sdtab ETAs and shrinkage are the VI posterior means, not modes.")
  }
  if (.ferx_vi_ofv_missing(fit)) {
    lines <- c(lines, strwrap(.ferx_vi_ofv_note(), width = width))
  }
  lines
}
