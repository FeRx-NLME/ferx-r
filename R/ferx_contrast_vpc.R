#' Visual predictive check of a treatment contrast
#'
#' A VPC of the \strong{difference} between each active arm and the reference
#' arm of the same study, at the same time - the quantity a model-based
#' meta-analysis (MBMA) or an active-comparator trial is about. Every other
#' ferx diagnostic works on arm-level \code{DV} /
#' \code{PRED}, and an arm-based model can fit every arm well while
#' describing the contrast poorly.
#'
#' Within each level of \code{within} (the study) and each \code{TIME}, the
#' function forms \code{DV[by != ref] - DV[by == ref]} on the observed data
#' and, identically, on every replicate of
#' \code{\link{ferx_simulate}(..., fit = fit)}. It returns the paired data
#' frame and the VPC of the difference; \code{plot()} draws it.
#'
#' @section Why it is sharper than the arm-level VPC:
#' The contrast is taken within a study, so anything the arms of one study
#' share cancels: the between-study random effect (an \code{omega} on the
#' baseline), and the study-by-visit placebo levels. What is left is the
#' structural drug effect plus the variability that differs between arms
#' (between-arm \code{kappa}, residual error). The between-study spread that
#' dominates an arm-level VPC is differenced away, so a misspecified dose- or
#' time-response is not hidden behind it.
#'
#' What cancels is whatever the arms of one \code{within} level share. That
#' is the between-study random effect only when the model's random-effect
#' unit (the data's \code{ID}) is the study, as with \code{within = "ID"} in
#' arm-level MBMA data. When each arm is its own \code{ID} and
#' \code{within} names a study column, the \code{omega} is between-arm and
#' stays in the contrast.
#'
#' @section Pairing:
#' Arms are matched within a level of \code{within} and a \code{TIME}. Each
#' active arm gets its own contrast against the reference arm; there is one
#' row of \code{$pairs} per study, time and active arm.
#' \itemize{
#'   \item A study, or a time within a study, with no reference arm
#'     contributes nothing. They are reported together in one warning that
#'     names the count and the studies.
#'   \item A study with more than one reference arm is refused: there is no
#'     single difference to form, and ferx does not pick one.
#'   \item An arm is identified by its \code{by} value, or by the \code{arm}
#'     column when given. Two arms of one study with the same \code{by} value
#'     at the same time (two placebo arms, two arms at one dose) are refused
#'     unless \code{arm} tells them apart.
#'   \item Arms pair only on an identical \code{TIME}. A design whose arms
#'     are observed at different clock times (the periods of a crossover)
#'     finds no pairs.
#'   \item The observation records must hold one endpoint. Data whose
#'     observation records span more than one \code{CMT} are refused, since
#'     pairing on study and time would difference one endpoint against
#'     another; pass \code{data =} holding one endpoint's observations.
#' }
#'
#' @section Scale:
#' \code{scale} says on which scale the difference is formed. It is an
#' argument, never inferred from the model, and it is a transform of the
#' \code{DV} column \strong{as the data file holds it}: the default contrast
#' is in the units the arm-level readout was reported in.
#' \describe{
#'   \item{\code{"observed"}}{(default) \code{DV_active - DV_ref} on the data
#'     file's scale: score points, percent, concentration. For a readout the
#'     data already hold transformed - a weighted logit \code{DV = logit(p)},
#'     a pre-logged \code{DV ~ log_additive(...)} - this is the contrast on
#'     that scale (a log odds ratio, a log ratio).}
#'   \item{\code{"log"}}{\code{log(DV_active) - log(DV_ref)}, the log ratio
#'     of a positive \code{DV}.}
#'   \item{\code{"logit"}}{\code{logit(DV_active) - logit(DV_ref)}, the log
#'     odds ratio of a proportion held in (0, 1).}
#'   \item{\code{"expit"}}{\code{expit(DV_active) - expit(DV_ref)}, the
#'     difference in proportions for data holding \code{DV = logit(p)}.}
#' }
#' The observed \code{DV} and every simulated \code{DV_SIM} are put on the
#' data file's scale first, then transformed the same way. They differ in
#' one case: under \code{log(DV) ~ ...} the engine logs \code{DV} at load and
#' simulates on the log scale, so \code{DV_SIM} is exponentiated back before
#' \code{scale} is applied (\code{"log"} then gives the contrast on the
#' model's own scale). A simulated value outside the transform's domain gives
#' a non-finite contrast, which is dropped from the quantiles and counted in
#' a warning.
#'
#' @param fit A \code{ferx_fit} result. Its \code{model_path} and
#'   \code{data_path} supply the model and the design unless \code{model} /
#'   \code{data} are given.
#' @param by Name of the column that identifies the treatment (e.g.
#'   \code{"DOSE"} or \code{"TRT"}).
#' @param ref The value of \code{by} that marks the reference arm (e.g.
#'   \code{0} for placebo). Must occur in the data.
#' @param within Name of the column contrasts are formed within - the study.
#'   Default \code{"ID"}.
#' @param arm Optional name of a column that identifies the arm within a
#'   study (e.g. \code{"OCC"} in data from \code{\link{ferx_mbma_data}}).
#'   Needed only when two arms of one study share a \code{by} value.
#' @param n_sim Number of simulated replicates. Default 500.
#' @param seed Random seed passed to \code{\link{ferx_simulate}}.
#' @param scale The scale the contrast is formed on: \code{"observed"}
#'   (default), \code{"log"}, \code{"logit"} or \code{"expit"}. See the
#'   Scale section.
#' @param probs Lower, middle and upper quantile of the simulated contrast.
#'   Default \code{c(0.05, 0.5, 0.95)}: the median and a 90 percent band.
#' @param model,data Optional model and data file paths, overriding
#'   \code{fit$model_path} / \code{fit$data_path}. The data are the design
#'   the simulation runs on and the observations the contrast is taken from.
#' @param ... For \code{ferx_contrast_vpc()}, further arguments to
#'   \code{\link{ferx_simulate}} (e.g. \code{match}); for \code{plot()},
#'   graphical parameters passed to \code{plot()}.
#'
#' @section Plot:
#' \code{plot()} draws the VPC of the contrast per time, in the style of a
#' VPC: the shaded band is the \code{probs} interval of the simulated median
#' contrast, the solid line its median, the points the observed contrasts
#' (one per study and arm) and the dashed line their median. It returns
#' \code{x} invisibly.
#'
#' @return An object of class \code{ferx_contrast_vpc}, a list with:
#'   \describe{
#'     \item{\code{pairs}}{One row per study, time and active arm:
#'       the \code{within} column, \code{TIME}, \code{ARM} (the arm's
#'       \code{arm} value, or its \code{by} value), the \code{by} column
#'       (the active arm's treatment), \code{DV_REF} and \code{DV_ACTIVE}
#'       (observed, untransformed), \code{CONTRAST} (observed, on
#'       \code{scale}), \code{SIM_LO}, \code{SIM_MED}, \code{SIM_HI} (the
#'       \code{probs} quantiles of that pair's simulated contrast) and
#'       \code{PCTL}, the fraction of simulated contrasts below the observed
#'       one - near 0 or 1 is a pair the model does not reproduce.}
#'     \item{\code{vpc}}{The VPC summary: one row per \code{by} value and
#'       \code{TIME}, with \code{N} (pairs), \code{OBS_MED} (median observed
#'       contrast) and \code{SIM_LO}, \code{SIM_MED}, \code{SIM_HI}, the
#'       \code{probs} quantiles of the median over pairs, across replicates -
#'       the band the observed median should fall in.}
#'     \item{\code{sim}}{The simulated contrasts, long: \code{SIM}, the
#'       \code{within} column, \code{TIME}, \code{ARM}, the \code{by} column
#'       and \code{CONTRAST}.}
#'     \item{\code{excluded}}{The \code{within} values with no reference arm
#'       at one or more times, as warned.}
#'     \item{\code{sim_exponentiated}}{\code{TRUE} when \code{DV_SIM} was
#'       exponentiated back to the data's scale (a \code{log(DV) ~ ...}
#'       readout; see the Scale section).}
#'     \item{\code{within}, \code{by}, \code{ref}, \code{scale},
#'       \code{probs}, \code{n_sim}}{The settings used.}
#'   }
#'
#' @examples
#' \donttest{
#' ex <- ferx_example("mbma_placebo")
#' fit <- ferx_fit(ex$model, ex$data)
#' cv <- ferx_contrast_vpc(fit, within = "ID", by = "DOSE", ref = 0,
#'                         n_sim = 200)
#' head(cv$pairs)
#' cv$vpc
#' plot(cv)
#' }
#'
#' @family simulation
#' @export
ferx_contrast_vpc <- function(fit, by, ref, within = "ID", arm = NULL,
                              n_sim = 500L, seed = 42L,
                              scale = c("observed", "log", "logit", "expit"),
                              probs = c(0.05, 0.5, 0.95),
                              model = NULL, data = NULL, ...) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit() or ferx_load_fit()).",
         call. = FALSE)
  }
  scale <- match.arg(scale)
  .ferx_cvpc_check_name(by, "by")
  .ferx_cvpc_check_name(within, "within")
  if (!is.null(arm)) .ferx_cvpc_check_name(arm, "arm")
  if (length(unique(c(by, within, arm, "TIME"))) != length(c(by, within, arm)) + 1L) {
    stop("`by`, `within`, `arm` and TIME must be four different columns.",
         call. = FALSE)
  }
  if (length(ref) != 1L || is.na(ref)) {
    stop("`ref` must be a single non-missing value of the `by` column.",
         call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  if (length(n_sim) != 1L || is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(probs) || length(probs) != 3L || anyNA(probs) ||
        any(probs <= 0 | probs >= 1) || is.unsorted(probs, strictly = TRUE)) {
    stop("`probs` must be three increasing probabilities in (0, 1), ",
         "e.g. c(0.05, 0.5, 0.95).", call. = FALSE)
  }

  model <- model %||% fit$model_path
  data <- data %||% fit$data_path
  if (is.null(model) || is.na(model) || !file.exists(model)) {
    stop("No usable model file: pass `model=` or refit so `fit$model_path` is set.",
         call. = FALSE)
  }
  if (is.null(data) || is.na(data) || !file.exists(data)) {
    stop("No usable data file: pass `data=` or refit so `fit$data_path` is set.",
         call. = FALSE)
  }

  obs <- .ferx_cvpc_obs_rows(data, model, c(by, within, arm))
  # Pairing is on (within, TIME) only: two endpoints in the observation rows
  # would be differenced against each other (PD minus PK) without a word.
  cmt <- if ("CMT" %in% names(obs)) sort(unique(obs$CMT[!is.na(obs$CMT)])) else NULL
  if (length(cmt) > 1L) {
    stop("The observation records span ", length(cmt), " compartments (CMT ",
         paste(cmt, collapse = ", "), "). ferx_contrast_vpc() pairs arms on ",
         within, " and TIME only, so it would difference one endpoint against ",
         "another. Pass `data =` holding the observation records of one ",
         "endpoint.", call. = FALSE)
  }
  pairing <- .ferx_cvpc_pair(obs, by = by, ref = ref, within = within, arm = arm)
  if (length(pairing$excluded_msg)) warning(pairing$excluded_msg, call. = FALSE)
  if (length(pairing$active) == 0L) {
    stop("No study has both a reference arm (", by, " == ", format(ref),
         ") and an active arm at the same TIME; there is no contrast to check.",
         call. = FALSE)
  }

  f <- .ferx_cvpc_transform(scale)
  y_obs <- suppressWarnings(as.numeric(obs$DV))
  c_obs <- f(y_obs[pairing$active]) - f(y_obs[pairing$ref])
  if (any(!is.finite(c_obs))) {
    bad <- sum(!is.finite(c_obs))
    stop(bad, " observed contrast", if (bad > 1L) "s are" else " is",
         " missing or non-finite on scale = \"", scale, "\" (a missing DV",
         if (scale == "observed") "" else ", or a DV outside the transform's domain",
         "). Drop those records, or choose another `scale`.", call. = FALSE)
  }

  sim <- ferx_simulate(model, data, n_sim = n_sim, seed = seed, fit = fit, ...)
  y_sim <- .ferx_cvpc_sim_matrix(sim, obs, n_sim)
  # Under `log(DV) ~ ...` the engine logs DV at load and simulates on the log
  # scale; the data file holds DV on the natural scale. Put DV_SIM back on the
  # data's scale so both sides of every contrast are on one scale.
  log_readout <- .ferx_cvpc_log_dv_readout(model)
  if (log_readout) y_sim <- exp(y_sim)
  c_sim <- f(y_sim[pairing$active, , drop = FALSE]) -
    f(y_sim[pairing$ref, , drop = FALSE])
  n_bad <- sum(!is.finite(c_sim))
  if (n_bad > 0L) {
    warning(n_bad, " of ", length(c_sim), " simulated contrasts are ",
            "non-finite on scale = \"", scale, "\" (a simulated DV outside ",
            "the transform's domain) and are left out of the quantiles.",
            call. = FALSE)
    c_sim[!is.finite(c_sim)] <- NA_real_
  }

  act <- obs[pairing$active, , drop = FALSE]
  arm_key <- if (is.null(arm)) act[[by]] else act[[arm]]
  q <- t(apply(c_sim, 1L, stats::quantile, probs = probs, na.rm = TRUE,
               names = FALSE))
  pairs <- data.frame(
    act[[within]], TIME = act$TIME, ARM = arm_key, act[[by]],
    DV_REF = y_obs[pairing$ref], DV_ACTIVE = y_obs[pairing$active],
    CONTRAST = c_obs,
    SIM_LO = q[, 1L], SIM_MED = q[, 2L], SIM_HI = q[, 3L],
    PCTL = rowMeans(c_sim < c_obs, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  names(pairs)[c(1L, 4L)] <- c(within, by)
  rownames(pairs) <- NULL

  vpc <- .ferx_cvpc_summary(pairs, c_sim, probs = probs)

  sim_long <- data.frame(
    SIM = rep(seq_len(n_sim), each = nrow(pairs)),
    pairs[[within]], TIME = pairs$TIME, ARM = pairs$ARM, pairs[[4L]],
    CONTRAST = as.vector(c_sim),
    stringsAsFactors = FALSE
  )
  names(sim_long)[c(2L, 5L)] <- names(pairs)[c(1L, 4L)]

  structure(
    list(pairs = pairs, vpc = vpc, sim = sim_long,
         excluded = pairing$excluded,
         sim_exponentiated = log_readout,
         within = within, by = by, ref = ref, scale = scale, probs = probs,
         n_sim = n_sim),
    class = "ferx_contrast_vpc"
  )
}

#' @export
print.ferx_contrast_vpc <- function(x, ...) {
  scale_txt <- switch(
    x$scale,
    observed = "the observed DV scale (active - reference)",
    log = "the log scale (log ratio, active / reference)",
    logit = "the logit scale (log odds ratio, active vs reference)",
    expit = "the proportion scale (expit(DV), active - reference)"
  )
  cat("<ferx_contrast_vpc>\n")
  cat("  contrast: each ", x$by, " vs ", x$by, " == ", format(x$ref),
      ", within ", x$within, ", on ", scale_txt, "\n", sep = "")
  n_study <- length(unique(x$pairs[[x$within]]))
  cat("  ", nrow(x$pairs), " pairs from ", n_study, " ", x$within,
      " level", if (n_study == 1L) "" else "s", "; ", x$n_sim,
      " replicates; band ", 100 * x$probs[1L], "-", 100 * x$probs[3L],
      "%\n", sep = "")
  if (length(x$excluded)) {
    cat("  no reference arm (excluded): ",
        paste(x$excluded, collapse = ", "), "\n", sep = "")
  }
  outside <- x$vpc$OBS_MED < x$vpc$SIM_LO | x$vpc$OBS_MED > x$vpc$SIM_HI
  cat("  observed median outside the band in ", sum(outside), " of ",
      nrow(x$vpc), " ", x$by, " x TIME cells\n", sep = "")
  invisible(x)
}

#' @rdname ferx_contrast_vpc
#' @param x A \code{ferx_contrast_vpc} object.
#' @param y Ignored.
#' @param facet \code{"by"} (default) draws one panel per active
#'   \code{by} value (e.g. per dose); \code{"none"} pools every active arm
#'   in one panel.
#' @export
plot.ferx_contrast_vpc <- function(x, y = NULL, facet = c("by", "none"), ...) {
  facet <- match.arg(facet)
  by_col <- names(x$pairs)[4L]
  vpc <- if (facet == "by") x$vpc else .ferx_cvpc_pooled_vpc(x)
  panels <- if (facet == "by") unique(vpc[[by_col]]) else "all"
  ylab <- switch(x$scale, observed = "Contrast (active - reference)",
                 log = "Log ratio (active / reference)",
                 logit = "Log odds ratio (active vs reference)",
                 expit = "Difference in proportion (active - reference)")
  n <- length(panels)
  nc <- ceiling(sqrt(n))
  old <- graphics::par(mfrow = c(ceiling(n / nc), nc),
                       mar = c(4, 4, 2.5, 1))
  on.exit(graphics::par(old), add = TRUE)
  for (p in panels) {
    v <- if (facet == "by") vpc[vpc[[by_col]] == p, , drop = FALSE] else vpc
    pts <- if (facet == "by") x$pairs[x$pairs[[by_col]] == p, , drop = FALSE] else x$pairs
    v <- v[order(v$TIME), , drop = FALSE]
    yr <- range(c(v$SIM_LO, v$SIM_HI, pts$CONTRAST, 0), finite = TRUE)
    graphics::plot(range(v$TIME), yr, type = "n", xlab = "TIME", ylab = ylab,
                   main = if (facet == "by") paste0(by_col, " = ", p) else
                     paste0("All active arms vs ", by_col, " = ", format(x$ref)),
                   ...)
    graphics::polygon(c(v$TIME, rev(v$TIME)), c(v$SIM_LO, rev(v$SIM_HI)),
                      col = "#4682B440",
                      border = NA)
    graphics::lines(v$TIME, v$SIM_MED, col = "steelblue", lwd = 2)
    graphics::abline(h = 0, col = "grey60", lty = 3)
    graphics::points(pts$TIME, pts$CONTRAST, pch = 16,
                     col = "#00000099")
    graphics::lines(v$TIME, v$OBS_MED, lty = 2)
  }
  invisible(x)
}

# -- Internals ------------------------------------------------------------

.ferx_cvpc_check_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be a single column name.", call. = FALSE)
  }
}

# The data rows ferx_simulate() returns a record for, in file order: not a
# dose (EVID 0, or no AMT where there is no EVID column) and not MDV = 1.
.ferx_cvpc_obs_rows <- function(data, model, cols) {
  d <- utils::read.csv(data, stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("NA", ".", ""), colClasses = "character")
  # Columns under the names the engine's reader gives them (a [data] map).
  d <- .ferx_apply_column_map(d, .ferx_model_column_map(model))
  # ID stays text, spelled as in the file (`001`), as ferx_simulate() returns
  # it; every other column is typed as read.csv() would type it.
  for (nm in setdiff(names(d), "ID")) {
    d[[nm]] <- utils::type.convert(d[[nm]], as.is = TRUE)
  }
  need <- unique(c("ID", "TIME", "DV", cols))
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop("Column", if (length(miss) > 1L) "s" else "", " not found in the data: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(d))
  if ("EVID" %in% names(d)) {
    keep <- keep & !is.na(d$EVID) & d$EVID == 0
  } else if ("AMT" %in% names(d)) {
    keep <- keep & (is.na(d$AMT) | d$AMT == 0)
  }
  if ("MDV" %in% names(d)) keep <- keep & (is.na(d$MDV) | d$MDV != 1)
  d <- d[keep, , drop = FALSE]
  rownames(d) <- NULL
  d
}

# Index the reference row for every active row, by (within, TIME).
.ferx_cvpc_pair <- function(obs, by, ref, within, arm) {
  b <- obs[[by]]
  if (anyNA(b)) {
    stop("The `by` column ", by, " is missing on ", sum(is.na(b)),
         " observation record(s).", call. = FALSE)
  }
  is_ref <- as.character(b) == as.character(ref)
  if (!any(is_ref)) {
    stop("No observation record has ", by, " == ", format(ref),
         "; `ref` must be a value of the `by` column (found: ",
         paste(utils::head(sort(unique(b)), 10L), collapse = ", "), ").",
         call. = FALSE)
  }
  w <- as.character(obs[[within]])
  arm_id <- as.character(if (is.null(arm)) b else obs[[arm]])

  # One reference arm per study: a second is refused, never guessed between.
  ref_arms <- tapply(arm_id[is_ref], w[is_ref], function(a) length(unique(a)))
  multi <- names(ref_arms)[ref_arms > 1L]
  if (length(multi)) {
    stop(length(multi), " ", within, " level", if (length(multi) > 1L) "s have" else " has",
         " more than one reference arm (", by, " == ", format(ref), "): ",
         paste(multi, collapse = ", "), ". A contrast needs exactly one ",
         "reference per study; drop or merge the extra reference arm.",
         call. = FALSE)
  }
  cell <- paste(w, obs$TIME, sep = "\r")
  dup <- duplicated(paste(cell, arm_id, sep = "\r"))
  if (any(dup)) {
    stop(sum(dup), " observation record(s) repeat an arm at the same ", within,
         " and TIME (first: ", within, " ", w[dup][1L], ", TIME ",
         obs$TIME[dup][1L], ", arm ", arm_id[dup][1L], "). ",
         if (is.null(arm)) {
           paste0("Two arms share a `by` value there; pass `arm =` naming ",
                  "the column that tells them apart.")
         } else "Each arm must have one record per TIME.", call. = FALSE)
  }

  ref_row <- match(cell, cell[is_ref])
  ref_row <- which(is_ref)[ref_row]
  active <- which(!is_ref & !is.na(ref_row))
  orphan <- !is_ref & is.na(ref_row)

  excluded <- character(0)
  msg <- character(0)
  if (any(orphan)) {
    has_ref <- unique(w[is_ref])
    no_ref_study <- unique(w[orphan & !(w %in% has_ref)])
    part_study <- unique(w[orphan & (w %in% has_ref)])
    excluded <- unique(w[orphan])
    n_cells <- length(unique(cell[orphan]))
    parts <- character(0)
    if (length(no_ref_study)) {
      parts <- c(parts, paste0(
        length(no_ref_study), " ", within, " level",
        if (length(no_ref_study) > 1L) "s have" else " has",
        " no reference arm at all (", paste(no_ref_study, collapse = ", "), ")"
      ))
    }
    if (length(part_study)) {
      parts <- c(parts, paste0(
        length(part_study), " ", if (length(no_ref_study)) "more " else "",
        "lack", if (length(part_study) == 1L) "s" else "",
        " one at some times (", paste(part_study, collapse = ", "), ")"
      ))
    }
    msg <- paste0(
      "ferx_contrast_vpc: ", n_cells, " ", within, " x TIME cell",
      if (n_cells > 1L) "s have" else " has", " no reference arm (", by,
      " == ", format(ref), ") and contribute", if (n_cells == 1L) "s" else "",
      " no contrast: ", paste(parts, collapse = "; "), "."
    )
  }
  list(active = active, ref = ref_row[active], excluded = excluded,
       excluded_msg = msg)
}

.ferx_cvpc_transform <- function(scale) {
  switch(scale,
         observed = identity,
         log = function(y) suppressWarnings(log(y)),
         logit = function(y) suppressWarnings(stats::qlogis(y)),
         expit = stats::plogis)
}

# TRUE when the model's [error_model] reads `log(DV) ~ ...`: the one readout
# the engine transforms itself, simulating DV on the log scale while the data
# hold it on the natural scale. `DV ~ log_additive(...)` keeps both on the
# data's (already log) scale and is FALSE here. One answer per model is
# enough: the engine refuses LTBS in a per-CMT (multi-endpoint) error model.
.ferx_cvpc_log_dv_readout <- function(model) {
  err <- .ferx_extract_blocks(model)[["error_model"]] %||% character(0)
  any(grepl("^(CMT\\s*=\\s*[0-9]+\\s*:\\s*)?log\\s*\\(\\s*DV\\s*\\)\\s*~",
            err, ignore.case = TRUE))
}

# DV_SIM as an (observation record) x (replicate) matrix. Each replicate holds
# the data's observation records in file order; that is checked here on ID and
# TIME rather than assumed.
.ferx_cvpc_sim_matrix <- function(sim, obs, n_sim) {
  n <- nrow(obs)
  cont <- is.na(sim$OBSERVED)
  sim <- sim[cont, , drop = FALSE]
  if (nrow(sim) != n * n_sim || is.unsorted(sim$SIM)) {
    stop("ferx_simulate() returned ", nrow(sim), " continuous records for ",
         n_sim, " replicates of ", n, " observation records; the contrast ",
         "cannot be paired by position. Records the simulation adds or drops ",
         "relative to the data file (binary or time-to-event rows, a data ",
         "selection the engine applies) are not supported.", call. = FALSE)
  }
  first <- sim[sim$SIM == sim$SIM[1L], , drop = FALSE]
  same <- as.character(first$ID) == as.character(obs$ID) &
    abs(as.numeric(first$TIME) - as.numeric(obs$TIME)) <= 1e-8 *
      pmax(1, abs(as.numeric(obs$TIME)))
  if (!all(same)) {
    i <- which(!same)[1L]
    stop("ferx_simulate()'s records do not line up with the data's ",
         "observation records (first mismatch at record ", i, ": data ID ",
         obs$ID[i], " TIME ", obs$TIME[i], ", simulated ID ", first$ID[i],
         " TIME ", first$TIME[i], ").", call. = FALSE)
  }
  matrix(sim$DV_SIM, nrow = n, ncol = n_sim)
}

# Per (by, TIME): the observed median contrast over pairs, and the quantiles
# across replicates of the simulated median over the same pairs.
.ferx_cvpc_summary <- function(pairs, c_sim, probs, groups = NULL) {
  by_col <- names(pairs)[4L]
  key <- if (is.null(groups)) {
    paste(pairs[[by_col]], pairs$TIME, sep = "\r")
  } else {
    paste(groups, pairs$TIME, sep = "\r")
  }
  ukey <- unique(key)
  rows <- lapply(ukey, function(k) {
    i <- which(key == k)
    med <- apply(c_sim[i, , drop = FALSE], 2L, stats::median, na.rm = TRUE)
    q <- stats::quantile(med, probs = probs, na.rm = TRUE, names = FALSE)
    data.frame(
      if (is.null(groups)) pairs[[by_col]][i[1L]] else groups[i[1L]],
      TIME = pairs$TIME[i[1L]], N = length(i),
      OBS_MED = stats::median(pairs$CONTRAST[i]),
      SIM_LO = q[1L], SIM_MED = q[2L], SIM_HI = q[3L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out)[1L] <- if (is.null(groups)) by_col else "GROUP"
  out <- out[order(out[[1L]], out$TIME), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Every active arm pooled per TIME, for plot(facet = "none").
.ferx_cvpc_pooled_vpc <- function(x) {
  c_sim <- matrix(x$sim$CONTRAST, nrow = nrow(x$pairs), ncol = x$n_sim)
  .ferx_cvpc_summary(x$pairs, c_sim, probs = x$probs,
                     groups = rep("all", nrow(x$pairs)))
}
