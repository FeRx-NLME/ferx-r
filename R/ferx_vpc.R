#' Visual predictive check of a fit's continuous observations
#'
#' Simulates \code{n_sim} replicates of the fit's design at its estimates and
#' hands them, with the fit's own observations, to \code{vpc::vpc()}. With
#' \code{pred_corr = TRUE} the result is a prediction-corrected VPC (pcVPC).
#' Binning, the prediction intervals and the plot are the \pkg{vpc}
#' package's; \code{ferx_vpc()} only decides which rows go in.
#'
#' @details
#' The observations are the continuous (Gaussian) rows of \code{fit$sdtab}: the
#' records the fit used, after any \code{[data_selection]} block and with its
#' \code{PRED}. The simulation (\code{\link{ferx_simulate}(..., fit = fit)})
#' covers the whole data file, so each replicate is matched to those
#' observations by compartment, subject, time and, for repeated times, the
#' order within them - never by position. Simulated records the fit did not use
#' (for example the rows a \code{[data_selection]} \code{ignore} rule dropped)
#' are left out, and every simulated row carries the \code{PRED} of the
#' observation it matches, which is what the pcVPC needs. If a replicate does
#' not match every observation exactly once the call is refused.
#'
#' \strong{Censoring.} An observation with \code{CENS = 1} (below the LLOQ) or
#' \code{CENS = -1} (above the ULOQ) carries the limit, not a measurement, in
#' \code{DV}; \code{ferx_vpc()} blanks it, so it is not plotted as an observed
#' value. \code{lloq} defaults to the \code{DV} of the \code{CENS = 1} rows and
#' \code{uloq} to that of the \code{CENS = -1} rows; \pkg{vpc} takes one limit,
#' so several values are refused unless you pass \code{lloq =} /
#' \code{uloq =}. It also takes only one of the two: a fit censored on both
#' sides is refused until you pass the limit to draw, and a limit you pass
#' replaces both defaults (the censored rows stay blanked either way).
#' \pkg{vpc} draws the limit. A prediction-corrected VPC with a limit is
#' refused: \pkg{vpc} (1.2.4) censors the simulated values at it and then fails
#' on them, so use \code{pred_corr = FALSE} for such a fit.
#'
#' \strong{Several endpoints.} A fit with more than one continuous endpoint
#' (\code{CMT} values in \code{fit$sdtab}) needs \code{cmt =}. A time-to-event
#' endpoint has its own VPC, \code{\link{ferx_vpc_tte}()}; a joint PK-TTE fit's
#' simulation needs \code{horizon =}, even though only the continuous endpoint
#' is shown here.
#'
#' @param fit A \code{ferx_fit} from \code{\link{ferx_fit}()} or
#'   \code{\link{ferx_load_fit}()}. Its model and data files are re-read for
#'   the simulation.
#' @param n_sim Number of simulated replicates.
#' @param seed Random seed for the simulation (a non-negative integer).
#' @param pred_corr \code{TRUE} for a prediction-corrected VPC.
#' @param cmt The continuous endpoint's compartment. Required when the fit has
#'   more than one, unless \code{stratify} includes \code{"CMT"}, which shows
#'   every continuous endpoint, one stratum each.
#' @param stratify Optional character vector of data columns to stratify by
#'   (read from the fit's data file, observation records), or \code{"CMT"}.
#'   Passed to \code{vpc::vpc()}.
#' @param lloq,uloq Lower / upper limit of quantification; at most one. With
#'   both \code{NULL}, each is taken from the fit's \code{CENS = 1} /
#'   \code{CENS = -1} rows.
#' @param horizon Administrative horizon for the simulation; required by a
#'   joint PK-TTE model (see \code{\link{ferx_simulate}()}).
#' @param ... Further arguments to \code{vpc::vpc()}, for example \code{bins},
#'   \code{n_bins}, \code{pi}, \code{ci}, \code{log_y},
#'   \code{pred_corr_lower_bnd}, \code{xlab} or \code{vpcdb = TRUE}. The
#'   arguments \code{ferx_vpc()} builds itself (\code{sim}, \code{obs},
#'   \code{sim_cols}, \code{obs_cols}, \code{psn_folder}, \code{software}) are
#'   refused.
#'
#' @return What \code{vpc::vpc()} returns: a ggplot object, or with
#'   \code{vpcdb = TRUE} its data (\code{obs} and \code{sim} with \pkg{vpc}'s
#'   column names, \code{TIME} as \code{idv}).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("vpc", quietly = TRUE)) {
#'   ex <- ferx_example("warfarin")
#'   fit <- ferx_fit(ex$model, ex$data, verbose = FALSE)
#'   ferx_vpc(fit, n_sim = 100)
#'   ferx_vpc(fit, n_sim = 100, pred_corr = TRUE)
#' }
#' }
#'
#' @seealso \code{\link{ferx_simulate}}, \code{\link{ferx_vpc_tte}},
#'   \code{\link{ferx_contrast_vpc}}
#' @family simulation
#' @export
ferx_vpc <- function(fit, n_sim = 200L, seed = 42L, pred_corr = FALSE,
                     cmt = NULL, stratify = NULL, lloq = NULL, uloq = NULL,
                     horizon = NULL, ...) {
  if (!.ferx_has_vpc()) {
    .ferx_vpc_stop(paste0(
      "ferx_vpc() needs the vpc package. Install it with ",
      "install.packages(\"vpc\")."
    ))
  }
  if (!inherits(fit, "ferx_fit")) {
    .ferx_vpc_stop(
      "`fit` must be a fit from ferx_fit() or ferx_load_fit()."
    )
  }
  dots <- list(...)
  .ferx_vpc_check_dots(dots)
  n_sim <- .ferx_vpc_count(n_sim, "n_sim", min = 1L)
  seed <- .ferx_vpc_count(seed, "seed", min = 0L)
  if (!is.logical(pred_corr) || length(pred_corr) != 1L || is.na(pred_corr)) {
    .ferx_vpc_stop("`pred_corr` must be TRUE or FALSE.")
  }
  lloq <- .ferx_vpc_limit(lloq, "lloq")
  uloq <- .ferx_vpc_limit(uloq, "uloq")
  if (!is.null(stratify) &&
      (!is.character(stratify) || length(stratify) == 0L ||
       anyNA(stratify) || anyDuplicated(stratify))) {
    .ferx_vpc_stop("`stratify` must be NULL or distinct column names.")
  }

  fit <- .ferx_fit_canonical(fit, "ferx_vpc")
  frames <- .ferx_vpc_frames(fit, n_sim = n_sim, seed = seed, cmt = cmt,
                             stratify = stratify, lloq = lloq, uloq = uloq,
                             horizon = horizon, pred_corr = pred_corr)

  args <- c(
    list(
      sim = frames$sim,
      obs = frames$obs,
      sim_cols = list(id = "ID", dv = "DV_SIM", idv = "TIME", sim = "SIM",
                      pred = "PRED"),
      obs_cols = list(id = "ID", dv = "DV", idv = "TIME", pred = "PRED"),
      pred_corr = pred_corr,
      stratify = stratify,
      lloq = frames$lloq,
      uloq = frames$uloq
    ),
    dots
  )
  do.call(vpc::vpc, args)
}

# The one place that asks whether vpc is installed, so a test can mock it
# (testthat cannot mock base::requireNamespace).
.ferx_has_vpc <- function() requireNamespace("vpc", quietly = TRUE)

.ferx_vpc_stop <- function(msg, ...) {
  stop(errorCondition(msg, ..., class = "ferx_vpc_error", call = NULL))
}

# vpc::vpc() arguments ferx_vpc() builds itself, and what to use instead.
.FERX_VPC_OWNED <- c(
  sim = "`n_sim =`, `seed =` and `horizon =`",
  obs = "`fit =`",
  psn_folder = "`fit =`",
  sim_cols = "nothing: ferx_vpc() names the columns itself",
  obs_cols = "nothing: ferx_vpc() names the columns itself",
  software = "nothing: the input is always a ferx fit"
)

.ferx_vpc_check_dots <- function(dots) {
  if (length(dots) == 0L) return(invisible())
  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    .ferx_vpc_stop("Arguments passed on to vpc::vpc() through `...` must be named.")
  }
  owned <- intersect(nms, names(.FERX_VPC_OWNED))
  if (length(owned) > 0L) {
    a <- owned[[1L]]
    .ferx_vpc_stop(paste0(
      "`", a, "` is built by ferx_vpc() from the fit and cannot be passed to ",
      "vpc::vpc() through `...`; use ", .FERX_VPC_OWNED[[a]], " instead. ",
      "To supply your own, call vpc::vpc() directly."
    ))
  }
  invisible()
}

.ferx_vpc_count <- function(x, name, min) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x == round(x) && x >= min
  if (!ok) {
    .ferx_vpc_stop(sprintf("`%s` must be a single whole number >= %d.", name, min))
  }
  as.integer(x)
}

.ferx_vpc_limit <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    .ferx_vpc_stop(sprintf("`%s` must be NULL or a single finite number.", name))
  }
  as.numeric(x)
}

.ferx_vpc_list <- function(x) {
  paste(format(sort(unique(x)), trim = TRUE, drop0trailing = TRUE),
        collapse = ", ")
}

# The observation and simulation frames vpc::vpc() gets, with the limits it
# draws. `obs` has ID, TIME, DV (censored rows blanked), PRED, CMT and the
# `stratify` columns; `sim` the same keys with SIM and DV_SIM, one block of
# rows per replicate in `obs` order.
.ferx_vpc_frames <- function(fit, n_sim, seed, cmt, stratify, lloq, uloq,
                             horizon, pred_corr = FALSE) {
  sdtab <- fit$sdtab
  blocks <- .fitrx_sdtab_blocks(fit)
  if (length(blocks$gaussian) == 0L) {
    .ferx_vpc_stop(paste0(
      "ferx_vpc() shows continuous observations, and this fit has none. ",
      "For a time-to-event endpoint use ferx_vpc_tte() or ferx_vpc_cif(); ",
      "for a categorical endpoint, vpc::vpc_cat()."
    ))
  }
  obs <- sdtab[blocks$gaussian, , drop = FALSE]
  rownames(obs) <- NULL

  # The endpoint: from the sdtab's CMT when the fit has several endpoints
  # (only then does it carry one), else from the simulation below.
  # `stratify = "CMT"` with no `cmt` keeps every continuous endpoint, one
  # stratum each; the key carries CMT, so each still matches its own rows.
  if ("CMT" %in% names(obs) && !(is.null(cmt) && "CMT" %in% stratify)) {
    tte_cmt <- if (length(blocks$tte) > 0L) sdtab$CMT[blocks$tte]
    cmt <- .ferx_vpc_pick_cmt(cmt, obs$CMT, tte_cmt)
    obs <- obs[obs$CMT == cmt, , drop = FALSE]
  }
  cens <- .ferx_vpc_censor(obs, lloq, uloq)
  obs <- cens$obs
  if (pred_corr && (!is.null(cens$lloq) || !is.null(cens$uloq))) {
    .ferx_vpc_refuse_pc_limit(cens$lloq, cens$uloq, lloq, uloq)
  }

  sim <- .ferx_vpc_simulate(fit, n_sim, seed, horizon)
  if (!is.null(stratify)) {
    sim <- .ferx_vpc_add_strata(sim, .ferx_fit_data_path(fit), stratify)
  }
  sim <- sim[is.na(sim$OBSERVED), , drop = FALSE]
  if (!"CMT" %in% names(obs)) {
    sim_cmt <- unique(sim$CMT)
    if (length(sim_cmt) != 1L) {
      .ferx_vpc_mismatch(sprintf(
        "the simulation has continuous records in %d compartments (%s) where the fit has one endpoint",
        length(sim_cmt), .ferx_vpc_list(sim_cmt)
      ))
    }
    cmt <- .ferx_vpc_pick_cmt(cmt, sim_cmt, NULL)
    obs$CMT <- rep(cmt, nrow(obs))
  }
  # The key carries CMT, so the other endpoints' rows simply go unmatched.
  idx <- .ferx_vpc_match(obs, sim, n_sim)
  keep <- c("SIM", "ID", "TIME", "CMT", "DV_SIM", stratify[stratify != "CMT"])
  sim_out <- sim[idx, keep, drop = FALSE]
  sim_out$PRED <- rep(obs$PRED, n_sim)
  rownames(sim_out) <- NULL
  first <- idx[seq_len(nrow(obs))]
  for (s in setdiff(stratify, "CMT")) obs[[s]] <- sim[[s]][first]

  obs_out <- obs[, c("ID", "TIME", "DV", "PRED", "CMT",
                     stratify[stratify != "CMT"]), drop = FALSE]
  list(obs = obs_out, sim = sim_out, lloq = cens$lloq, uloq = cens$uloq)
}

# `cmt` checked against the continuous endpoints `cmts` (the TTE ones,
# `tte_cmt`, named in the refusal when `cmt` is one of them).
.ferx_vpc_pick_cmt <- function(cmt, cmts, tte_cmt) {
  cmts <- sort(unique(cmts))
  if (is.null(cmt)) {
    if (length(cmts) > 1L) {
      .ferx_vpc_stop(paste0(
        "This fit has ", length(cmts), " continuous endpoints (CMT ",
        .ferx_vpc_list(cmts), "); pass `cmt =` to choose the one to show."
      ))
    }
    return(cmts[[1L]])
  }
  if (!is.numeric(cmt) || length(cmt) != 1L || !is.finite(cmt) ||
      !cmt %in% cmts) {
    shown <- if (is.atomic(cmt) && length(cmt) == 1L) format(cmt) else "this value"
    msg <- paste0(
      "`cmt = ", shown, "` is not a continuous endpoint of this fit; ",
      if (length(cmts) == 1L) "its continuous endpoint is CMT " else
        "its continuous endpoints are CMT ",
      .ferx_vpc_list(cmts), "."
    )
    if (is.numeric(cmt) && length(cmt) == 1L && cmt %in% tte_cmt) {
      msg <- paste0(msg, " CMT ", format(cmt), " is a time-to-event endpoint: ",
                    "use ferx_vpc_tte() for its VPC.")
    }
    .ferx_vpc_stop(msg)
  }
  cmt
}

.ferx_vpc_simulate <- function(fit, n_sim, seed, horizon) {
  tryCatch(
    ferx_simulate(fit$model_path, .ferx_fit_data_path(fit), n_sim = n_sim,
                  seed = seed, fit = fit, horizon = horizon),
    error = function(e) {
      if (is.null(horizon) && identical(e$code, "E_SIMULATE_HORIZON")) {
        .ferx_vpc_stop(paste0(
          "This model's simulation needs an administrative horizon (it has a ",
          "drug-driven time-to-event endpoint): pass `horizon =`, the end of ",
          "follow-up, to ferx_vpc(). The VPC still shows the continuous ",
          "endpoint only."
        ), code = "E_SIMULATE_HORIZON")
      }
      stop(e)
    }
  )
}

# Key of each row within its group: compartment, subject, time and the row's
# occurrence among rows sharing those three (in table order). The ID is the
# only text field and sits between numbers, so a key reads back one way.
# Rows sharing (CMT, ID, TIME) are exchangeable draws, so pairing them by
# occurrence is right in distribution even where a selection dropped one.
.ferx_vpc_key <- function(cmt, id, time) {
  k <- paste(cmt, as.character(id), sprintf("%.17g", as.numeric(time)),
             sep = "|")
  occ <- stats::ave(seq_along(k), k, FUN = seq_along)
  paste(k, occ, sep = "#")
}

# For each replicate, in `obs` order, the row of `sim` that simulates each
# observation. Refuses unless every replicate matches every observation once.
#
# Every replicate simulates the same records in the same order, so only the
# first is keyed. A replicate with the same CMT, ID and TIME at every position
# takes its pairing; any other is keyed on its own. Keying every row instead
# cost ~5 s at 200 x 10,000 observations, more than the simulation (#796).
.ferx_vpc_match <- function(obs, sim, n_sim) {
  reps <- sort(unique(sim$SIM))
  if (length(reps) != n_sim) {
    .ferx_vpc_mismatch(sprintf(
      "the simulation returned %d replicates where %d were asked for",
      length(reps), n_sim
    ))
  }
  n_obs <- nrow(obs)
  k_obs <- .ferx_vpc_key(obs$CMT, obs$ID, obs$TIME)
  rows <- unname(split(seq_len(nrow(sim)), factor(sim$SIM, levels = reps)))
  pair <- function(j) {
    r <- rows[[j]]
    m <- match(k_obs, .ferx_vpc_key(sim$CMT[r], sim$ID[r], sim$TIME[r]))
    if (anyNA(m)) {
      .ferx_vpc_mismatch(sprintf(
        "replicate %s matches %d of the fit's %d observations",
        format(reps[[j]]), sum(!is.na(m)), n_obs
      ))
    }
    m
  }
  m1 <- pair(1L)
  first <- rows[[1L]]
  same <- lengths(rows) == length(first)
  if (any(same)) {
    pos <- unlist(rows[same], use.names = FALSE)
    ref <- rep(first, sum(same))
    eq <- (sim$CMT[pos] == sim$CMT[ref] & sim$ID[pos] == sim$ID[ref] &
             sim$TIME[pos] == sim$TIME[ref]) %in% TRUE
    same[same] <- colSums(matrix(!eq, nrow = length(first))) == 0L
  }
  idx <- lapply(seq_along(rows), function(j) {
    rows[[j]][if (same[[j]]) m1 else pair(j)]
  })
  unlist(idx, use.names = FALSE)
}

.ferx_vpc_mismatch <- function(what) {
  .ferx_vpc_stop(paste0(
    "ferx_vpc() could not pair the simulation with the fit's observations by ",
    "(CMT, ID, TIME, occurrence): ", what, ". This is a ferx bug, not a ",
    "problem with your data; please report it at ",
    "https://github.com/FeRx-NLME/ferx-r/issues."
  ))
}

# vpc 1.2.4 sets the simulated values past a limit to NA when it
# prediction-corrects and then takes their quantiles without na.rm, so a
# pcVPC with a limit fails inside vpc whenever a replicate crosses it.
.ferx_vpc_refuse_pc_limit <- function(lloq, uloq, lloq_given, uloq_given) {
  lim <- c(if (!is.null(lloq)) sprintf("lloq = %s", format(lloq)),
           if (!is.null(uloq)) sprintf("uloq = %s", format(uloq)))
  from <- if (is.null(lloq_given) && is.null(uloq_given)) {
    " (from the fit's censored observations, CENS = 1 / -1)"
  } else {
    ""
  }
  .ferx_vpc_stop(paste0(
    "A prediction-corrected VPC cannot use a quantification limit (",
    paste(lim, collapse = ", "), from, "): vpc::vpc() censors the ",
    "simulated values at it and then fails on them. Use pred_corr = FALSE ",
    "for a VPC of this fit; the censored observations stay blanked."
  ))
}

# Blank the censored observations and, when neither limit is given, take each
# from them. vpc::vpc() takes one limit, so two are refused: a limit you pass
# replaces both defaults.
.ferx_vpc_censor <- function(obs, lloq, uloq) {
  if (!is.null(lloq) && !is.null(uloq)) {
    .ferx_vpc_stop(paste0(
      "vpc::vpc() takes one quantification limit; pass `lloq =` or `uloq =` ",
      "to ferx_vpc(), not both."
    ))
  }
  if (!"CENS" %in% names(obs)) {
    return(list(obs = obs, lloq = lloq, uloq = uloq))
  }
  cens <- suppressWarnings(as.numeric(obs$CENS))
  below <- cens %in% 1
  above <- cens %in% -1
  given <- !is.null(lloq) || !is.null(uloq)
  pick <- function(rows, name, flag) {
    if (!any(rows)) return(NULL)
    v <- unique(obs$DV[rows])
    if (length(v) > 1L) {
      .ferx_vpc_stop(paste0(
        "The fit's CENS = ", flag, " observations carry ", length(v),
        " different ", toupper(name), " values (", .ferx_vpc_list(v),
        "), and vpc::vpc() takes a single `", name, "`. Pass `", name,
        " =` to ferx_vpc()."
      ))
    }
    v
  }
  if (!given) {
    lloq <- pick(below, "lloq", "1")
    uloq <- pick(above, "uloq", "-1")
    if (!is.null(lloq) && !is.null(uloq)) {
      .ferx_vpc_stop(paste0(
        "This fit has observations censored on both sides (lloq = ",
        .ferx_vpc_list(lloq), " from CENS = 1, uloq = ", .ferx_vpc_list(uloq),
        " from CENS = -1), and vpc::vpc() takes one quantification limit. ",
        "Pass `lloq =` or `uloq =` to ferx_vpc() to choose the one it draws; ",
        "the censored observations on both sides stay blanked."
      ))
    }
  }
  obs$DV[below | above] <- NA_real_
  list(obs = obs, lloq = lloq, uloq = uloq)
}

# The `stratify` columns of the fit's data, attached to every simulated row.
# The simulation has one row per observation record of the data (EVID = 0,
# MDV not 1: the engine does not simulate an EVID = 0, MDV = 1 record) per
# replicate, in data order; that is checked (ID, and TIME on the continuous
# rows) before a column is attached, and refused when it does not hold.
.ferx_vpc_add_strata <- function(sim, data_path, stratify) {
  data <- utils::read.csv(data_path, colClasses = "character",
                          check.names = FALSE, na.strings = c(".", ""))
  cols <- setdiff(stratify, "CMT")
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    .ferx_vpc_stop(paste0(
      "`stratify` names ", paste0("`", missing, "`", collapse = ", "),
      ", not a column of the fit's data. Its columns are: ",
      paste(names(data), collapse = ", "), " (or \"CMT\")."
    ))
  }
  reserved <- intersect(cols, c("ID", "TIME", "DV", "DV_SIM", "SIM", "PRED",
                                "IPRED", "OBSERVED", "DRAW"))
  if (length(reserved) > 0L) {
    .ferx_vpc_stop(paste0("`stratify` cannot use the column `", reserved[[1L]],
                          "`, which the VPC itself is built from."))
  }
  if (length(cols) == 0L) return(sim)
  num <- function(col) suppressWarnings(as.numeric(data[[col]]))
  rec <- rep(TRUE, nrow(data))
  if ("EVID" %in% names(data)) rec <- rec & num("EVID") %in% 0
  if ("MDV" %in% names(data)) rec <- rec & !num("MDV") %in% 1
  data <- data[rec, , drop = FALSE]
  n <- nrow(data)
  reps <- split(seq_len(nrow(sim)), sim$SIM)
  cont <- is.na(sim$OBSERVED)
  lined_up <- all(vapply(reps, function(i) {
    length(i) == n && all(sim$ID[i] == data$ID) &&
      all(!cont[i] | as.numeric(sim$TIME[i]) == as.numeric(data$TIME))
  }, logical(1)))
  if (!lined_up) {
    .ferx_vpc_stop(paste0(
      "ferx_vpc() could not line the fit's data records up with the ",
      "simulation to read `stratify` from them. Stratify on \"CMT\", or ",
      "call vpc::vpc() directly with your own columns."
    ))
  }
  pos <- unlist(lapply(reps, seq_along), use.names = FALSE)
  ord <- unlist(reps, use.names = FALSE)
  for (s in cols) {
    v <- utils::type.convert(data[[s]], as.is = TRUE)
    out <- v[rep(NA_integer_, nrow(sim))]
    out[ord] <- v[pos]
    sim[[s]] <- out
  }
  sim
}
