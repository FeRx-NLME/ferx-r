# Xpose interoperability ------------------------------------------------------
#
# `ferx_xpose()` turns a fitted `ferx_fit` into a ready-to-use Xpose object,
# in memory, with no NONMEM table files written to disk. Two backends are
# supported:
#
#   * "xpose"  - the modern tidyverse package (CRAN `xpose`), an `xpose_data`
#                list with a nested `data` tibble + per-column `index`.
#   * "xpose4" - the classic S4 package (`xpose4`), an `xpose.data` object.
#
# Both backends are fed from one observation-grain data frame plus a role map
# (which column is DV, IPRED, a covariate, an eta, ...). The role map is built
# from the `fit` object alone; see `.ferx_xpose_frame()`. That builder is pure
# and has no Xpose dependency, so it is unit-tested directly. The thin backend
# constructors are guarded by `requireNamespace()` and only run when the user
# has the relevant package installed.

#' Build an Xpose object from a ferx fit
#'
#' Produces an in-memory Xpose object so that all downstream Xpose diagnostics
#' (goodness-of-fit, covariate, and parameter plots) work out-of-the-box on a
#' `ferx` fit, without writing NONMEM-style table files to disk.
#'
#' The mapping from `fit` to the standard NONMEM table columns is:
#' \itemize{
#'   \item \strong{sdtab} from `fit$sdtab`: `ID`, `TIME`, `DV`, `PRED`,
#'     `IPRED`, `CWRES`, `IWRES`, and (when present) `CMT`, `OCC`, `CENS`,
#'     `TAFD`, `TAD`. `RES = DV - PRED` and `IRES = DV - IPRED` are derived;
#'     `WRES` is set to `NA` (ferx does not compute the FO-weighted residual).
#'   \item \strong{patab}: individual parameter values (`fit$individual_estimates`)
#'     and empirical-Bayes etas (`fit$ebe_etas`), repeated for every observation
#'     row of each subject.
#'   \item \strong{cotab}/\strong{catab}: covariates, split into continuous and
#'     categorical using `fit$covariate_types` (override with the `continuous`
#'     and `categorical` arguments). Covariate values are taken from `fit$sdtab`
#'     when echoed there, otherwise carried forward (LOCF) from `fit$covtab`.
#' }
#'
#' @section Limitations:
#' Only the table data is populated, not the estimation-iteration trace. Xpose
#' functions that read the NONMEM `.ext` / `.grd` files — notably
#' [xpose::prm_vs_iteration()] and [xpose::grd_vs_iteration()] — therefore do
#' not work on the returned object: ferx records only a scalar OFV and gradient
#' norm per iteration, not the per-parameter value/gradient trajectory those
#' plots need. For an OFV-over-iterations view use [ferx_plot_trace()] (requires
#' `optimizer_trace = TRUE` at fit time).
#'
#' @param fit A `ferx_fit` object returned by [ferx_fit()].
#' @param backend Which Xpose package to target: `"xpose"` (the modern
#'   tidyverse package, the default) or `"xpose4"` (the classic S4 package).
#' @param continuous,categorical Optional character vectors naming covariates to
#'   treat as continuous / categorical, overriding the `fit$covariate_types`
#'   classification. Names not present among the fit's covariates are ignored
#'   with a warning.
#' @param runno Run number recorded on the resulting Xpose object (cosmetic).
#'
#' @return For `backend = "xpose"`, an `xpose_data` object. For
#'   `backend = "xpose4"`, an `xpose.data` (S4) object.
#'
#' @examples
#' \dontrun{
#' fit  <- ferx_fit("warfarin.ferx", data = "warfarin.csv")
#' xpdb <- ferx_xpose(fit)
#' xpose::dv_vs_ipred(xpdb)
#'
#' xpdb4 <- ferx_xpose(fit, backend = "xpose4")
#' xpose4::basic.gof(xpdb4)
#' }
#' @export
ferx_xpose <- function(fit,
                       backend = c("xpose", "xpose4"),
                       continuous = NULL,
                       categorical = NULL,
                       runno = 1L) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit()).", call. = FALSE)
  }
  backend <- match.arg(backend)

  tbl <- .ferx_xpose_frame(fit, continuous = continuous, categorical = categorical)

  switch(
    backend,
    xpose  = .ferx_xpose_new(tbl, fit, runno = runno),
    xpose4 = .ferx_xpose_xpose4(tbl, fit, runno = runno)
  )
}

# --- table + role builder (pure; no Xpose dependency) ------------------------

# Returns a list:
#   data   : observation-grain data.frame with NONMEM-standard column names
#   roles  : named list mapping each role -> column name(s)
#   cont   : continuous covariate column names
#   cat    : categorical covariate column names
#   params : individual-parameter column names
#   etas   : eta column names
.ferx_xpose_frame <- function(fit, continuous = NULL, categorical = NULL) {
  sdtab <- fit$sdtab
  if (is.null(sdtab) || !is.data.frame(sdtab) || nrow(sdtab) == 0L) {
    stop("`fit$sdtab` is empty; cannot build an Xpose object.", call. = FALSE)
  }
  if (!all(c("ID", "TIME") %in% names(sdtab))) {
    stop("`fit$sdtab` must contain ID and TIME columns.", call. = FALSE)
  }

  df <- sdtab

  # Derived residuals. RES/IRES are cheap; WRES is unavailable in ferx (FO
  # weighting is not computed) and is emitted as NA so Xpose still finds the
  # column its default templates reference.
  if (all(c("DV", "PRED") %in% names(df)) && is.null(df[["RES"]])) {
    df[["RES"]] <- df[["DV"]] - df[["PRED"]]
  }
  if (all(c("DV", "IPRED") %in% names(df)) && is.null(df[["IRES"]])) {
    df[["IRES"]] <- df[["DV"]] - df[["IPRED"]]
  }
  if (is.null(df[["WRES"]])) {
    df[["WRES"]] <- NA_real_
  }

  # Parameters + etas, one value per subject, repeated down the observation
  # rows by an ID join.
  params <- .ferx_join_by_id(df, fit$individual_estimates, "individual_estimates")
  df <- params$data
  param_cols <- params$added

  etas <- .ferx_join_by_id(df, fit$ebe_etas, "ebe_etas")
  df <- etas$data
  eta_cols <- etas$added

  # Covariate classification.
  cov_split <- .ferx_xpose_covariates(fit, continuous, categorical)

  # Covariate values at observation grain: prefer columns already echoed into
  # sdtab (via [output]); otherwise LOCF-carry them from covtab.
  cov_added <- .ferx_attach_covariates(df, fit, cov_split$all)
  df <- cov_added$data
  cont_cols <- intersect(cov_split$continuous, names(df))
  cat_cols  <- intersect(cov_split$categorical, names(df))

  present <- function(x) x[x %in% names(df)]
  roles <- list(
    id      = "ID",
    idv     = "TIME",
    dv      = present("DV"),
    pred    = present(c("PRED", "IPRED")),
    res     = present(c("RES", "IRES", "WRES", "CWRES", "IWRES")),
    occ     = present("OCC"),
    param   = param_cols,
    eta     = eta_cols,
    contcov = cont_cols,
    catcov  = cat_cols
  )

  list(
    data   = df,
    roles  = roles,
    cont   = cont_cols,
    cat    = cat_cols,
    params = param_cols,
    etas   = eta_cols
  )
}

# Repeat one-row-per-subject columns down `df` by an ID join. `src` is a data
# frame whose first column is ID (e.g. fit$ebe_etas). Returns updated data and
# the names of the columns that were added.
.ferx_join_by_id <- function(df, src, what) {
  if (is.null(src) || !is.data.frame(src) || nrow(src) == 0L) {
    return(list(data = df, added = character(0)))
  }
  if (!"ID" %in% names(src)) {
    warning(sprintf("`fit$%s` has no ID column; skipping.", what), call. = FALSE)
    return(list(data = df, added = character(0)))
  }
  add_cols <- setdiff(names(src), "ID")
  # Avoid clobbering an existing column (e.g. a param already echoed in sdtab).
  add_cols <- setdiff(add_cols, names(df))
  if (length(add_cols) == 0L) {
    return(list(data = df, added = character(0)))
  }
  idx <- match(df[["ID"]], src[["ID"]])
  for (col in add_cols) {
    df[[col]] <- src[[col]][idx]
  }
  list(data = df, added = add_cols)
}

# Resolve which covariates are continuous vs categorical. Starts from
# `fit$covariate_types` (named "continuous"/"categorical") and applies the
# user's overrides.
.ferx_xpose_covariates <- function(fit, continuous, categorical) {
  types <- fit$covariate_types
  cont <- character(0)
  cat  <- character(0)
  if (!is.null(types) && length(types) > 0L) {
    nms <- names(types)
    cont <- nms[types == "continuous"]
    cat  <- nms[types == "categorical"]
  }

  known <- names(fit$covariate_types)
  warn_unknown <- function(x, kind) {
    if (is.null(x)) return(invisible())
    bad <- setdiff(x, known)
    if (length(bad) > 0L && !is.null(known)) {
      warning(sprintf("%s covariate(s) not declared in the model, ignored: %s",
                      kind, paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  warn_unknown(continuous, "continuous")
  warn_unknown(categorical, "categorical")

  # Keep only overrides that name a declared covariate; unknown names were
  # warned about above and are dropped here so they don't leak into the tables.
  if (!is.null(known)) {
    continuous  <- intersect(continuous, known)
    categorical <- intersect(categorical, known)
  }

  if (!is.null(continuous)) {
    cont <- union(setdiff(cont, categorical %||% character(0)), continuous)
    cat  <- setdiff(cat, continuous)
  }
  if (!is.null(categorical)) {
    cat  <- union(setdiff(cat, continuous %||% character(0)), categorical)
    cont <- setdiff(cont, categorical)
  }

  list(continuous = unique(cont), categorical = unique(cat),
       all = unique(c(cont, cat)))
}

# Make covariate columns available at observation grain. If a covariate is
# already a column of `df` (echoed into sdtab via [output]) it is left as-is.
# Otherwise it is carried forward (LOCF) from `fit$covtab`, which has one row
# per dataset record (doses + observations).
.ferx_attach_covariates <- function(df, fit, cov_names) {
  if (length(cov_names) == 0L) {
    return(list(data = df, added = character(0)))
  }
  missing <- setdiff(cov_names, names(df))
  covtab <- fit$covtab
  if (length(missing) > 0L) {
    if (is.null(covtab) || !is.data.frame(covtab) ||
        !all(c("ID", "TIME") %in% names(covtab))) {
      warning(sprintf("covariate(s) not in sdtab and no usable covtab, dropped: %s",
                      paste(missing, collapse = ", ")), call. = FALSE)
    } else {
      for (cov in missing) {
        if (!cov %in% names(covtab)) {
          warning(sprintf("covariate '%s' not found in covtab, dropped.", cov),
                  call. = FALSE)
          next
        }
        df[[cov]] <- .ferx_locf_join(df[["ID"]], df[["TIME"]],
                                     covtab[["ID"]], covtab[["TIME"]], covtab[[cov]])
      }
    }
  }
  list(data = df, added = intersect(cov_names, names(df)))
}

# Last-observation-carried-forward join: for each (id, t) in the target, pick
# the source value at the latest source record with the same ID and time <= t.
.ferx_locf_join <- function(id, t, src_id, src_t, src_val) {
  out <- rep(NA_real_, length(id))
  for (i in seq_along(id)) {
    sel <- which(src_id == id[i] & src_t <= t[i])
    if (length(sel) > 0L) {
      out[i] <- src_val[[sel[which.max(src_t[sel])]]]
    } else {
      # No prior record: fall back to the earliest record for that ID.
      first <- which(src_id == id[i])
      if (length(first) > 0L) out[i] <- src_val[[first[which.min(src_t[first])]]]
    }
  }
  out
}

# --- backend: modern xpose ---------------------------------------------------

.ferx_xpose_new <- function(tbl, fit, runno = 1L) {
  if (!requireNamespace("xpose", quietly = TRUE)) {
    stop("Package 'xpose' is required for backend = \"xpose\". ",
         "Install it or use backend = \"xpose4\".", call. = FALSE)
  }
  df <- tbl$data

  # Per-column role index expected by xpose. Column `type` drives which plot
  # templates pick up which variable.
  type_of <- .ferx_xpose_new_types(tbl)
  index <- data.frame(
    table = "sdtab",
    col   = names(df),
    type  = type_of,
    label = NA_character_,
    units = NA_character_,
    stringsAsFactors = FALSE
  )

  data_tbl <- data.frame(
    problem  = 1L,
    simtab   = FALSE,
    stringsAsFactors = FALSE
  )
  data_tbl$index    <- list(index)
  data_tbl$data     <- list(df)
  data_tbl$modified <- FALSE

  summary <- .ferx_xpose_new_summary(fit, runno)

  structure(
    list(
      code     = NULL,
      summary  = summary,
      data     = data_tbl,
      files    = NULL,
      gg_theme = xpose::theme_readable(),
      xp_theme = xpose::theme_xp_default(),
      options  = list(dir = NULL, quiet = TRUE, manual_import = NULL),
      software = "nonmem"
    ),
    class = c("xpose_data", "uneval")
  )
}

# Map ferx/NONMEM column names to xpose column `type` strings.
.ferx_xpose_new_types <- function(tbl) {
  roles <- tbl$roles
  cols  <- names(tbl$data)
  type  <- rep("na", length(cols))
  set <- function(role, value) type[cols %in% roles[[role]]] <<- value
  set("id", "id")
  set("idv", "idv")
  set("dv", "dv")
  set("occ", "occ")
  set("pred", "pred")
  set("res", "res")
  set("param", "param")
  set("eta", "eta")
  set("contcov", "contcov")
  set("catcov", "catcov")
  # xpose treats individual predictions as their own column type (`ipred`),
  # distinct from population `pred`; templates like dv_vs_ipred() require it.
  type[cols == "IPRED"] <- "ipred"
  type
}

.ferx_xpose_new_summary <- function(fit, runno) {
  rows <- function(label, value, ...) {
    data.frame(problem = 1L, subprob = 0L, descr = label,
               label = label, value = as.character(value),
               stringsAsFactors = FALSE)
  }
  do.call(rbind, list(
    rows("software", "nonmem"),
    rows("runno", runno),
    rows("ofv", if (!is.null(fit$ofv)) fit$ofv else NA),
    rows("nobs", nrow(fit$sdtab)),
    rows("nind", length(unique(fit$sdtab$ID)))
  ))
}

# --- backend: classic xpose4 -------------------------------------------------

.ferx_xpose_xpose4 <- function(tbl, fit, runno = 1L) {
  if (!requireNamespace("xpose4", quietly = TRUE)) {
    stop("Package 'xpose4' is required for backend = \"xpose4\". ",
         "Install it or use backend = \"xpose\".", call. = FALSE)
  }
  df <- tbl$data

  # xpose4 detects categorical covariates partly by factor-ness; coerce the
  # declared categorical columns so they are not treated as continuous.
  for (col in tbl$cat) {
    if (col %in% names(df)) df[[col]] <- as.factor(df[[col]])
  }

  # Build an empty xpose.data shell (its prototype seeds default Prefs whose
  # variable definitions already expect the NONMEM-standard names we use), then
  # assign Data through the exported replacement method.
  xpdb <- methods::new("xpose.data", Runno = as.character(runno), Data = NULL)
  xpose4::Data(xpdb) <- df

  # Point the variable definitions at the covariate / parameter columns so the
  # covariate- and parameter-plot families discover them.
  vardef <- xpdb@Prefs@Xvardef
  vardef$covariates <- c(tbl$cont, tbl$cat)
  vardef$parms      <- c(tbl$params, tbl$etas)
  if ("TAD" %in% names(df)) vardef$tad <- "TAD"
  xpdb@Prefs@Xvardef <- vardef

  xpdb
}
