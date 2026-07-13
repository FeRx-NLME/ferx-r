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
#'     `IPRED`, `CWRES`, `IWRES`, and (when present) `NPDE`, `NPD`, `CMT`, `OCC`,
#'     `CENS`, `TAFD`, `TAD`. `RES = DV - PRED` and `IRES = DV - IPRED` are
#'     derived; `WRES` is set to `NA` (ferx does not compute the FO-weighted
#'     residual). `NPDE`/`NPD` (present when the fit ran with `npde_nsim > 0`)
#'     are mapped to the residual role, so Xpose residual plots can use them.
#'   \item \strong{patab}: individual parameter values (`fit$individual_estimates`)
#'     and empirical-Bayes etas (`fit$ebe_etas`), repeated for every observation
#'     row of each subject.
#'   \item \strong{cotab}/\strong{catab}: covariates, split into continuous and
#'     categorical using `fit$covariate_types` (override with the `continuous`
#'     and `categorical` arguments). Covariate values are taken from `fit$sdtab`
#'     when echoed there, otherwise carried forward (LOCF) from `fit$covtab`.
#' }
#'
#' @section Estimation-iteration trace:
#' When the fit was run with `optimizer_trace = TRUE` and `iterations = TRUE`
#' (the default), the per-iteration parameter and gradient trajectories are
#' populated into the `$files` slot of the returned `"xpose"` object as
#' synthetic NONMEM `.ext` / `.grd` tables, so
#' [xpose::prm_vs_iteration()] (parameter value vs iteration) and
#' [xpose::grd_vs_iteration()] (gradient vs iteration) work out-of-the-box. The
#' `.ext` table carries one column per parameter (using the fit's declared
#' parameter names, e.g. `TVCL`, `ETA_CL`) plus `OBJ`; the `.grd` table carries
#' one `GRD(n)` column per parameter and is only emitted for gradient-based
#' methods (a derivative-free trace has no gradient to plot). When no
#' per-parameter trace is present - the fit was run without
#' `optimizer_trace = TRUE`, or predates ferx recording per-parameter values -
#' the slot is left empty: the iteration plots then raise xpose's usual
#' "no files" message, while the table (goodness-of-fit / covariate) plots are
#' unaffected. For an OFV-over-iterations view use [plot.ferx_fit()]. Only the
#' `"xpose"` backend populates this slot; `"xpose4"` ignores `iterations`.
#'
#' @param fit A `ferx_fit` object returned by [ferx_fit()].
#' @param backend Which Xpose package to target: `"xpose"` (the modern
#'   tidyverse package, the default) or `"xpose4"` (the classic S4 package).
#' @param continuous,categorical Optional character vectors naming covariates to
#'   treat as continuous / categorical, overriding the `fit$covariate_types`
#'   classification. Names not present among the fit's covariates are ignored
#'   with a warning.
#' @param runno Run number recorded on the resulting Xpose object (cosmetic).
#' @param iterations Logical; when `TRUE` (the default) and the fit carries a
#'   per-parameter optimizer trace (from `optimizer_trace = TRUE`), populate the
#'   `$files` slot so [xpose::prm_vs_iteration()] and
#'   [xpose::grd_vs_iteration()] work. Ignored by the `"xpose4"` backend and a
#'   no-op when no trace is present.
#'
#' @return For `backend = "xpose"`, an `xpose_data` object (with a populated
#'   `$files` slot when an optimizer trace is present, see the
#'   "Estimation-iteration trace" section). For `backend = "xpose4"`, an
#'   `xpose.data` (S4) object.
#'
#' @examples
#' \dontrun{
#' fit  <- ferx_fit("warfarin.ferx", data = "warfarin.csv",
#'                  optimizer_trace = TRUE)
#' xpdb <- ferx_xpose(fit)
#' xpose::dv_vs_ipred(xpdb)
#' xpose::prm_vs_iteration(xpdb)   # parameter trajectories (needs the trace)
#' xpose::grd_vs_iteration(xpdb)   # gradient trajectories (gradient methods)
#'
#' xpdb4 <- ferx_xpose(fit, backend = "xpose4")
#' xpose4::basic.gof(xpdb4)
#' }
#' @export
ferx_xpose <- function(fit,
                       backend = c("xpose", "xpose4"),
                       continuous = NULL,
                       categorical = NULL,
                       runno = 1L,
                       iterations = TRUE) {
  if (!inherits(fit, "ferx_fit")) {
    stop("`fit` must be a ferx_fit object (from ferx_fit()).", call. = FALSE)
  }
  if (!is.logical(iterations) || length(iterations) != 1L || is.na(iterations)) {
    stop("`iterations` must be TRUE or FALSE.", call. = FALSE)
  }
  backend <- match.arg(backend)

  tbl <- .ferx_xpose_frame(fit, continuous = continuous, categorical = categorical)

  switch(
    backend,
    xpose  = .ferx_xpose_new(tbl, fit, runno = runno, iterations = iterations),
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

  # Warn at build time (not deep inside a plot call) when an sdtab is missing a
  # column the standard GOF templates need.
  for (need in c("PRED", "IPRED", "CWRES", "IWRES")) {
    if (!need %in% names(df)) {
      warning(sprintf(
        "sdtab has no %s column; Xpose plots that use it will not work.", need),
        call. = FALSE)
    }
  }

  nm <- names(df)
  present <- function(x) x[x %in% nm]
  roles <- list(
    id      = "ID",
    idv     = "TIME",
    dv      = present("DV"),
    pred    = present(c("PRED", "IPRED")),
    res     = present(c("RES", "IRES", "WRES", "CWRES", "IWRES", "NPDE", "NPD")),
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
  if (length(add_cols) == 0L) {
    return(list(data = df, added = character(0)))
  }
  # A param/eta may already be echoed into sdtab (per-observation grain, e.g. a
  # covariate-driven CL via [output]). The per-subject value is the canonical
  # one for the param/eta role, so prefer it (overwrite) and flag the collision.
  collide <- intersect(add_cols, names(df))
  if (length(collide) > 0L) {
    warning(sprintf(
      "`fit$%s` column(s) also present in sdtab; using the per-subject value: %s",
      what, paste(collide, collapse = ", ")), call. = FALSE)
  }
  # sdtab IDs are numeric (Rust writes them as f64) while ebe_etas /
  # individual_estimates IDs are character; match in character space so the
  # join works for non-numeric or zero-padded IDs (cf. diagnostics.R).
  idx <- match(as.character(df[["ID"]]), as.character(src[["ID"]]))
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
    # Normalise the type labels; an NA or unexpected value must not silently
    # drop a covariate from both sets - warn and treat it as unclassified.
    tt  <- tolower(trimws(as.character(types)))
    cont <- nms[!is.na(tt) & tt == "continuous"]
    cat  <- nms[!is.na(tt) & tt == "categorical"]
    bad_type <- nms[is.na(tt) | !(tt %in% c("continuous", "categorical"))]
    bad_type <- bad_type[!is.na(bad_type)]
    if (length(bad_type) > 0L) {
      warning(sprintf(
        "covariate(s) with unrecognized type (not continuous/categorical), omitted: %s",
        paste(bad_type, collapse = ", ")), call. = FALSE)
    }
  }

  # Empty character vector (not NULL) when the fit declares no covariates, so an
  # override that names an undeclared covariate is always warned about and
  # dropped - not silently leaked into the tables.
  known <- names(fit$covariate_types) %||% character(0)
  warn_unknown <- function(x, kind) {
    if (is.null(x)) return(invisible())
    bad <- setdiff(x, known)
    if (length(bad) > 0L) {
      warning(sprintf("%s covariate(s) not declared in the model, ignored: %s",
                      kind, paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  warn_unknown(continuous, "continuous")
  warn_unknown(categorical, "categorical")

  # Keep only overrides that name a declared covariate.
  continuous  <- intersect(continuous, known)
  categorical <- intersect(categorical, known)

  if (length(continuous) > 0L) {
    cont <- union(setdiff(cont, categorical), continuous)
    cat  <- setdiff(cat, continuous)
  }
  if (length(categorical) > 0L) {
    cat  <- union(setdiff(cat, continuous), categorical)
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
# When several source records tie at that latest time, the last one in source
# order wins (true LOCF). Observations before any source record for that ID fall
# back to the earliest record. Vectorised per ID group (no per-row scan).
.ferx_locf_join <- function(id, t, src_id, src_t, src_val) {
  id     <- as.character(id)       # sdtab IDs are numeric, covtab IDs character
  src_id <- as.character(src_id)
  out <- rep(NA_real_, length(id))
  src_by_id <- split(seq_along(src_id), src_id)
  for (g in unique(id)) {
    rows <- which(id == g)
    sidx <- src_by_id[[g]]
    if (is.null(sidx) || length(sidx) == 0L) next
    # Stable order by time, ties broken by original position so the last tied
    # record sorts to the right; findInterval then returns it for an equal time.
    ord <- sidx[order(src_t[sidx], sidx)]
    st  <- src_t[ord]
    j   <- findInterval(t[rows], st)   # largest index with st <= t; 0 if before all
    j[j == 0L] <- 1L                   # fallback to earliest record for that ID
    out[rows] <- src_val[ord][j]
  }
  out
}

# --- backend: modern xpose ---------------------------------------------------

.ferx_xpose_new <- function(tbl, fit, runno = 1L, iterations = TRUE) {
  if (!requireNamespace("xpose", quietly = TRUE)) {
    stop("Package 'xpose' is required for backend = \"xpose\". ",
         "Install it or use backend = \"xpose4\".", call. = FALSE)
  }
  df <- tbl$data

  # Estimation-iteration trace -> synthetic .ext / .grd tables (NULL when the
  # fit carries no per-parameter trace, leaving the files slot empty so the
  # iteration plots raise xpose's own "no files" message).
  files <- if (isTRUE(iterations)) .ferx_xpose_iteration_files(fit) else NULL

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

  summary <- .ferx_xpose_new_summary(fit, runno, files = files)

  structure(
    list(
      code     = NULL,
      summary  = summary,
      data     = data_tbl,
      files    = files,
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

.ferx_xpose_new_summary <- function(fit, runno, files = NULL) {
  rows <- function(label, value, ...) {
    data.frame(problem = 1L, subprob = 0L, descr = label,
               label = label, value = as.character(value),
               stringsAsFactors = FALSE)
  }
  out <- list(
    rows("software", "nonmem"),
    # xpose's title-template engine looks up the key "run" (not "runno").
    rows("run", runno),
    rows("ofv", if (!is.null(fit$ofv)) fit$ofv else NA),
    rows("nobs", nrow(fit$sdtab)),
    rows("nind", length(unique(fit$sdtab$ID)))
  )
  # The iteration plots' title/subtitle/caption templates reference the
  # `method`, `runtime`, `term` and `dir` keys; supply them (only when a trace
  # is present) so prm_vs_iteration() / grd_vs_iteration() render without
  # xpose's "not part of the available keywords" warnings.
  if (!is.null(files)) {
    method <- files$method[files$extension == "ext"][1]
    out <- c(out, list(
      rows("method", method %||% NA),
      rows("runtime", .ferx_trace_runtime_string(fit)),
      rows("term", if (isTRUE(fit$converged))
        "OPTIMIZATION COMPLETED" else "OPTIMIZATION NOT COMPLETED"),
      rows("dir", "")
    ))
  }
  do.call(rbind, out)
}

# --- estimation-iteration trace -> synthetic .ext / .grd file tables ---------

# Build the xpose `$files` tibble (an `.ext` table and, for gradient-based
# methods, a `.grd` table) from a fit's per-parameter optimizer trace. Returns
# NULL - so the caller leaves the files slot empty - when there is no trace, or
# the trace predates ferx recording per-parameter values (only scalar
# ofv/grad_norm, no `val:*` columns).
.ferx_xpose_iteration_files <- function(fit) {
  tr <- tryCatch(ferx_trace(fit), error = function(e) NULL)
  if (is.null(tr) || !is.data.frame(tr) || nrow(tr) == 0L) return(NULL)

  val_cols <- grep("^val:", names(tr), value = TRUE)
  if (length(val_cols) == 0L) return(NULL)

  iter   <- suppressWarnings(as.numeric(tr[["iter"]]))
  method <- .ferx_trace_method_label(tr)

  # .ext table: ITERATION, one column per parameter (declared name, e.g. TVCL /
  # ETA_CL), then OBJ. prm_vs_iteration() drops parameters non-varying across
  # ITERATION, tidies, and facets by variable.
  ext <- data.frame(ITERATION = iter, check.names = FALSE,
                    stringsAsFactors = FALSE)
  for (col in val_cols) {
    ext[[sub("^val:", "", col)]] <- suppressWarnings(as.numeric(tr[[col]]))
  }
  ext[["OBJ"]] <- suppressWarnings(as.numeric(tr[["ofv"]]))

  data_list <- list(ext)
  name_v    <- "ferx.ext"
  ext_v     <- "ext"

  # .grd table: ITERATION + one GRD(n) column per parameter. Emitted only when
  # at least one gradient is finite - a derivative-free trace (e.g. SAEM, or a
  # BOBYQA eval) writes every grad:* column NA and has no gradient to plot.
  # grd_vs_iteration() keys the facet on the digits in each column name, so the
  # NONMEM-style GRD(n) naming (not the parameter name) is required.
  grad_cols <- grep("^grad:", names(tr), value = TRUE)
  grad_num  <- lapply(grad_cols, function(col) suppressWarnings(as.numeric(tr[[col]])))
  has_grad  <- length(grad_cols) > 0L &&
    any(vapply(grad_num, function(x) any(is.finite(x)), logical(1)))
  if (has_grad) {
    grd <- data.frame(ITERATION = iter, check.names = FALSE,
                      stringsAsFactors = FALSE)
    for (i in seq_along(grad_cols)) {
      grd[[sprintf("GRD(%d)", i)]] <- grad_num[[i]]
    }
    data_list <- c(data_list, list(grd))
    name_v    <- c(name_v, "ferx.grd")
    ext_v     <- c(ext_v, "grd")
  }

  files <- data.frame(name = name_v, extension = ext_v, problem = 1L,
                      subprob = 1L, method = method %||% NA_character_,
                      modified = FALSE, stringsAsFactors = FALSE)
  # `data` is a list-column holding one table per row; assign after the frame
  # is built so it is not recycled column-wise.
  files$data <- data_list
  files[c("name", "extension", "problem", "subprob", "method", "data",
          "modified")]
}

# The estimation-method label recorded in the files tibble and the plot
# subtitle: the last non-empty method the trace visited (e.g. "focei" for the
# polish phase of a gn_hybrid fit). NA when the trace has no method column.
.ferx_trace_method_label <- function(tr) {
  m <- tr[["method"]]
  if (is.null(m)) return(NA_character_)
  m <- m[!is.na(m) & nzchar(as.character(m))]
  if (length(m) == 0L) return(NA_character_)
  as.character(m[length(m)])
}

# Format the fit's wall-clock runtime (max `wall_ms` in the trace) as
# "HH:MM:SS" for the iteration-plot subtitle. "na" when unavailable, mirroring
# xpose's own placeholder. Reads via ferx_trace() (not fit$trace directly) so a
# fit carrying only a trace_path - not an in-memory trace - is still covered.
.ferx_trace_runtime_string <- function(fit) {
  tr <- tryCatch(ferx_trace(fit), error = function(e) NULL)
  ms <- if (is.data.frame(tr)) suppressWarnings(as.numeric(tr[["wall_ms"]])) else NULL
  ms <- ms[is.finite(ms)]
  if (length(ms) == 0L) return("na")
  secs <- as.integer(max(ms) %/% 1000)
  sprintf("%02d:%02d:%02d", secs %/% 3600L, (secs %% 3600L) %/% 60L, secs %% 60L)
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
