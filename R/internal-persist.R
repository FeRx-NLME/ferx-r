# Reading the engine's optional values back on a fit: `.fitrx_unwrap_opt_*`.
# ferx_load_fit() reads bundles with ferx-core's loader.

# A scalar that may be `null`: the JSON writers turn NaN into `null`, and a VI
# fit under `vi_final_ofv = none` has a NaN OFV by design (#508), which
# `as.numeric(NULL)` would read as `numeric(0)`.
.fitrx_num_or_nan <- function(x) {
  v <- as.numeric(unlist(x, use.names = FALSE))
  if (length(v) == 0L) NaN else v[[1L]]
}

# Read back a JSON array that may contain nulls - a Rust `Vec<Option<T>>` such
# as the per-kappa `weight =` source text and its typical value. `unlist()`
# would silently *drop* the nulls and shift every later element onto the wrong
# kappa, so map them to NA element-wise instead. NULL when absent or empty.
.fitrx_unwrap_nullable_vec <- function(x, na) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  if (!is.list(x)) return(c(x, na)[seq_along(x)])
  vapply(x, function(el) {
    if (is.null(el) || length(el) == 0L) na else c(el, na)[[1L]]
  }, na)
}

.fitrx_unwrap_opt_chr_vec <- function(x) .fitrx_unwrap_nullable_vec(x, NA_character_)

.fitrx_unwrap_opt_lgl <- function(x) {
  if (is.null(x)) return(NA)
  v <- suppressWarnings(as.logical(unlist(x, use.names = FALSE)))
  if (length(v) != 1L || is.na(v)) return(NA)
  v
}

.fitrx_unwrap_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.character(x)
  if (length(v) == 0L || is.na(v) || !nzchar(v)) return(NULL)
  v
}

# Fields derived from other fields already present on `result` (correlation
# matrix, tidy estimates table, eta-covariate correlations). Shared between
# ferx_fit() and ferx_load_fit() - the two fit-construction entry points -
# so they can't drift on how a derived field is computed. Requires
# `result$data_path` to already be set (normalised, in ferx_fit()'s case).
.ferx_populate_derived_fields <- function(result) {
  result$cor_matrix <- .ferx_compute_cor_matrix(result$cov_matrix, result$cov_fixed)
  result$estimates  <- .ferx_compute_estimates(result)
  result$eta_cov    <- .ferx_compute_eta_cov(result$ebe_etas, result$data_path,
                                           result$model_path)
  result
}

# -- Parameter priors (ferx-core #254) ---------------------------------------
#
# `PriorSummary` on the wire is a plain array of objects with eight required
# fields (no serde defaults), so the write path emits every column and the read
# path tolerates a missing one only by filling NA. The pair is written *only*
# for a priored fit, matching ferx-core's own writer: absent means "this bundle
# predates priors", and the loader then reconstructs `ofv_data = ofv`, which is
# right for every such file because no prior could have been applied.

# -- Exact doubles in the bundle --
#
# A .fitrx round trip must give back the doubles it was handed. Two things
# stood in the way. jsonlite's `digits = NA` and `write.table()` /
# `write.csv()` all write a double with 15 significant digits, which is not
# enough to identify it. And R's own number parser (`as.numeric()`,
# `read.csv()`) is not correctly rounded: on a build whose long double is a
# plain double (aarch64) it reads up to four in five 17-digit strings an ulp
# or more away, depending on magnitude. So `fit.json` is written with 17 significant digits and read by
# jsonlite, whose parser is correctly rounded; the CSV entries get the
# shortest text that reads back as the same double, and are read back
# through jsonlite's parser too.

# Correctly rounded parse of decimal text, as `as.numeric()` would read it:
# "" and "NA" are NA, "NaN" is NaN, "Inf" / "-Inf" (and Rust's "inf" /
# "-inf") are infinite. Text that is not a plain decimal number falls back to
# `as.numeric()`, so an unexpected bundle still loads as it did before.
.fitrx_parse_doubles <- function(txt) {
  s <- trimws(as.character(txt))
  out <- rep(NA_real_, length(s))
  out[s %in% c("Inf", "+Inf", "inf", "+inf")] <- Inf
  out[s %in% c("-Inf", "-inf")] <- -Inf
  out[s %in% c("NaN", "nan")] <- NaN
  # JSON's number grammar exactly, so the batch below never fails on one odd
  # token; text outside it ("1.", "007") takes the per-element fallback.
  num <- !is.na(s) &
    grepl("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", s)
  if (any(num)) {
    out[num] <- tryCatch(
      as.numeric(jsonlite::parse_json(
        paste0("[", paste(s[num], collapse = ","), "]"),
        simplifyVector = TRUE
      )),
      error = function(e) suppressWarnings(as.numeric(s[num]))
    )
  }
  odd <- !is.na(s) & nzchar(s) & s != "NA" & !num &
    !(s %in% c("Inf", "+Inf", "inf", "+inf", "-Inf", "-inf", "NaN", "nan"))
  if (any(odd)) out[odd] <- suppressWarnings(as.numeric(s[odd]))
  out
}

# The shortest decimal text (15, 16 or 17 significant digits) that
# `.fitrx_parse_doubles()` reads back as exactly `x`. NA and NaN become NA
# (an empty cell, as `write.table(na = "")` wrote them), +/-Inf "Inf" /
# "-Inf". 17 digits always suffice for a double.
.fitrx_double_text <- function(x) {
  out <- rep(NA_character_, length(x))
  inf <- is.infinite(x)
  out[inf] <- ifelse(x[inf] > 0, "Inf", "-Inf")
  fin <- is.finite(x)
  xf <- x[fin]
  txt <- sprintf("%.15g", xf)
  for (d in 16:17) {
    bad <- .fitrx_parse_doubles(txt) != xf
    if (!any(bad)) break
    txt[bad] <- sprintf(paste0("%.", d, "g"), xf[bad])
  }
  out[fin] <- txt
  out
}

# `read.csv()`, with every column it types as double re-parsed by
# `.fitrx_parse_doubles()`. read.csv() still decides the column types, so a
# loaded table has the shape it always had - except the `id_cols`, which are
# read as text, verbatim: inferring a type turns the IDs `001` and `1.0` into
# `1`, and `as.character()` afterwards only puts the type back (#468). They
# come from a pass with no NA strings, since no writer puts an NA in an
# ID column and read.csv() reads `NA` as one even when quoted (#475).
#
# With `integer_cols` given, the table is one whose live columns are double
# but for these (and its text columns): any other column read.csv() types as
# integer - every value whole, as `1`, `100000` - or as logical - every value
# missing - is read as double too, so the table is identical() to the live
# fit's (#417). NULL leaves the typing to read.csv().
.fitrx_read_csv_exact <- function(path, id_cols = character(),
                                  integer_cols = NULL) {
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  id_cols <- intersect(id_cols, header)
  col_classes <- if (length(id_cols) > 0L) {
    stats::setNames(rep("character", length(id_cols)), id_cols)
  } else {
    NA
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses = col_classes)
  as_dbl <- if (is.null(integer_cols)) {
    is.double
  } else {
    function(v) is.numeric(v) || (is.logical(v) && all(is.na(v)))
  }
  dbl <- which(vapply(df, as_dbl, logical(1)) &
                 !(names(df) %in% c(id_cols, integer_cols)))
  if (length(dbl) > 0L || length(id_cols) > 0L) {
    raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = "character", na.strings = character())
    for (j in dbl) df[[j]] <- .fitrx_parse_doubles(raw[[j]])
    for (j in id_cols) df[[j]] <- raw[[j]]
  }
  df
}

# -- Warning details (ferx-core #2268) --
#
# `fit$warning_details` holds the structured warnings that carry a `details`
# payload, one list per entry with ferx-core's `WarningEntry` fields
# (`severity`, `category`, `message`, `source_method` when set, `details`):
# exactly the entries ferx-core writes to a `.fitrx` `fit.json:warning_details`,
# so an R bundle and an engine bundle carry them alike. Every warning's text
# stays in `fit$warnings`; `ferx_get_warnings(as_df = TRUE)` puts each entry's
# `details` beside its row. NULL when no warning carries details.

# Entries from the engine's JSON strings, parsed the one way both the live fit
# and a reloaded one take, so the two are `identical()`.
.ferx_parse_warning_details <- function(json) {
  json <- as.character(unlist(json %||% character(), use.names = FALSE))
  if (length(json) == 0L) return(NULL)
  lapply(json, jsonlite::parse_json)
}

# -- The sdtab's row blocks --

# The subject (1-based position in fit order) of every sdtab row, from the
# per-subject row counts: each subject's `n_obs` Gaussian rows in subject
# order, then each subject's `n_obs_discrete` discrete rows in subject order,
# then each subject's `n_obs_tte` time-to-event rows (ferx-core #901), as
# `io::output::sdtab` writes them. Not from the ID number, which runs two
# adjacent subjects with the same number together and has no rows for a
# subject without observations (#540). NULL when the counts do not add up to
# `n_rows`.
.fitrx_sdtab_subject <- function(n_obs, n_obs_discrete, n_rows, n_obs_tte = 0L) {
  n_obs <- as.integer(n_obs)
  n_obs_discrete <- as.integer(n_obs_discrete)
  n_obs_tte <- as.integer(n_obs_tte)
  # A scalar is the same count for every subject (0: no TTE rows); any other
  # length is not one count per subject.
  if (length(n_obs_tte) == 1L) n_obs_tte <- rep_len(n_obs_tte, length(n_obs))
  if (anyNA(n_obs) || anyNA(n_obs_discrete) || anyNA(n_obs_tte) ||
      length(n_obs_tte) != length(n_obs) ||
      sum(n_obs) + sum(n_obs_discrete) + sum(n_obs_tte) != n_rows) {
    return(NULL)
  }
  pos <- seq_along(n_obs)
  list(gaussian = rep.int(pos, n_obs), discrete = rep.int(pos, n_obs_discrete),
       tte = rep.int(pos, n_obs_tte))
}

# The sdtab columns ferx-core creates for time-to-event rows only (#901).
.FITRX_TTE_COLS <- c("COXSNELL", "MARTINGALE", "DEVIANCE", "PIT")

# The one place that splits `sdtab` (default `fit$sdtab`) into its Gaussian,
# discrete (ferx-core #900) and time-to-event (#901) row blocks (#696). Used
# by `ferx_xpose()`, `ferx_save_fit()` and `ferx_calc_npde()`.
#
# Returns `gaussian`, `discrete` and `tte` (integer row indices into `sdtab`,
# increasing) and `subject`: NULL, or a list with the same three names giving
# the subject (1-based position in fit order) of each of those rows, with an
# `n_subjects` attribute.
#
# One fallback order:
# 1. The `fit$subject_summary` counts (`.fitrx_sdtab_subject()`), a missing
#    `n_obs_discrete` / `n_obs_tte` column counting as 0, when they add up
#    to the sdtab's rows. `subject` is set.
# 2. Otherwise the columns, and `subject` is NULL: the TTE rows are those
#    with a COXSNELL value, in a table carrying all of `.FITRX_TTE_COLS`; the
#    discrete rows the others without an N_OBS value (ferx-core writes N_OBS,
#    the subject's Gaussian count, on every Gaussian row only); the rest are
#    Gaussian.
.fitrx_sdtab_blocks <- function(fit, sdtab = fit$sdtab) {
  n_rows <- if (is.null(sdtab)) 0L else nrow(sdtab)
  summ <- fit$subject_summary
  count <- function(col) {
    v <- summ[[col]]
    if (is.null(v)) 0L else v
  }
  pos <- if (is.data.frame(summ) && nrow(summ) > 0L && !is.null(summ$n_obs)) {
    .fitrx_sdtab_subject(summ$n_obs, count("n_obs_discrete"), n_rows,
                         count("n_obs_tte"))
  }
  if (!is.null(pos)) {
    n_gauss <- length(pos$gaussian)
    n_disc <- length(pos$discrete)
    return(list(
      gaussian = seq_len(n_gauss),
      discrete = n_gauss + seq_len(n_disc),
      tte = n_gauss + n_disc + seq_along(pos$tte),
      subject = structure(pos, n_subjects = nrow(summ))
    ))
  }
  tte <- if (n_rows > 0L && all(.FITRX_TTE_COLS %in% names(sdtab))) {
    !is.na(sdtab$COXSNELL)
  } else {
    rep(FALSE, n_rows)
  }
  disc <- if (n_rows > 0L && "N_OBS" %in% names(sdtab)) {
    is.na(sdtab$N_OBS) & !tte
  } else {
    rep(FALSE, n_rows)
  }
  list(gaussian = which(!(tte | disc)), discrete = which(disc),
       tte = which(tte), subject = NULL)
}
