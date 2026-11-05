# NULL-or-value coercion helpers: `.fitrx_opt_*` shape R values for the JSON
# ferx_save_fit() writes, `.fitrx_unwrap_opt_*` read the engine's optional
# values back on a fit. ferx_load_fit() reads bundles with ferx-core's loader.

# A scalar that may be `null`: the JSON writers turn NaN into `null`, and a VI
# fit under `vi_final_ofv = none` has a NaN OFV by design (#508), which
# `as.numeric(NULL)` would read as `numeric(0)`.
.fitrx_num_or_nan <- function(x) {
  v <- as.numeric(unlist(x, use.names = FALSE))
  if (length(v) == 0L) NaN else v[[1L]]
}

.fitrx_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.numeric(x)
  if (is.na(v)) return(NULL)
  v
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

# A Rust `Option<bool>`. Both directions treat NA as absent: a tri-state flag
# reaches R as NA when the engine recorded no verdict, and writing NA back
# would claim a verdict of "false" to any reader that coerces it.
.fitrx_opt_lgl <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.logical(x))
  if (length(v) != 1L || is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_lgl <- function(x) {
  if (is.null(x)) return(NA)
  v <- suppressWarnings(as.logical(unlist(x, use.names = FALSE)))
  if (length(v) != 1L || is.na(v)) return(NA)
  v
}

.fitrx_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.character(x)
  if (is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_chr <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.character(x)
  if (length(v) == 0L || is.na(v) || !nzchar(v)) return(NULL)
  v
}

.fitrx_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(x)
  if (length(v) == 0L) return(NULL)
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

.fitrx_prior_summary_to_wire <- function(ps) {
  if (!is.data.frame(ps) || nrow(ps) == 0L) return(NULL)
  lapply(seq_len(nrow(ps)), function(i) {
    list(
      name               = as.character(ps$name[[i]]),
      prior_value        = as.numeric(ps$prior_value[[i]]),
      estimate           = as.numeric(ps$estimate[[i]]),
      shift_in_prior_sds = as.numeric(ps$shift_in_prior_sds[[i]]),
      penalty            = as.numeric(ps$penalty[[i]]),
      family             = as.character(ps$family[[i]]),
      prior_lower_95     = as.numeric(ps$prior_lower_95[[i]]),
      prior_upper_95     = as.numeric(ps$prior_upper_95[[i]])
    )
  })
}

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

# A text cell as ferx-core's `.fitrx` writer spells it (`csv_escape()` in
# `io/fitrx.rs`): quoted, with `"` doubled, only when it holds a comma, a quote
# or a newline; verbatim otherwise. NA stays NA (an empty cell). R's reader
# and ferx-core's both take a quoted comma or quote, so a subject ID such as
# `Smith, 2019` round-trips (#475). A newline makes the trip through R only:
# ferx-core's reader is line-based (`parse_csv_row`).
.fitrx_csv_escape <- function(x) {
  x <- as.character(x)
  hit <- !is.na(x) & grepl("[,\"\n]", x)
  x[hit] <- paste0("\"", gsub("\"", "\"\"", x[hit], fixed = TRUE), "\"")
  x
}

# Write a data frame as a bundle CSV with every double column in
# `.fitrx_double_text()` form and every text column in `.fitrx_csv_escape()`
# form, the quoting ferx-core's own writer uses. A number is never quoted.
.fitrx_write_csv_exact <- function(df, path) {
  for (j in which(vapply(df, function(v) is.character(v) || is.factor(v),
                         logical(1)))) {
    df[[j]] <- .fitrx_csv_escape(df[[j]])
  }
  for (j in which(vapply(df, is.double, logical(1)))) {
    df[[j]] <- .fitrx_double_text(df[[j]])
  }
  utils::write.table(df, path, row.names = FALSE, sep = ",", na = "",
                     quote = FALSE)
}

# `read.csv()`, with every column it types as double re-parsed by
# `.fitrx_parse_doubles()`. read.csv() still decides the column types, so a
# loaded table has the shape it always had - except the `id_cols`, which are
# read as text, verbatim: inferring a type turns the IDs `001` and `1.0` into
# `1`, and `as.character()` afterwards only puts the type back (#468). They
# come from a pass with no NA strings, since the writers never put an NA in an
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

# The bundle entry of one level's kappa estimates, as ferx-core names it
# (`level_kappa_entry` in io/fitrx.rs).
.fitrx_level_kappa_entry <- function(level) paste0("ebes_kappa_", level, ".csv")

# `fit$iov_occasion` (the `settings = list(iov_occasion = ...)` spelling the
# glue writes: "column", "dose" or "time(24, 48.5)") to and from the `.fitrx`
# key ferx-core writes for `FitResult::iov_occasion` (#1783): "column",
# "per_dose" or {"time_windows": [24, 48.5]}. NULL both ways when the fit
# records no rule, so `ferx_sir()` / `ferx_covariance()` on a reloaded fit
# derive its occasions as the fit did (#512).
#
# A wire value of any other shape (a newer ferx-core's variant, a hand-edited
# key) drops the field with a `ferx_iov_occasion_dropped_warning` naming the
# value, and the rest of the bundle loads (#524): the fit is then one that
# records no rule, as a bundle saved before ferx-core #1783, and the post-hoc
# steps fall back to the model file's rule.
.fitrx_iov_occasion_to_wire <- function(rule) {
  rule <- .fitrx_opt_chr(rule)
  if (is.null(rule)) return(NULL)
  if (identical(rule, "column")) return("column")
  if (identical(rule, "dose")) return("per_dose")
  edges <- sub("^time\\((.*)\\)$", "\\1", rule)
  if (identical(edges, rule)) {
    stop("ferx_save_fit: `fit$iov_occasion` is \"", rule, "\", not \"column\", ",
         "\"dose\" or \"time(...)\".", call. = FALSE)
  }
  edges <- suppressWarnings(as.numeric(strsplit(edges, ",", fixed = TRUE)[[1]]))
  if (!length(edges) || anyNA(edges)) {
    stop("ferx_save_fit: `fit$iov_occasion` \"", rule, "\" has a breakpoint ",
         "that is not a number.", call. = FALSE)
  }
  # A list, so `auto_unbox` keeps a one-edge rule an array.
  list(time_windows = as.list(edges))
}

# `fit$scoring_settings` / `fit$sir_settings` (#511, #472) on the wire, in
# ferx-core's own layout (`ScoringSettingsWire` / `SirSettingsWire`), so the
# engine reads an R bundle's records and R reads an engine bundle's: the
# top-level `scoring_settings` block keyed as the R list is, and `sir.settings`
# with the `scoring` half flattened beside the SIR keys. Enums are their
# `[fit_options]` tokens on both sides; NA (`ode_stiff_abort_after` off) is
# written as JSON null. The record is first read through the engine's own
# decoder, so one edited into something the engine would refuse is refused
# here, by name, rather than written into a bundle nothing loads.
.fitrx_settings_to_wire <- function(fit, kind) {
  record <- fit[[kind]]
  if (is.null(record)) return(NULL)
  record <- ferx_rust_settings_record(kind, record, "ferx_save_fit", paste0("fit$", kind))
  if (identical(kind, "sir_settings")) {
    # A multivariate-normal proposal records `df = Inf`, which JSON cannot
    # hold; ferx-core's wire spells it `"normal"` (ferx-core#1819), and both
    # its reader and `ferx_load_fit()` read that back as `Inf`.
    if (is.infinite(record$df)) record$df <- "normal"
    record <- c(record[names(record) != "scoring"], record$scoring)
  }
  record
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

# Whether an entry's warning is among `messages` (a fit's flat warnings). A
# chained fit's flat message carries the stage as a `[METHOD] ` prefix that the
# entry keeps apart in `source_method`.
.ferx_warning_detail_in <- function(entry, messages) {
  msg <- entry$message
  if (!is.character(msg) || length(msg) != 1L) return(FALSE)
  msg %in% messages ||
    (is.character(entry$source_method) &&
       paste0("[", entry$source_method, "] ", msg) %in% messages)
}

# The entries of `old` whose warning the fit still carries in `messages`, then
# those of `new` (a post-hoc step's) not already there. ferx_covariance() and
# ferx_sir() drop a stale warning from the flat vector; its details go with it.
.ferx_merge_warning_details <- function(old, new, messages) {
  keep <- Filter(function(e) .ferx_warning_detail_in(e, messages), old %||% list())
  for (e in new %||% list()) {
    if (!any(vapply(keep, identical, logical(1L), e))) keep <- c(keep, list(e))
  }
  if (length(keep) == 0L) NULL else keep
}

# Each entry as JSON, in R's spelling. On save they are written verbatim; on
# load they go through the engine's `WarningEntry` and come back in serde's.
# `details` is an untyped `serde_json::Value`, which keeps `50` an integer and
# `50.0` a float, so a whole-number double is written with its decimal point
# (`always_decimal`): otherwise it would reload as an integer, and the loaded
# list would not be `identical()` to the fit's.
.fitrx_warning_details_json <- function(entries) {
  vapply(entries, function(e) {
    as.character(jsonlite::toJSON(e, auto_unbox = TRUE, digits = I(17),
                                  always_decimal = TRUE, null = "null", na = "null"))
  }, character(1L))
}

# `fit.json:warning_details` for the fit: its entries whose warning it still
# carries, each written verbatim as the object it is. NULL (key absent, as
# ferx-core leaves it) when there are none.
.fitrx_warning_details_to_wire <- function(fit) {
  entries <- .ferx_merge_warning_details(fit$warning_details, NULL,
                                         as.character(fit$warnings %||% character()))
  if (is.null(entries)) return(NULL)
  lapply(.fitrx_warning_details_json(entries), function(j) structure(j, class = "json"))
}
