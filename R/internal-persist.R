# Write-path (ferx_save_fit) and read-path (ferx_load_fit) scalar/vector
# NULL-or-value coercion helpers. Kept side by side (rather than split across
# ferx_save_fit.R / ferx_load_fit.R) since they are mirror images of each
# other and a change to one direction's coercion rule usually needs checking
# against the other. Not merged into single functions: the two directions
# have different semantics (write coerces R values for JSON serialisation,
# read coerces JSON-decoded values back, including its own NA/length quirks),
# so the small differences between e.g. .fitrx_opt_num and
# .fitrx_unwrap_opt_num are intentional, not accidental drift.

.fitrx_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- as.numeric(x)
  if (is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_num <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
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

.fitrx_unwrap_nullable_num_vec <- function(x) .fitrx_unwrap_nullable_vec(x, NA_real_)

# A float vector the engine writes with non-finite entries as `null` (its
# `serde_nan` adapter), such as kappa shrinkage, whose slot reads NaN when no
# observation informs it (ferx-core #1911). `null` comes back as NaN, as the
# engine reads it, so a round trip returns the live fit's value.
.fitrx_unwrap_nan_num_vec <- function(x) .fitrx_unwrap_nullable_vec(x, NaN)

.fitrx_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0L) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) return(NULL)
  v
}

.fitrx_unwrap_opt_int <- function(x) {
  if (is.null(x)) return(NULL)
  v <- suppressWarnings(as.integer(x))
  if (length(v) == 0L || is.na(v)) return(NULL)
  v
}

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

.fitrx_unwrap_opt_num_vec <- function(x) {
  if (is.null(x)) return(NULL)
  v <- as.numeric(unlist(x, use.names = FALSE))
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

.fitrx_prior_summary_from_wire <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  num <- function(el, k) {
    v <- suppressWarnings(as.numeric(el[[k]] %||% NA_real_))
    if (length(v) != 1L) NA_real_ else v
  }
  chr <- function(el, k) {
    v <- as.character(el[[k]] %||% NA_character_)
    if (length(v) != 1L) NA_character_ else v
  }
  data.frame(
    name               = vapply(x, chr, character(1L), "name"),
    prior_value        = vapply(x, num, numeric(1L), "prior_value"),
    estimate           = vapply(x, num, numeric(1L), "estimate"),
    shift_in_prior_sds = vapply(x, num, numeric(1L), "shift_in_prior_sds"),
    penalty            = vapply(x, num, numeric(1L), "penalty"),
    family             = vapply(x, chr, character(1L), "family"),
    prior_lower_95     = vapply(x, num, numeric(1L), "prior_lower_95"),
    prior_upper_95     = vapply(x, num, numeric(1L), "prior_upper_95"),
    stringsAsFactors   = FALSE
  )
}

# Recover the objective's data / prior split from a bundle that does not carry
# it, returning list(ofv_data =, ofv_prior =).
#
# ferx-core's own loader reads a missing split as "this file predates priors",
# and for a bundle *it* wrote that is sound: once the feature existed, its
# writer always emitted the split for a priored fit. It is NOT sound for a
# bundle this package wrote. `ferx_save_fit()` shipped before ferx-r #366 and
# emitted neither half, while R could already fit a priored model - so an
# R-written bundle of a priored fit carries a penalized `ofv` and nothing to
# say so. Reading that as the data half relabels the penalized objective as the
# likelihood, breaks the `aic == ofv_data + 2k` invariant the stored AIC was
# computed under, and hands `ferx_sir()` back the #366 double count.
#
# The split is recoverable from what such a bundle does carry: the engine
# computes `aic = ofv_data + 2 * n_parameters` for every fit, priored or not
# (ferx-core api/fit.rs), and both `aic` and `n_parameters` are on the wire. So
# `ofv_data = aic - 2 * n_parameters` and `ofv_prior = ofv - ofv_data`.
#
# Guarded, because the identity is only as good as the two fields it reads. A
# penalty is a sum of squares, so a recovered prior half below zero means the
# identity does not hold for this file (a hand-edited bundle, a writer that
# computed AIC differently) and the unpriored reading is restored. A recovered
# half within rounding distance of zero is taken as exactly zero, which is what
# the overwhelmingly common unpriored bundle should report.
.fitrx_recover_ofv_split <- function(w) {
  # A null OFV (a VI fit's NaN, #508) gives a NaN data half, as on the fit.
  ofv <- suppressWarnings(.fitrx_num_or_nan(w$ofv))
  unpriored <- list(ofv_data = ofv, ofv_prior = 0, recovered = FALSE)
  if (length(ofv) != 1L || !is.finite(ofv)) return(unpriored)

  aic <- suppressWarnings(as.numeric(w$aic %||% NA_real_))
  k <- suppressWarnings(as.numeric(w$n_parameters %||% NA_real_))
  if (length(aic) != 1L || !is.finite(aic) ||
        length(k) != 1L || !is.finite(k) || k < 0) {
    return(unpriored)
  }

  ofv_data <- aic - 2 * k
  ofv_prior <- ofv - ofv_data
  if (!is.finite(ofv_prior)) return(unpriored)
  # The JSON carries full f64 precision, so the identity round-trips to within
  # rounding; anything below this is noise, not a prior.
  tol <- 1e-9 * max(1, abs(ofv))
  if (ofv_prior < -tol) return(unpriored)
  if (ofv_prior <= tol) return(unpriored)

  list(ofv_data = ofv_data, ofv_prior = ofv_prior, recovered = TRUE)
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

# The per-level kappa estimates of a nested-level fit (ferx-core #1082 PR4):
# one data frame per entry of fit.json's `iov.levels`, read from its
# `ebes_kappa_<LEVEL>.csv`, in the live fit's shape - ID (text, verbatim),
# GROUP (integer), the level's label columns and its kappas (double). An empty
# kappa cell is a NaN estimate, as ferx-core writes and reads it. NULL when
# fit.json names no level. Whatever ferx-core's own loader refuses is refused
# here too, in its words (`parse_level_kappa_csv` and the `iov` checks in
# io/fitrx.rs): a level naming a kappa that `kappa_levels` puts elsewhere, a
# missing CSV, a header that is not the level's, a row with the wrong number
# of fields, an ID that is not a subject of the fit (`subject_ids`), and a
# subject whose GROUPs do not run 1..n.
.fitrx_read_level_kappas <- function(staging, levels, kappa_names = character(),
                                     kappa_levels = character(),
                                     subject_ids = NULL) {
  if (length(levels) == 0L) return(NULL)
  out <- list()
  for (lv in levels) {
    name <- as.character(lv$name)
    columns <- as.character(unlist(lv$columns, use.names = FALSE))
    kappas <- as.character(unlist(lv$kappa_names, use.names = FALSE))
    for (k in kappas) {
      at <- unname(kappa_levels[match(k, kappa_names)])
      if (length(at) != 1L || is.na(at) || !identical(at, name)) {
        stop(sprintf(
          "iov.levels names kappa `%s` at level `%s`, but iov.kappa_levels puts it at %s",
          k, name, if (length(at) != 1L || is.na(at)) "no level" else sprintf("`%s`", at)
        ))
      }
    }
    entry <- .fitrx_level_kappa_entry(name)
    path <- file.path(staging, entry)
    if (!file.exists(path)) {
      stop(sprintf("fit.json names level `%s`, but the bundle has no %s", name, entry))
    }
    expected <- c("ID", "GROUP", columns, kappas)
    # read.csv() pads a short row and folds a long one into the next, so the
    # field count is checked on the raw lines first.
    n_fields <- utils::count.fields(path, sep = ",", quote = "\"", blank.lines.skip = TRUE)
    bad <- which(n_fields != length(expected))
    if (length(bad) > 0L && bad[1L] > 1L) {
      stop(sprintf("%s: row %d has %d fields, expected %d",
                   entry, bad[1L] - 1L, n_fields[bad[1L]], length(expected)))
    }
    df <- .fitrx_read_csv_exact(path, id_cols = "ID")
    if (!identical(names(df), expected)) {
      stop(sprintf("%s: header %s does not match fit.json's level (expected %s)",
                   entry, paste(names(df), collapse = ","),
                   paste(expected, collapse = ",")))
    }
    if (!is.null(subject_ids)) {
      unknown <- setdiff(df$ID, subject_ids)
      if (length(unknown) > 0L) {
        stop(sprintf("%s: unknown ID \"%s\"", entry, unknown[1L]))
      }
    }
    df$GROUP <- as.integer(df$GROUP)
    for (id in unique(df$ID)) {
      g <- df$GROUP[df$ID == id]
      if (anyNA(g) || !identical(sort(g), seq_along(g))) {
        stop(sprintf("%s: subject \"%s\": GROUP must run 1..n without gaps or repeats",
                     entry, id))
      }
    }
    for (col in columns) df[[col]] <- as.double(df[[col]])
    for (col in kappas) {
      v <- as.double(df[[col]])
      v[is.na(v)] <- NaN
      df[[col]] <- v
    }
    out[[name]] <- df
  }
  out
}

# The bundle entry of one level's kappa estimates, as ferx-core names it
# (`level_kappa_entry` in io/fitrx.rs).
.fitrx_level_kappa_entry <- function(level) paste0("ebes_kappa_", level, ".csv")

# `predictions.csv`'s subject IDs (text) as the number ferx-core puts in
# `sdtab$ID`: the ID parsed as Rust parses an f64, else the subject's
# position, 1-based, among ALL subjects (`io::output::sdtab`). A subject with
# no observations has no rows here but still takes a position, so positions
# come from `subject_ids` (the bundle's ebes.csv IDs, every subject in order):
# each block takes the next subject with its text. Without them, or when a
# block's text is not found, a block takes its own index. A bundle written
# before #475 holds the number already, which parses to itself.
#
# `row_subject`, when given, is every row's position, from the per-subject row
# counts (`.fitrx_sdtab_subject()`), and replaces that walk: the walk reads a
# block from the text, so two adjacent subjects with the same text are one
# block (#540). The walk is for a bundle whose counts are missing or do not
# fit its rows.
.fitrx_sdtab_id_number <- function(ids, subject_ids = NULL, row_subject = NULL) {
  if (is.null(row_subject)) {
    blk <- .ferx_subject_index(ids)
    pos <- seq_len(if (length(blk)) max(blk) else 0L)
    if (length(subject_ids) > 0L) {
      first <- ids[!duplicated(blk)]
      walked <- integer(length(first))
      j <- 0L
      for (k in seq_along(first)) {
        rest <- subject_ids[seq.int(j + 1L, length.out = length(subject_ids) - j)]
        hit <- match(first[k], rest)
        if (is.na(hit)) break
        j <- j + hit
        walked[k] <- j
      }
      if (all(walked > 0L)) pos <- walked
    }
    row_subject <- pos[blk]
  }
  out <- as.numeric(row_subject)
  s <- tolower(ids)
  dec <- grepl("^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)(e[-+]?[0-9]+)?$", s)
  out[dec] <- .fitrx_parse_doubles(sub("^\\+", "", s[dec]))
  inf <- grepl("^[-+]?inf(inity)?$", s)
  out[inf] <- ifelse(startsWith(s[inf], "-"), -Inf, Inf)
  out[grepl("^[-+]?nan$", s)] <- NaN
  out
}

# `fit$subject_summary` from ebes.csv's per-subject columns, as the live fit
# carries it. A bundle without n_obs_discrete has no discrete rows (0); one
# written before #554 may have empty n_obs (NA). NULL when ebes.csv has no
# ofv_contribution / n_obs.
.fitrx_read_subject_summary <- function(ebes) {
  if (!all(c("ID", "ofv_contribution", "n_obs") %in% names(ebes))) return(NULL)
  n <- nrow(ebes)
  if (n == 0L) return(NULL)
  disc <- if ("n_obs_discrete" %in% names(ebes)) ebes$n_obs_discrete else rep(0L, n)
  data.frame(
    ID = as.character(ebes$ID),
    ofv_contribution = as.numeric(ebes$ofv_contribution),
    n_obs = as.integer(ebes$n_obs),
    n_obs_discrete = as.integer(disc),
    stringsAsFactors = FALSE
  )
}

# The bundle's sdtab: predictions.csv's Gaussian rows, then
# predictions_discrete.csv's discrete rows when present, in the live fit's
# order. Each row's ID goes back to the engine's number for its subject, the
# subject read from the per-subject row counts when they fit the rows and the
# rows' IDs (`.fitrx_sdtab_subject()`, #540), else walked from the text.
.fitrx_read_sdtab <- function(staging, subject_ids, summ) {
  # Every sdtab column is double on the live fit, ID included (#417).
  sdtab <- .fitrx_read_csv_exact(file.path(staging, "predictions.csv"),
                                 id_cols = "ID", integer_cols = character())
  n_gauss <- nrow(sdtab)
  disc_path <- file.path(staging, "predictions_discrete.csv")
  if (file.exists(disc_path)) {
    sdtab <- rbind(sdtab, .fitrx_read_csv_exact(disc_path, id_cols = "ID",
                                                 integer_cols = character()))
  }
  if (is.null(sdtab$ID)) return(sdtab)
  rows <- if (!is.null(summ) && length(subject_ids) == nrow(summ)) {
    .fitrx_sdtab_subject(summ$n_obs, summ$n_obs_discrete, nrow(sdtab))
  }
  row_subject <- NULL
  if (!is.null(rows) && length(rows$gaussian) == n_gauss) {
    pos <- c(rows$gaussian, rows$discrete)
    if (identical(sdtab$ID, subject_ids[pos])) row_subject <- pos
  }
  sdtab$ID <- .fitrx_sdtab_id_number(sdtab$ID, subject_ids, row_subject)
  sdtab
}

# `fit$iov_occasion` (the `settings = list(iov_occasion = ...)` spelling the
# glue writes: "column", "dose" or "time(24, 48.5)") to and from the `.fitrx`
# key ferx-core writes for `FitResult::iov_occasion` (#1783): "column",
# "per_dose" or {"time_windows": [24, 48.5]}. NULL both ways when the fit
# records no rule, so `ferx_sir()` / `ferx_covariance()` on a reloaded fit
# derive its occasions as the fit did (#512).
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

.fitrx_iov_occasion_from_wire <- function(w) {
  if (is.null(w)) return(NULL)
  if (is.character(w) && length(w) == 1L) {
    if (identical(w, "column")) return("column")
    if (identical(w, "per_dose")) return("dose")
  }
  if (is.list(w) && identical(names(w), "time_windows")) {
    edges <- as.numeric(unlist(w$time_windows, use.names = FALSE))
    # The shortest spelling that reads back to the same double, as the glue
    # writes it: 120.1, not 120.09999999999999; 100000, not 1e+05 (Rust's
    # `{}` never writes an exponent).
    spell <- vapply(edges, function(e) {
      for (d in 15:17) {
        s <- format(e, digits = d, trim = TRUE, scientific = FALSE)
        if (as.numeric(s) == e) return(s)
      }
      s
    }, character(1))
    return(sprintf("time(%s)", paste(spell, collapse = ", ")))
  }
  stop("ferx_load_fit: the bundle's `iov_occasion` is not \"column\", ",
       "\"per_dose\" or {\"time_windows\": [...]}.", call. = FALSE)
}

# A JSON-carried fit field (`fit$reader_settings`, `fit$population_fingerprint`,
# #462) back from the wire, where `ferx_save_fit()` wrote it verbatim and
# `read_json(simplifyVector = FALSE)` read it as nested lists. jsonlite spells
# numbers its own way (`[0, 1]` for serde's `[0.0, 1.0]`), so the lists are
# written back losslessly - arrays stay arrays (lists never unbox), scalars
# unbox, 17 significant digits identify every double - and the engine
# re-serialises them in its own spelling, which is the string the fit carried
# (#526 review 2). NULL when the bundle carries none.
.fitrx_json_from_wire <- function(w, field) {
  if (is.null(w)) return(NULL)
  json <- as.character(jsonlite::toJSON(w, auto_unbox = TRUE, null = "null",
                                        digits = I(17)))
  tryCatch(
    ferx_rust_fit_json_canonical(field, json),
    error = function(e) {
      stop("ferx_load_fit: the bundle's ", conditionMessage(e), call. = FALSE)
    }
  )
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
    # hold: it would be written as null, which neither ferx-core nor
    # `ferx_load_fit()` reads as a number, and the bundle would not load. Such
    # a fit is saved without the record and reloads as one that has none
    # (ferx-core#1819), which the caller is told.
    if (!is.finite(record$df)) {
      warning(
        "ferx_save_fit: the fit's SIR used a normal proposal (`sir_df = Inf`), ",
        "which a .fitrx bundle cannot record yet (ferx-core#1819), so its SIR ",
        "settings are not saved. On the loaded fit, ferx_sir() runs with the ",
        "defaults unless you pass the settings, `sir_df = Inf` included.",
        call. = FALSE
      )
      return(NULL)
    }
    record <- c(record[names(record) != "scoring"], record$scoring)
  }
  record
}

# The two records back from the wire, through the engine's decoder, so the
# loaded list is `identical()` to the one the fit carried: `read_json()` reads
# a whole-number double as an integer and JSON null as NULL, and the decoder
# writes each field back in the type `ferx_fit()` gives it. `sir.settings`
# comes back un-flattened. NULL when the bundle carries none (one saved before
# the records existed), which keeps today's behaviour.
#
# Migration (#535): a key the record lacks is filled by the decoder as
# ferx-core's wire fills it (the glue's `SCORING_KEYS` / `SIR_KEYS` table). A
# key this ferx does not know - a bundle written by a newer ferx-core, whose
# record gained a setting - drops that record with a
# `ferx_settings_record_dropped_warning` naming it, and the rest of the bundle
# loads: the fit is then one that records none, as a bundle saved before the
# records existed. ferx-core's own reader ignores such a key and keeps the
# record; R does not, since re-scoring under the record without a setting the
# writer used would not repeat the writer's run.
.fitrx_settings_from_wire <- function(w, kind) {
  if (is.null(w)) return(NULL)
  is_sir <- identical(kind, "sir_settings")
  unknown <- setdiff(names(w), .fitrx_settings_wire_keys(kind))
  if (length(unknown) > 0L) {
    block <- if (is_sir) "sir.settings" else "scoring_settings"
    warning(warningCondition(paste0(
      "ferx_load_fit: the bundle's `", block, "` record names ",
      if (length(unknown) == 1L) "a setting" else "settings",
      " this ferx does not know (", paste0("\"", unknown, "\"", collapse = ", "),
      "), so fit$", kind, " is not loaded; it was likely written by a newer ",
      "ferx-core. The rest of the fit loads, and ferx_covariance() / ferx_sir() ",
      "run as on a fit that records no ", kind, ". Update ferx to keep it."),
      class = "ferx_settings_record_dropped_warning", keys = unknown,
      record = kind))
    return(NULL)
  }
  # A normal-proposal record (`df = Inf`) is written by ferx-core's serde as
  # `"df": null` (ferx-core#1819). It is read as no record, as
  # `ferx_save_fit()` writes it, so the bundle still loads.
  if (is_sir && "df" %in% names(w) && is.null(w$df)) return(NULL)
  w <- lapply(w, function(v) if (is.null(v)) NA else v)
  tryCatch(
    ferx_rust_settings_record(
      if (is_sir) "sir_settings_wire" else kind, w, "ferx_load_fit",
      if (is_sir) "sir.settings" else "scoring_settings"
    ),
    error = function(e) stop(conditionMessage(e), call. = FALSE)
  )
}

# The names a record may carry on the wire: the engine's default record's, with
# `sir.settings`' scoring half flattened beside the SIR keys.
.fitrx_settings_wire_keys <- function(kind) {
  if (!identical(kind, "sir_settings")) return(names(.ferx_default_settings(kind)))
  c(setdiff(names(.ferx_default_settings("sir_settings")), "scoring"),
    names(.ferx_default_settings("scoring_settings")))
}
