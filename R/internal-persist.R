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
  ofv <- suppressWarnings(as.numeric(w$ofv %||% NA_real_))
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

# Write a data frame as a bundle CSV with every double column in
# `.fitrx_double_text()` form. `quote = TRUE` keeps write.csv()'s behaviour
# for the entries that used it: the columns that were text are quoted, and a
# number never is.
.fitrx_write_csv_exact <- function(df, path, quote = FALSE) {
  text_cols <- which(vapply(df, function(v) is.character(v) || is.factor(v),
                            logical(1)))
  for (j in which(vapply(df, is.double, logical(1)))) {
    df[[j]] <- .fitrx_double_text(df[[j]])
  }
  utils::write.table(
    df, path,
    row.names = FALSE, sep = ",", na = "",
    quote = if (isTRUE(quote)) text_cols else FALSE,
    qmethod = "double"
  )
}

# `read.csv()`, with every column it types as double re-parsed by
# `.fitrx_parse_doubles()`. read.csv() still decides the column types, so a
# loaded table has the shape it always had - except the `id_cols`, which are
# read as text, verbatim: inferring a type turns the IDs `001` and `1.0` into
# `1`, and `as.character()` afterwards only puts the type back (#468).
.fitrx_read_csv_exact <- function(path, id_cols = character()) {
  header <- names(utils::read.csv(path, nrows = 0L, check.names = FALSE))
  id_cols <- intersect(id_cols, header)
  col_classes <- if (length(id_cols) > 0L) {
    stats::setNames(rep("character", length(id_cols)), id_cols)
  } else {
    NA
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses = col_classes)
  dbl <- which(vapply(df, is.double, logical(1)))
  if (length(dbl) > 0L) {
    raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = "character", na.strings = character())
    for (j in dbl) df[[j]] <- .fitrx_parse_doubles(raw[[j]])
  }
  df
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
.fitrx_settings_from_wire <- function(w, kind) {
  if (is.null(w)) return(NULL)
  is_sir <- identical(kind, "sir_settings")
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
