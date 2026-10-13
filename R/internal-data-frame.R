# -- A data.frame wherever `data` is accepted (#419) --------------------------
#
# The engine reads a NONMEM CSV from a path, and stays that way: a data.frame
# is written once, by `.ferx_data_arg()`, to a CSV the engine then reads like
# any other. The file is content-addressed (`ferx-data/data-<md5>.csv` under
# the session's `tempdir()`), so writing the same frame twice gives one file,
# and the same bytes always land at the same name.
#
# A fit read from such a file carries the file's text as `fit$data_source`
# (as `fit$model_source` carries the model's), set by `.ferx_fit_from_raw()`.
# Every step that re-reads a fit's data finds it through
# `.ferx_fit_data_path()`, which writes the text back when the file has gone -
# after `tempdir()` was cleaned, or in a later session - so the bytes, and the
# `data_hash` the engine checks them against, are the ones the fit read.

# The session's cache of written data files.
.ferx_data_cache_dir <- function() file.path(tempdir(), "ferx-data")

# `data` as the engine takes it: a path stays a path (each caller keeps its
# own checks on it), NULL stays NULL, and a data.frame - a tibble, a
# `ferx_mbma_data()` table, a `ferx_data` (its retained rows) - becomes the
# path of the CSV written from it. Anything else is refused here, naming
# `arg`, before a caller's `file.exists()` would fail on it.
.ferx_data_arg <- function(data, arg = "data") {
  if (is.null(data) || is.character(data)) return(data)
  if (is.data.frame(data)) return(.ferx_write_data_text(.ferx_data_csv_text(data, arg)))
  stop(sprintf(
    "`%s` must be a path to a NONMEM-format CSV file or a data.frame, not %s.",
    arg, paste(class(data), collapse = "/")
  ), call. = FALSE)
}

# The CSV text the engine reads for `df`: a header and one line per row,
# "\n"-terminated, UTF-8.
#  - A double is written as the shortest decimal that reads back as exactly
#    that double (`.fitrx_double_text()`, at most 17 significant digits), so
#    a fit of the frame equals a fit of the CSV it was read from.
#  - A missing value (NA, NaN) is `.`, the NONMEM spelling the reader takes as
#    missing in every column.
#  - A factor is its labels, a logical 1/0, a Date its `format()`, a
#    date-time its `format()` in UTC (so the bytes, and the fit's
#    `data_hash`, do not depend on the session's time zone), a
#    `bit64::integer64` its digits. Text is written verbatim - an ID `007`
#    stays `007` - and quoted, with `"` doubled, when it holds a comma, a
#    quote or a line break. (The engine reads the text `NA` as missing in
#    every column but ID, as it does in any CSV.)
#  - Every column is written, under its own name, in its own order: the
#    engine decides which columns it reads.
# A column of any other type (a list, a matrix, complex) is refused, naming
# it, and so are duplicate or empty column names.
.ferx_data_csv_text <- function(df, arg = "data") {
  df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
  nms <- names(df)
  if (length(nms) == 0L) {
    stop(sprintf("`%s` has no columns.", arg), call. = FALSE)
  }
  if (anyNA(nms) || !all(nzchar(nms))) {
    stop(sprintf("`%s` has a column without a name.", arg), call. = FALSE)
  }
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0L) {
    stop(sprintf(
      "`%s` has more than one column named %s; the engine would read only one of them.",
      arg, paste0("`", dup, "`", collapse = ", ")
    ), call. = FALSE)
  }
  bad <- !vapply(df, .ferx_data_column_ok, logical(1))
  if (any(bad)) {
    stop(sprintf(
      "`%s` column(s) %s cannot be written as CSV: give each a number, text, factor, logical or date.",
      arg, paste0("`", nms[bad], "`", collapse = ", ")
    ), call. = FALSE)
  }
  cells <- lapply(df, .ferx_data_column_text)
  lines <- c(
    paste(.ferx_data_csv_quote(enc2utf8(nms)), collapse = ","),
    if (nrow(df) > 0L) do.call(paste, c(unname(cells), sep = ","))
  )
  text <- paste0(paste(lines, collapse = "\n"), "\n")
  Encoding(text) <- "UTF-8"
  text
}

# A `bit64::integer64` is stored as a double holding the integer's bits, so
# it is written through bit64's own methods; without bit64 it is refused
# rather than written as the bits' meaningless double value.
.ferx_data_column_ok <- function(v) {
  if (inherits(v, "integer64")) return(requireNamespace("bit64", quietly = TRUE))
  # Checked before is.list(): a POSIXlt is a list.
  if (inherits(v, c("Date", "POSIXt"))) return(TRUE)
  is.null(dim(v)) && !is.list(v) &&
    (is.numeric(v) || is.character(v) || is.factor(v) || is.logical(v))
}

.ferx_data_column_text <- function(v) {
  out <- if (inherits(v, "integer64")) {
    as.character(v)
  } else if (inherits(v, "POSIXt")) {
    .ferx_data_csv_quote(format(v, tz = "UTC", usetz = TRUE))
  } else if (is.factor(v) || inherits(v, "Date")) {
    .ferx_data_csv_quote(enc2utf8(as.character(if (is.factor(v)) v else format(v))))
  } else if (is.logical(v)) {
    ifelse(v, "1", "0")
  } else if (is.integer(v)) {
    as.character(v)
  } else if (is.double(v)) {
    .fitrx_double_text(v)
  } else {
    .ferx_data_csv_quote(enc2utf8(v))
  }
  out[is.na(v) | is.na(out)] <- "."
  out
}

# Quote a text cell when the CSV reader would otherwise split or misread it.
.ferx_data_csv_quote <- function(x) {
  hit <- !is.na(x) & grepl("[,\"\r\n]", x)
  x[hit] <- paste0("\"", gsub("\"", "\"\"", x[hit], fixed = TRUE), "\"")
  x
}

# Write `text` to the session cache under its md5 and return the path. A file
# already there with that md5 is reused; any other is (re)written, through a
# temporary name so a reader never sees half a file.
.ferx_write_data_text <- function(text) {
  dir <- .ferx_data_cache_dir()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  bytes <- charToRaw(text)
  tmp <- tempfile("partial-", tmpdir = dir, fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(bytes, tmp)
  md5 <- unname(tools::md5sum(tmp))
  path <- file.path(dir, paste0("data-", md5, ".csv"))
  if (!(file.exists(path) && identical(unname(tools::md5sum(path)), md5))) {
    if (!file.rename(tmp, path)) {
      stop("Could not write the data.frame to ", path, call. = FALSE)
    }
  }
  normalizePath(path)
}

# The text of a session-cache data file, for `fit$data_source`; NULL for any
# other path, so a fit of a user's own CSV does not copy the file into the fit.
.ferx_data_cache_text <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) ||
      !file.exists(path)) {
    return(NULL)
  }
  dir <- normalizePath(.ferx_data_cache_dir(), mustWork = FALSE)
  if (!identical(normalizePath(dirname(path), mustWork = FALSE), dir) ||
      is.null(.ferx_data_cache_md5(path))) {
    return(NULL)
  }
  .ferx_read_data_text(path)
}

# The md5 a data file written from a data.frame carries in its name
# (`ferx-data/data-<md5>.csv`, in any session's tempdir()); NULL for any other
# path. `ferx_load_fit()` uses it to tell a bundle's data that came from a
# data.frame - which exists nowhere else - from a copy of the user's own CSV.
.ferx_data_cache_md5 <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) ||
      basename(dirname(path)) != "ferx-data") {
    return(NULL)
  }
  m <- regmatches(basename(path),
                  regexec("^data-([0-9a-f]{32})\\.csv$", basename(path)))[[1L]]
  if (length(m) == 2L) m[[2L]] else NULL
}

.ferx_read_data_text <- function(path) {
  text <- rawToChar(readBin(path, "raw", file.size(path)))
  Encoding(text) <- "UTF-8"
  text
}

# `data` left out of a call that takes `fit =` (ferx_simulate(), ferx_predict(),
# ferx_predict_survival()): the data that fit was fitted to, when it can still
# be found - its file, or its `data_source` written back - so a fit of a
# data.frame simulates its own design; otherwise, as before, the model file's
# `[data]` block.
.ferx_default_data <- function(model, fit = NULL) {
  if (inherits(fit, "ferx_fit")) {
    path <- .ferx_fit_data_path(fit)
    if (!is.null(path) && length(path) == 1L && !is.na(path) && nzchar(path) &&
        file.exists(path)) {
      return(path)
    }
  }
  .ferx_model_data_path(model)
}

# The data file a fit was read from, for a step that reads it again. The
# recorded path when it is still there; otherwise, for a fit that carries its
# data (`fit$data_source`), the path that text is written back to. Any other
# fit gets its recorded path back unchanged, and the caller reports it missing
# as before.
.ferx_fit_data_path <- function(fit) {
  path <- fit$data_path
  if (!is.null(path) && length(path) == 1L && !is.na(path) && nzchar(path) &&
      file.exists(path)) {
    return(path)
  }
  src <- fit$data_source
  if (is.character(src) && length(src) == 1L && !is.na(src)) {
    return(.ferx_write_data_text(src))
  }
  path
}
