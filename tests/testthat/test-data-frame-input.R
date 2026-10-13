# A data.frame wherever `data` is accepted (#419). The frame is written to a
# CSV in the session's temporary directory and fitted from there; these tests
# pin what that CSV holds, that a fit of it equals the fit of the CSV the frame
# came from, and that a fit keeps working once the file is gone.

# Local aliases avoid the ::: operator (undesirable_operator_linter).
.df_csv_text    <- getFromNamespace(".ferx_data_csv_text",    "ferx")
.data_arg       <- getFromNamespace(".ferx_data_arg",         "ferx")
.fit_data_path  <- getFromNamespace(".ferx_fit_data_path",    "ferx")
.data_cache_dir <- getFromNamespace(".ferx_data_cache_dir",   "ferx")
.parse_doubles  <- getFromNamespace(".fitrx_parse_doubles",   "ferx")

# Every cell of a written CSV as text, exactly as written: no type guessing,
# no NA strings.
.read_cells <- function(path) {
  utils::read.csv(path, colClasses = "character", na.strings = character(),
                  check.names = FALSE)
}

# Remove every file a data.frame was written to, as a later session (or a
# cleaned tempdir()) would find it.
.forget_written_data <- function() {
  unlink(.data_cache_dir(), recursive = TRUE)
}

# -- The CSV text -------------------------------------------------------------

test_that("a double is written as text that reads back as exactly that double", {
  x <- c(0.1, 1 / 3, 2 / 3, pi * 1e-300, 2^-1074, .Machine$double.xmax,
         1e22, 123456789.123456789, -0.30000000000000004, 5.3653, 100)
  path <- .data_arg(data.frame(ID = seq_along(x), DV = x))
  cells <- .read_cells(path)
  expect_identical(.parse_doubles(cells$DV), x)
  # The shortest such text: a value 15 digits hold is not padded to 17.
  expect_identical(cells$DV[c(1L, 10L, 11L)], c("0.1", "5.3653", "100"))
  expect_identical(cells$DV[c(2L, 8L)], c("0.3333333333333333", "123456789.12345679"))
})

test_that("missing values are `.`, factors their labels, logicals 1/0", {
  df <- data.frame(
    ID   = 1:4,
    DV   = c(1.5, NA, NaN, 2),
    AMT  = c(NA_integer_, 100L, NA_integer_, NA_integer_),
    SEX  = factor(c("F", "M", NA, "F"), levels = c("M", "F")),
    FLAG = c(TRUE, FALSE, NA, TRUE),
    NOTE = c("a", NA, "b", "."),
    stringsAsFactors = FALSE
  )
  cells <- .read_cells(.data_arg(df))
  expect_identical(cells$DV,   c("1.5", ".", ".", "2"))
  expect_identical(cells$AMT,  c(".", "100", ".", "."))
  expect_identical(cells$SEX,  c("F", "M", ".", "F"))
  expect_identical(cells$FLAG, c("1", "0", ".", "1"))
  expect_identical(cells$NOTE, c("a", ".", "b", "."))
})

test_that("text IDs are written verbatim, quoted only where the reader needs it", {
  ids <- c("007", "NA", "Smith, 2019", "say \"hi\"", "1.0")
  path <- .data_arg(data.frame(ID = ids, TIME = 0, stringsAsFactors = FALSE))
  expect_identical(.read_cells(path)$ID, ids)
  lines <- readLines(path)
  expect_identical(lines[2:3], c("007,0", "NA,0"))
  expect_identical(lines[4:5], c("\"Smith, 2019\",0", "\"say \"\"hi\"\"\",0"))
})

test_that("every column is written, under its own name, in its own order", {
  df <- data.frame(ZZ = 1, ID = 1, `My col` = 2, RESGRP = 3, TIME = 0,
                   check.names = FALSE)
  path <- .data_arg(df)
  expect_identical(readLines(path)[1L], "ZZ,ID,My col,RESGRP,TIME")
  expect_identical(names(.read_cells(path)), names(df))
})

test_that("a tibble-like subclass and a ferx_data are written as their rows", {
  df <- data.frame(ID = 1:2, DV = c(1, 2))
  sub <- structure(df, class = c("tbl_df", "tbl", "data.frame"))
  expect_identical(.data_arg(sub), .data_arg(df))
  sel <- ferx_apply_selection(df, ignore = "DV > 1.5")
  expect_identical(.read_cells(.data_arg(sel))$ID, "1")
})

test_that("a column that is not a vector of numbers, text or dates is refused, by name", {
  df <- data.frame(ID = 1:2)
  df$L <- list(1, 2)
  expect_error(.data_arg(df), "`L`")
  expect_error(.data_arg(list(ID = 1)), "path to a NONMEM-format CSV file or a data.frame")
  expect_error(.data_arg(42), "not numeric")
})

test_that("the same frame lands at one content-addressed file holding exactly its text", {
  df <- data.frame(ID = 1:3, DV = c(0.1, 0.2, 0.3))
  a <- .data_arg(df)
  b <- .data_arg(df)
  expect_identical(a, b)
  expect_match(basename(a), "^data-[0-9a-f]{32}\\.csv$")
  expect_identical(rawToChar(readBin(a, "raw", file.size(a))), .df_csv_text(df))
  expect_false(identical(.data_arg(transform(df, DV = DV + 1)), a))
})

test_that("a path and NULL pass through untouched", {
  expect_identical(.data_arg("some/file.csv"), "some/file.csv")
  expect_null(.data_arg(NULL))
})

# -- Fits of a data.frame equal fits of the CSV it came from ------------------

# The data.frame a user gets from the bundled CSV: `.` read as missing, so the
# columns are numbers and the writer's double text is what the engine reads.
.read_example <- function(name) {
  utils::read.csv(ferx_example(name)$data, na.strings = ".",
                  stringsAsFactors = FALSE, check.names = FALSE)
}

.fit_quiet <- function(model, data, ...) {
  ferx_fit(model, data, verbose = FALSE, threads = 1L, ...)
}

.expect_same_fit <- function(a, b) {
  expect_identical(a$ofv, b$ofv)
  expect_identical(a$theta, b$theta)
  expect_identical(a$omega, b$omega)
  expect_identical(a$sigma, b$sigma)
  expect_identical(a$sdtab$IPRED, b$sdtab$IPRED)
}

test_that("a warfarin fit of the data.frame is identical to the fit of its CSV", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  settings <- list(maxiter = 30L)
  from_csv <- .fit_quiet(ex$model, ex$data, settings = settings)
  from_df  <- .fit_quiet(ex$model, .read_example("warfarin"), settings = settings)
  .expect_same_fit(from_df, from_csv)
  expect_null(from_csv$data_source)
  expect_type(from_df$data_source, "character")
})

test_that("an mbma_placebo fit of the data.frame is identical to the fit of its CSV", {
  skip_on_cran()
  ex <- ferx_example("mbma_placebo")
  from_csv <- .fit_quiet(ex$model, ex$data, covariance = FALSE)
  from_df  <- .fit_quiet(ex$model, .read_example("mbma_placebo"),
                         covariance = FALSE)
  .expect_same_fit(from_df, from_csv)
  expect_identical(from_df$theta_levels$value, from_csv$theta_levels$value)
})

test_that("values that need 17 significant digits reach the engine exactly", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  df <- .read_example("warfarin")
  # One ulp off each observation: 15 digits no longer hold these values.
  df$DV <- df$DV * (1 + .Machine$double.eps)
  # The reference CSV is written by hand at 17 digits, which any correctly
  # rounding reader takes back exactly.
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  txt <- df
  txt$DV <- ifelse(is.na(df$DV), ".", sprintf("%.17g", df$DV))
  txt$AMT <- ifelse(is.na(df$AMT), ".", sprintf("%.17g", df$AMT))
  utils::write.table(txt, csv, sep = ",", quote = FALSE, row.names = FALSE)
  settings <- list(maxiter = 5L)
  .expect_same_fit(.fit_quiet(ex$model, df, settings = settings),
                   .fit_quiet(ex$model, csv, settings = settings))
})

test_that("text subject IDs reach the fit exactly", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  df <- .read_example("warfarin")
  ids <- unique(df$ID)
  text <- sprintf("%03d", ids)
  text[1:3] <- c("NA", "Smith, 2019", "say \"hi\"")
  df$ID <- text[match(df$ID, ids)]
  fit <- .fit_quiet(ex$model, df, covariance = FALSE,
                    settings = list(maxiter = 2L))
  expect_identical(as.character(fit$subject_ids), text)
})

test_that("ferx_model() keeps a data.frame and ferx_fit() fits it", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  df <- .read_example("warfarin")
  m <- ferx_model(df, ex$model)
  expect_identical(m$data, df)
  expect_output(print(m), "<data.frame: \\d+ rows, 8 columns>")
  settings <- list(maxiter = 5L)
  .expect_same_fit(ferx_fit(m, verbose = FALSE, threads = 1L, settings = settings),
                   .fit_quiet(ex$model, ex$data, settings = settings))
})

test_that("predictions and simulations of a data.frame equal those of its CSV", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  df <- .read_example("warfarin")
  fit <- warfarin_fit()
  expect_identical(ferx_predict(ex$model, df, fit = fit),
                   ferx_predict(ex$model, ex$data, fit = fit))
  expect_identical(ferx_simulate(ex$model, df, n_sim = 2L, seed = 3L, fit = fit),
                   ferx_simulate(ex$model, ex$data, n_sim = 2L, seed = 3L, fit = fit))
  expect_identical(ferx_get_columns(df), ferx_get_columns(ex$data))
})

# -- The fit keeps its data ---------------------------------------------------

test_that("post-hoc steps work after the written file is gone", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  fit <- .fit_quiet(ex$model, .read_example("warfarin"), covariance = TRUE,
                    settings = list(maxiter = 30L))
  written <- fit$data_path
  expect_true(file.exists(written))
  .forget_written_data()
  expect_false(file.exists(written))

  # Written back from fit$data_source, byte for byte, so the engine's
  # data_hash check passes.
  path <- .fit_data_path(fit)
  expect_true(file.exists(path))
  expect_identical(rawToChar(readBin(path, "raw", file.size(path))),
                   fit$data_source)

  .forget_written_data()
  cov <- ferx_covariance(fit)
  expect_true(is.matrix(cov$cov_matrix))
  .forget_written_data()
  sir <- ferx_sir(fit, sir_samples = 200L, sir_resamples = 50L, sir_seed = 1L)
  expect_false(is.null(sir$sir_ci_theta))
  .forget_written_data()
  npde <- ferx_calc_npde(fit, nsim = 20L, seed = 1L)
  expect_true("NPDE" %in% names(npde$sdtab))
})

test_that("ferx_save_fit() bundles a data.frame fit's data and ferx_load_fit() restores it", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  fit <- .fit_quiet(ex$model, .read_example("warfarin"), covariance = TRUE,
                    settings = list(maxiter = 30L))
  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out), add = TRUE)
  .forget_written_data()
  # include_data = FALSE, the default: the data goes in anyway.
  ferx_save_fit(fit, out)
  expect_true("data.csv" %in% utils::unzip(out, list = TRUE)$Name)
  .forget_written_data()

  loaded <- ferx_load_fit(out)
  expect_identical(loaded$data_source, fit$data_source)
  unlink(loaded$data_path)
  cov <- ferx_covariance(loaded)
  expect_true(is.matrix(cov$cov_matrix))
})

test_that("a CSV fit without include_data bundles no data, as before", {
  skip_on_cran()
  fit <- warfarin_fit()
  out <- tempfile(fileext = ".fitrx")
  on.exit(unlink(out), add = TRUE)
  ferx_save_fit(fit, out)
  expect_false("data.csv" %in% utils::unzip(out, list = TRUE)$Name)
  expect_null(ferx_load_fit(out)$data_source)
})

test_that("a fit without data_source keeps reporting a missing file", {
  fake <- list(data_path = file.path(tempdir(), "no-such-data.csv"))
  expect_identical(.fit_data_path(fake), fake$data_path)
})

test_that("ferx_collect() gives a background data.frame fit a file in this session", {
  text <- .df_csv_text(data.frame(ID = 1, TIME = 0, DV = 1))
  fake <- structure(
    list(data_path = file.path(tempdir(), "gone-worker", "data.csv"),
         data_source = text),
    class = "ferx_fit"
  )
  local_mocked_bindings(.ferx_collect_callr = function(...) fake,
                        .package = "ferx")
  job <- structure(list(backend = "callr", user_wanted_trace = TRUE),
                   class = "ferx_job")
  fit <- ferx_collect(job)
  expect_true(file.exists(fit$data_path))
  expect_identical(rawToChar(readBin(fit$data_path, "raw", file.size(fit$data_path))),
                   text)
})
