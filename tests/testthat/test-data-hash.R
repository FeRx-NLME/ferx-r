# fit$data_hash names the bytes the fit read, and ferx_save_fit(include_data =
# TRUE) bundles only those bytes (#585) - the data-side twin of #577 / #584.

dh_bytes <- function(path) readBin(path, "raw", n = file.size(path))

# The bytes of one entry of a .fitrx bundle.
fitrx_entry_bytes <- function(path, entry) {
  staging <- tempfile("fitrx_entry_")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::unzip(path, files = entry, exdir = staging, junkpaths = TRUE)
  dh_bytes(file.path(staging, entry))
}

# The warfarin data with one DV changed: same shape, different rows.
dh_edited_csv <- function(csv) {
  lines <- readLines(csv)
  hdr <- strsplit(lines[1], ",", fixed = TRUE)[[1]]
  dv <- match("DV", hdr)
  row <- which(vapply(strsplit(lines[-1], ",", fixed = TRUE),
                      function(f) suppressWarnings(!is.na(as.numeric(f[dv])) &&
                                                     as.numeric(f[dv]) > 0),
                      logical(1)))[1] + 1L
  f <- strsplit(lines[row], ",", fixed = TRUE)[[1]]
  f[dv] <- format(as.numeric(f[dv]) * 2)
  lines[row] <- paste(f, collapse = ",")
  out <- tempfile(fileext = ".csv")
  writeLines(lines, out)
  out
}

dh_fit <- function(model, data) {
  suppressWarnings(ferx_fit(model, data, method = "focei", covariance = FALSE,
                            verbose = FALSE, settings = list(maxiter = 2L)))
}

test_that("an edit while the glue reads the data does not reach fit$data_hash (#585)", {
  # The glue used to read the population from the path and, after fit()
  # returned, hash the path a second time, so an edit made during the fit was
  # recorded as the data that was fitted. The data is served through a FIFO:
  # the moment the glue opens it, the server renames an edited copy onto its
  # path - so any later open gets the edited file - and only then sends the
  # original rows. No timing is involved.
  skip_on_cran()
  skip_on_os("windows")
  skip_if(!nzchar(Sys.which("mkfifo")), "no mkfifo")
  ex <- ferx_example("warfarin")
  original <- dh_bytes(ex$data)
  edited <- dh_edited_csv(ex$data)
  fifo <- tempfile(fileext = ".csv")
  expect_identical(system2("mkfifo", fifo), 0L)
  plain <- tempfile(tmpdir = dirname(fifo), fileext = ".csv")
  file.copy(edited, plain)
  server <- callr::r_bg(function(fifo, plain, csv) {
    rows <- readBin(csv, "raw", n = file.size(csv))
    con <- file(fifo, "wb")  # blocks until the glue opens the FIFO
    stopifnot(file.rename(plain, fifo))
    writeBin(rows, con)
    close(con)
  }, args = list(fifo, plain, ex$data))
  withr::defer(server$kill())
  fit <- dh_fit(ex$model, fifo)
  server$wait(5000)
  expect_identical(server$get_exit_status(), 0L)
  # The edit did land during the fit ...
  expect_identical(dh_bytes(fifo), dh_bytes(edited))
  # ... and the hash is the one of the rows that were fitted.
  expect_identical(fit$data_hash, ferx_rust_sha256_raw(original))
  expect_false(identical(fit$data_hash, ferx_rust_sha256_raw(dh_bytes(edited))))
})

test_that("ferx_save_fit(include_data = TRUE) refuses a data file edited since the fit (#585)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  data <- withr::local_tempfile(fileext = ".csv")
  file.copy(ex$data, data)
  fitted <- dh_bytes(data)
  fit <- dh_fit(ex$model, data)
  expect_identical(fit$data_hash, ferx_rust_sha256_raw(fitted))

  # Unchanged: bundled, and the bundle's data.csv is the fitted bytes.
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path, include_data = TRUE)
  expect_identical(fitrx_entry_bytes(path, "data.csv"), fitted)

  # Edited after the fit: refused, classed, naming the file, nothing written.
  file.copy(dh_edited_csv(ex$data), data, overwrite = TRUE)
  refused <- withr::local_tempfile(fileext = ".fitrx")
  err <- expect_error(ferx_save_fit(fit, refused, include_data = TRUE),
                      class = "ferx_data_edited_error")
  expect_match(conditionMessage(err), basename(data), fixed = TRUE)
  expect_match(conditionMessage(err), "no longer matches fit$data_hash",
               fixed = TRUE)
  expect_identical(err$data_path, fit$data_path)
  expect_false(file.exists(refused))
  # Without the data there is nothing to check.
  expect_no_error(ferx_save_fit(fit, refused))
  expect_false("data.csv" %in% utils::unzip(refused, list = TRUE)$Name)

  # A fit with no recorded hash (an older ferx) bundles the file as it is.
  nohash <- fit
  nohash$data_hash <- NA_character_
  ferx_save_fit(nohash, path, include_data = TRUE)
  expect_identical(fitrx_entry_bytes(path, "data.csv"), dh_bytes(data))
})

test_that("a data.frame fit bundles its own rows even if the cached file changed (#419, #585)", {
  # The data.frame is written to a session-cache CSV the user never edits, so
  # the file normally hashes to fit$data_hash. Should it not, the fit's own
  # `data_source` text - the bytes that were fitted - is bundled instead.
  skip_on_cran()
  ex <- ferx_example("warfarin")
  df <- utils::read.csv(ex$data)
  fit <- dh_fit(ex$model, df)
  src <- charToRaw(fit$data_source)
  expect_identical(fit$data_hash, ferx_rust_sha256_raw(src))
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  expect_identical(fitrx_entry_bytes(path, "data.csv"), src)

  cat("garbage\n", file = fit$data_path, append = TRUE)
  withr::defer(unlink(fit$data_path))
  ferx_save_fit(fit, path)
  expect_identical(fitrx_entry_bytes(path, "data.csv"), src)
})

test_that("the staged copy of the data is owner-only, in R's session tempdir, and removed (#585 review)", {
  # The glue parses the population from a copy of the bytes it read. The
  # dataset may be owner-only, so the copy must not be world-readable in a
  # shared temporary directory, nor outlive the read.
  bytes <- as.raw(c(0x49, 0x44, 0x2c, 0xff, 0x0a))  # not valid UTF-8
  probe <- ferx_rust_stage_data_probe(bytes)
  expect_identical(probe$bytes, bytes)
  expect_identical(normalizePath(dirname(probe$path)), normalizePath(tempdir()))
  expect_false(file.exists(probe$path))
  if (.Platform$OS.type == "unix") {
    expect_identical(probe$mode, strtoi("600", 8L))
  }
})

test_that("a fit leaves no staged data copy behind (#585 review)", {
  skip_on_cran()
  ex <- ferx_example("warfarin")
  before <- list.files(tempdir(), pattern = "^ferx-data-")
  fit <- dh_fit(ex$model, ex$data)
  expect_identical(list.files(tempdir(), pattern = "^ferx-data-"), before)
})

test_that("a data.frame fit saves from data_source when its cache file is unreadable (#585 review)", {
  skip_on_cran()
  skip_on_os("windows")
  ex <- ferx_example("warfarin")
  fit <- dh_fit(ex$model, utils::read.csv(ex$data))
  src <- charToRaw(fit$data_source)
  Sys.chmod(fit$data_path, "000")
  withr::defer({
    Sys.chmod(fit$data_path, "644")
    unlink(fit$data_path)
  })
  skip_if(file.access(fit$data_path, 4L) == 0L, "running as root")
  path <- withr::local_tempfile(fileext = ".fitrx")
  expect_no_error(ferx_save_fit(fit, path))
  expect_identical(fitrx_entry_bytes(path, "data.csv"), src)
})
