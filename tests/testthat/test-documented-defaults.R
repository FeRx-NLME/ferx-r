# The defaults ?ferx_fit documents for its `settings` keys must be the
# engine's own (#510). The engine's default scoring / SIR records
# (ferx-core `ScoringSettings::default()` / `SirSettings::default()`, both
# built from `FitOptions::default()`) are the reference, so a default the
# engine changes fails here until the help page follows.

dd_rd_text <- function() {
  src <- testthat::test_path("..", "..", "man", "ferx_fit.Rd")
  lines <- if (file.exists(src)) {
    readLines(src, warn = FALSE)
  } else {
    db <- tools::Rd_db("ferx")
    strsplit(paste(as.character(db[["ferx_fit.Rd"]]), collapse = ""), "\n")[[1]]
  }
  gsub("\\s+", " ", paste(lines, collapse = " "))
}

# The text of the `\item{\code{key}}{...}` entry, up to the next item.
dd_item <- function(rd, key) {
  start <- regexpr(paste0("\\item{\\code{", key, "}}{"), rd, fixed = TRUE)
  if (start < 0L) stop("no ?ferx_fit item for `", key, "`")
  body <- substring(rd, start + attr(start, "match.length"))
  end <- regexpr("\\item{", body, fixed = TRUE)
  if (end > 0L) substring(body, 1L, end - 1L) else body
}

# The documented default of `key`, typed: `\code{"auto"} (default)` and
# `(default \code{"packed"})` give a string, `TRUE` / `FALSE` a logical, and
# anything else a number. The earliest of the two spellings wins. The first
# needs a bare `(default)`, so that `as for \code{ode_reltol} (default
# \code{1e-6})` reads as 1e-6 rather than as the name of another key.
dd_documented_default <- function(rd, key) {
  item <- dd_item(rd, key)
  patterns <- c(
    "\\(default (?:\\\\code\\{)?([^)};,]+)",
    "\\\\code\\{([^}]+)\\} \\(default\\)"
  )
  hits <- lapply(patterns, function(p) regexec(p, item, perl = TRUE)[[1]])
  pos <- vapply(hits, function(h) as.numeric(h[1L]), numeric(1))
  if (all(pos < 0)) stop("?ferx_fit item `", key, "` states no default")
  best <- hits[[which(pos == min(pos[pos > 0]))[1L]]]
  value <- trimws(substring(item, best[2L], best[2L] + attr(best, "match.length")[2L] - 1L))
  if (grepl('^".*"$', value)) {
    gsub('"', "", value, fixed = TRUE)
  } else if (value %in% c("TRUE", "FALSE")) {
    as.logical(value)
  } else {
    as.numeric(value)
  }
}

dd_expect_default <- function(rd, key, engine) {
  documented <- dd_documented_default(rd, key)
  if (is.numeric(documented)) engine <- as.numeric(engine)
  expect_identical(documented, engine, info = paste0("settings key `", key, "`"))
}

test_that("?ferx_fit documents ode_method's engine default, \"auto\" (#510)", {
  rd <- dd_rd_text()
  engine <- ferx:::.ferx_default_settings("scoring_settings")$ode_method
  expect_identical(engine, "auto")
  expect_identical(dd_documented_default(rd, "ode_method"), engine)
})

test_that("?ferx_fit's scoring-setting defaults are the engine's (#510)", {
  rd <- dd_rd_text()
  engine <- ferx:::.ferx_default_settings("scoring_settings")
  for (key in c("ode_method", "ode_reltol", "ode_abstol", "ode_max_steps",
                "inner_maxiter", "inner_tol", "inner_restarts",
                "inner_optimizer", "ebe_warm_start")) {
    dd_expect_default(rd, key, engine[[key]])
  }
})

test_that("?ferx_fit's SIR defaults are the engine's (#510)", {
  rd <- dd_rd_text()
  engine <- ferx:::.ferx_default_settings("sir_settings")
  for (field in c("samples", "resamples", "df", "keep_samples", "scale")) {
    dd_expect_default(rd, paste0("sir_", field), engine[[field]])
  }
})

test_that("the README gives the engine's default ode_method (#510)", {
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md is not shipped with the package")
  text <- readLines(readme, warn = FALSE)
  engine <- ferx:::.ferx_default_settings("scoring_settings")$ode_method
  ode_line <- grep("ode_method", text, value = TRUE, fixed = TRUE)
  expect_length(ode_line, 1L)
  expect_true(grepl(sprintf('`ode_method = "%s"` (default)', engine), ode_line, fixed = TRUE))
  expect_false(any(grepl("RK45 (default)", text, fixed = TRUE)))
})
