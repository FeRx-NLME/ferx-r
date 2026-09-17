# src/Makevars must let CARGO_FEATURES be replaced from the environment, so a
# constrained machine can build without `nn` (Deep Compartment Models / neural
# ODEs, the heaviest part of the dependency graph) without editing a tracked
# file (#288). A plain `=` assignment beats the environment - `R CMD INSTALL`
# does not run make with -e - so the check below runs make against a probe
# makefile that includes Makevars and prints the resolved value.

make_program <- function() {
  cmd <- Sys.getenv("MAKE", "make")
  # MAKE may carry flags (e.g. "make -j4"); only the program is needed here.
  Sys.which(strsplit(cmd, " ", fixed = TRUE)[[1]][1])
}

# testthat runs with the working directory at tests/testthat; the package
# source (and therefore src/Makevars) is absent from an installed package and
# from R CMD check's test directory.
makevars_path <- function() {
  p <- file.path("..", "..", "src", "Makevars")
  if (file.exists(p)) normalizePath(p, winslash = "/") else NA_character_
}

# Resolve one variable the way make itself would, under the environment `env`
# (a named character vector; `NA` means "unset this one"). The variables are
# set on this process and restored afterwards rather than passed to
# `system2(env=)`, which is not portable to Windows and does not quote values
# containing spaces.
#
# Every caller states the full environment it needs, including the variables it
# needs *absent*: whatever the person running the suite exported would otherwise
# leak in and decide the answer. `CARGO_FEATURES` set in the caller's shell -
# the very thing this file is about, and documented in the README - would make
# the default case read back that value, and an inherited `MAKEFLAGS=-e` makes
# the environment beat a makefile `=` too, so the override case would pass
# against the old, broken assignment.
resolve_make_var <- function(makevars, var, env = character()) {
  dir <- tempfile("ferx-makevars-probe")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  probe <- file.path(dir, "probe.mk")
  writeLines(
    c(
      sprintf("ferx-print-var:\n\t@echo $(%s)", var),
      sprintf("include %s", makevars)
    ),
    probe
  )
  if (length(env)) {
    before <- Sys.getenv(names(env), unset = NA_character_, names = TRUE)
    wanted <- !is.na(env)
    if (any(wanted)) do.call(Sys.setenv, as.list(env[wanted]))
    if (any(!wanted)) Sys.unsetenv(names(env)[!wanted])
    on.exit({
      had <- !is.na(before)
      if (any(had)) do.call(Sys.setenv, as.list(before[had]))
      if (any(!had)) Sys.unsetenv(names(before)[!had])
    }, add = TRUE)
  }
  out <- suppressWarnings(system2(
    make_program(),
    c("-s", "-f", shQuote(normalizePath(probe, winslash = "/")), "ferx-print-var"),
    stdout = TRUE, stderr = FALSE
  ))
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    return(NA_character_)
  }
  trimws(paste(out, collapse = " "))
}

# The environment every probe starts from: the three variables that can change
# what the makefile resolves to are unset, and a case that wants one of them
# names it. Spelling the absences out is what keeps the suite's verdict
# independent of the shell it was launched from.
clean_env <- function(...) {
  base <- c(
    CARGO_FEATURES = NA_character_,
    CARGO_PROFILE_RELEASE_LTO = NA_character_,
    MAKEFLAGS = NA_character_
  )
  named <- c(...)
  if (length(named)) base[names(named)] <- named
  base
}

skip_unless_make_probe_works <- function() {
  mk <- makevars_path()
  skip_if(is.na(mk), "package source src/Makevars not available")
  skip_if(!nzchar(make_program()), "no make on PATH")
  # Nothing here should fail merely because make could not parse the probe.
  skip_if(
    is.na(resolve_make_var(mk, "STATLIB", env = clean_env())),
    "make could not evaluate src/Makevars"
  )
  mk
}

test_that("CARGO_FEATURES defaults to the shipped feature set", {
  mk <- skip_unless_make_probe_works()
  expect_equal(
    resolve_make_var(mk, "CARGO_FEATURES", env = clean_env()),
    "--no-default-features --features ci,nn,survival"
  )
})

test_that("CARGO_FEATURES can be overridden from the environment without -e", {
  mk <- skip_unless_make_probe_works()
  without_nn <- "--no-default-features --features ci,survival"
  expect_equal(
    resolve_make_var(
      mk, "CARGO_FEATURES",
      # MAKEFLAGS stays unset: under -e the environment beats a makefile `=`
      # too, so this case would pass against the assignment it exists to
      # reject.
      env = clean_env(CARGO_FEATURES = without_nn)
    ),
    without_nn
  )
})

test_that("CARGO_PROFILE_RELEASE_LTO stays thin whatever the environment says", {
  mk <- skip_unless_make_probe_works()
  expect_equal(
    resolve_make_var(
      mk, "CARGO_PROFILE_RELEASE_LTO",
      env = clean_env(CARGO_PROFILE_RELEASE_LTO = "fat", MAKEFLAGS = "-e")
    ),
    "thin"
  )
})
