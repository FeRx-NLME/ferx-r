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

# Resolve one variable the way make itself would, with `env` (a named character
# vector) added to the environment of the make process. The variables are set
# on this process and restored afterwards rather than passed to `system2(env=)`,
# which is not portable to Windows and does not quote values containing spaces.
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
    set <- Sys.getenv(names(env), unset = NA_character_, names = TRUE)
    do.call(Sys.setenv, as.list(env))
    on.exit({
      keep <- !is.na(set)
      if (any(keep)) do.call(Sys.setenv, as.list(set[keep]))
      if (any(!keep)) Sys.unsetenv(names(set)[!keep])
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

skip_unless_make_probe_works <- function() {
  mk <- makevars_path()
  skip_if(is.na(mk), "package source src/Makevars not available")
  skip_if(!nzchar(make_program()), "no make on PATH")
  # Nothing here should fail merely because make could not parse the probe.
  skip_if(
    is.na(resolve_make_var(mk, "STATLIB")),
    "make could not evaluate src/Makevars"
  )
  mk
}

test_that("CARGO_FEATURES defaults to the shipped feature set", {
  mk <- skip_unless_make_probe_works()
  expect_equal(
    resolve_make_var(mk, "CARGO_FEATURES"),
    "--no-default-features --features ci,nn,survival"
  )
})

test_that("CARGO_FEATURES can be overridden from the environment without -e", {
  mk <- skip_unless_make_probe_works()
  without_nn <- "--no-default-features --features ci,survival"
  expect_equal(
    resolve_make_var(
      mk, "CARGO_FEATURES",
      env = c(CARGO_FEATURES = without_nn)
    ),
    without_nn
  )
})

test_that("CARGO_PROFILE_RELEASE_LTO stays thin whatever the environment says", {
  mk <- skip_unless_make_probe_works()
  expect_equal(
    resolve_make_var(
      mk, "CARGO_PROFILE_RELEASE_LTO",
      env = c(CARGO_PROFILE_RELEASE_LTO = "fat", MAKEFLAGS = "-e")
    ),
    "thin"
  )
})
