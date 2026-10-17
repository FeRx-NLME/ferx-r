# .fitrx CSV doubles written exactly by both writers (ferx-core #1847, #1631).
#
# ferx_save_fit() has written every CSV double as the shortest text that reads
# back as the same value since #415; ferx-core's own writer wrote them at six
# decimals until #1847 and now writes its own shortest round-trip spelling
# (`1e-7`, `0.0`). The oracle for the engine's half is ferx-core's own
# `load_fit` / `save_fit` (`fitrx_engine_trip()`), never ferx_load_fit()
# reading what ferx_save_fit() wrote. ferx_covariance() is held to the inline
# covariance step bit for bit, which #1847 extends to every engine: the warm-
# started EBEs of gn / SAEM are no longer replaced by a cold solve.

nested_model <- test_path("fixtures", "nested_levels_3.ferx")
nested_data <- test_path("fixtures", "nested_levels_3.csv")

quiet <- function(x) suppressWarnings(suppressMessages(x))

# Every standard error and the covariance matrix, unnamed.
cov_parts <- function(fit) {
  lapply(fit[c("se_theta", "se_omega", "se_sigma", "se_kappa", "cov_matrix")], unname)
}

# The lines of one entry of the bundle at `path`.
bundle_entry <- function(path, entry) {
  dir <- withr::local_tempdir()
  utils::unzip(path, files = entry, exdir = dir)
  readLines(file.path(dir, entry))
}

# A copy of the bundle at `path` with every decimal cell of its fit-output CSVs
# (EBEs, kappas, predictions, conditional distribution) rewritten at six
# decimals, as ferx-core's writer spelled them before #1847. The bundled data
# is left alone: its hash is checked against the fit's.
six_decimal_bundle <- function(path, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  utils::unzip(path, exdir = dir)
  csvs <- list.files(dir, pattern = "^(ebes.*|predictions|conddist)\\.csv$",
                     full.names = TRUE)
  for (f in csvs) {
    l <- readLines(f)
    l[-1] <- vapply(strsplit(l[-1], ",", fixed = TRUE), function(r) {
      dec <- grepl("^-?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?$|^-?[0-9]+[eE][-+]?[0-9]+$", r)
      r[dec] <- sprintf("%.6f", as.numeric(r[dec]))
      paste(r, collapse = ",")
    }, "")
    writeLines(l, f)
  }
  out <- withr::local_tempfile(fileext = ".fitrx", .local_envir = env)
  withr::with_dir(dir, utils::zip(out, list.files(), flags = "-q"))
  out
}

nested_cov_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- quiet(ferx_fit(nested_model, nested_data, method = "focei",
                             covariance = TRUE, settings = list(maxiter = 0L)))
    }
    fit
  }
})

warfarin_gn_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("warfarin")
      fit <<- quiet(ferx_fit(ex$model, ex$data, method = "gn", covariance = TRUE,
                             verbose = FALSE))
    }
    fit
  }
})

test_that("a nested fit's bundle reads back exactly through ferx-core's writer", {
  fit <- nested_cov_fit()
  expect_identical(fit$covariance_status, "computed")
  trip <- quiet(fitrx_engine_trip(fit))
  # ferx-core's writer spells an EBE and a kappa in full now, not at six
  # decimals.
  ebes <- bundle_entry(trip$core_path, "ebes.csv")
  eta_cells <- vapply(strsplit(ebes[-1], ",", fixed = TRUE), `[`, "", 2L)
  expect_true(any(nchar(sub("^-?[0-9]*\\.", "", eta_cells)) > 6L))
  arm <- bundle_entry(trip$core_path, "ebes_kappa_ARM.csv")
  kappa_cells <- vapply(strsplit(arm[-1], ",", fixed = TRUE), `[`, "", 5L)
  expect_true(any(nchar(sub("^-?[0-9]*\\.", "", kappa_cells)) > 6L))

  for (path in c(trip$r_path, trip$core_path)) {
    back <- quiet(ferx_load_fit(path))
    expect_identical(back$ebe_kappas_by_level, fit$ebe_kappas_by_level)
    expect_identical(back$ebe_kappas, fit$ebe_kappas)
    expect_identical(back$ebe_etas, fit$ebe_etas)
    expect_identical(cov_parts(quiet(ferx_covariance(back))), cov_parts(fit))
  }
  expect_identical(cov_parts(quiet(ferx_covariance(fit))), cov_parts(fit))
})

test_that("ferx_covariance() repeats a gn fit's inline step bit for bit, live and reloaded", {
  # Gauss-Newton reports the EBEs of its own warm-started solves, about 1e-7
  # off a cold one; before #1847 run_covariance restarted from the cold solve.
  fit <- warfarin_gn_fit()
  expect_identical(fit$covariance_status, "computed")
  expect_identical(cov_parts(quiet(ferx_covariance(fit))), cov_parts(fit))
  trip <- quiet(fitrx_engine_trip(fit))
  for (path in c(trip$r_path, trip$core_path)) {
    back <- quiet(ferx_load_fit(path))
    expect_identical(back$ebe_etas, fit$ebe_etas)
    expect_identical(back$sdtab$IPRED, fit$sdtab$IPRED)
    expect_identical(cov_parts(quiet(ferx_covariance(back))), cov_parts(fit))
  }
})

test_that("conddist.csv crosses both writers: COND_MODE is checked by value", {
  # ferx-core compared COND_MODE with the EBE as six-decimal text, so it
  # refused every R-written conddist.csv whose EBEs were not six-decimal
  # numbers (#1847); fitrx_engine_trip() loads the R bundle with core's loader.
  fit <- warfarin_saem_conddist_fit()
  # Without a computed step the identity below would compare NULL with NULL.
  expect_identical(fit$covariance_status, "computed")
  trip <- quiet(fitrx_engine_trip(fit))
  expect_true("conddist.csv" %in% utils::unzip(trip$core_path, list = TRUE)$Name)
  back <- quiet(ferx_load_fit(trip$core_path))
  expect_identical(back$cond_dist$data, fit$cond_dist$data)
  expect_identical(back$ebe_etas, fit$ebe_etas)
  expect_identical(cov_parts(quiet(ferx_covariance(back))), cov_parts(fit))
})

test_that("a bundle written at six decimals before #1847 still loads", {
  fit <- nested_cov_fit()
  f <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(fit, f, include_data = TRUE))
  old <- six_decimal_bundle(f)
  expect_match(bundle_entry(old, "ebes.csv")[2], "^1,-?[0-9]+\\.[0-9]{6},")
  # ferx-core's loader takes it too.
  expect_no_error(quiet(ferx:::ferx_rust_fitrx_engine_resave(
    old, withr::local_tempfile(fileext = ".fitrx"))))
  back <- quiet(ferx_load_fit(old))
  expect_equal(back$ebe_kappas_by_level, fit$ebe_kappas_by_level, tolerance = 1e-5)
  expect_false(identical(back$ebe_etas, fit$ebe_etas))
  # Six-decimal EBEs within 1e-6 of the cold solve run from the cold solve,
  # which is what the inline FOCEI step started from.
  expect_identical(cov_parts(quiet(ferx_covariance(back))), cov_parts(fit))

  # gn's inline step started from its warm EBEs, which the bundle rounded
  # away: the reloaded step is close, not exact.
  gn <- warfarin_gn_fit()
  g <- withr::local_tempfile(fileext = ".fitrx")
  quiet(ferx_save_fit(gn, g, include_data = TRUE))
  gn_back <- quiet(ferx_load_fit(six_decimal_bundle(g)))
  gn_cov <- quiet(ferx_covariance(gn_back))
  expect_equal(gn_cov$se_theta, gn$se_theta, tolerance = 1e-6)
  expect_equal(gn_cov$se_omega, gn$se_omega, tolerance = 1e-6)
  # ...and not the warm path: from the rounded EBEs the step runs cold.
  expect_false(identical(cov_parts(gn_cov), cov_parts(gn)))
})
