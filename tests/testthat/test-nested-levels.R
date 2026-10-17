# Nested variability levels (ferx-core #1082 PR3 / PR4) on ferx-core's exact3
# anchor (`nonmem_anchor/nested_levels_3.*`, copied to fixtures/): 12 sponsors
# (the subject, ID), two studies each, two arms per study, ARM numbered 1, 2
# inside each study (nested by key). The model's initial estimates are
# NONMEM's converged ones, so `maxiter = 0` evaluates the fit at one point.
#
# The oracles are R-side and independent of the glue's own conversion: the
# per-level estimates against the per-cell `ebe_kappas` and the sdtab columns,
# shrinkage against its definition over each level's groups, the printed KAPPA
# rows against ferx-core's own loader of the saved bundle, and the simulated
# draws against a degenerate design where only the study kappa varies.

nested_model <- test_path("fixtures", "nested_levels_3.ferx")
nested_data <- test_path("fixtures", "nested_levels_3.csv")

nested_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      fit <<- suppressWarnings(ferx_fit(nested_model, nested_data,
                                        method = "focei", covariance = FALSE,
                                        settings = list(maxiter = 0L)))
    }
    fit
  }
})

test_that("a nested fit reports each kappa's level and each level's estimates once per group", {
  fit <- nested_fit()
  expect_identical(fit$kappa_levels, c(K_ST = "STUDY", K_ARM = "ARM"))
  by_level <- fit$ebe_kappas_by_level
  expect_named(by_level, c("STUDY", "ARM"))
  expect_named(by_level$STUDY, c("ID", "GROUP", "STUDY", "K_ST"))
  expect_named(by_level$ARM, c("ID", "GROUP", "STUDY", "ARM", "K_ARM"))
  expect_type(by_level$STUDY$ID, "character")
  expect_type(by_level$STUDY$GROUP, "integer")
  expect_type(by_level$ARM$ARM, "double")
  # 24 studies and 48 arms, each once: not one study row per arm.
  expect_identical(nrow(by_level$STUDY), 24L)
  expect_identical(nrow(by_level$ARM), 48L)
  expect_identical(anyDuplicated(by_level$STUDY[c("ID", "STUDY")]), 0L)
  expect_identical(anyDuplicated(by_level$ARM[c("ID", "STUDY", "ARM")]), 0L)
  expect_identical(by_level$ARM$GROUP, rep(1:4, 12L))

  # Each group's estimate is the one the per-cell `ebe_kappas` repeats: the
  # cells come in first-record order, as the arms do.
  cells <- fit$ebe_kappas
  expect_identical(nrow(cells), 48L)
  expect_identical(by_level$ARM$K_ARM, cells$K_ARM)
  st <- merge(by_level$ARM[c("ID", "STUDY", "GROUP")], by_level$STUDY[c("ID", "STUDY", "K_ST")],
              by = c("ID", "STUDY"), sort = FALSE)
  st <- st[order(match(st$ID, unique(cells$ID)), st$GROUP), ]
  expect_identical(st$K_ST, cells$K_ST)
})

test_that("nested kappa shrinkage is pooled over the groups of the kappa's own level", {
  fit <- nested_fit()
  for (k in c("K_ST", "K_ARM")) {
    lv <- fit$kappa_levels[[k]]
    v <- fit$ebe_kappas_by_level[[lv]][[k]]
    want <- 1 - sqrt(mean(v^2)) / sqrt(fit$omega_iov[k, k])
    expect_equal(unname(fit$shrinkage_kappa[[k]]), want, tolerance = 1e-12, info = k)
  }
})

test_that("the sdtab carries each row's own study and arm kappa", {
  fit <- nested_fit()
  data <- utils::read.csv(nested_data)
  sd <- fit$sdtab
  expect_true(all(c("K_ST", "K_ARM") %in% names(sd)))
  expect_identical(nrow(sd), nrow(data))
  expect_identical(as.numeric(sd$ID), as.numeric(data$ID))
  key <- function(...) paste(..., sep = "/")
  st <- fit$ebe_kappas_by_level$STUDY
  arm <- fit$ebe_kappas_by_level$ARM
  expect_identical(sd$K_ST, st$K_ST[match(key(data$ID, data$STUDY), key(st$ID, st$STUDY))])
  expect_identical(sd$K_ARM, arm$K_ARM[match(key(data$ID, data$STUDY, data$ARM),
                                             key(arm$ID, arm$STUDY, arm$ARM))])
})

test_that("print() labels a nested fit's kappa rows and shrinkage with their level", {
  fit <- nested_fit()
  kappa_rows <- getFromNamespace(".ferx_kappa_rows", "ferx")
  fitrx_kappa_rows <- getFromNamespace("ferx_rust_fitrx_kappa_rows", "ferx")
  rows <- kappa_rows(fit, ascii = FALSE)
  expect_true(any(grepl("K_ST at STUDY", rows, fixed = TRUE)))
  expect_true(any(grepl("K_ARM at ARM", rows, fixed = TRUE)))
  # ferx-core's own console rows, from its loader of the R-written bundle:
  # the level reaches the engine through iov.kappa_levels.
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  expect_identical(rows, strsplit(fitrx_kappa_rows(f), "\n", fixed = TRUE)[[1L]])

  out <- utils::capture.output(print(fit))
  sh <- out[grep("SHRINKAGE", out, fixed = TRUE) + 2L]
  expect_match(sh, "K_ST at STUDY: ", fixed = TRUE)
  expect_match(sh, "K_ARM at ARM: ", fixed = TRUE)
})

test_that("a one-level IOV fit has no level fields and prints no level", {
  ex <- ferx_example("warfarin_iov")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, method = "foce",
                                   covariance = FALSE, settings = list(maxiter = 0L)))
  expect_null(fit$kappa_levels)
  expect_null(fit$ebe_kappas_by_level)
  expect_false(any(grepl(" at ", getFromNamespace(".ferx_kappa_rows", "ferx")(fit), fixed = TRUE)))
})

test_that("a nested fit's levels survive ferx_save_fit() / ferx_load_fit()", {
  fit <- nested_fit()
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  entries <- utils::unzip(f, list = TRUE)$Name
  expect_true(all(c("ebes_kappa_STUDY.csv", "ebes_kappa_ARM.csv") %in% entries))
  dir <- withr::local_tempdir()
  utils::unzip(f, files = "fit.json", exdir = dir)
  json <- jsonlite::read_json(file.path(dir, "fit.json"))
  expect_identical(unlist(json$iov$kappa_levels), c("STUDY", "ARM"))
  expect_identical(vapply(json$iov$levels, `[[`, "", "name"), c("STUDY", "ARM"))
  expect_identical(unlist(json$iov$levels[[2]]$columns), c("STUDY", "ARM"))
  expect_identical(unlist(json$iov$levels[[2]]$kappa_names), "K_ARM")

  back <- suppressWarnings(ferx_load_fit(f))
  expect_identical(back$kappa_levels, fit$kappa_levels)
  expect_identical(back$ebe_kappas_by_level, fit$ebe_kappas_by_level)
  expect_identical(getFromNamespace(".ferx_kappa_rows", "ferx")(back, ascii = FALSE),
                   getFromNamespace(".ferx_kappa_rows", "ferx")(fit, ascii = FALSE))
})

test_that("a bundle naming a level without its CSV is refused, as ferx-core refuses it", {
  fit <- nested_fit()
  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  dir <- withr::local_tempdir()
  utils::unzip(f, exdir = dir)
  unlink(file.path(dir, "ebes_kappa_ARM.csv"))
  g <- withr::local_tempfile(fileext = ".fitrx")
  withr::with_dir(dir, utils::zip(g, list.files(), flags = "-q"))
  expect_error(suppressWarnings(ferx_load_fit(g)),
               "names level `ARM`, but the bundle has no ebes_kappa_ARM.csv", fixed = TRUE)
})

test_that("ferx_simulate() draws a nested model's study kappa once per study", {
  # Degenerate design: no subject eta, no arm kappa, so EMAX varies only by
  # the study kappa. It is then one value on every record of a study, across
  # both of its arms, and differs between studies.
  txt <- readLines(nested_model)
  txt <- sub("omega ETA_SP ~ 0.0699388", "omega ETA_SP ~ 0.0", txt, fixed = TRUE)
  txt <- sub("kappa K_ARM ~ 3.75478 at ARM weight = NARM", "kappa K_ARM ~ 0.0 at ARM", txt,
             fixed = TRUE)
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(txt, model)
  data <- utils::read.csv(nested_data)
  sim <- suppressWarnings(ferx_simulate(model, nested_data, n_sim = 2L, seed = 7L))
  for (s in 1:2) {
    one <- sim[sim$SIM == s, , drop = FALSE]
    expect_identical(nrow(one), nrow(data))
    expect_identical(as.numeric(one$TIME), as.numeric(data$TIME))
    f <- data$DOSE / (33.4048 + data$DOSE) * data$TIME / (data$TIME + 3.20150)
    emax <- (one$IPRED - 10.1818) / f
    by_study <- split(emax, paste(data$ID, data$STUDY))
    spread <- vapply(by_study, function(v) diff(range(v)) / mean(v), 0)
    expect_lt(max(spread), 1e-10)
    expect_identical(length(unique(signif(vapply(by_study, mean, 0), 8))), 24L)
  }
  # Seeded: the same draws again.
  expect_identical(sim, suppressWarnings(ferx_simulate(model, nested_data, n_sim = 2L,
                                                       seed = 7L)))
})

test_that("ferx_simulate(fit = ) simulates a nested fit", {
  fit <- nested_fit()
  sim <- suppressWarnings(ferx_simulate(nested_model, nested_data, n_sim = 1L, seed = 3L,
                                        fit = fit))
  expect_identical(nrow(sim), 240L)
  expect_true(all(is.finite(sim$IPRED)))
})
