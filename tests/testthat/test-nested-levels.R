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

# The fit's bundle with `edit(dir)` applied to its unpacked entries, re-zipped.
tampered_bundle <- function(fit, edit, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".fitrx", .local_envir = env)
  suppressWarnings(ferx_save_fit(fit, f))
  dir <- withr::local_tempdir(.local_envir = env)
  utils::unzip(f, exdir = dir)
  edit(dir)
  g <- withr::local_tempfile(fileext = ".fitrx", .local_envir = env)
  withr::with_dir(dir, utils::zip(g, list.files(), flags = "-q"))
  g
}

edit_lines <- function(dir, entry, fn) {
  path <- file.path(dir, entry)
  writeLines(fn(readLines(path)), path)
}

edit_fit_json <- function(dir, fn) {
  path <- file.path(dir, "fit.json")
  json <- jsonlite::read_json(path)
  jsonlite::write_json(fn(json), path, auto_unbox = TRUE, digits = NA, null = "null")
}

test_that("ferx_load_fit() refuses a malformed nested bundle as ferx-core's loader does", {
  fit <- nested_fit()
  refuses <- function(edit, msg) {
    expect_error(suppressWarnings(ferx_load_fit(tampered_bundle(fit, edit))), msg,
                 fixed = TRUE)
  }
  refuses(function(d) unlink(file.path(d, "ebes_kappa_ARM.csv")),
          "names level `ARM`, but the bundle has no ebes_kappa_ARM.csv")
  # #587 review, finding 3: a header that is not the level's.
  refuses(function(d) edit_lines(d, "ebes_kappa_ARM.csv",
                                 function(l) { l[1] <- sub("K_ARM", "K_X", l[1]); l }),
          "ebes_kappa_ARM.csv: header ID,GROUP,STUDY,ARM,K_X does not match fit.json's level")
  # Finding 4: what core's loader refuses beyond that.
  refuses(function(d) edit_fit_json(d, function(j) {
            j$iov$kappa_levels <- list("STUDY"); j
          }),
          "iov.kappa_levels (1) does not match iov.kappa_names (2)")
  refuses(function(d) edit_fit_json(d, function(j) {
            j$iov$kappa_levels <- list("STUDY", "STUDY"); j
          }),
          "iov.levels names kappa `K_ARM` at level `ARM`, but iov.kappa_levels puts it at `STUDY`")
  refuses(function(d) edit_lines(d, "ebes_kappa_ARM.csv",
                                 function(l) { l[2] <- paste0(l[2], ",1"); l }),
          "ebes_kappa_ARM.csv: row 1 has 6 fields, expected 5")
  refuses(function(d) edit_lines(d, "ebes_kappa_ARM.csv",
                                 function(l) { l[2] <- sub("^[^,]*", "999", l[2]); l }),
          "ebes_kappa_ARM.csv: unknown ID \"999\"")
  refuses(function(d) edit_lines(d, "ebes_kappa_ARM.csv",
                                 function(l) { l[3] <- sub("^([^,]*),2,", "\\1,3,", l[3]); l }),
          "ebes_kappa_ARM.csv: subject \"1\": GROUP must run 1..n without gaps or repeats")
})

test_that("an empty kappa cell in a level's CSV loads as a NaN estimate", {
  # ferx-core writes a NaN estimate as an empty cell and reads it back as NaN;
  # R's reader would make it NA (#587 review, finding 3).
  fit <- nested_fit()
  g <- tampered_bundle(fit, function(d) edit_lines(d, "ebes_kappa_ARM.csv", function(l) {
    l[2] <- sub("[^,]*$", "", l[2])
    l
  }))
  back <- suppressWarnings(ferx_load_fit(g))
  k <- back$ebe_kappas_by_level$ARM$K_ARM
  expect_true(is.nan(k[1]))
  expect_identical(k[-1], fit$ebe_kappas_by_level$ARM$K_ARM[-1])
})

test_that("two kappas at one level keep their columns apart, live and through .fitrx", {
  # One kappa per level cannot show a swapped pair (#587 review, finding 5):
  # a second study kappa on the baseline, with its own variance.
  txt <- readLines(nested_model)
  txt <- sub("kappa K_ST ~ 0.0904362 at STUDY",
             "kappa K_ST ~ 0.0904362 at STUDY\n  kappa K_E0 ~ 0.01 at STUDY", txt, fixed = TRUE)
  txt <- sub("EMAX = TVEMAX * exp(ETA_SP + K_ST + K_ARM)",
             "EMAX = TVEMAX * exp(ETA_SP + K_ST + K_ARM)\n  E0 = TVE0 * exp(K_E0)", txt,
             fixed = TRUE)
  txt <- sub("y = TVE0 + EMAX", "y = E0 + EMAX", txt, fixed = TRUE)
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(txt, model)
  fit <- suppressWarnings(ferx_fit(model, nested_data, method = "focei", covariance = FALSE,
                                   settings = list(maxiter = 0L)))
  expect_identical(fit$kappa_levels, c(K_ST = "STUDY", K_E0 = "STUDY", K_ARM = "ARM"))
  st <- fit$ebe_kappas_by_level$STUDY
  expect_named(st, c("ID", "GROUP", "STUDY", "K_ST", "K_E0"))
  # Each column against the per-cell table, by name: the first arm of each
  # study carries the study's kappas.
  cells <- fit$ebe_kappas
  first <- cells[c(TRUE, FALSE), ]
  expect_identical(st$K_ST, first$K_ST)
  expect_identical(st$K_E0, first$K_E0)
  expect_false(isTRUE(all.equal(st$K_ST, st$K_E0)))

  f <- withr::local_tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, f))
  back <- suppressWarnings(ferx_load_fit(f))
  expect_identical(back$ebe_kappas_by_level, fit$ebe_kappas_by_level)
  fitrx_kappa_rows <- getFromNamespace("ferx_rust_fitrx_kappa_rows", "ferx")
  expect_identical(getFromNamespace(".ferx_kappa_rows", "ferx")(fit, ascii = FALSE),
                   strsplit(fitrx_kappa_rows(f), "\n", fixed = TRUE)[[1L]])
})

test_that("the post-hoc steps run on a nested fit", {
  # NEWS says simulate-with-uncertainty and the NPDE reference simulate a
  # nested model (#587 review, finding 6); SIR, covariance and predict ride
  # the same fit.
  fit <- suppressWarnings(ferx_fit(nested_model, nested_data, method = "focei",
                                   covariance = TRUE, settings = list(maxiter = 0L)))
  expect_identical(fit$covariance_status, "computed")
  pred <- suppressWarnings(ferx_predict(nested_model, nested_data, fit = fit))
  expect_identical(nrow(pred), 240L)
  expect_true(all(is.finite(pred$PRED)))
  unc <- suppressWarnings(ferx_simulate_with_uncertainty(nested_model, nested_data, fit,
                                                         n_uncertainty_draws = 3L,
                                                         seed = 1L))
  expect_identical(nrow(unc), 3L * 240L)
  expect_true(all(is.finite(unc$IPRED)))
  npde <- suppressWarnings(ferx_calc_npde(fit, nsim = 50L, seed = 1L))
  expect_identical(nrow(npde$sdtab), 240L)
  expect_true(all(is.finite(npde$sdtab$NPDE)))
  cov <- suppressWarnings(ferx_covariance(fit))
  expect_true(all(is.finite(cov$se_kappa)))
  sir <- suppressWarnings(ferx_sir(fit, sir_samples = 200L, sir_resamples = 50L,
                                   sir_seed = 1L))
  expect_identical(dim(sir$sir_ci_kappa), c(2L, 2L))
})

# The bundled nested_levels example with a `theta PL[cols]` level block added
# (ferx-core #2034 / #2325). On the baseline it is an offset per level; in the
# exponent beside the kappas it is read like them.
nested_with_block <- function(cols, in_exponent = FALSE) {
  ex <- ferx_example("nested_levels")
  t <- readLines(ex$model)
  t <- append(t, sprintf("  theta PL[%s](0.0, -50.0, 50.0)", cols), grep("^  sigma ADD", t) - 1L)
  t <- if (in_exponent) {
    sub("exp(ETA_SP + K_ST + K_ARM)", "exp(ETA_SP + K_ST + K_ARM + PL)", t, fixed = TRUE)
  } else {
    sub("y = TVE0 +", "y = TVE0 + PL +", t, fixed = TRUE)
  }
  model <- withr::local_tempfile(fileext = ".ferx", .local_envir = parent.frame())
  writeLines(t, model)
  list(model = model, data = ex$data)
}

test_that("a theta level block fits beside kappas at nested levels (ferx-core #2034)", {
  skip_on_cran()
  # Each block's levels lie inside one study or one arm, which ferx-core
  # refused outright with E_NESTED_LEVEL_THETA_BLOCK_UNSUPPORTED before #2325.
  for (cols in c("STUDYID", "ARMID", "STUDYID, TIME")) {
    m <- nested_with_block(cols)
    fit <- suppressWarnings(ferx_fit(m$model, m$data, covariance = FALSE, verbose = FALSE))
    expect_true(fit$converged, label = cols)
    expect_true(is.finite(fit$ofv), label = cols)
    expect_true(nrow(fit$theta_levels) > 0L, label = cols)
  }
  # A block read like the kappas it lies inside absorbs them: refused with the
  # binder's own code, naming what it absorbs, not the lifted blanket refusal.
  m <- nested_with_block("ARMID", in_exponent = TRUE)
  e <- tryCatch(suppressWarnings(ferx_fit(m$model, m$data, covariance = FALSE,
                                          verbose = FALSE)),
                error = function(e) e)
  expect_s3_class(e, "error")
  expect_identical(e$code, "E_THETA_LEVEL_BINDING")
  expect_match(conditionMessage(e), "absorb `ETA_SP`, `K_ST` and `K_ARM`", fixed = TRUE)
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
