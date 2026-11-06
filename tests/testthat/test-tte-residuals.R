# Time-to-event residual rows in fit$sdtab (ferx-core #901), and what
# ferx_save_fit() does with them (ferx-core #911: core keeps them in
# predictions_tte.csv, the glue in sdtab_endpoint_rows.csv).

tte_cols <- c("COXSNELL", "MARTINGALE", "DEVIANCE", "PIT")

joint_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("pktte_joint")
      fit <<- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                        verbose = FALSE))
    }
    fit
  }
})

test_that("a joint PK-TTE fit's sdtab carries one residual row per TTE record, last", {
  fit <- joint_fit()
  sdtab <- fit$sdtab
  summ <- fit$subject_summary
  expect_true(all(tte_cols %in% names(sdtab)))
  expect_identical(summ$n_obs_tte, rep(1L, nrow(summ)))
  expect_identical(nrow(sdtab),
                   sum(summ$n_obs) + sum(summ$n_obs_discrete) + sum(summ$n_obs_tte))

  n_tte <- sum(summ$n_obs_tte)
  tte <- sdtab[nrow(sdtab) - rev(seq_len(n_tte)) + 1L, , drop = FALSE]
  gauss <- sdtab[seq_len(nrow(sdtab) - n_tte), , drop = FALSE]
  # Subject order, on the endpoint's CMT, no concentration prediction.
  expect_identical(tte$ID, as.numeric(summ$ID))
  expect_true(all(tte$CMT == 3))
  expect_true(all(is.na(tte$PRED)))
  expect_true(all(is.na(as.matrix(gauss[tte_cols]))))
  # The residuals' definitions, as ferx-core documents them.
  expect_equal(tte$MARTINGALE, tte$DV - tte$COXSNELL)
  expect_equal(tte$PIT, exp(-tte$COXSNELL))
  expect_identical(sign(tte$DEVIANCE), sign(tte$MARTINGALE))
})

test_that("a fit without a TTE endpoint has no TTE columns and n_obs_tte = 0", {
  ex <- ferx_example("warfarin")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                   verbose = FALSE))
  expect_false(any(tte_cols %in% names(fit$sdtab)))
  expect_identical(fit$subject_summary$n_obs_tte,
                   integer(nrow(fit$subject_summary)))
})

test_that("ferx_save_fit() keeps the TTE rows through a reload (#689, #711)", {
  fit <- joint_fit()
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  loaded <- suppressWarnings(ferx_load_fit(path))
  # Core restores the rows from its own entry (ferx-core #911) and the glue's
  # entry carries them too: each row counts once.
  expect_identical(loaded$sdtab, fit$sdtab)
  expect_identical(loaded$subject_summary, fit$subject_summary)
  expect_gt(sum(loaded$subject_summary$n_obs_tte), 0L)
})

test_that("a bundle without the glue's endpoint entry rebuilds its TTE rows from core's", {
  skip_if(!nzchar(Sys.which("zip")), "needs a zip binary")
  fit <- joint_fit()
  path <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, path)
  # The TTE rows only in core's predictions_tte.csv (ferx-core #911), as in a
  # bundle written by ferx-core alone. This bundle keeps the R package's wide
  # predictions.csv, so it already has CMT and the TTE columns: the branch
  # that adds them for core's narrow predictions.csv is not reached here.
  status <- utils::zip(path, "sdtab_endpoint_rows.csv", flags = "-qd")
  skip_if(status != 0L, "zip -d failed")
  expect_false("sdtab_endpoint_rows.csv" %in% utils::unzip(path, list = TRUE)$Name)

  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(loaded$subject_summary$n_obs_tte, fit$subject_summary$n_obs_tte)
  expect_identical(nrow(loaded$sdtab), nrow(fit$sdtab))
  n_tte <- sum(fit$subject_summary$n_obs_tte)
  tte <- seq.int(nrow(fit$sdtab) - n_tte + 1L, nrow(fit$sdtab))
  # Every column of the TTE block, EBE_OFV included, as core's sdtab fills it.
  for (col in names(fit$sdtab)) {
    expect_identical(loaded$sdtab[[col]][tte], fit$sdtab[[col]][tte], info = col)
  }
})

test_that("a bundle ferx-core wrote alone loads its TTE rows, with CMT and the TTE columns", {
  skip_if(!nzchar(Sys.which("zip")), "needs a zip binary")
  fit <- joint_fit()
  saved <- tempfile(fileext = ".fitrx")
  ferx_save_fit(fit, saved)
  # Core's narrow predictions.csv: no CMT and none of the TTE columns, so the
  # glue adds both before it appends core's TTE rows.
  path <- fitrx_core_only_bundle(saved)
  expect_false(any(c("CMT", tte_cols) %in%
                     names(utils::read.csv(unz(path, "predictions.csv")))))

  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(loaded$subject_summary$n_obs_tte, fit$subject_summary$n_obs_tte)
  expect_identical(nrow(loaded$sdtab), nrow(fit$sdtab))
  expect_true(all(c("CMT", tte_cols, "EBE_OFV") %in% names(loaded$sdtab)))
  n_tte <- sum(fit$subject_summary$n_obs_tte)
  tte <- seq.int(nrow(fit$sdtab) - n_tte + 1L, nrow(fit$sdtab))
  # Every column the bundle carries, on the TTE rows, as the live fit has it.
  for (col in names(loaded$sdtab)) {
    expect_identical(loaded$sdtab[[col]][tte], fit$sdtab[[col]][tte], info = col)
  }
})

test_that("a non-TTE fit keeps an [output] column called PIT through save/load", {
  # The four TTE column names are dropped only with the TTE rows; on a fit
  # without them, a user column that shares a name is the user's.
  ex <- ferx_example("warfarin")
  model <- tempfile(fileext = ".ferx")
  src <- readLines(ex$model)
  src <- sub("^  KA = TVKA \\* exp\\(ETA_KA\\)$",
             "  KA = TVKA * exp(ETA_KA)\n  PIT = CL / V", src)
  writeLines(c(src, "", "[output]", "  PIT"), model)
  fit <- suppressWarnings(ferx_fit(model, ex$data, covariance = FALSE,
                                   verbose = FALSE))
  expect_true("PIT" %in% names(fit$sdtab))
  path <- tempfile(fileext = ".fitrx")
  suppressWarnings(ferx_save_fit(fit, path))
  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(loaded$sdtab$PIT, fit$sdtab$PIT)
})

test_that("ferx_calc_npde() on a joint PK-TTE fit fills the Gaussian rows, NA on TTE rows", {
  fit <- joint_fit()
  out <- suppressWarnings(ferx_calc_npde(fit, nsim = 50L, seed = 1L))
  n_tte <- sum(fit$subject_summary$n_obs_tte)
  tte <- seq.int(nrow(out$sdtab) - n_tte + 1L, nrow(out$sdtab))
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
  expect_true(all(is.na(out$sdtab$NPDE[tte])))
  expect_true(all(is.finite(out$sdtab$NPDE[-tte])))
  expect_true(all(is.finite(out$sdtab$NPD[-tte])))
})
