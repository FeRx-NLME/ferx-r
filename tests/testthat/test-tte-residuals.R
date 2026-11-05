# Time-to-event residual rows in fit$sdtab (ferx-core #901), and what
# ferx_save_fit() does with them (ferx-core #911: a bundle keeps none).

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
  expect_no_warning(ferx_save_fit(fit, path),
                    class = "ferx_tte_rows_dropped_warning")
  loaded <- suppressWarnings(ferx_load_fit(path))
  expect_identical(loaded$sdtab, fit$sdtab)
  expect_identical(loaded$subject_summary, fit$subject_summary)
  expect_gt(sum(loaded$subject_summary$n_obs_tte), 0L)
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
  dropped <- 0L
  withCallingHandlers(
    ferx_save_fit(fit, path),
    ferx_tte_rows_dropped_warning = function(w) dropped <<- dropped + 1L,
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_identical(dropped, 0L)
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
