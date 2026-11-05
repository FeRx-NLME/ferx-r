# .fitrx_sdtab_blocks(): the one sdtab block split behind ferx_xpose(),
# ferx_save_fit() and ferx_calc_npde() (#696). Synthetic fits, since no
# bundled example has Gaussian, binary and TTE rows together (#689).

blocks_fn <- getFromNamespace(".fitrx_sdtab_blocks", "ferx")
tte_cols <- c("COXSNELL", "MARTINGALE", "DEVIANCE", "PIT")

# Two subjects as ferx-core writes them: Gaussian rows (N_OBS set), then
# discrete rows, then TTE rows (COXSNELL etc. set), each block in subject
# order. Subject 1: 2 Gaussian, 1 discrete, 1 TTE; subject 2: 1, 2, 1.
mixed_fit <- function(n_obs = c(2L, 1L), n_disc = c(1L, 2L), n_tte = c(1L, 1L)) {
  pos <- c(rep.int(1:2, n_obs), rep.int(1:2, n_disc), rep.int(1:2, n_tte))
  kind <- rep(c("g", "d", "t"), c(sum(n_obs), sum(n_disc), sum(n_tte)))
  sdtab <- data.frame(ID = c(11, 12)[pos], TIME = seq_along(pos), DV = 0,
                      CMT = 1, PRED = 1, IPRED = 1, IWRES = 0,
                      N_OBS = ifelse(kind == "g", n_obs[pos], NA))
  for (col in tte_cols) sdtab[[col]] <- ifelse(kind == "t", 0.5, NA_real_)
  list(sdtab = sdtab,
       subject_summary = data.frame(ID = c("11", "12"), ofv_contribution = 1,
                                    n_obs = n_obs, n_obs_discrete = n_disc,
                                    n_obs_tte = n_tte))
}

expect_blocks <- function(b, gaussian, discrete, tte) {
  expect_identical(b$gaussian, as.integer(gaussian))
  expect_identical(b$discrete, as.integer(discrete))
  expect_identical(b$tte, as.integer(tte))
}

test_that("a live Gaussian + binary + TTE fit splits by its counts", {
  b <- blocks_fn(mixed_fit())
  expect_blocks(b, 1:3, 4:6, 7:8)
  expect_identical(b$subject$gaussian, c(1L, 1L, 2L))
  expect_identical(b$subject$discrete, c(1L, 2L, 2L))
  expect_identical(b$subject$tte, c(1L, 2L))
  expect_identical(attr(b$subject, "n_subjects"), 2L)
})

test_that("without a subject_summary the columns give the same blocks", {
  fit <- mixed_fit()
  fit$subject_summary <- NULL
  b <- blocks_fn(fit)
  expect_blocks(b, 1:3, 4:6, 7:8)
  expect_null(b$subject)
})

test_that("counts that do not add up fall back to the columns", {
  # Missing n_obs_tte counts as 0, which leaves the two TTE rows uncounted.
  fit <- mixed_fit()
  fit$subject_summary$n_obs_tte <- NULL
  b <- blocks_fn(fit)
  expect_blocks(b, 1:3, 4:6, 7:8)
  expect_null(b$subject)
  # A count that is NA (a pre-#554 bundle's n_obs) does the same.
  fit <- mixed_fit()
  fit$subject_summary$n_obs[1] <- NA
  expect_null(blocks_fn(fit)$subject)
  expect_blocks(blocks_fn(fit), 1:3, 4:6, 7:8)
})

test_that("the counts win when they and the columns disagree", {
  # A TTE row without a COXSNELL value: by the columns it would be a discrete
  # row (no N_OBS), and a Gaussian row with no N_OBS a discrete one too.
  fit <- mixed_fit()
  fit$sdtab$COXSNELL[8] <- NA
  fit$sdtab$N_OBS[3] <- NA
  b <- blocks_fn(fit)
  expect_blocks(b, 1:3, 4:6, 7:8)
  expect_identical(b$subject$tte, c(1L, 2L))
  fit$subject_summary <- NULL
  expect_blocks(blocks_fn(fit), 1:2, c(3:6, 8L), 7L)
})

test_that("an old fit without n_obs_discrete / n_obs_tte uses n_obs alone", {
  fit <- mixed_fit(n_disc = c(0L, 0L), n_tte = c(0L, 0L))
  fit$sdtab <- fit$sdtab[, setdiff(names(fit$sdtab), tte_cols)]
  fit$subject_summary <- fit$subject_summary[, c("ID", "ofv_contribution", "n_obs")]
  b <- blocks_fn(fit)
  expect_blocks(b, 1:3, integer(), integer())
  expect_identical(b$subject$gaussian, c(1L, 1L, 2L))
})

test_that("a reloaded .fitrx (TTE rows dropped, n_obs_tte = 0) has no TTE block", {
  fit <- mixed_fit()
  fit$sdtab <- fit$sdtab[1:6, setdiff(names(fit$sdtab), tte_cols)]
  fit$subject_summary$n_obs_tte <- c(0L, 0L)
  b <- blocks_fn(fit)
  expect_blocks(b, 1:3, 4:6, integer())
  expect_identical(b$subject$discrete, c(1L, 2L, 2L))
})

test_that("a binary-only fit is all discrete, by counts or by N_OBS", {
  fit <- mixed_fit(n_obs = c(0L, 0L), n_tte = c(0L, 0L))
  fit$sdtab <- fit$sdtab[, setdiff(names(fit$sdtab), tte_cols)]
  expect_blocks(blocks_fn(fit), integer(), 1:3, integer())
  fit$subject_summary <- NULL
  expect_blocks(blocks_fn(fit), integer(), 1:3, integer())
})

test_that("a TTE-only fit is all TTE, by counts or by COXSNELL", {
  fit <- mixed_fit(n_obs = c(0L, 0L), n_disc = c(0L, 0L))
  b <- blocks_fn(fit)
  expect_blocks(b, integer(), integer(), 1:2)
  expect_identical(b$subject$tte, 1:2)
  fit$subject_summary <- NULL
  expect_blocks(blocks_fn(fit), integer(), integer(), 1:2)
})

test_that("without counts, the TTE columns need all four, and N_OBS is optional", {
  fit <- mixed_fit()
  fit$subject_summary <- NULL
  # An [output] column called PIT alone is not a TTE table.
  only_pit <- fit
  only_pit$sdtab <- fit$sdtab[, setdiff(names(fit$sdtab), tte_cols[1:3])]
  expect_blocks(blocks_fn(only_pit), 1:3, 4:8, integer())
  # Without N_OBS only the TTE rows can be told apart.
  no_nobs <- fit
  no_nobs$sdtab$N_OBS <- NULL
  expect_blocks(blocks_fn(no_nobs), 1:6, integer(), 7:8)
})

test_that("an empty or missing sdtab has no rows in any block", {
  fit <- mixed_fit()
  fit$subject_summary <- NULL
  fit$sdtab <- fit$sdtab[0, ]
  expect_blocks(blocks_fn(fit), integer(), integer(), integer())
  fit$sdtab <- NULL
  expect_blocks(blocks_fn(fit), integer(), integer(), integer())
})

test_that(".ferx_attach_npde leaves the TTE rows NA wherever they are", {
  attach_npde <- getFromNamespace(".ferx_attach_npde", "ferx")
  sdtab <- data.frame(ID = c(1, 1, 2, 2), TIME = c(1, 2, 1, 2))
  npde_tbl <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 2),
                         NPDE = c(0.1, 0.2, 0.3), NPD = c(1, 2, 3))
  out <- attach_npde(sdtab, npde_tbl, tte_rows = 3L)
  expect_identical(out$NPDE, c(0.1, 0.2, NA, 0.3))
  expect_identical(out$NPD, c(1, 2, NA, 3))
  # Leaving out the TTE rows does not skip the ID/TIME check on the rest.
  npde_tbl$TIME[3] <- 1
  expect_error(attach_npde(sdtab, npde_tbl, tte_rows = 3L), "do not line up")
})

# Each caller on the synthetic Gaussian + binary + TTE fit, the one geometry
# no bundled example covers (#689): with the counts, and with the columns.
for (case in c("counts", "columns")) {
  test_that(sprintf("the three callers split a binary + TTE fit (%s)", case), {
    fit <- mixed_fit()
    if (case == "columns") fit$subject_summary <- NULL
    fit$subject_ids <- c("11", "12")

    # ferx_xpose(): the Gaussian rows only, each on its own subject.
    gauss_rows <- getFromNamespace(".ferx_xpose_gaussian_rows", "ferx")
    expect_message(g <- gauss_rows(fit, fit$sdtab),
                   "3 discrete-endpoint and 2 time-to-event sdtab rows")
    expect_identical(g$data$TIME, 1:3)
    expect_false(any(tte_cols %in% names(g$data)))
    if (case == "counts") expect_identical(as.integer(g$subject), c(1L, 1L, 2L))

    # ferx_save_fit(): the Gaussian rows in predictions.csv; with the counts the
    # discrete rows in core's layout and the others in full. Without the counts
    # the rows cannot be placed, so it refuses.
    write_preds <- getFromNamespace(".fitrx_write_predictions_csv", "ferx")
    staging <- tempfile()
    dir.create(staging)
    if (case == "counts") {
      expect_identical(write_preds(fit, staging),
                       c("predictions.csv", "predictions_discrete.csv",
                         "sdtab_endpoint_rows.csv"))
      preds <- utils::read.csv(file.path(staging, "predictions.csv"))
      disc <- utils::read.csv(file.path(staging, "predictions_discrete.csv"))
      rest <- utils::read.csv(file.path(staging, "sdtab_endpoint_rows.csv"))
      expect_identical(preds$TIME, 1:3)
      expect_identical(disc$TIME, 4:6)
      expect_identical(disc$ID, c(11L, 12L, 12L))
      expect_identical(rest$TIME, 4:8)
      expect_identical(rest$SUBJECT, c(1L, 2L, 2L, 1L, 2L))
    } else {
      expect_error(write_preds(fit, staging), "cannot be placed", fixed = TRUE)
    }

    # ferx_calc_npde(): the engine's NPDE table has the Gaussian rows only
    # (ferx-core `compute_npde_npd` loops over `subject.observations`), so
    # this is a Gaussian + TTE fit; a fit with discrete rows is still refused
    # (#708).
    attach_npde <- getFromNamespace(".ferx_attach_npde", "ferx")
    pk_tte <- mixed_fit(n_disc = c(0L, 0L))
    if (case == "columns") pk_tte$subject_summary <- NULL
    npde_tbl <- cbind(pk_tte$sdtab[1:3, c("ID", "TIME")], NPDE = 1:3 / 10, NPD = 1:3)
    out <- attach_npde(pk_tte$sdtab, npde_tbl, tte_rows = blocks_fn(pk_tte)$tte)
    expect_identical(out$NPDE, c(1:3 / 10, NA, NA))
  })
}
