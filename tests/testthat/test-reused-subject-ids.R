# Regression tests for #252: a dataset that reuses a subject ID in a
# non-contiguous block carries two distinct subjects sharing one textual ID
# (ferx-core, like NONMEM, splits subjects sequentially). The per-subject joins
# in xpose / save_fit / eta_cov / cov_screen must key on subject order, not on
# the raw ID, or the second block silently inherits the first subject's values.

test_that(".ferx_subject_index breaks on every ID change", {
  expect_identical(
    .ferx_subject_index(c("1", "1", "2", "2", "2", "1")),
    c(1L, 1L, 2L, 2L, 2L, 3L)
  )
  expect_identical(.ferx_subject_index(character(0)), integer(0))
  expect_identical(.ferx_subject_index(c(12, 12, 5, 12)), c(1L, 1L, 2L, 3L))
  expect_identical(.ferx_n_subjects(c("12", "12", "5", "12")), 3L)
})

test_that("xpose join gives a reused ID's second block its own subject value", {
  # sdtab: subject A (ID 12, 3 obs), B (ID 5, 2 obs), C (ID 12 reused, 1 obs).
  df  <- data.frame(ID = c(12, 12, 12, 5, 5, 12), TIME = 1:6)
  src <- data.frame(ID = c("12", "5", "12"), ETA_CL = c(0.1, 0.2, 0.9))
  res <- .ferx_join_by_id(df, src, "ebe_etas")
  expect_identical(res$data$ETA_CL, c(0.1, 0.1, 0.1, 0.2, 0.2, 0.9))
})

test_that("xpose join still matches by ID when counts don't line up", {
  # src covers a subset (2 subjects) but sdtab shows 3 -> fall back to ID match.
  df  <- data.frame(ID = c(1, 1, 2, 3), TIME = 1:4)
  src <- data.frame(ID = c("1", "2"), P = c(7, 8))
  res <- .ferx_join_by_id(df, src, "individual_estimates")
  expect_identical(res$data$P, c(7, 7, 8, NA_real_))
})

test_that("save_fit per-subject summary keeps both blocks of a reused ID", {
  fit <- list(
    sdtab = data.frame(
      ID      = c(12, 12, 12, 5, 5, 12),
      EBE_OFV = c(1, 1, 1, 2, 2, 9),
      N_OBS   = c(3, 3, 3, 2, 2, 1)
    ),
    ebe_etas = data.frame(ID = c("12", "5", "12"), ETA_CL = c(0.1, 0.2, 0.9))
  )
  ps <- .fitrx_per_subject_ofv_nobs(fit)
  expect_identical(ps$id,    c("12", "5", "12"))
  expect_identical(ps$n_obs, c(3L, 2L, 1L))
  expect_identical(ps$ofv,   c(1, 2, 9))

  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  .fitrx_write_ebes_csv(fit, tf)
  eb <- utils::read.csv(tf)
  # The two ID-12 rows must carry their own n_obs (3 and 1), not both = 3,
  # otherwise ebes.csv disagrees with predictions.csv and the loader rejects it.
  expect_identical(as.integer(eb$n_obs), c(3L, 2L, 1L))
})

test_that("eta_cov treats a reused ID as two subjects with distinct covariates", {
  data <- data.frame(
    ID   = c(1, 1, 1, 2, 2, 1, 1),
    TIME = c(0, 1, 2, 0, 1, 0, 1),
    DV   = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
    WT   = c(70, 70, 70, 80, 80, 50, 50)
  )
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(data, csv, row.names = FALSE)
  ebe <- data.frame(ID = c("1", "2", "1"), ETA_CL = c(0.1, 0.5, -0.4))

  res <- .ferx_compute_eta_cov(ebe, csv)
  expect_true(is.data.frame(res) && nrow(res) == 1L)
  # WT is genuinely constant within each of the 3 subjects (70/80/50), so it is
  # NOT flagged non-constant (which a split-by-ID would do, seeing 70 and 50
  # under one ID), and the correlation runs over all three distinct subjects.
  expect_true(is.finite(res$r))
})

test_that("cov_screen treats a reused ID as two subjects", {
  fit <- list(
    covtab = data.frame(
      ID   = c(1, 1, 1, 2, 2, 1, 1),
      TIME = c(0, 1, 2, 0, 1, 0, 1),
      WT   = c(70, 70, 70, 80, 80, 50, 50)
    ),
    covariate_types      = c(WT = "continuous"),
    ebe_etas             = data.frame(ID = c("1", "2", "1"),
                                      ETA_CL = c(0.1, 0.5, -0.4)),
    individual_estimates = data.frame(ID = c("1", "2", "1"),
                                      CL = c(1.1, 1.5, 0.6))
  )
  out <- ferx_cov_screen(fit, threshold = 0)
  expect_true(is.data.frame(out) && nrow(out) == 1L)
  expect_identical(out$parameter, "CL")
  expect_identical(out$covariate, "WT")
  # ETA and EBE associations both run over the three distinct subjects.
  expect_true(is.finite(out$eta) && is.finite(out$ebe))
})

test_that("xpose covariate LOCF does not leak across a reused ID's blocks", {
  # covtab: subject A (ID 12, WT 70), B (ID 5, WT 80), C (ID 12 reused, WT 50).
  covtab <- data.frame(
    ID   = c("12", "12", "5", "12"),
    TIME = c(0, 12, 0, 0),
    WT   = c(70, 70, 80, 50)
  )
  df <- data.frame(
    ID   = c(12, 12, 12, 5, 5, 12),
    TIME = c(0, 6, 12, 0, 6, 0)
  )
  out <- .ferx_attach_covariates(df, list(covtab = covtab), "WT")
  # Keying on the raw ID would put both ID-12 blocks in one LOCF group, so
  # subject A's rows would pick up subject C's WT = 50.
  expect_identical(out$data$WT, c(70, 70, 70, 80, 80, 50))
})
