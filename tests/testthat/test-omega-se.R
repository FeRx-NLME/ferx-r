# Unit tests for internal omega-SE lookup and IMPMAP trace reconstruction.
# Local aliases avoid the ::: operator (undesirable_operator_linter).
omega_se_at <- getFromNamespace(".omega_se_at", "ferx")
reconstruct_impmap_trace <- getFromNamespace(".reconstruct_impmap_trace", "ferx")

test_that(".omega_se_at returns NA when se_omega is NULL", {
  expect_identical(omega_se_at(NULL, 2L, 1L, 1L), NA_real_)
})

test_that(".omega_se_at reads diagonal-only se_omega by position", {
  se <- c(0.11, 0.22)  # length n_eta -> diagonal-only
  expect_equal(omega_se_at(se, 2L, 1L, 1L), 0.11)
  expect_equal(omega_se_at(se, 2L, 2L, 2L), 0.22)
})

test_that(".omega_se_at returns NA for off-diagonal of a diagonal-only se_omega", {
  se <- c(0.11, 0.22)
  expect_identical(omega_se_at(se, 2L, 2L, 1L), NA_real_)
})

test_that(".omega_se_at returns NA for diagonal index beyond se_omega length", {
  # n_eta = 2 but only one SE present -> (2,2) is out of range.
  expect_identical(omega_se_at(c(0.11), 2L, 2L, 2L), NA_real_)
})

test_that(".omega_se_at reads a 2x2 full lower-triangle (column-major)", {
  # Layout for n_eta=2: [(1,1), (2,1), (2,2)]
  se <- c(0.10, 0.20, 0.30)
  expect_equal(omega_se_at(se, 2L, 1L, 1L), 0.10)  # (1,1)
  expect_equal(omega_se_at(se, 2L, 2L, 1L), 0.20)  # (2,1) off-diagonal
  expect_equal(omega_se_at(se, 2L, 2L, 2L), 0.30)  # (2,2)
})

test_that(".omega_se_at is symmetric in (i, j)", {
  se <- c(0.10, 0.20, 0.30)
  expect_equal(omega_se_at(se, 2L, 1L, 2L), omega_se_at(se, 2L, 2L, 1L))
})

test_that(".omega_se_at reads a 3x3 full lower-triangle (column-major)", {
  # Layout for n_eta=3: [(1,1),(2,1),(3,1),(2,2),(3,2),(3,3)]
  se <- c(1, 2, 3, 4, 5, 6)
  expect_equal(omega_se_at(se, 3L, 1L, 1L), 1)
  expect_equal(omega_se_at(se, 3L, 3L, 1L), 3)  # (3,1)
  expect_equal(omega_se_at(se, 3L, 2L, 2L), 4)  # (2,2)
  expect_equal(omega_se_at(se, 3L, 3L, 2L), 5)  # (3,2)
  expect_equal(omega_se_at(se, 3L, 3L, 3L), 6)  # (3,3)
})

test_that(".reconstruct_impmap_trace rebuilds a data.frame column-major", {
  tr <- list(
    flat      = c(1, 2, 3, 4, 5, 6),  # 3 rows x 2 cols, column-major
    n_rows    = 3L,
    n_cols    = 2L,
    col_names = c("CL", "V")
  )
  df <- reconstruct_impmap_trace(tr)
  expect_s3_class(df, "data.frame")
  expect_identical(dim(df), c(3L, 2L))
  expect_identical(colnames(df), c("CL", "V"))
  expect_equal(df$CL, c(1, 2, 3))
  expect_equal(df$V, c(4, 5, 6))
})
