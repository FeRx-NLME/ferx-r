# print.ferx_fit()'s `--- Correlations ---` blocks, under the OMEGA rows and a
# block_kappa's KAPPA rows, report the pairs ferx-core's offdiag_corrs reports
# (ferx-core #2231), with its correlations; R only formats the lines (#547).

offdiag_corrs <- getFromNamespace(".ferx_offdiag_corrs", "ferx")

test_that("pairs come in row-major lower-triangle order, structural zeros skipped", {
  m <- diag(c(0.04, 0.09, 0.16, 0.25))
  m[2, 1] <- m[1, 2] <- 0.012
  m[4, 1] <- m[1, 4] <- -0.02
  m[4, 3] <- m[3, 4] <- 0.03
  m[3, 2] <- m[2, 3] <- 1e-16 # at the threshold: not reported
  got <- offdiag_corrs(m)
  expect_identical(got$i, c(2L, 4L, 4L))
  expect_identical(got$j, c(1L, 1L, 3L))
  expect_identical(got$cov, c(0.012, -0.02, 0.03))
  # No param_corr matrix: the cov / sqrt(var_i * var_j) fallback.
  expect_false(any(got$param_corr))
  expect_equal(got$corr, c(0.012 / (0.2 * 0.3), -0.02 / (0.2 * 0.5), 0.03 / (0.4 * 0.5)))
})

test_that("a param_corr matrix supplies each correlation and says so", {
  m <- matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
  pc <- matrix(c(1, 0.42, 0.42, 1), 2, 2)
  got <- offdiag_corrs(m, pc)
  expect_identical(got$corr, 0.42)
  expect_true(got$param_corr)
})

test_that("a non-positive variance gives the fallback correlation 0", {
  m <- matrix(c(0, 0.01, 0.01, 0.09), 2, 2)
  got <- offdiag_corrs(m)
  expect_identical(got$corr, 0)
  expect_false(got$param_corr)
})

test_that("a diagonal or 1 x 1 matrix reports no pairs", {
  expect_identical(nrow(offdiag_corrs(diag(c(0.04, 0.09)))), 0L)
  expect_identical(nrow(offdiag_corrs(0.04)), 0L)
})

test_that("the glue refuses a matrix it cannot read", {
  glue <- getFromNamespace("ferx_rust_offdiag_corrs", "ferx")
  expect_error(glue(c(0.04, 0, 0), 2L, NULL), "2 x 2")
  expect_error(glue(c(0.04, 0, 0, 0.09), 2L, c(1, 0.4, 0.4, 1)), "numeric matrix")
})

test_that("print: block_kappa correlations use the engine's pairs, R's labels", {
  iov <- diag(c(0.05, 0.04, 0.03))
  iov[1, 2] <- iov[2, 1] <- 0.01
  iov[3, 2] <- iov[2, 3] <- -0.006
  base <- list(
    omega = matrix(0.10, 1, 1), omega_iov = iov,
    kappa_names = c("KAPPA_CL", "KAPPA_V", "KAPPA_KA"),
    eta_param_types = "log_normal"
  )
  cov_lines <- function(fit) {
    out <- utils::capture.output(print(fit))
    out[grepl(" : cov = ", out, fixed = TRUE)]
  }
  expect_identical(cov_lines(do.call(make_fake_fit, base)), c(
    sprintf("  KAPPA_V ~ KAPPA_CL : cov = 0.010000  (corr = %.4f)", 0.01 / sqrt(0.05 * 0.04)),
    sprintf("  KAPPA_KA ~ KAPPA_V : cov = -0.006000  (corr = %.4f)", -0.006 / sqrt(0.04 * 0.03))
  ))
  pc <- diag(3)
  pc[2, 1] <- pc[1, 2] <- 0.31
  pc[3, 2] <- pc[2, 3] <- -0.17
  fit <- do.call(make_fake_fit, c(base, list(omega_iov_param_corr = pc)))
  expect_identical(cov_lines(fit), c(
    "  KAPPA_V ~ KAPPA_CL : cov = 0.010000  (param corr = 0.3100)",
    "  KAPPA_KA ~ KAPPA_V : cov = -0.006000  (param corr = -0.1700)"
  ))
})
