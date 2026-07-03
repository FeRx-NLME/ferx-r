# fake_cond_dist() mirrors the shape fit_result_to_list() (Rust) and
# ferx_load_fit() (persist) both produce for fit$cond_dist -- these tests
# exercise the R-side accessor/print with a mocked posterior so they run in
# milliseconds. The real SAEM round-trip is covered in test-ferx_save_fit.R.
fake_cond_dist <- function() {
  list(
    data = data.frame(
      ID = c("1", "1", "2", "2"),
      ETA = c("ETA_CL", "ETA_V", "ETA_CL", "ETA_V"),
      COND_MEAN = c(0.05, -0.02, -0.03, 0.01),
      COND_SD = c(0.12, 0.15, 0.11, 0.14),
      COND_MODE = c(0.04, -0.01, -0.02, 0.02),
      stringsAsFactors = FALSE
    ),
    shrinkage = c(0.1, 0.05),
    nsamp = 200L,
    burnin = 20L
  )
}

test_that("ferx_conddist returns the long-format data frame with attributes", {
  fit <- make_fake_fit(cond_dist = fake_cond_dist(), eta_names = c("ETA_CL", "ETA_V"))
  cd <- ferx_conddist(fit)

  expect_s3_class(cd, "ferx_conddist")
  expect_equal(nrow(cd), 4L)
  expect_equal(names(cd), c("ID", "ETA", "COND_MEAN", "COND_SD", "COND_MODE"))
  expect_equal(attr(cd, "nsamp"), 200L)
  expect_equal(attr(cd, "burnin"), 20L)
  expect_equal(unname(attr(cd, "shrinkage")), c(0.1, 0.05))
  expect_equal(names(attr(cd, "shrinkage")), c("ETA_CL", "ETA_V"))
})

test_that("ferx_conddist errors when the fit has no cond_dist", {
  fit <- make_fake_fit()
  expect_error(ferx_conddist(fit), "no conditional-distribution results")
})

test_that("ferx_conddist rejects non-ferx_fit input", {
  expect_error(ferx_conddist(list()), "ferx_fit object")
  expect_error(ferx_conddist(42), "ferx_fit object")
})

test_that("print.ferx_conddist shows shrinkage and the data frame", {
  fit <- make_fake_fit(cond_dist = fake_cond_dist(), eta_names = c("ETA_CL", "ETA_V"))
  cd <- ferx_conddist(fit)

  out <- capture.output(print(cd))
  expect_true(any(grepl("SAEM conditional distribution", out)))
  expect_true(any(grepl("200 draws retained, 20 burn-in", out)))
  expect_true(any(grepl("ETA_CL", out)))
  expect_true(any(grepl("COND_MEAN", out)))
})

test_that("ferx_conddist falls back to an unnamed shrinkage vector on a length mismatch", {
  # eta_names/shrinkage length disagreement shouldn't happen via ferx_fit()/
  # ferx_load_fit(), but a hand-built fit object (as this file's own
  # make_fake_fit() calls construct) shouldn't hit setNames()'s raw error
  # either (#248 review).
  cd <- fake_cond_dist()
  cd$shrinkage <- 0.1
  fit <- make_fake_fit(cond_dist = cd, eta_names = c("ETA_CL", "ETA_V"))

  out <- ferx_conddist(fit)
  expect_equal(unname(attr(out, "shrinkage")), 0.1)
  expect_null(names(attr(out, "shrinkage")))
})

test_that("print.ferx_conddist reports unknown provenance when nsamp/burnin are NA (loaded fit)", {
  cd <- fake_cond_dist()
  cd$nsamp <- NA_integer_
  cd$burnin <- NA_integer_
  fit <- make_fake_fit(cond_dist = cd, eta_names = c("ETA_CL", "ETA_V"))

  out <- capture.output(print(ferx_conddist(fit)))
  expect_true(any(grepl("unknown", out)))
  expect_false(any(grepl("NA draws", out)))
})

test_that("ferx_conddist survives a real SAEM conddist fit end-to-end", {
  skip_on_cran()
  fit <- warfarin_saem_conddist_fit()
  skip_if(is.null(fit$cond_dist), "no cond_dist result in this fit")

  cd <- ferx_conddist(fit)
  expect_s3_class(cd, "ferx_conddist")
  expect_equal(sort(unique(cd$ID)), sort(unique(as.character(fit$ebe_etas$ID))))
  expect_true(all(c("ID", "ETA", "COND_MEAN", "COND_SD", "COND_MODE") %in% names(cd)))
  expect_equal(length(attr(cd, "shrinkage")), length(fit$eta_names))
})
