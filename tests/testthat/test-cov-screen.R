# ferx_cov_screen(): quick covariate screen against individual parameters and
# ETAs, built on fit$covtab. See issue #128. These are synthetic unit tests -
# they assemble a minimal fit-shaped list so the correlation/aggregation logic
# is exercised deterministically without a real estimation. A live-fit
# integration check lives in test-covtab.R.

# Build a fit-shaped list with a known covariate->parameter structure:
#   WT  (continuous)  drives ETA_CL exactly
#   SEX (categorical) shifts ETA_V1 by group
#   AGE (continuous)  is moderate noise, tunable via `age_r`
make_fit <- function(n = 40L, age_r = 0.1, seed = 42L) {
  set.seed(seed)
  ids  <- as.character(seq_len(n))
  z    <- function(v) (v - mean(v)) / stats::sd(v)
  wt_c <- seq(50, 100, length.out = n)
  sex  <- rep(c("M", "F"), length.out = n)

  eta_cl <- z(wt_c) * 0.3                                  # r(WT, ETA_CL) = 1
  eta_v1 <- ifelse(sex == "F", 0.4, -0.4) + rnorm(n, 0, 0.05)
  eta_ka <- rnorm(n, 0, 0.3)                               # unrelated to all
  age    <- z(wt_c) * age_r + rnorm(n, 0, 1)               # weak WT-linked noise

  ebe_etas <- data.frame(ID = ids, ETA_CL = eta_cl, ETA_V1 = eta_v1,
                         ETA_KA = eta_ka, stringsAsFactors = FALSE)
  indiv <- data.frame(ID = ids, CL = 5 * exp(eta_cl), V1 = 50 * exp(eta_v1),
                      KA = 1.2 * exp(eta_ka), stringsAsFactors = FALSE)

  # covtab: 3 rows/subject, WT time-varying around wt_c so the median recovers it
  covtab <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(ID = ids[i], TIME = c(0, 12, 24), EVID = c(1L, 0L, 0L),
               WT = wt_c[i] + c(-2, 0, 2), SEX = sex[i], AGE = age[i],
               stringsAsFactors = FALSE)
  }))

  list(covtab = covtab,
       covariate_types = c(WT = "continuous", SEX = "categorical",
                           AGE = "continuous"),
       ebe_etas = ebe_etas,
       individual_estimates = indiv)
}

test_that("returns the documented columns, ordered by descending strength", {
  res <- ferx_cov_screen(make_fit())
  expect_s3_class(res, "data.frame")
  expect_identical(names(res), c("parameter", "covariate", "type", "ebe", "eta"))
  strength <- pmax(abs(res$ebe), abs(res$eta), na.rm = TRUE)
  expect_false(is.unsorted(rev(strength)))   # non-increasing
})

test_that("surfaces the strong covariate-parameter pairs", {
  res <- ferx_cov_screen(make_fit())
  # WT drives CL (continuous, near-perfect), SEX shifts V1 (categorical)
  wt_cl <- res[res$covariate == "WT" & res$parameter == "CL", ]
  expect_equal(nrow(wt_cl), 1L)
  expect_gt(abs(wt_cl$eta), 0.9)
  expect_true(any(res$covariate == "SEX" & res$parameter == "V1"))
})

test_that("ETA names map to their individual parameter for the EBE column", {
  res <- ferx_cov_screen(make_fit())
  # parameter labels come from individual_estimates (CL/V1/KA), not ETA_*
  expect_true(all(res$parameter %in% c("CL", "V1", "KA")))
  expect_false(any(grepl("^ETA", res$parameter)))
})

test_that("median aggregation collapses time-varying covariate to one value/subject", {
  # WT is given as wt_c +/- 2 across 3 rows; the median is exactly wt_c, so
  # the WT-CL correlation stays ~1. A non-aggregating impl would mismatch rows.
  res <- ferx_cov_screen(make_fit())
  expect_equal(abs(res$eta[res$covariate == "WT" & res$parameter == "CL"]),
               1, tolerance = 1e-6)
})

test_that("categorical association is a correlation ratio in [0, 1]", {
  res <- ferx_cov_screen(make_fit())
  cat_rows <- res[res$type == "categorical", ]
  expect_true(all(cat_rows$eta >= 0 & cat_rows$eta <= 1, na.rm = TRUE))
})

test_that("threshold filters weak associations", {
  fit <- make_fit()
  res_lo <- ferx_cov_screen(fit, threshold = 0.2)
  res_hi <- ferx_cov_screen(fit, threshold = 0.9)
  # Every reported row clears its threshold...
  expect_true(all(pmax(abs(res_lo$ebe), abs(res_lo$eta), na.rm = TRUE) >= 0.2))
  expect_true(all(pmax(abs(res_hi$ebe), abs(res_hi$eta), na.rm = TRUE) >= 0.9))
  # ...and raising the threshold never adds rows.
  expect_lte(nrow(res_hi), nrow(res_lo))

  # Constant covariates have no variance -> nothing to report -> empty + message.
  flat <- fit
  flat$covtab$WT  <- 70
  flat$covtab$AGE <- 40
  flat$covtab$SEX <- "M"
  expect_message(empty <- ferx_cov_screen(flat), "No covariate correlations")
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
})

test_that("returns NULL with a message when the model has no covariate table", {
  fit <- make_fit()
  fit$covtab <- NULL
  expect_message(out <- ferx_cov_screen(fit), "No covariate table")
  expect_null(out)
})

test_that("returns NULL with a message when the model declares no ETAs", {
  fit <- make_fit()
  fit$ebe_etas <- NULL
  expect_message(out <- ferx_cov_screen(fit), "no .*inter-individual")
  expect_null(out)
})

test_that("falls back to inferring covariate types when covariate_types is absent", {
  fit <- make_fit()
  fit$covariate_types <- NULL  # force inference from column storage mode
  res <- ferx_cov_screen(fit)
  expect_s3_class(res, "data.frame")
  # WT is numeric -> inferred continuous; SEX is character -> categorical
  expect_identical(unique(res$type[res$covariate == "WT"]), "continuous")
  expect_identical(unique(res$type[res$covariate == "SEX"]), "categorical")
})

# ---- edge / guard paths -------------------------------------------------

test_that("returns NULL when declared covariates are absent from covtab", {
  fit <- make_fit()
  fit$covtab <- fit$covtab[, c("ID", "TIME", "EVID")]  # drop covariate columns
  expect_message(out <- ferx_cov_screen(fit), "No declared covariates")
  expect_null(out)
})

test_that("returns NULL when ebe_etas has no ETA columns", {
  fit <- make_fit()
  fit$ebe_etas <- data.frame(ID = fit$ebe_etas$ID)  # ID only, no etas
  expect_message(out <- ferx_cov_screen(fit), "No ETA columns")
  expect_null(out)
})

test_that("works without individual_estimates: ebe is NA, labels fall back to ETA names", {
  fit <- make_fit()
  fit$individual_estimates <- NULL
  res <- ferx_cov_screen(fit)
  expect_s3_class(res, "data.frame")
  expect_true(all(is.na(res$ebe)))                 # no EBE source -> NA column
  expect_true(all(res$parameter %in% c("CL", "V1", "KA")))  # stripped ETA names
  expect_true(any(!is.na(res$eta)))                # ETA correlations still run
})

test_that("fewer than 3 subjects yields NA associations (empty screen)", {
  fit <- make_fit(n = 2L)
  expect_message(out <- ferx_cov_screen(fit), "No covariate correlations")
  expect_equal(nrow(out), 0L)
})

test_that("categorical association is dropped when the parameter response is constant", {
  set.seed(1)
  n <- 6L
  ids <- as.character(seq_len(n))
  sex <- rep(c("M", "F"), length.out = n)
  fit <- list(
    covtab = data.frame(ID = rep(ids, each = 2L), TIME = rep(c(0, 1), n),
                        EVID = 0L, SEX = rep(sex, each = 2L),
                        stringsAsFactors = FALSE),
    covariate_types = c(SEX = "categorical"),
    ebe_etas = data.frame(ID = ids, ETA_CL = rnorm(n), ETA_V = rep(0.1, n)),
    individual_estimates = data.frame(ID = ids, CL = rnorm(n), V = rep(50, n))
  )
  res <- ferx_cov_screen(fit, threshold = 0)
  expect_false("V" %in% res$parameter)  # constant V -> correlation ratio undefined
  expect_true("CL" %in% res$parameter)  # varying CL is still screened
})

test_that("a subject whose categorical covariate is entirely missing aggregates to NA", {
  n <- 6L
  ids <- as.character(seq_len(n))
  sex <- rep(c("M", "F"), length.out = n)
  covtab <- data.frame(ID = rep(ids, each = 2L), TIME = rep(c(0, 1), n),
                       EVID = 0L, SEX = rep(sex, each = 2L),
                       stringsAsFactors = FALSE)
  covtab$SEX[covtab$ID == "1"] <- NA  # subject 1 fully missing -> mode_fn -> NA
  fit <- list(
    covtab = covtab,
    covariate_types = c(SEX = "categorical"),
    ebe_etas = data.frame(ID = ids, ETA_CL = as.numeric(seq_len(n))),
    individual_estimates = data.frame(ID = ids, CL = as.numeric(seq_len(n)))
  )
  expect_s3_class(ferx_cov_screen(fit, threshold = 0), "data.frame")
})
