# ferx_gam_screen() -- GAM covariate pre-screening.
#
# Synthetic unit tests that assemble a minimal fit-shaped list so the
# OLS/spline/AIC logic is exercised deterministically without a real
# estimation run. A live-fit integration check at the bottom uses the
# cached two_cpt_oral_cov fit (same pattern as test-ferx_cov_screen.R).

# ---- synthetic fit builder --------------------------------------------------
# Known ground truth:
#   ETA_CL ~ 0.40 * (WT - 70) / 70  (continuous, linear signal)
#   ETA_V  ~ 0.35 * (SEX - 0.5)     (categorical signal)
#   CRCL   ~ noise                   (should not rank first for either ETA)

make_gam_fit <- function(n = 60L, seed = 1L) {
  set.seed(seed)
  ids  <- as.character(seq_len(n))
  wt   <- seq(50, 90, length.out = n)
  sex  <- rep(c(0, 1), length.out = n)
  crcl <- 80 + rnorm(n, 0, 20)

  eta_cl <- 0.40 * (wt - 70) / 70 + rnorm(n, 0, 0.15)
  eta_v  <- 0.35 * (sex - 0.5)    + rnorm(n, 0, 0.20)

  ebe_etas <- data.frame(
    ID     = ids,
    ETA_CL = eta_cl,
    ETA_V  = eta_v,
    stringsAsFactors = FALSE
  )

  # covtab: two rows per subject so median aggregation is tested.
  covtab <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      ID   = ids[i],
      TIME = c(0, 24),
      EVID = c(1L, 0L),
      WT   = wt[i]   + c(-1, 1),
      CRCL = crcl[i] + c(0, 0),
      SEX  = sex[i],
      stringsAsFactors = FALSE
    )
  }))

  list(
    ebe_etas        = ebe_etas,
    covtab          = covtab,
    covariate_types = c(WT = "continuous", CRCL = "continuous",
                        SEX = "categorical"),
    shrinkage_eta   = c(ETA_CL = 0.10, ETA_V = 0.12)
  )
}


# ---- column contract --------------------------------------------------------

test_that("returns the documented columns in the right types", {
  res <- ferx_gam_screen(make_gam_fit())
  expect_s3_class(res, "data.frame")
  expected_cols <- c("eta_name", "covariate", "delta_aic", "best_form",
                     "aic", "aic_null", "r_squared", "shrinkage")
  expect_identical(names(res), expected_cols)
  expect_type(res$eta_name,  "character")
  expect_type(res$covariate, "character")
  expect_type(res$delta_aic, "double")
  expect_type(res$r_squared, "double")
  expect_type(res$shrinkage, "double")
})


# ---- ranking ----------------------------------------------------------------

test_that("WT ranks first for ETA_CL and SEX ranks first for ETA_V", {
  res <- ferx_gam_screen(make_gam_fit())

  cl_rows <- res[res$eta_name == "ETA_CL", ]
  cl_rows <- cl_rows[order(-cl_rows$delta_aic), ]
  expect_equal(cl_rows$covariate[1L], "WT")

  v_rows <- res[res$eta_name == "ETA_V", ]
  v_rows <- v_rows[order(-v_rows$delta_aic), ]
  expect_equal(v_rows$covariate[1L], "SEX")
})

test_that("within each ETA rows are ordered by delta_aic descending", {
  res <- ferx_gam_screen(make_gam_fit())
  for (eta in unique(res$eta_name)) {
    sub <- res[res$eta_name == eta, ]
    expect_false(is.unsorted(rev(sub$delta_aic)))
  }
})


# ---- AIC formula ------------------------------------------------------------

test_that("delta_aic is aic_null minus aic (not the reverse)", {
  res <- ferx_gam_screen(make_gam_fit())
  expect_equal(res$delta_aic, res$aic_null - res$aic, tolerance = 1e-10)
})

test_that("r_squared is in [0, 1]", {
  res <- ferx_gam_screen(make_gam_fit())
  expect_true(all(res$r_squared >= -1e-10 & res$r_squared <= 1 + 1e-10))
})


# ---- best_form values -------------------------------------------------------

test_that("best_form is one of the documented values", {
  res <- ferx_gam_screen(make_gam_fit())
  allowed <- c("Linear", "Categorical",
                paste0("Spline(df=", c(2L, 3L), ")"))
  expect_true(all(res$best_form %in% allowed))
})

test_that("SEX always gets Categorical form", {
  res <- ferx_gam_screen(make_gam_fit())
  sex_rows <- res[res$covariate == "SEX", ]
  expect_true(all(sex_rows$best_form == "Categorical"))
})


# ---- shrinkage warning ------------------------------------------------------

test_that("emits a warning when an ETA's shrinkage exceeds the threshold", {
  fit <- make_gam_fit()
  fit$shrinkage_eta <- c(ETA_CL = 0.50, ETA_V = 0.10)
  expect_warning(
    ferx_gam_screen(fit, shrinkage_warn = 0.30),
    "ETA_CL"
  )
})

test_that("no warning when shrinkage is below the threshold", {
  fit <- make_gam_fit()
  fit$shrinkage_eta <- c(ETA_CL = 0.10, ETA_V = 0.12)
  expect_no_warning(ferx_gam_screen(fit, shrinkage_warn = 0.30))
})


# ---- guard paths ------------------------------------------------------------

test_that("returns NULL with a message when fit has no ebe_etas", {
  fit <- make_gam_fit()
  fit$ebe_etas <- NULL
  expect_message(out <- ferx_gam_screen(fit), "No ETAs")
  expect_null(out)
})

test_that("returns NULL with a message when fit has no covtab", {
  fit <- make_gam_fit()
  fit$covtab <- NULL
  expect_message(out <- ferx_gam_screen(fit), "No covariate table")
  expect_null(out)
})

test_that("returns NULL when no ETA columns remain after filtering", {
  fit <- make_gam_fit()
  expect_message(
    out <- ferx_gam_screen(fit, etas = "ETA_DOES_NOT_EXIST"),
    "No ETA columns"
  )
  expect_null(out)
})

test_that("returns empty data frame when no covariates match after filtering", {
  fit <- make_gam_fit()
  expect_message(
    out <- ferx_gam_screen(fit, covariates = "NOSUCHCOV"),
    "No declared covariates"
  )
  expect_null(out)
})

test_that("returns empty data frame when declared covariates absent from covtab", {
  fit <- make_gam_fit()
  fit$covtab <- fit$covtab[, c("ID", "TIME", "EVID")]
  expect_message(out <- ferx_gam_screen(fit), "No declared covariates")
  expect_null(out)
})


# ---- subsetting arguments ---------------------------------------------------

test_that("etas argument restricts which ETAs are screened", {
  res <- ferx_gam_screen(make_gam_fit(), etas = "ETA_CL")
  expect_equal(unique(res$eta_name), "ETA_CL")
})

test_that("covariates argument restricts which covariates are screened", {
  res <- ferx_gam_screen(make_gam_fit(), covariates = "WT")
  expect_equal(unique(res$covariate), "WT")
})

test_that("include_linear = FALSE still returns results via splines", {
  res <- ferx_gam_screen(make_gam_fit(), include_linear = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0L)
  expect_false(any(res$best_form == "Linear"))
})

test_that("custom spline_df is respected", {
  res <- ferx_gam_screen(make_gam_fit(), spline_df = 4L, include_linear = FALSE)
  non_cat <- res[res$covariate != "SEX", ]
  if (nrow(non_cat) > 0L) {
    expect_true(all(grepl("Spline", non_cat$best_form)))
  }
})


# ---- covariate_types fallback -----------------------------------------------

test_that("infers covariate types from storage mode when covariate_types is NULL", {
  fit <- make_gam_fit()
  fit$covariate_types <- NULL
  res <- ferx_gam_screen(fit)
  expect_s3_class(res, "data.frame")
  # WT and CRCL are numeric -> continuous; SEX is numeric too here (0/1),
  # so it is inferred as continuous. The test just verifies no crash.
  expect_true(nrow(res) > 0L)
})

test_that("shrinkage column is NA when fit has no shrinkage_eta", {
  fit <- make_gam_fit()
  fit$shrinkage_eta <- NULL
  res <- ferx_gam_screen(fit)
  expect_true(all(is.na(res$shrinkage)))
})


# ---- integration: real two_cpt_oral_cov fit ---------------------------------

cov_fit_gam <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex  <- ferx_example("two_cpt_oral_cov")
      fit <<- ferx_fit(
        ex$model, ex$data,
        method    = "focei",
        verbose   = FALSE,
        covariance = FALSE,
        settings  = list(maxiter = 5L)
      )
    }
    fit
  }
})

test_that("runs on a real two_cpt_oral_cov fit and returns a valid table", {
  fit <- cov_fit_gam()
  res <- ferx_gam_screen(fit)

  expect_s3_class(res, "data.frame")
  expected_cols <- c("eta_name", "covariate", "delta_aic", "best_form",
                     "aic", "aic_null", "r_squared", "shrinkage")
  expect_identical(names(res), expected_cols)
  expect_setequal(unique(res$covariate), c("WT", "CRCL"))
  # delta_aic = aic_null - aic everywhere.
  expect_equal(res$delta_aic, res$aic_null - res$aic, tolerance = 1e-10)
})
