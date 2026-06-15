# Tests for the Xpose interop layer. These exercise the pure table/role
# builder (`.ferx_xpose_frame()` and helpers), which has no Xpose dependency.
# Backend object construction is only smoke-tested for its guard behaviour,
# since xpose / xpose4 are optional Suggests and may be absent.

# Minimal fake ferx_fit with covariates, etas, and individual parameters.
fake_fit <- function(with_covtab = TRUE, echo_wt = FALSE) {
  sdtab <- data.frame(
    ID    = c(1, 1, 2, 2),
    TIME  = c(1, 4, 1, 4),
    DV    = c(10, 8, 12, 9),
    PRED  = c(11, 7, 13, 8),
    IPRED = c(10.5, 7.5, 12.5, 8.5),
    CWRES = c(0.1, -0.2, 0.3, -0.1),
    IWRES = c(0.2, -0.1, 0.2, -0.2),
    TAD   = c(1, 4, 1, 4),
    stringsAsFactors = FALSE
  )
  if (echo_wt) sdtab$WT <- c(70, 70, 80, 80)

  fit <- list(
    sdtab = sdtab,
    individual_estimates = data.frame(ID = c(1, 2), CL = c(1.1, 1.4),
                                      V = c(30, 35)),
    ebe_etas = data.frame(ID = c(1, 2), ETA_CL = c(0.05, -0.03),
                          ETA_V = c(-0.1, 0.2)),
    covariate_types = c(WT = "continuous", SEX = "categorical"),
    ofv = 123.45,
    theta = c(CL = 1.0, V = 32),
    omega = matrix(c(0.04, 0, 0, 0.09), 2, 2)
  )
  if (with_covtab) {
    fit$covtab <- data.frame(
      ID   = c(1, 1, 2, 2),
      TIME = c(0, 1, 0, 1),
      EVID = c(1, 0, 1, 0),
      WT   = c(70, 70, 80, 80),
      SEX  = c(0, 0, 1, 1)
    )
  }
  class(fit) <- "ferx_fit"
  fit
}

test_that("frame derives RES/IRES and NA WRES", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  expect_equal(df$RES, df$DV - df$PRED)
  expect_equal(df$IRES, df$DV - df$IPRED)
  expect_true(all(is.na(df$WRES)))
})

test_that("etas and individual parameters are joined per subject (repeated by row)", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  expect_true(all(c("CL", "V", "ETA_CL", "ETA_V") %in% names(df)))
  # subject 1 rows get subject 1's values
  expect_equal(df$CL[df$ID == 1], c(1.1, 1.1))
  expect_equal(df$ETA_V[df$ID == 2], c(0.2, 0.2))
  expect_setequal(tbl$params, c("CL", "V"))
  expect_setequal(tbl$etas, c("ETA_CL", "ETA_V"))
})

test_that("covariates split into continuous / categorical via covariate_types", {
  tbl <- .ferx_xpose_frame(fake_fit())
  expect_setequal(tbl$cont, "WT")
  expect_setequal(tbl$cat, "SEX")
  expect_true(all(c("WT", "SEX") %in% names(tbl$data)))
})

test_that("continuous/categorical args override classification", {
  tbl <- .ferx_xpose_frame(fake_fit(), categorical = "WT")
  expect_false("WT" %in% tbl$cont)
  expect_true("WT" %in% tbl$cat)
})

test_that("unknown covariate override warns", {
  expect_warning(
    .ferx_xpose_frame(fake_fit(), continuous = "NOPE"),
    "not declared"
  )
})

test_that("covariate values come from covtab via LOCF when not echoed in sdtab", {
  tbl <- .ferx_xpose_frame(fake_fit())
  df <- tbl$data
  # WT is time-constant per subject in the fixture
  expect_equal(df$WT[df$ID == 1], c(70, 70))
  expect_equal(df$WT[df$ID == 2], c(80, 80))
})

test_that("covariate already echoed in sdtab is used directly (no covtab needed)", {
  # SEX (categorical, declared) has neither an echoed column nor a covtab here,
  # so it is dropped with a warning; this test only asserts the echoed WT.
  tbl <- suppressWarnings(
    .ferx_xpose_frame(fake_fit(with_covtab = FALSE, echo_wt = TRUE)))
  expect_true("WT" %in% names(tbl$data))
  expect_equal(tbl$data$WT, c(70, 70, 80, 80))
})

test_that("LOCF join carries the latest record at or before each obs time", {
  id  <- c(1, 1)
  t   <- c(2, 6)
  sid <- c(1, 1, 1)
  st  <- c(0, 3, 5)
  sv  <- c(10, 20, 30)
  expect_equal(.ferx_locf_join(id, t, sid, st, sv), c(10, 30))
})

test_that("roles map points at the expected columns", {
  roles <- .ferx_xpose_frame(fake_fit())$roles
  expect_equal(roles$id, "ID")
  expect_equal(roles$idv, "TIME")
  expect_equal(roles$dv, "DV")
  expect_setequal(roles$pred, c("PRED", "IPRED"))
  expect_setequal(roles$res, c("RES", "IRES", "WRES", "CWRES", "IWRES"))
  expect_setequal(roles$contcov, "WT")
  expect_setequal(roles$catcov, "SEX")
})

test_that("ferx_xpose rejects non-fit input", {
  expect_error(ferx_xpose(list()), "ferx_fit")
})

test_that("ferx_xpose errors clearly when the backend package is absent", {
  skip_if(requireNamespace("xpose", quietly = TRUE),
          "xpose installed; guard path not exercised")
  expect_error(ferx_xpose(fake_fit(), backend = "xpose"), "xpose")
})

# --- backend integration (only when the packages are installed) --------------

test_that("backend 'xpose' builds an xpose_data with correct column roles", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose")
  expect_s3_class(xpdb, "xpose_data")
  idx <- xpdb$data$index[[1]]
  role <- function(col) idx$type[idx$col == col]
  expect_equal(role("ID"), "id")
  expect_equal(role("TIME"), "idv")
  expect_equal(role("DV"), "dv")
  expect_equal(role("PRED"), "pred")
  expect_equal(role("IPRED"), "ipred")   # distinct from pred in xpose
  expect_equal(role("CWRES"), "res")
  expect_equal(role("CL"), "param")
  expect_equal(role("ETA_CL"), "eta")
  expect_equal(role("WT"), "contcov")
  expect_equal(role("SEX"), "catcov")
})

test_that("backend 'xpose' GOF templates render", {
  skip_if_not_installed("xpose")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(print(xpose::dv_vs_ipred(xpdb))))
  expect_no_error(suppressWarnings(print(xpose::res_vs_idv(xpdb, res = "CWRES"))))
})

test_that("backend 'xpose4' builds an xpose.data with covariate/param vardefs", {
  skip_if_not_installed("xpose4")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose4")
  expect_s4_class(xpdb, "xpose.data")
  expect_setequal(xpdb@Prefs@Xvardef$covariates, c("WT", "SEX"))
  expect_setequal(xpdb@Prefs@Xvardef$parms, c("CL", "V", "ETA_CL", "ETA_V"))
  # categorical covariate coerced to factor in the data slot
  expect_s3_class(xpdb@Data$SEX, "factor")
})

test_that("backend 'xpose4' default plot renders", {
  skip_if_not_installed("xpose4")
  xpdb <- ferx_xpose(fake_fit(), backend = "xpose4")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(print(
    xpose4::xpose.plot.default("PRED", "DV", xpdb)
  )))
})
