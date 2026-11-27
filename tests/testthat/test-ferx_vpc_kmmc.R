# ferx_vpc_kmmc(): the Kaplan-Meier mean covariate VPC (ferx-core #1886 PR3,
# #2955). Every number is core's; these tests check that the glue carries it
# into R unchanged, against the at-risk mean computed by hand (core checks it
# against survival::survfit's risk sets), and the shapes, settings and
# refusals.

kmmc_data <- function() {
  ex <- ferx_example("tte_weibull")
  d <- utils::read.csv(ex$data)
  d$WT <- 50 + 3 * (d$ID %% 10)
  d$SEX <- d$ID %% 2
  d
}

kmmc_quiet <- function(...) suppressWarnings(ferx_vpc_kmmc(...))

test_that("ferx_vpc_kmmc returns the observed curve and the band", {
  ex <- ferx_example("tte_weibull")
  v <- kmmc_quiet(ex$model, kmmc_data(), covariate = "WT", n_sim = 30L, seed = 1L)
  expect_s3_class(v, "ferx_vpc_kmmc")
  expect_identical(v$covariate, "WT")
  expect_true(is.na(v$level))
  expect_identical(v$min_at_risk, 10L)
  expect_identical(v$n_sim, 30L)
  expect_identical(v$cmt, 2L)
  expect_identical(v$censoring, "observed")
  expect_named(v$observed, c("STRATUM", "time", "n_risk", "n_event",
                             "n_censor", "mean"))
  expect_named(v$band, c("STRATUM", "time", "obs", "lo", "med", "hi", "n_used"))
  b <- v$band[!is.na(v$band$lo), ]
  expect_true(all(b$lo <= b$med & b$med <= b$hi))
  # At time 0 everyone is at risk everywhere: the band has width 0 there.
  d <- kmmc_data()
  base <- mean(d$WT[!duplicated(d$ID)])
  expect_equal(unlist(v$band[1, c("obs", "lo", "med", "hi")]),
               c(obs = base, lo = base, med = base, hi = base))
  expect_output(print(v), "Mean WT at risk; at least 10 at risk")
})

test_that("the observed curve is the at-risk mean, row for row", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  v <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 1L)
  r <- d[d$CMT == 2, ]
  # mean[j] is over the subjects with time >= time[j] (core's risk set).
  expect_equal(v$observed$time, sort(unique(r$TIME)))
  expect_equal(v$observed$n_risk,
               vapply(v$observed$time, function(t) sum(r$TIME >= t), integer(1)))
  expect_equal(v$observed$mean,
               vapply(v$observed$time, function(t) mean(r$WT[r$TIME >= t]),
                      numeric(1)),
               tolerance = 1e-12)
  # level = v: the share at risk at exactly that value.
  s <- kmmc_quiet(ex$model, d, covariate = "SEX", level = 1, n_sim = 20L,
                  seed = 1L)
  expect_equal(s$level, 1)
  expect_equal(s$observed$mean,
               vapply(s$observed$time, function(t) mean(r$SEX[r$TIME >= t] == 1),
                      numeric(1)),
               tolerance = 1e-12)
  expect_output(print(s), "Share at risk with SEX = 1")
})

test_that("min_at_risk, the seed, strata and a fit reach the engine", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  a <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L)
  b <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L)
  c <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 4L)
  expect_identical(a$band, b$band)
  expect_false(identical(a$band$med, c$band$med))
  m25 <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L,
                    min_at_risk = 25L)
  expect_identical(m25$min_at_risk, 25L)
  # obs is undefined once fewer than min_at_risk remain.
  expect_gt(sum(is.na(m25$band$obs)), sum(is.na(a$band$obs)))

  s <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, strata = "SEX",
                  min_at_risk = 5L)
  expect_identical(levels(s$strata$STRATUM), c("SEX=0", "SEX=1"))

  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                   verbose = FALSE))
  vf <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L,
                   fit = fit)
  expect_identical(vf$observed, a$observed)
  expect_false(identical(vf$band$med, a$band$med))
})

test_that("the observed records are cut at the horizon (#2955)", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  full <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L)
  h <- kmmc_quiet(ex$model, d, covariate = "WT", n_sim = 20L, seed = 3L,
                  horizon = 20)
  expect_equal(h$horizon, 20)
  expect_true(all(h$observed$time <= 20))
  keep <- full$observed$time < 20
  expect_equal(h$observed$mean[seq_len(sum(keep))], full$observed$mean[keep])
})

test_that("each refusal comes back as E_VPC_TTE_INPUT", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  refused <- function(expr, pattern) {
    e <- tryCatch(expr, error = function(e) e)
    expect_s3_class(e, "ferx_engine_error")
    expect_identical(e$code, "E_VPC_TTE_INPUT")
    expect_match(conditionMessage(e), pattern)
  }
  refused(ferx_vpc_kmmc(ex$model, d, covariate = "AGE", n_sim = 20L),
          "has no covariate `AGE`")
  refused(ferx_vpc_kmmc(ex$model, d, covariate = "WT", level = 999,
                        n_sim = 20L), "`level` asked for")
  refused(ferx_vpc_kmmc(ex$model, d, covariate = "WT", strata = "WT",
                        n_sim = 20L), "makes the mean constant")
  # A repeated-events endpoint names the R function that draws its VPC.
  m <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVLAMBDA(0.05, 0.001, 10.0)",
               "  omega ETA_LAMBDA ~ 0.09", "[event_model]", "  cmt    = 2",
               "  type   = rtte", "  family = exponential",
               "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"), m)
  dr <- data.frame(ID = rep(1:3, each = 2), TIME = c(3, 30, 5, 30, 8, 30),
                   DV = c(1, 0, 1, 0, 1, 0), EVID = 0L, AMT = 0, CMT = 2L,
                   MDV = 0L, WT = 70)
  e <- tryCatch(ferx_vpc_kmmc(m, dr, covariate = "WT", n_sim = 20L),
                error = function(e) e)
  expect_identical(e$code, "E_VPC_TTE_INPUT")
  expect_match(conditionMessage(e), "`ferx_vpc_tte()`", fixed = TRUE)
  expect_false(grepl("`vpc_rtte`", conditionMessage(e), fixed = TRUE))

  expect_error(ferx_vpc_kmmc(ex$model, d, n_sim = 20L),
               "`covariate` must be a single covariate name")
  expect_error(ferx_vpc_kmmc(ex$model, d, covariate = "WT", level = NA),
               "`level` must be NULL or a single finite number")
  expect_error(ferx_vpc_kmmc(ex$model, d, covariate = "WT", min_at_risk = 0),
               "`min_at_risk` must be a single whole number >= 1")
})

test_that("a time-varying covariate warns, with guidance", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  # Two records for subject 1: a dose row with a different WT first.
  extra <- d[d$ID == 1, ][1, ]
  extra$TIME <- 0
  extra$EVID <- 1L
  extra$AMT <- 0
  extra$MDV <- 1L
  extra$WT <- 99
  d2 <- rbind(extra, d)
  d2 <- d2[order(d2$ID, d2$TIME, -d2$EVID), ]
  expect_warning(v <- ferx_vpc_kmmc(ex$model, d2, covariate = "WT", n_sim = 20L),
                 "W_VPC_KMMC_TIME_VARYING_COVARIATE")
  expect_true(any(grepl("W_VPC_KMMC_TIME_VARYING_COVARIATE: The curve uses",
                        capture.output(print(v)), fixed = TRUE)))
})

test_that("the KMMC band is drawn left-continuous, the row at H kept (#754 review, 2)", {
  bins <- getFromNamespace(".ferx_vpc_kmmc_bins", "ferx")
  b <- data.frame(time = c(10, 0, 5, 20, 25), lo = c(3, 1, 2, 4, 5),
                  med = c(3, 1, 2, 4, 5), hi = c(3, 1, 2, 4, 5))
  out <- bins(b, horizon = 20)
  # Sorted; each row's interval is (previous grid time, its own time].
  expect_identical(out$time, c(0, 5, 10, 20))
  expect_identical(out$bin_min, c(0, 0, 5, 10))
  expect_identical(out$lo, c(1, 2, 3, 4))
  # An undefined row is dropped, and the next row still starts at its time.
  b$lo[b$time == 5] <- NA
  out2 <- bins(b, horizon = 20)
  expect_identical(out2$time, c(0, 10, 20))
  expect_identical(out2$bin_min, c(0, 5, 10))
})

test_that("plot() draws one panel per stratum", {
  ex <- ferx_example("tte_weibull")
  v <- kmmc_quiet(ex$model, kmmc_data(), covariate = "WT", n_sim = 20L,
                  strata = "SEX", min_at_risk = 5L)
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot(v), v)
})

test_that("an empty stratum warns, with guidance, and its panel says so (ferx-core #2947)", {
  ex <- ferx_example("tte_weibull")
  d <- kmmc_data()
  n_subjects <- length(unique(d$ID))
  # A floor above the subject count: no risk set ever reaches it.
  expect_warning(
    v <- ferx_vpc_kmmc(ex$model, d, covariate = "WT", n_sim = 20L, seed = 1L,
                       min_at_risk = n_subjects + 1L),
    "W_VPC_EMPTY_STRATUM"
  )
  w <- grep("^W_VPC_EMPTY_STRATUM", v$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "fewer than `min_at_risk", fixed = TRUE)
  expect_true(all(is.na(v$band$lo)))
  expect_match(ferx:::.ferx_vpc_warning_guidance(w), "lower min_at_risk", fixed = TRUE)
  expect_true(any(grepl("W_VPC_EMPTY_STRATUM: This stratum's panel has no simulated band",
                        capture.output(print(v)), fixed = TRUE)))
  # The panel carries the note instead of the observed curve alone.
  notes <- 0L
  testthat::local_mocked_bindings(
    .ferx_vpc_no_band_note = function() notes <<- notes + 1L,
    .package = "ferx"
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(v)
  expect_identical(notes, 1L)
})

test_that("the empty-stratum guidance follows the engine's reason (ferx-core #2947)", {
  g <- function(why) ferx:::.ferx_vpc_warning_guidance(
    paste0("W_VPC_EMPTY_STRATUM: stratum `all` has no band at any grid time.", why)
  )
  expect_match(g(" It has 5 subjects, fewer than `min_at_risk = 10`, so ..."),
               "lower min_at_risk", fixed = TRUE)
  expect_match(g(" Every grid time is past the last time of every replicate curve; ..."),
               "grid times inside the simulated follow-up", fixed = TRUE)
  expect_match(g(" At most 3 of the 20 replicate curves are defined at any grid time ..."),
               "raise n_sim", fixed = TRUE)
})
