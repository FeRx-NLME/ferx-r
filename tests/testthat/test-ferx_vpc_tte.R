# ferx_vpc_tte(): Kaplan-Meier VPC for a single-event TTE endpoint (ferx-core
# #1886). Every number is core's; these tests check that the glue carries it
# into R unchanged, and that each refusal arrives as E_VPC_TTE_INPUT.

vpc_tte_quiet <- function(...) suppressWarnings(ferx_vpc_tte(...))

# A one-endpoint exponential model and a data frame of TTE records, for the
# refusals the bundled examples cannot reach.
vpc_tte_model <- function(extra = character()) {
  m <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.05, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "",
    "[event_model]",
    "  cmt    = 2",
    extra,
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), m)
  m
}

vpc_tte_data <- function(n = 30L) {
  set.seed(7)
  t <- round(stats::rexp(n, 0.05), 2)
  data.frame(ID = seq_len(n), TIME = pmin(t, 30), DV = as.integer(t < 30),
             EVID = 0L, AMT = 0, CMT = 2L, MDV = 0L)
}

expect_vpc_tte_refusal <- function(expr, pattern) {
  e <- tryCatch(expr, error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_VPC_TTE_INPUT")
  expect_match(conditionMessage(e), pattern)
}

test_that("ferx_vpc_tte returns the observed curve and the band per stratum", {
  ex <- ferx_example("tte_exponential")
  v <- vpc_tte_quiet(ex$model, ex$data, n_sim = 50L, seed = 1L)
  expect_s3_class(v, "ferx_vpc_tte")
  expect_identical(v$cmt, 2L)
  expect_identical(v$n_sim, 50L)
  expect_identical(v$censoring, "observed")
  expect_equal(v$probs, c(0.05, 0.5, 0.95))
  # H defaults to the latest time recorded on the endpoint.
  d <- utils::read.csv(ex$data)
  expect_equal(v$horizon, max(d$TIME[d$CMT == 2]))

  expect_named(v$observed, c("STRATUM", "time", "n_risk", "n_event",
                             "n_censor", "surv", "cumhaz"))
  expect_named(v$band, c("STRATUM", "time", "obs", "lo", "med", "hi",
                         "n_used"))
  expect_identical(levels(v$strata$STRATUM), "all")
  expect_identical(v$strata$n_subjects, 30L)
  expect_true(is.na(v$strata$value))
  # The band is ordered wherever it is defined, and starts at S(0) = 1.
  b <- v$band[!is.na(v$band$lo), ]
  expect_true(all(b$lo <= b$med & b$med <= b$hi))
  expect_equal(v$band$time[1], 0)
  expect_equal(unlist(v$band[1, c("obs", "lo", "med", "hi")]),
               c(obs = 1, lo = 1, med = 1, hi = 1))
  expect_true(all(v$band$n_used <= 50L))
  # The automatic grid holds 101 evenly spaced points up to H.
  expect_true(all(vapply(v$horizon * (0:100) / 100, function(t)
    any(abs(v$band$time - t) < 1e-9), logical(1))))
  expect_output(print(v), "30 subjects, 19 events")
})

test_that("the observed curve matches survival::survfit row for row", {
  skip_if_not_installed("survival")
  ex <- ferx_example("tte_exponential")
  d <- utils::read.csv(ex$data)
  d <- d[d$CMT == 2, ]
  d$ARM <- d$ID %% 2L
  v <- vpc_tte_quiet(ex$model, d, n_sim = 20L, strata = "ARM")
  expect_identical(levels(v$strata$STRATUM), c("ARM=0", "ARM=1"))
  expect_equal(v$strata$value, c(0, 1))
  for (arm in 0:1) {
    sf <- survival::survfit(survival::Surv(TIME, DV) ~ 1,
                            data = d[d$ARM == arm, ])
    o <- v$observed[v$observed$STRATUM == paste0("ARM=", arm), ]
    expect_equal(o$time, sf$time)
    expect_equal(o$n_risk, as.integer(sf$n.risk))
    expect_equal(o$n_event, as.integer(sf$n.event))
    expect_equal(o$n_censor, as.integer(sf$n.censor))
    expect_equal(o$surv, sf$surv, tolerance = 1e-14)
    expect_equal(o$cumhaz, sf$cumhaz, tolerance = 1e-14)
  }
})

test_that("ties follow survfit: events before censorings at one time", {
  skip_if_not_installed("survival")
  d <- data.frame(ID = 1:8, TIME = c(2, 2, 2, 5, 5, 5, 9, 9),
                  DV = c(1, 1, 0, 1, 0, 0, 1, 0), EVID = 0L, AMT = 0,
                  CMT = 2L, MDV = 0L)
  v <- vpc_tte_quiet(vpc_tte_model(), d, n_sim = 20L)
  sf <- survival::survfit(survival::Surv(TIME, DV) ~ 1, data = d)
  expect_equal(v$observed$surv, sf$surv, tolerance = 1e-14)
  expect_equal(v$observed$cumhaz, sf$cumhaz, tolerance = 1e-14)
  expect_equal(v$observed$n_risk, as.integer(sf$n.risk))
})

test_that("the seed, censoring, horizon and grid reach the engine", {
  ex <- ferx_example("tte_exponential")
  a <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, seed = 3L)
  b <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, seed = 3L)
  c <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, seed = 4L)
  expect_identical(a$band, b$band)
  expect_false(identical(a$band$med, c$band$med))

  # Every censoring in this example is at H, so censoring the replicates
  # like the data is censoring them at H.
  h <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, seed = 3L,
                     censoring = "horizon")
  expect_identical(h$censoring, "horizon")
  expect_identical(h$band, a$band)
  # With dropouts before H the two differ; the observed curve does not.
  m <- vpc_tte_model()
  dd <- vpc_tte_data()
  early <- dd$DV == 0L | seq_len(nrow(dd)) %% 3L == 0L
  dd$TIME[early] <- round(dd$TIME[early] / 2, 2)
  dd$DV[early] <- 0L
  ho <- vpc_tte_quiet(m, dd, n_sim = 30L, seed = 3L)
  hh <- vpc_tte_quiet(m, dd, n_sim = 30L, seed = 3L, censoring = "horizon")
  expect_false(identical(ho$band, hh$band))
  expect_identical(ho$observed, hh$observed)

  g <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, grid = c(20, 0, 10, 10),
                     horizon = 40)
  expect_equal(g$horizon, 40)
  expect_equal(g$band$time, c(0, 10, 20))

  s <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, grid = "all_steps")
  expect_gt(nrow(s$band), nrow(a$band))
})

test_that("by_dose stratifies on the first dose", {
  ex <- ferx_example("tte_exponential")
  v <- vpc_tte_quiet(ex$model, ex$data, n_sim = 20L, by_dose = TRUE)
  expect_identical(levels(v$strata$STRATUM), "dose=0")
  expect_equal(v$strata$value, 0)
  expect_error(ferx_vpc_tte(ex$model, ex$data, strata = "ARM", by_dose = TRUE),
               "not both")
})

test_that("a fit's parameters drive the replicates", {
  ex <- ferx_example("tte_exponential")
  fit <- suppressWarnings(suppressMessages(
    ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  ))
  # data defaults to the fit's.
  vf <- vpc_tte_quiet(ex$model, n_sim = 30L, seed = 2L, fit = fit)
  v0 <- vpc_tte_quiet(ex$model, ex$data, n_sim = 30L, seed = 2L)
  expect_identical(vf$observed, v0$observed)
  expect_false(identical(vf$band$med, v0$band$med))
})

test_that("a competing-risks model needs cmt, and takes it", {
  cr <- ferx_example("tte_competing_risks")
  expect_vpc_tte_refusal(ferx_vpc_tte(cr$model, cr$data, n_sim = 20L),
                         "set `cmt`")
  v <- vpc_tte_quiet(cr$model, cr$data, n_sim = 20L, cmt = 3L)
  expect_identical(v$cmt, 3L)
  expect_vpc_tte_refusal(ferx_vpc_tte(cr$model, cr$data, n_sim = 20L, cmt = 1L),
                         "not a TTE endpoint")
})

test_that("each refusal comes back as E_VPC_TTE_INPUT", {
  m <- vpc_tte_model()
  d <- vpc_tte_data()

  expect_vpc_tte_refusal(ferx_vpc_tte(m, d, n_sim = 19L), "too few")
  expect_vpc_tte_refusal(ferx_vpc_tte(m, d, n_sim = 20L, probs = c(0.5, 0.05, 0.95)),
                         "strictly")
  expect_vpc_tte_refusal(ferx_vpc_tte(m, d, n_sim = 20L, grid = c(0, -1)),
                         "not a finite, non-negative")

  # Repeated events.
  rtte <- vpc_tte_model("  type   = rtte")
  dr <- data.frame(ID = rep(1:3, each = 2), TIME = c(3, 30, 5, 30, 8, 30),
                   DV = c(1, 0, 1, 0, 1, 0), EVID = 0L, AMT = 0, CMT = 2L,
                   MDV = 0L)
  expect_vpc_tte_refusal(ferx_vpc_tte(rtte, dr, n_sim = 20L), "repeated-events")

  # Left truncation (delayed entry through TENTRY).
  dl <- d
  dl$TENTRY <- 0
  dl$TENTRY[1] <- min(1, dl$TIME[1] / 2)
  expect_vpc_tte_refusal(ferx_vpc_tte(m, dl, n_sim = 20L), "Left truncation")

  # Interval censoring: a DV = 0 row followed by its DV = 2 right bound.
  di <- rbind(d, data.frame(ID = 31L, TIME = c(4, 9), DV = c(0L, 2L), EVID = 0L,
                            AMT = 0, CMT = 2L, MDV = 0L))
  expect_vpc_tte_refusal(ferx_vpc_tte(m, di, n_sim = 20L), "interval-censored")

  # Two records on the endpoint for one subject of a single-event model.
  d2 <- rbind(d, data.frame(ID = 1L, TIME = 40, DV = 1L, EVID = 0L, AMT = 0,
                            CMT = 2L, MDV = 0L))
  expect_vpc_tte_refusal(ferx_vpc_tte(m, d2[order(d2$ID, d2$TIME), ], n_sim = 20L),
                         "records on cmt 2")

  # A stratification covariate the data lacks, and one with too many values.
  expect_vpc_tte_refusal(ferx_vpc_tte(m, d, n_sim = 20L, strata = "ARM"),
                         "no covariate `ARM`")
  d11 <- d
  d11$SITE <- d11$ID %% 11L
  expect_vpc_tte_refusal(ferx_vpc_tte(m, d11, n_sim = 20L, strata = "SITE"),
                         "11 distinct values")
})

test_that("arguments are checked before the engine is called", {
  ex <- ferx_example("tte_exponential")
  expect_error(ferx_vpc_tte(ex$model, ex$data, n_sim = 2.5), "whole number")
  expect_error(ferx_vpc_tte(ex$model, ex$data, censoring = "none"))
  expect_error(ferx_vpc_tte(ex$model, ex$data, probs = c(0.1, 0.9)), "three")
  expect_error(ferx_vpc_tte(ex$model, ex$data, grid = "fine"), "grid")
  expect_error(ferx_vpc_tte(ex$model, ex$data, horizon = -1), "horizon")
  expect_error(ferx_vpc_tte(ex$model, ex$data, strata = c("A", "B")), "strata")
})

test_that("plot() draws with base graphics and with vpc", {
  ex <- ferx_example("tte_exponential")
  d <- utils::read.csv(ex$data)
  d$ARM <- d$ID %% 2L
  v <- vpc_tte_quiet(ex$model, d, n_sim = 20L, strata = "ARM")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(v, engine = "base"), v)

  skip_if_not_installed("vpc")
  p <- plot(v, engine = "vpc")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  # One facet per stratum, labelled as the engine labels it.
  expect_identical(as.character(built$layout$layout$STRATUM),
                   c("ARM=0", "ARM=1"))
  # The ribbon is the band: bin_min/bin_max are consecutive grid times.
  db <- ferx:::.ferx_vpc_tte_vpcdb(v)
  s0 <- db$sim_km[db$sim_km$strat == "ARM=0", ]
  b0 <- v$band[v$band$STRATUM == "ARM=0" & !is.na(v$band$lo), ]
  expect_equal(s0$bin_min, b0$time)
  expect_equal(s0$qmin, b0$lo)
  expect_equal(s0$qmed, b0$med)
  expect_equal(s0$qmax, b0$hi)
  expect_equal(utils::head(s0$bin_max, -1), utils::tail(s0$bin_min, -1))
})
