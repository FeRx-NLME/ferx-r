# ferx_vpc_tte() on a repeated-events (RTTE) endpoint: core's vpc_rtte()
# (ferx-core #1886 PR2, #2908). Every number is core's; these tests check that
# the glue carries it into R unchanged: the observed curves against
# survival::survfit and a hand-computed mean cumulative function (the
# references core validates against), the shapes, the settings and the
# refusals.

rtte_model <- function() {
  m <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "",
    "[event_model]",
    "  cmt    = 2",
    "  type   = rtte",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), m)
  m
}

# An observed RTTE dataset drawn from the model itself: each subject's events
# on [0, 30], then a censoring record at 30.
rtte_data <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      design <- data.frame(ID = 1:60, TIME = 30, DV = 0, EVID = 0L, AMT = 0,
                           CMT = 2L, MDV = 0L)
      sim <- suppressWarnings(ferx_simulate(rtte_model(), design, n_sim = 1L,
                                            seed = 2L, horizon = 30))
      cache <<- data.frame(ID = as.integer(sim$ID), TIME = sim$TIME,
                           DV = sim$OBSERVED, EVID = 0L, AMT = 0, CMT = 2L,
                           MDV = 0L)
    }
    cache
  }
})

rtte_vpc <- function(...) suppressWarnings(ferx_vpc_tte(rtte_model(), ...))

test_that("an RTTE endpoint returns the RTTE VPC, in its shape", {
  d <- rtte_data()
  expect_true(any(table(d$ID[d$DV == 1]) > 1L))
  v <- rtte_vpc(d, n_sim = 30L, seed = 1L)
  expect_s3_class(v, "ferx_vpc_rtte")
  expect_identical(v$cmt, 2L)
  expect_identical(v$n_sim, 30L)
  expect_identical(v$censoring, "observed")
  expect_equal(v$horizon, 30)
  expect_named(v$mcf, c("STRATUM", "time", "n_risk", "n_event", "mcf"))
  for (b in c("mcf_band", "absolute_band", "gap_band")) {
    expect_true(all(c("STRATUM", "time", "obs", "lo", "med", "hi", "n_used")
                    %in% names(v[[b]])), info = b)
  }
  expect_named(v$absolute, c("STRATUM", "order", "time", "n_risk", "n_event",
                             "n_censor", "surv"))
  expect_named(v$gap, names(v$absolute))
  # Orders 1..K on absolute time and 2..K on gap time, one row each.
  K <- v$max_order
  expect_gte(K, 2L)
  o <- v$orders
  expect_identical(o$order[o$clock == "absolute"], seq_len(K))
  expect_identical(o$order[o$clock == "gap"], seq(2L, K))
  expect_true(all(o$n_subjects[o$clock == "absolute"] == 60L))
  expect_identical(sort(unique(v$absolute_band$order)), seq_len(K))
  # The bands are ordered where they are defined.
  for (b in c("mcf_band", "absolute_band", "gap_band")) {
    bb <- v[[b]][!is.na(v[[b]]$lo), ]
    expect_true(all(bb$lo <= bb$med & bb$med <= bb$hi), info = b)
    expect_true(all(bb$n_used <= 30L), info = b)
  }
  expect_output(print(v), "60 subjects, ")
  expect_output(print(v), "event orders 1-")
})

test_that("the observed curves match survfit and the MCF by hand", {
  skip_if_not_installed("survival")
  d <- rtte_data()
  v <- rtte_vpc(d, n_sim = 20L, seed = 1L)
  ends <- tapply(d$TIME, d$ID, max)
  ev <- split(d$TIME[d$DV == 1], d$ID[d$DV == 1])
  events_of <- function(id) sort(ev[[as.character(id)]] %||% numeric())
  ids <- as.integer(names(ends))

  # Absolute time, order k: the k-th event, else the follow-up end censored.
  for (k in seq_len(v$max_order)) {
    tk <- vapply(ids, function(i) {
      e <- events_of(i)
      if (length(e) >= k) e[k] else ends[[as.character(i)]]
    }, numeric(1))
    dk <- vapply(ids, function(i) length(events_of(i)) >= k, logical(1))
    sf <- survival::survfit(survival::Surv(tk, dk) ~ 1)
    o <- v$absolute[v$absolute$order == k, ]
    expect_equal(o$time, sf$time, info = paste("absolute", k))
    expect_equal(o$n_risk, as.integer(sf$n.risk), info = paste("absolute", k))
    expect_equal(o$surv, sf$surv, tolerance = 1e-12, info = paste("absolute", k))
  }
  # Gap time, order k >= 2: subjects with k - 1 events, from the (k-1)-th.
  for (k in seq(2L, v$max_order)) {
    keep <- ids[vapply(ids, function(i) length(events_of(i)) >= k - 1L, logical(1))]
    gk <- vapply(keep, function(i) {
      e <- events_of(i)
      (if (length(e) >= k) e[k] else ends[[as.character(i)]]) - e[k - 1L]
    }, numeric(1))
    dk <- vapply(keep, function(i) length(events_of(i)) >= k, logical(1))
    sf <- survival::survfit(survival::Surv(gk, dk) ~ 1)
    o <- v$gap[v$gap$order == k, ]
    expect_identical(v$orders$n_subjects[v$orders$clock == "gap" &
                                           v$orders$order == k], length(keep))
    expect_equal(o$time, sf$time, info = paste("gap", k))
    expect_equal(o$surv, sf$surv, tolerance = 1e-12, info = paste("gap", k))
  }
  # Mean cumulative function: sum over event times of dN(s) / Y(s), with
  # Y(s) the subjects whose follow-up end is at or after s. Core has one row per
  # distinct event time or follow-up end (a row with no event leaves it flat).
  et <- sort(unique(c(d$TIME[d$DV == 1], ends)))
  dn <- vapply(et, function(t) sum(d$TIME[d$DV == 1] == t), numeric(1))
  y <- vapply(et, function(t) sum(ends >= t), numeric(1))
  expect_equal(v$mcf$time, et)
  expect_equal(v$mcf$n_event, as.integer(dn))
  expect_equal(v$mcf$n_risk, as.integer(y))
  expect_equal(v$mcf$mcf, cumsum(dn / y), tolerance = 1e-12)
})

test_that("seed, censoring, max_order, horizon and strata reach the engine", {
  d <- rtte_data()
  a <- rtte_vpc(d, n_sim = 20L, seed = 3L)
  b <- rtte_vpc(d, n_sim = 20L, seed = 3L)
  c <- rtte_vpc(d, n_sim = 20L, seed = 4L)
  expect_identical(a$mcf_band, b$mcf_band)
  expect_false(identical(a$mcf_band$med, c$mcf_band$med))

  h <- rtte_vpc(d, n_sim = 20L, seed = 3L, censoring = "horizon")
  expect_identical(h$censoring, "horizon")
  expect_identical(h$mcf, a$mcf)

  k2 <- rtte_vpc(d, n_sim = 20L, seed = 3L, max_order = 2)
  expect_identical(k2$max_order, 2L)
  expect_identical(sort(unique(k2$absolute$order)), 1:2)
  expect_identical(unique(k2$gap$order), 2L)

  # The observed streams are cut at the horizon (#2908): nothing after it.
  cut <- rtte_vpc(d, n_sim = 20L, seed = 3L, horizon = 15)
  expect_equal(cut$horizon, 15)
  expect_true(all(cut$mcf$time <= 15))
  # Every follow-up now ends at 15, which adds a flat row there.
  expect_identical(max(cut$mcf$time), 15)
  before <- a$mcf$time < 15
  expect_equal(cut$mcf$mcf[seq_len(sum(before))], a$mcf$mcf[before])

  ds <- d
  ds$ARM <- ds$ID %% 2L
  s <- rtte_vpc(ds, n_sim = 20L, seed = 3L, strata = "ARM")
  expect_identical(levels(s$strata$STRATUM), c("ARM=0", "ARM=1"))
  expect_identical(s$strata$n_subjects, c(30L, 30L))
  expect_setequal(as.character(unique(s$mcf$STRATUM)), c("ARM=0", "ARM=1"))
})

test_that("a fit's parameters drive the RTTE replicates", {
  d <- rtte_data()
  m <- rtte_model()
  fit <- suppressWarnings(ferx_fit(m, d, covariance = FALSE, verbose = FALSE))
  vf <- suppressWarnings(ferx_vpc_tte(m, n_sim = 20L, seed = 2L, fit = fit))
  v0 <- rtte_vpc(d, n_sim = 20L, seed = 2L)
  expect_s3_class(vf, "ferx_vpc_rtte")
  expect_identical(vf$mcf, v0$mcf)
  expect_false(identical(vf$mcf_band$med, v0$mcf_band$med))
})

test_that("a subject ending on an event warns, with guidance", {
  d <- rtte_data()
  # Drop the closing censoring record of the first subject with an event: its
  # last record is then an event.
  id <- d$ID[d$DV == 1][1L]
  last <- max(which(d$ID == id))
  expect_identical(d$DV[last], 0)
  d <- d[-last, ]
  expect_identical(d$DV[max(which(d$ID == id))], 1)
  expect_warning(v <- ferx_vpc_tte(rtte_model(), d, n_sim = 20L, seed = 1L),
                 "W_VPC_RTTE_ENDS_ON_EVENT")
  w <- grep("^W_VPC_RTTE_ENDS_ON_EVENT", v$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "^W_VPC_RTTE_ENDS_ON_EVENT: 1 of 60 subjects on cmt 2")
  out <- capture.output(print(v))
  expect_true(any(grepl("W_VPC_RTTE_ENDS_ON_EVENT: A subject whose last record",
                        out, fixed = TRUE)))
})

test_that("max_order is refused on a single-event endpoint, and checked in R", {
  ex <- ferx_example("tte_exponential")
  expect_error(ferx_vpc_tte(ex$model, ex$data, n_sim = 20L, max_order = 2),
               "`max_order` applies to a repeated-events")
  expect_error(rtte_vpc(rtte_data(), n_sim = 20L, max_order = 0),
               "`max_order` must be a single whole number >= 1")
  e <- tryCatch(rtte_vpc(rtte_data(), n_sim = 19L), error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_VPC_TTE_INPUT")
})

test_that("plot() draws the MCF and the curves by order", {
  v <- rtte_vpc(rtte_data(), n_sim = 20L, seed = 1L, max_order = 2)
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot(v), v)
  expect_identical(plot(v, which = "mcf"), v)
  expect_identical(plot(v, which = "gap", orders = 2), v)
  expect_error(plot(v, which = "gap", orders = 9), "Nothing to draw")
})
