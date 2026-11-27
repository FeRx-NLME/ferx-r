# ferx_vpc_cif(): the cumulative-incidence VPC for competing risks (ferx-core
# #1887). Every number is core's; these tests check that the glue carries it
# into R unchanged, against survival::survfit's multi-state fit, and the
# shapes, settings and refusals.

vpc_cif_quiet <- function(...) suppressWarnings(ferx_vpc_cif(...))

# Competing-risks data (a path or a data frame), one row per subject: its time
# and its cause (0 censored, else the cmt of the cause that occurred).
cif_subjects <- function(data) {
  d <- if (is.character(data)) utils::read.csv(data) else data
  d <- d[d$EVID == 0 & d$CMT %in% c(2L, 3L), ]
  ids <- unique(d$ID)
  data.frame(
    ID = ids,
    time = vapply(ids, function(i) d$TIME[d$ID == i][1L], numeric(1)),
    cause = vapply(ids, function(i) {
      hit <- d$CMT[d$ID == i & d$DV == 1]
      if (length(hit)) hit else 0L
    }, integer(1))
  )
}

# The bundled data with censorings part-way through follow-up. The bundled
# file censors only at t = 14, after the last event, and on such data the
# Aalen-Johansen estimate is the plain proportion cumsum(n_event_k) / n, so
# an oracle on it cannot tell the two apart (#778 review, finding 2). Here the
# events of every fifth subject before t = 10 become censorings, and an ARM
# covariate gives two strata.
cif_censored_data <- function() {
  d <- utils::read.csv(ferx_example("tte_competing_risks")$data)
  d$DV[d$ID %% 5 == 0 & d$DV == 1 & d$TIME < 10] <- 0
  d$ARM <- d$ID %% 2
  d
}

# One stratum's observed rows against survival::survfit's multi-state fit of
# `s` (cif_subjects() rows): times exactly, risk sets and per-cause event
# counts identically, the cumulative incidences and the survival to 1e-12.
expect_cif_survfit <- function(o, s) {
  status <- factor(s$cause, levels = c(0L, 2L, 3L))
  sf <- survival::survfit(survival::Surv(s$time, status) ~ 1)
  # survfit keeps a row per distinct time, event or censoring, as core does.
  expect_equal(o$time, sf$time, tolerance = 0)
  expect_identical(o$n_risk, as.integer(sf$n.risk[, 1L]))
  for (k in c("2", "3")) {
    expect_equal(o[[paste0("cif_", k)]], unname(sf$pstate[, match(k, sf$states)]),
                 tolerance = 1e-12)
    expect_identical(o[[paste0("n_event_", k)]],
                     as.integer(sf$n.event[, match(k, sf$states)]))
  }
  expect_equal(o$surv, unname(sf$pstate[, match("(s0)", sf$states)]),
               tolerance = 1e-12)
  # And the all-cause survival is the plain Kaplan-Meier curve.
  km <- survival::survfit(survival::Surv(s$time, s$cause > 0) ~ 1)
  expect_equal(o$surv, km$surv, tolerance = 1e-12)
}

expect_vpc_cif_refusal <- function(expr, pattern) {
  e <- tryCatch(expr, error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_VPC_TTE_INPUT")
  expect_match(conditionMessage(e), pattern)
  invisible(e)
}

test_that("ferx_vpc_cif returns the observed curves and a band per cause", {
  ex <- ferx_example("tte_competing_risks")
  v <- vpc_cif_quiet(ex$model, ex$data, n_sim = 30L, seed = 1L)
  expect_s3_class(v, "ferx_vpc_cif")
  expect_identical(v$cmts, c(2L, 3L))
  expect_identical(v$n_sim, 30L)
  expect_identical(v$censoring, "observed")
  expect_equal(v$probs, c(0.05, 0.5, 0.95))
  d <- utils::read.csv(ex$data)
  expect_equal(v$horizon, max(d$TIME))
  expect_named(v$observed, c("STRATUM", "time", "n_risk", "n_event", "n_censor",
                             "surv", "cif_2", "cif_3", "n_event_2", "n_event_3"))
  expect_named(v$band, c("STRATUM", "cmt", "time", "obs", "lo", "med", "hi",
                         "n_used"))
  expect_named(v$survival_band, c("STRATUM", "time", "obs", "lo", "med", "hi",
                                  "n_used"))
  expect_identical(levels(v$strata$STRATUM), "all")
  expect_identical(v$strata$n_subjects, length(unique(d$ID)))
  # Each cause's band is on the same grid as the survival band.
  for (k in v$cmts) {
    expect_identical(v$band$time[v$band$cmt == k], v$survival_band$time)
  }
  b <- v$band[!is.na(v$band$lo), ]
  expect_true(all(b$lo <= b$med & b$med <= b$hi))
  # The observed identities: the causes' events make up every event, and the
  # incidences and the survival add up to 1.
  o <- v$observed
  expect_identical(o$n_event_2 + o$n_event_3, o$n_event)
  expect_equal(o$cif_2 + o$cif_3 + o$surv, rep(1, nrow(o)), tolerance = 1e-12)
  # At time 0 nothing has happened yet, in the data or any replicate.
  expect_equal(v$band$obs[v$band$time == 0], c(0, 0))
  expect_equal(v$survival_band$obs[v$survival_band$time == 0], 1)
})

test_that("the observed CIF equals survival::survfit's multi-state fit", {
  skip_if_not_installed("survival")
  ex <- ferx_example("tte_competing_risks")
  # The bundled data as shipped.
  v <- vpc_cif_quiet(ex$model, ex$data, n_sim = 20L, seed = 1L)
  expect_cif_survfit(v$observed, cif_subjects(ex$data))

  # Censorings before the last event: censored subjects leave the risk set
  # early, so the estimate is no longer the plain proportion.
  d <- cif_censored_data()
  s <- cif_subjects(d)
  expect_gt(sum(s$cause == 0L & s$time < 10), 5L)
  v <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 1L)
  expect_cif_survfit(v$observed, s)
  naive <- cumsum(v$observed$n_event_2) / nrow(s)
  expect_gt(max(abs(naive - v$observed$cif_2)), 0.05)

  # Per stratum: each stratum's curves are its own subjects' survfit fit.
  vs <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 1L, strata = "ARM")
  arm <- d$ARM[match(s$ID, d$ID)]
  for (a in c(0, 1)) {
    expect_cif_survfit(
      vs$observed[vs$observed$STRATUM == paste0("ARM=", a), , drop = FALSE],
      s[arm == a, , drop = FALSE]
    )
  }

  # A horizon below the follow-up: the data are read on [0, H], a record after
  # H being a censoring at H.
  h <- 8
  vh <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 1L, horizon = h)
  sh <- s
  sh$cause[sh$time > h] <- 0L
  sh$time <- pmin(sh$time, h)
  expect_cif_survfit(vh$observed, sh)
})

test_that("strata, horizon, censoring, probs and grid reach the engine", {
  ex <- ferx_example("tte_competing_risks")
  d <- utils::read.csv(ex$data)
  d$ARM <- d$ID %% 2
  v <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 3L, strata = "ARM",
                     censoring = "horizon", horizon = 10,
                     probs = c(0.1, 0.5, 0.9), grid = c(0, 2, 4, 8))
  expect_identical(levels(v$strata$STRATUM), c("ARM=0", "ARM=1"))
  expect_equal(v$strata$value, c(0, 1))
  expect_identical(sum(v$strata$n_subjects), length(unique(d$ID)))
  expect_identical(v$censoring, "horizon")
  expect_equal(v$horizon, 10)
  expect_equal(v$probs, c(0.1, 0.5, 0.9))
  expect_equal(unique(v$survival_band$time), c(0, 2, 4, 8))
  # The data are read on [0, H] only.
  expect_true(all(v$observed$time <= 10))
  # The seed reproduces the band; another seed moves it.
  v2 <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 3L, strata = "ARM",
                      censoring = "horizon", horizon = 10,
                      probs = c(0.1, 0.5, 0.9), grid = c(0, 2, 4, 8))
  expect_identical(v2$band, v$band)
  v3 <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 4L, strata = "ARM",
                      censoring = "horizon", horizon = 10,
                      probs = c(0.1, 0.5, 0.9), grid = c(0, 2, 4, 8))
  expect_false(identical(v3$band$med, v$band$med))
  expect_identical(v3$observed, v$observed)
})

test_that("one cause is ferx_vpc_tte() upside down", {
  ex <- ferx_example("tte_exponential")
  cv <- vpc_cif_quiet(ex$model, ex$data, n_sim = 30L, seed = 5L)
  tv <- suppressWarnings(ferx_vpc_tte(ex$model, ex$data, n_sim = 30L, seed = 5L))
  expect_identical(cv$cmts, tv$cmt)
  expect_equal(cv$observed$time, tv$observed$time)
  expect_equal(cv$observed$surv, tv$observed$surv, tolerance = 1e-12)
  expect_equal(cv$observed[[paste0("cif_", tv$cmt)]], 1 - tv$observed$surv,
               tolerance = 1e-12)
  # The same replicates (core keys the censoring draws alike), so the CIF band
  # is the Kaplan-Meier band reflected.
  expect_equal(cv$band$time, tv$band$time)
  expect_equal(cv$band$lo, 1 - tv$band$hi, tolerance = 1e-12)
  expect_equal(cv$band$med, 1 - tv$band$med, tolerance = 1e-12)
  expect_equal(cv$survival_band$med, tv$band$med, tolerance = 1e-12)
})

test_that("fit = simulates at the fitted parameters", {
  ex <- ferx_example("tte_competing_risks")
  fit <- suppressWarnings(suppressMessages(
    ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE)
  ))
  vf <- vpc_cif_quiet(ex$model, n_sim = 20L, seed = 2L, fit = fit)
  v0 <- vpc_cif_quiet(ex$model, ex$data, n_sim = 20L, seed = 2L)
  expect_identical(vf$observed, v0$observed)
  expect_false(identical(vf$band$med, v0$band$med))
})

test_that("each refusal comes back as E_VPC_TTE_INPUT", {
  ex <- ferx_example("tte_competing_risks")
  d <- utils::read.csv(ex$data)

  expect_vpc_cif_refusal(ferx_vpc_cif(ex$model, d, n_sim = 19L), "too few")
  # A subject with a record on one cause only.
  d1 <- d[!(d$ID == d$ID[1L] & d$CMT == 3L), ]
  expect_vpc_cif_refusal(ferx_vpc_cif(ex$model, d1, n_sim = 20L),
                         "no record on cmt 3")
  # Events on both causes.
  d2 <- d
  d2$DV[d2$ID == d2$ID[1L]] <- 1
  expect_vpc_cif_refusal(ferx_vpc_cif(ex$model, d2, n_sim = 20L),
                         "an event on cmt 2 and on cmt 3")
  # Cause records at different times.
  d3 <- d
  r <- which(d3$ID == d3$ID[1L] & d3$CMT == 3L)
  d3$TIME[r] <- d3$TIME[r] / 2
  expect_vpc_cif_refusal(ferx_vpc_cif(ex$model, d3, n_sim = 20L),
                         "Every cause is followed to the same time")
  # No TTE endpoint at all.
  w <- ferx_example("warfarin")
  expect_vpc_cif_refusal(ferx_vpc_cif(w$model, w$data, n_sim = 20L),
                         "no TTE endpoint")
})

test_that("a repeated-events cause is refused, naming ferx_vpc_tte()", {
  m <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVLAMBDA(0.1, 0.001, 10.0)",
    "  omega ETA_LAMBDA ~ 0.09",
    "[event_model]",
    "  cmt    = 2",
    "  type   = rtte",
    "  family = exponential",
    "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"
  ), m)
  d <- data.frame(ID = 1:30, TIME = 30, DV = 0, EVID = 0L, AMT = 0, CMT = 2L,
                  MDV = 0L)
  e <- expect_vpc_cif_refusal(ferx_vpc_cif(m, d, n_sim = 20L),
                              "repeated-events")
  expect_match(conditionMessage(e), "`ferx_vpc_tte()`", fixed = TRUE)
  expect_false(grepl("`vpc_rtte`", conditionMessage(e), fixed = TRUE))
})

test_that("the glue-side argument checks refuse before the engine", {
  ex <- ferx_example("tte_competing_risks")
  expect_error(ferx_vpc_cif(ex$model, ex$data, n_sim = 0L), "`n_sim`")
  expect_error(ferx_vpc_cif(ex$model, ex$data, strata = "ARM", by_dose = TRUE),
               "not both")
  expect_error(ferx_vpc_cif(ex$model, ex$data, horizon = -1), "`horizon`")
  expect_error(ferx_vpc_cif(ex$model, ex$data, grid = "every"), "`grid`")
  expect_error(ferx_vpc_cif(ex$model, ex$data, probs = c(0.1, 0.9)), "`probs`")
  expect_false("cmt" %in% names(formals(ferx_vpc_cif)))
  # The glue's own check, which the R checks above keep out of reach, names
  # the R function it serves (#778 review, finding 4).
  e <- tryCatch(
    ferx:::ferx_rust_vpc_cif(normalizePath(ex$model), normalizePath(ex$data),
                             20L, 1L, c(0.05, 0.5, 0.95), "bogus", numeric(0),
                             "observed", "none", "", -1),
    error = function(e) e
  )
  expect_match(conditionMessage(e), "internal: ferx_vpc_cif got grid `bogus`",
               fixed = TRUE)
})

test_that("print() lists the settings and the events per cause", {
  ex <- ferx_example("tte_competing_risks")
  v <- vpc_cif_quiet(ex$model, ex$data, n_sim = 20L, seed = 1L)
  s <- cif_subjects(ex$data)
  out <- capture.output(print(v))
  expect_match(out[1L], "Cumulative-incidence VPC", fixed = TRUE)
  expect_true(any(grepl("causes cmt 2, 3", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("all: %d subjects, events %d on cmt 2, %d on cmt 3",
                                nrow(s), sum(s$cause == 2L), sum(s$cause == 3L)),
                        out, fixed = TRUE)))
})

test_that("plot() draws a panel per stratum and cause, survival optional", {
  ex <- ferx_example("tte_competing_risks")
  d <- utils::read.csv(ex$data)
  d$ARM <- d$ID %% 2
  v <- vpc_cif_quiet(ex$model, d, n_sim = 20L, seed = 1L, strata = "ARM")
  panels <- 0L
  testthat::local_mocked_bindings(
    .ferx_vpc_step_panel = function(...) panels <<- panels + 1L,
    .package = "ferx"
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(v), v)
  expect_identical(panels, 4L)
  panels <- 0L
  plot(v, survival = TRUE)
  expect_identical(panels, 6L)
  panels <- 0L
  plot(v, cmts = 3)
  expect_identical(panels, 2L)
  expect_error(plot(v, cmts = 4), "`cmts`")
})

test_that("a panel with no band says so", {
  ex <- ferx_example("tte_competing_risks")
  # Every grid time past the horizon: no replicate is under follow-up there.
  v <- vpc_cif_quiet(ex$model, ex$data, n_sim = 20L, seed = 1L,
                     grid = c(100, 200))
  expect_true(all(is.na(v$band$lo)))
  expect_true(all(is.na(v$survival_band$lo)))
  # Core's vpc_cif raises no W_VPC_EMPTY_STRATUM at this pin (#778 review,
  # finding 1), so print() has no reason to show. When this starts failing,
  # core raises it: give CIF panels the shared note back and drop the CIF one.
  expect_false(any(grepl("^W_VPC_EMPTY_STRATUM", v$warnings)))
  # So the panel's note says the reason itself and cites no warning.
  label <- ferx:::.ferx_vpc_cif_no_band_label
  expect_false(grepl("W_VPC_EMPTY_STRATUM", label, fixed = TRUE))
  expect_false(grepl("$warnings", label, fixed = TRUE))
  expect_match(label, "too few replicates under follow-up", fixed = TRUE)
  notes <- character()
  testthat::local_mocked_bindings(
    .ferx_vpc_no_band_note = function(label = "shared") notes <<- c(notes, label),
    .package = "ferx"
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  plot(v, survival = TRUE)
  # Two cause panels and the survival panel, each with the CIF note.
  expect_identical(notes, rep(label, 3L))
})

test_that("the vpc_tte / vpc_kmmc refusal for several causes names ferx_vpc_cif()", {
  ex <- ferx_example("tte_competing_risks")
  for (call in list(
    quote(ferx_vpc_tte(ex$model, ex$data, n_sim = 20L)),
    quote(ferx_vpc_kmmc(ex$model, ex$data, covariate = "ID", n_sim = 20L))
  )) {
    e <- expect_vpc_cif_refusal(eval(call), "set `cmt`")
    expect_match(conditionMessage(e), "`ferx_vpc_cif()`", fixed = TRUE)
    expect_false(grepl("`vpc_cif`", conditionMessage(e), fixed = TRUE))
  }
})
