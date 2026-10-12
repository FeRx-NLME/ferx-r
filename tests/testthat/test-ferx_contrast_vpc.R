# ferx_contrast_vpc() (#297): the VPC of the within-study treatment contrast.

mbma_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      ex <- ferx_example("mbma_placebo")
      fit <<- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
    }
    fit
  }
})

mbma_rows <- function() {
  utils::read.csv(ferx_example("mbma_placebo")$data, stringsAsFactors = FALSE)
}

write_rows <- function(d) {
  p <- tempfile(fileext = ".csv")
  utils::write.csv(d, p, row.names = FALSE)
  p
}

# Every warning a call raises, and its value.
collect_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("(a) pairs every active arm with its own study's reference arm at the same time", {
  fit <- mbma_fit()
  d <- mbma_rows()
  # Two studies, rows reversed so the reference arm comes last: pairing is by
  # (study, TIME), not by position in the file.
  small <- d[d$ID %in% c(1, 2), ]
  small <- small[rev(seq_len(nrow(small))), ]
  cv <- ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 20L,
                          data = write_rows(small))

  act <- small[small$DOSE != 0, ]
  ref <- small[small$DOSE == 0, c("ID", "TIME", "DV")]
  names(ref)[3L] <- "DV_REF"
  want <- merge(act, ref, by = c("ID", "TIME"))
  want$CONTRAST <- want$DV - want$DV_REF
  got <- merge(cv$pairs, want, by = c("ID", "TIME", "DOSE"),
               suffixes = c("", ".want"))
  expect_identical(nrow(got), nrow(act))
  expect_identical(nrow(cv$pairs), nrow(act))
  expect_equal(got$CONTRAST, got$CONTRAST.want, tolerance = 1e-12)
  expect_equal(got$DV_REF, got$DV_REF.want, tolerance = 1e-12)
  expect_equal(got$DV_ACTIVE, got$DV, tolerance = 1e-12)
  # The arm is named by its `by` value when no `arm` column is given.
  expect_identical(as.character(cv$pairs$ARM), as.character(cv$pairs$DOSE))
  # At week 0 the drug effect is 0: every simulated contrast median is
  # centred near 0, while later weeks are pulled below it by EMAX.
  expect_lt(max(abs(cv$pairs$SIM_MED[cv$pairs$TIME == 0])), 1)
  expect_true(all(cv$pairs$SIM_LO <= cv$pairs$SIM_MED &
                    cv$pairs$SIM_MED <= cv$pairs$SIM_HI))
})

test_that("(a) `arm =` tells apart two arms sharing a `by` value", {
  fit <- mbma_fit()
  d <- mbma_rows()
  small <- d[d$ID == 1, ]
  # Study 1's 200 mg arm relabelled 50 mg: two arms at one dose.
  small$DOSE[small$DOSE == 200] <- 50
  expect_error(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 5L,
                      data = write_rows(small)),
    "repeat an arm at the same ID and TIME.*pass `arm =`"
  )
  cv <- ferx_contrast_vpc(fit, by = "DOSE", ref = 0, arm = "OCC", n_sim = 5L,
                          data = write_rows(small))
  expect_identical(nrow(cv$pairs), sum(small$DOSE != 0))
  expect_setequal(cv$pairs$ARM, c(2L, 3L, 4L))
})

test_that("a study with two reference arms is refused, not guessed between", {
  fit <- mbma_fit()
  d <- mbma_rows()
  small <- d[d$ID == 1, ]
  small$DOSE[small$OCC == 2] <- 0
  # Without `arm` the two placebo arms are one `by` value at one time...
  expect_error(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 5L,
                      data = write_rows(small)),
    "pass `arm =`"
  )
  # ... with it, they are two reference arms.
  err <- expect_error(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, arm = "OCC", n_sim = 5L,
                      data = write_rows(small)),
    "1 ID level has more than one reference arm \\(DOSE == 0\\): 1\\."
  )
  expect_match(conditionMessage(err),
               "exactly one reference per study; drop or merge", fixed = FALSE)
})

test_that("`ref` must be a value of `by`; argument checks", {
  fit <- mbma_fit()
  expect_error(ferx_contrast_vpc(fit, by = "DOSE", ref = 7, n_sim = 5L),
               "No observation record has DOSE == 7")
  expect_error(ferx_contrast_vpc(fit, by = "NOPE", ref = 0, n_sim = 5L),
               "Column not found in the data: NOPE")
  expect_error(ferx_contrast_vpc(fit, by = "ID", ref = 0, within = "ID"),
               "four different columns")
  expect_error(ferx_contrast_vpc(fit, by = "DOSE", ref = 0, probs = c(0.5, 0.05, 0.95)),
               "three increasing probabilities")
  expect_error(ferx_contrast_vpc(fit, by = "DOSE", ref = 0, scale = "sqrt"),
               "should be one of")
  expect_error(ferx_contrast_vpc(list(), by = "DOSE", ref = 0), "ferx_fit object")
})

test_that("(b) study-level random effects cancel: a 100x study omega widens the arm VPC, not the contrast", {
  ex <- ferx_example("mbma_placebo")
  fit <- mbma_fit()
  big <- fit
  big$omega <- fit$omega * 100   # ETA_E0 SD x 10
  n_sim <- 500L

  cv0 <- ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = n_sim, seed = 3L)
  cv1 <- ferx_contrast_vpc(big, by = "DOSE", ref = 0, n_sim = n_sim, seed = 3L)
  w0 <- cv0$pairs$SIM_HI - cv0$pairs$SIM_LO
  w1 <- cv1$pairs$SIM_HI - cv1$pairs$SIM_LO
  # Measured at 594151d: the simulated contrasts agree to 3.6e-14 (the eta
  # cancels up to rounding), band-width ratio exactly 1 on every pair.
  expect_lt(max(abs(cv0$sim$CONTRAST - cv1$sim$CONTRAST)), 1e-9)
  expect_equal(w1 / w0, rep(1, length(w0)), tolerance = 1e-9)

  # The arm-level 90% band of the same simulation: measured width ratio
  # median 8.69, range 6.65 - 9.75.
  arm_width <- function(f) {
    s <- ferx_simulate(ex$model, ex$data, n_sim = n_sim, seed = 3L, fit = f)
    m <- matrix(s$DV_SIM, ncol = n_sim)
    apply(m, 1L, function(x) diff(stats::quantile(x, c(0.05, 0.95))))
  }
  r_arm <- arm_width(big) / arm_width(fit)
  expect_gt(min(r_arm), 5)
  expect_gt(stats::median(r_arm), 7)
})

test_that("(c) single-arm studies and missing reference times: one warning, excluded", {
  fit <- mbma_fit()
  d <- mbma_rows()
  # Study 2: a single active arm (placebo and the other arms dropped).
  # Study 3: active arms only. Study 1: no placebo record at week 2.
  first_active_2 <- min(d$OCC[d$ID == 2 & d$DOSE != 0])
  keep <- !(d$ID == 2 & d$OCC != first_active_2) &
    !(d$ID == 3 & d$DOSE == 0) &
    !(d$ID == 1 & d$DOSE == 0 & d$TIME == 2)
  dd <- d[keep, ]
  n_cells <- length(unique(paste(dd$ID, dd$TIME)[dd$ID %in% c(2, 3)])) + 1L

  res <- collect_warnings(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 20L,
                      data = write_rows(dd))
  )
  expect_length(res$warnings, 1L)
  w <- res$warnings
  expect_match(w, paste0(n_cells, " ID x TIME cells have no reference arm ",
                         "\\(DOSE == 0\\) and contribute no contrast"))
  expect_match(w, "2 ID levels have no reference arm at all \\(2, 3\\)")
  expect_match(w, "1 more lacks one at some times \\(1\\)")

  cv <- res$value
  expect_setequal(cv$excluded, c("1", "2", "3"))
  expect_true(any(grepl("no reference arm (excluded): 1, 2, 3",
                        capture.output(print(cv)), fixed = TRUE)))
  expect_false(any(cv$pairs$ID %in% c(2, 3)))
  expect_false(any(cv$pairs$ID == 1 & cv$pairs$TIME == 2))
  # Study 1's other weeks still contribute.
  expect_setequal(unique(cv$pairs$TIME[cv$pairs$ID == 1]),
                  setdiff(unique(d$TIME[d$ID == 1]), 2))
})

test_that("(c) no reference arm anywhere is an error, not an empty result", {
  fit <- mbma_fit()
  d <- mbma_rows()
  expect_error(
    suppressWarnings(ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 5L,
                                       data = write_rows(d[d$DOSE != 0, ]))),
    "No observation record has DOSE == 0"
  )
  only_ref_t0 <- d[!(d$DOSE == 0 & d$TIME > 0) & !(d$DOSE != 0 & d$TIME == 0), ]
  expect_error(
    suppressWarnings(ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 5L,
                                       data = write_rows(only_ref_t0))),
    "No study has both a reference arm \\(DOSE == 0\\) and an active arm"
  )
})

test_that("(d) one row per study, time and active arm; the VPC table per dose and time", {
  fit <- mbma_fit()
  d <- mbma_rows()
  cv <- ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 50L)
  expect_s3_class(cv, "ferx_contrast_vpc")
  expect_identical(nrow(cv$pairs), sum(d$DOSE != 0))
  expect_false(anyDuplicated(cv$pairs[, c("ID", "TIME", "ARM")]) > 0L)
  expect_identical(
    names(cv$pairs),
    c("ID", "TIME", "ARM", "DOSE", "DV_REF", "DV_ACTIVE", "CONTRAST",
      "SIM_LO", "SIM_MED", "SIM_HI", "PCTL")
  )
  expect_true(all(cv$pairs$PCTL >= 0 & cv$pairs$PCTL <= 1))
  expect_identical(nrow(cv$sim), 50L * nrow(cv$pairs))
  expect_length(cv$excluded, 0L)
  expect_false(cv$sim_exponentiated)

  cells <- unique(d[d$DOSE != 0, c("DOSE", "TIME")])
  expect_identical(nrow(cv$vpc), nrow(cells))
  expect_identical(sum(cv$vpc$N), nrow(cv$pairs))
  i <- which(cv$vpc$DOSE == 50 & cv$vpc$TIME == 12)
  expect_equal(cv$vpc$OBS_MED[i],
               stats::median(cv$pairs$CONTRAST[cv$pairs$DOSE == 50 &
                                                 cv$pairs$TIME == 12]))
  expect_true(all(cv$vpc$SIM_LO <= cv$vpc$SIM_MED &
                    cv$vpc$SIM_MED <= cv$vpc$SIM_HI))

  # The band is the quantiles, across replicates, of the median over the
  # cell's pairs - recomputed here from the returned simulations.
  in_cell <- cv$sim$DOSE == 50 & cv$sim$TIME == 12
  med <- tapply(cv$sim$CONTRAST[in_cell], cv$sim$SIM[in_cell], stats::median)
  expect_equal(unname(unlist(cv$vpc[i, c("SIM_LO", "SIM_MED", "SIM_HI")])),
               unname(stats::quantile(med, c(0.05, 0.5, 0.95))))
  # PCTL is the fraction of a pair's simulated contrasts below the observed.
  j <- which.max(abs(cv$pairs$PCTL - 0.5))
  s <- cv$sim$CONTRAST[cv$sim$ID == cv$pairs$ID[j] & cv$sim$TIME == cv$pairs$TIME[j] &
                         cv$sim$ARM == cv$pairs$ARM[j]]
  expect_length(s, 50L)
  expect_equal(cv$pairs$PCTL[j], mean(s < cv$pairs$CONTRAST[j]))
  expect_equal(cv$pairs$SIM_MED[j], stats::median(s))
})

test_that("an arm repeated at one time, and an observed contrast off the scale's domain, are refused", {
  fit <- mbma_fit()
  d <- mbma_rows()
  small <- d[d$ID == 1, ]
  expect_error(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, arm = "OCC", n_sim = 5L,
                      data = write_rows(rbind(small, small[2L, ]))),
    "1 observation record\\(s\\) repeat an arm at the same ID and TIME \\(first: ID 1, TIME 2, arm 1\\)\\. Each arm must have one record per TIME\\."
  )
  small$DV[small$DOSE == 25 & small$TIME == 8] <- -1
  expect_error(
    ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 5L, scale = "log",
                      data = write_rows(small)),
    "1 observed contrast is missing or non-finite on scale = \"log\" \\(a missing DV, or a DV outside the transform's domain\\)\\. Drop those records, or choose another `scale`\\."
  )
})

test_that("the observed records are the rows ferx_simulate() returns: no doses, no MDV = 1", {
  rows <- function(d) {
    p <- tempfile(fileext = ".csv")
    utils::write.csv(d, p, row.names = FALSE, na = ".")
    ferx:::.ferx_cvpc_obs_rows(p, NULL, character(0))
  }
  # EVID without MDV: the dose records go on EVID alone.
  d <- data.frame(ID = 1, TIME = 0:3, DV = c(NA, 1, 2, 3), EVID = c(1, 0, 0, 0),
                  AMT = c(100, NA, NA, NA))
  expect_identical(rows(d)$TIME, 1:3)
  # No EVID column: a record with an AMT is a dose.
  expect_identical(rows(d[, c("ID", "TIME", "DV", "AMT")])$TIME, 1:3)
  # MDV = 1 is excluded even on an EVID = 0 record.
  d$MDV <- c(1, 0, 1, 0)
  expect_identical(rows(d)$TIME, c(1L, 3L))
})

test_that("simulated records are paired by position only after checking they line up", {
  obs <- data.frame(ID = c(1, 1, 2), TIME = c(0, 1, 0), DV = 1:3)
  sim <- data.frame(SIM = rep(1:2, each = 3L), ID = c("1", "1", "2"),
                    TIME = c(0, 1, 0), DV_SIM = 1:6, OBSERVED = NA)
  m <- ferx:::.ferx_cvpc_sim_matrix(sim, obs, 2L)
  expect_identical(dim(m), c(3L, 2L))
  expect_identical(m[, 2L], 4:6)
  swapped <- sim
  swapped$ID[1:3] <- c("1", "2", "1")
  expect_error(ferx:::.ferx_cvpc_sim_matrix(swapped, obs, 2L),
               "do not line up.*first mismatch at record 2: data ID 1 TIME 1, simulated ID 2")
  expect_error(ferx:::.ferx_cvpc_sim_matrix(sim[-6L, ], obs, 2L),
               "returned 5 continuous records for 2 replicates of 3 observation records")
})

test_that("print() and plot() describe and draw the contrast", {
  fit <- mbma_fit()
  cv <- ferx_contrast_vpc(fit, by = "DOSE", ref = 0, n_sim = 20L)
  out <- capture.output(print(cv))
  expect_match(out[2L], "contrast: each DOSE vs DOSE == 0, within ID, on the observed DV scale",
               fixed = TRUE)
  expect_match(out[3L], "64 pairs from 6 ID levels; 20 replicates; band 5-95%",
               fixed = TRUE)
  n_out <- sum(cv$vpc$OBS_MED < cv$vpc$SIM_LO | cv$vpc$OBS_MED > cv$vpc$SIM_HI)
  expect_identical(out[4L], paste0("  observed median outside the band in ",
                                   n_out, " of 24 DOSE x TIME cells"))
  expect_length(out, 4L)
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(cv), cv)
  expect_identical(plot(cv, facet = "none"), cv)
})

test_that("a log(DV) readout: DV_SIM is put back on the data's scale before the contrast", {
  # warfarin_ltbs simulates on the log scale while the data hold DV on the
  # natural scale. Paired subjects (same sampling times) stand in for the
  # two arms of a study.
  ex <- ferx_example("warfarin_ltbs")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                   verbose = FALSE))
  d <- utils::read.csv(ex$data, na.strings = c(".", "NA"))
  d$STUDY <- ceiling(d$ID / 2)
  d$TRT <- (d$ID + 1) %% 2
  p <- tempfile(fileext = ".csv")
  utils::write.csv(d, p, row.names = FALSE, na = ".")
  covered <- function(cv) {
    mean(cv$pairs$CONTRAST >= cv$pairs$SIM_LO & cv$pairs$CONTRAST <= cv$pairs$SIM_HI)
  }

  obs <- ferx_contrast_vpc(fit, by = "TRT", ref = 0, within = "STUDY",
                           n_sim = 200L, data = p)
  expect_true(obs$sim_exponentiated)
  expect_equal(obs$pairs$CONTRAST, obs$pairs$DV_ACTIVE - obs$pairs$DV_REF)
  # Measured at 594151d: 52 of 55 observed contrasts (0.945) inside the 90%
  # band; without the back-transform 0.29 on this scale.
  expect_gt(covered(obs), 0.85)

  lg <- ferx_contrast_vpc(fit, by = "TRT", ref = 0, within = "STUDY",
                          n_sim = 200L, data = p, scale = "log")
  expect_equal(lg$pairs$CONTRAST, log(lg$pairs$DV_ACTIVE / lg$pairs$DV_REF))
  # Measured: 0.945; without the back-transform 0.745.
  expect_gt(covered(lg), 0.85)
  expect_match(capture.output(print(lg))[2L], "on the log scale (log ratio",
               fixed = TRUE)
})

test_that("only `log(DV) ~` is a readout the engine simulates off the data's scale", {
  detect <- function(lines) {
    m <- tempfile(fileext = ".ferx")
    writeLines(c("[error_model]", lines), m)
    ferx:::.ferx_cvpc_log_dv_readout(m)
  }
  expect_true(detect("  log(DV) ~ additive(ADD_LOG)"))
  expect_true(detect("  LOG( DV )~additive(ADD_LOG)"))
  expect_true(detect("  CMT=2: log(DV) ~ additive(ADD_LOG)"))
  expect_false(detect("  DV ~ log_additive(ADD_LOG)"))
  expect_false(detect("  DV ~ additive(ADD_ERR) weight = SE"))
  expect_false(detect("  # log(DV) ~ additive(ADD_LOG)"))
})

test_that("scale = 'expit' and 'logit' transform both sides alike", {
  f <- ferx:::.ferx_cvpc_transform
  expect_equal(f("expit")(0), 0.5)
  expect_equal(f("logit")(0.5), 0)
  expect_identical(f("observed")(c(1, -2)), c(1, -2))
  expect_true(is.nan(f("log")(-1)))
})
