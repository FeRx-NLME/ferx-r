
# ---- header from test-npde.R ----
# Tests for ferx_calc_npde(): post-hoc simulation-based NPDE/NPD from a fit.
#
# The argument-guard tests are pure R and run without the compiled engine: each
# error is raised before the FFI call. The end-to-end tests need the real
# backend and the bundled warfarin example (shared `warfarin_fit()` helper).

# --- argument validation (no engine needed) ----------------------------------


# A minimal structurally-valid fit that passes validate_fit_for_params and the
# sdtab gate but carries no model/data path, so the nsim/seed/file guards fire
# before any FFI call. Its bytes are never parsed: the guards come first, or the
# test mocks `.ferx_fit_canonical()`.
.stub_fit <- function() {
  list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2), sigma = 0.05,
       sdtab = data.frame(ID = 1, TIME = 1, DV = 5), bundle = as.raw(1))
}





# --- merge helper (pure R) ---------------------------------------------------

.attach <- getFromNamespace(".ferx_attach_npde", "ferx")




# --- end-to-end against the real engine --------------------------------------








test_that("ferx_calc_npde requires a fit with theta/omega/sigma", {
  expect_error(ferx_calc_npde(list()), "theta, omega, and sigma")
})
test_that("ferx_calc_npde rejects non-positive or non-scalar nsim", {
  expect_error(ferx_calc_npde(.stub_fit(), nsim = 0L),  "positive integer")
  expect_error(ferx_calc_npde(.stub_fit(), nsim = -5L), "positive integer")
  expect_error(ferx_calc_npde(.stub_fit(), nsim = c(10L, 20L)), "positive integer")
})
test_that("ferx_calc_npde rejects a non-integer / NA seed", {
  expect_error(ferx_calc_npde(.stub_fit(), seed = c(1L, 2L)), "single integer or NULL")
  expect_error(ferx_calc_npde(.stub_fit(), seed = NA_integer_), "single integer or NULL")
})
test_that("ferx_calc_npde rejects a negative seed (would silently become the default)", {
  expect_error(ferx_calc_npde(.stub_fit(), seed = -5L), "non-negative integer")
})
test_that("ferx_calc_npde errors when no model/data path is available", {
  testthat::local_mocked_bindings(.ferx_fit_canonical = function(fit, who) fit)
  expect_error(ferx_calc_npde(.stub_fit()), "model file")
  expect_error(ferx_calc_npde(.stub_fit(), model = "no_such.ferx"), "model file")
})
test_that(".ferx_attach_npde copies NPDE/NPD positionally and keeps existing cols", {
  sd <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1), DV = c(5, 6, 7))
  np <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1),
                   NPDE = c(0.1, 0.2, 0.3), NPD = c(0.4, 0.5, 0.6))
  out <- .attach(sd, np)
  expect_equal(out$NPDE, c(0.1, 0.2, 0.3))
  expect_equal(out$NPD,  c(0.4, 0.5, 0.6))
  expect_true(all(c("DV", "NPDE", "NPD") %in% names(out)))
})
test_that(".ferx_attach_npde errors (not silently joins) when row counts differ", {
  sd <- data.frame(ID = c(1, 9), TIME = c(1, 9))
  np <- data.frame(ID = 1, TIME = 1, NPDE = 0.1, NPD = 0.2)
  expect_error(.attach(sd, np), "cannot align")
})
test_that(".ferx_attach_npde errors when rows do not line up by ID/TIME", {
  sd <- data.frame(ID = c(2, 1), TIME = c(1, 2))
  np <- data.frame(ID = c(1, 2), TIME = c(2, 1), NPDE = c(0.9, 0.1), NPD = c(0.8, 0.2))
  expect_error(.attach(sd, np), "do not line up")
})
test_that("ferx_calc_npde returns the fit with NPDE/NPD in sdtab, one row per obs", {
  ex  <- ferx_example("warfarin")
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 1L,
                   model = ex$model, data = ex$data)
  expect_s3_class(fit, "ferx_fit")
  expect_true(all(c("NPDE", "NPD") %in% names(fit$sdtab)))

  dat   <- read.csv(ex$data)
  n_obs <- sum(dat$EVID == 0, na.rm = TRUE)
  expect_equal(nrow(fit$sdtab), n_obs)
})
test_that("ferx_calc_npde NPD is finite and roughly mean-zero on a sensible fit", {
  ex  <- ferx_example("warfarin")
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 500L, seed = 1L,
                   model = ex$model, data = ex$data)
  expect_true(all(is.finite(fit$sdtab$NPD)))
  expect_lt(abs(mean(fit$sdtab$NPD)), 0.3)
})
test_that("ferx_calc_npde is reproducible for a fixed seed", {
  ex <- ferx_example("warfarin")
  a  <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  b  <- ferx_calc_npde(warfarin_fit(), nsim = 200L, seed = 42L,
                  model = ex$model, data = ex$data)
  expect_equal(a$sdtab$NPDE, b$sdtab$NPDE)
  expect_equal(a$sdtab$NPD,  b$sdtab$NPD)
})
test_that("ferx_calc_npde uses fit$model_path / fit$data_path by default", {
  # The cached fit records the paths it was run on, so no explicit model/data
  # is needed.
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 100L, seed = 1L)
  expect_true(all(c("NPDE", "NPD") %in% names(fit$sdtab)))
})
test_that("ferx_calc_npde errors when the fit carries no sdtab", {
  bad <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2),
              sigma = 0.05, bundle = as.raw(1))  # no sdtab
  testthat::local_mocked_bindings(.ferx_fit_canonical = function(fit, who) fit)
  expect_error(ferx_calc_npde(bad), "sdtab` is empty")
})
test_that("ferx_calc_npde surfaces a clean error when the engine returns NULL", {
  # The FFI raises on failure (#385), so a NULL table is not expected; should one
  # arrive anyway, ferx_calc_npde must say so rather than crash in the alignment step.
  fit <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2),
              sigma = 0.05, sdtab = data.frame(ID = 1, TIME = 1, DV = 5),
              model_path = tempfile(fileext = ".ferx"),
              data_path  = tempfile(fileext = ".csv"), bundle = as.raw(1))
  file.create(fit$model_path, fit$data_path)
  testthat::local_mocked_bindings(ferx_rust_npde_from_fit = function(...) NULL,
                                  .ferx_fit_canonical = function(fit, who) fit)
  expect_error(ferx_calc_npde(fit), "engine returned no NPDE table")
})

# --- #416: the fit's own data selection -----------------------------------

test_that("ferx_calc_npde on an `ignore =` fit scores the fit's rows (#416)", {
  # Before #416 the re-read kept the two records the fit dropped and the
  # alignment check refused: "110 row(s) but fit$sdtab has 108".
  fit <- rs_fit(ignore = rs_ignore)
  out <- ferx_calc_npde(fit, nsim = 50L, seed = 1L)
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
  expect_true(all(is.finite(out$sdtab$NPD)))
})

test_that("ferx_calc_npde refuses a legacy fit with a record-only `ignore =` (#416)", {
  fit <- rs_legacy(rs_fit(ignore = rs_ignore))
  e <- rs_steps(fit)$npde
  expect_s3_class(e, "error")
  msg <- conditionMessage(e)
  expect_match(msg, "ferx_calc_npde: this fit predates", fixed = TRUE)
  expect_match(msg, "`ignore: EVID == 0 && DV < 1.0`", fixed = TRUE)
  expect_match(msg, "settings = list(npde_nsim = 20)", fixed = TRUE)
  # Not the alignment check's cause, which names the wrong one.
  expect_no_match(msg, "differ from the fit", fixed = TRUE)
})

test_that("ferx_calc_npde runs on a legacy fit whose selection the model file states (#416)", {
  fit <- rs_legacy(warfarin_sel_fit())
  out <- rs_steps(fit)$npde
  expect_false(inherits(out, "condition"))
})

test_that("ferx_calc_npde runs on a legacy fit whose record-only clause fired nothing (#416)", {
  fit <- rs_legacy(rs_fit(ignore = "DV < -1"))
  out <- rs_steps(fit)$npde
  expect_false(inherits(out, "condition"))
})

# --- #526 review: what npde reads, and what its legacy guard compares --------

test_that("ferx_calc_npde refuses a legacy fit whose model file was edited (#526)", {
  # npde has no hash check of its own: without the guard's verified read it
  # re-read the data with the dropped late dose put back, and the alignment
  # check, which sees only observation rows, let it through.
  fit <- rs_edited_model(rs_legacy(rs_late_dose_fit()))
  expect_identical(fit$exclusions$fired_ignore, "ignore: EVID == 1 && TIME > 400")
  e <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 5L), error = function(e) e)
  expect_s3_class(e, "error")
  msg <- conditionMessage(e)
  expect_match(msg, "ferx_calc_npde: model hash mismatch", fixed = TRUE)
  expect_no_match(msg, "predates the record", fixed = TRUE)
  # Coded since ferx-core #2027 (#2451), the same code as an edited data file.
  expect_identical(e$code, "E_FIT_INPUT_CHANGED")
})

test_that("ferx_calc_npde on a `model =` override keeps that file's reader settings (#526)", {
  # The override reads `DV` from `CONC` through its own `[data]` block; the
  # fit's selection is replayed on it. Replaying the whole record would read
  # with the fit's (empty) column map and find no DV; dropping the selection
  # would put the two records back (110 rows, not 108).
  fit <- rs_fit(ignore = rs_ignore)
  ex <- ferx_example("warfarin")
  rows <- utils::read.csv(ex$data)
  names(rows)[names(rows) == "DV"] <- "CONC"
  data <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  model <- withr::local_tempfile(fileext = ".ferx")
  writeLines(c(readLines(ex$model), "", "[data]",
               paste0("  path = ", gsub("/+", "/", normalizePath(data))),
               "  dv = CONC"), model)
  out <- ferx_calc_npde(fit, nsim = 20L, seed = 5L, model = model, data = data)
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
})

test_that("ferx_calc_npde's legacy guard compares the `model =` file it reads (#526)", {
  fit <- rs_legacy(rs_fit(ignore = rs_ignore))
  ex <- ferx_example("warfarin")
  # A file that states the fit's clause: nothing unstated, it runs.
  stating <- withr::local_tempfile(fileext = ".ferx")
  writeLines(c(readLines(ex$model), "", "[data_selection]",
               "  ignore = EVID == 0 && DV < 1.0"), stating)
  out <- ferx_calc_npde(fit, nsim = 20L, seed = 5L, model = stating)
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
  # A moved `fit$model_path` is not what is compared: the `model =` file is,
  # and it does not state the clause.
  fit$model_path <- file.path(tempdir(), "moved-away.ferx")
  e <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 5L, model = ex$model),
                error = function(e) e)
  expect_match(conditionMessage(e), "`ignore: EVID == 0 && DV < 1.0`", fixed = TRUE)
})

test_that("ferx_calc_npde names a `data =` override it cannot check (#526)", {
  fit <- rs_legacy(rs_fit(ignore = rs_ignore))
  data <- withr::local_tempfile(fileext = ".csv")
  file.copy(ferx_example("warfarin")$data, data)
  e <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 5L, data = data),
                error = function(e) e)
  msg <- conditionMessage(e)
  expect_match(msg, "ferx_calc_npde: this fit predates", fixed = TRUE)
  expect_match(msg, paste0("`data = \"", data, "\"` cannot be checked"), fixed = TRUE)
})

test_that("ferx_calc_npde refuses a legacy fit read with a `settings =` iov_column (#526)", {
  fit <- rs_legacy(rs_iov_fit())
  e <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 5L), error = function(e) e)
  expect_match(conditionMessage(e), '`settings = list(iov_column = "VISIT")`',
               fixed = TRUE)
})

# --- ferx-core #1814: the fit's population fingerprint -----------------------

test_that("ferx_calc_npde(fit) refuses data that are not the fit's (ferx-core #1814)", {
  ex <- ferx_example("warfarin")
  data <- withr::local_tempfile(fileext = ".csv")
  file.copy(ex$data, data)
  fit <- suppressWarnings(ferx_fit(ex$model, data, method = "gn",
                                   covariance = FALSE, verbose = FALSE))
  expect_false(is.null(fit$population_fingerprint))
  # Premise: on the fit's own data it runs, and says nothing.
  expect_no_warning(out <- ferx_calc_npde(fit, nsim = 20L, seed = 1L))
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))

  # One observation edited in place: the same rows, IDs and times, so the
  # sdtab alignment check alone let it through.
  rows <- utils::read.csv(data, na.strings = ".")
  obs <- which(rows$EVID == 0)[1L]
  rows$DV[obs] <- rows$DV[obs] * 2
  utils::write.csv(rows, data, row.names = FALSE, quote = FALSE, na = ".")
  e <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 1L), error = function(e) e)
  expect_s3_class(e, "ferx_engine_error")
  expect_identical(e$code, "E_POPULATION_MISMATCH")
  expect_match(conditionMessage(e), "npde error: ", fixed = TRUE)

  # `data =` naming the fit's own file is still the fit's data, checked too.
  e2 <- tryCatch(ferx_calc_npde(fit, nsim = 20L, seed = 1L, data = data),
                 error = function(e) e)
  expect_s3_class(e2, "ferx_engine_error")
  expect_identical(e2$code, "E_POPULATION_MISMATCH")

  # Another file is not the fit's population by construction: not checked,
  # so the same edited rows under another name run.
  other <- withr::local_tempfile(fileext = ".csv")
  file.copy(data, other)
  out2 <- ferx_calc_npde(fit, nsim = 20L, seed = 1L, data = other)
  expect_identical(nrow(out2$sdtab), nrow(fit$sdtab))
})

test_that("ferx_calc_npde(fit) verifies a level-block fit's population after binding it (ferx-core #1814)", {
  # `verify` compares the population as the model scores it, so a theta level
  # block needs the population bound with the fit's layout first. Checked with
  # the fingerprint on (no `model =`, no `data =`): verifying before the bind
  # refuses the fit's own data.
  ex <- ferx_example("mbma_placebo")
  fit <- suppressWarnings(ferx_fit(ex$model, ex$data, verbose = FALSE))
  expect_false(is.null(fit$population_fingerprint))
  expect_gt(nrow(fit$theta_levels), 0L)
  expect_no_warning(out <- ferx_calc_npde(fit, nsim = 50L, seed = 1L))
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
  expect_false(anyNA(out$sdtab$NPD))
})

test_that("ferx_calc_npde(fit) raises the fingerprint check's note as a warning (#559)", {
  # A fingerprint of another scheme cannot be compared: the engine checks
  # nothing, returns its note in the glue's `npde_notes` attribute, and R raises
  # it. NPDE is still computed: the note is not a refusal.
  fit <- foreign_scheme_fit(warfarin_fit())
  w <- character(0)
  out <- withCallingHandlers(
    ferx_calc_npde(fit, nsim = 20L, seed = 1L),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(w, paste0("ferx_calc_npde: ", foreign_scheme_note))
  expect_identical(nrow(out$sdtab), nrow(fit$sdtab))
  expect_true(all(c("NPDE", "NPD") %in% names(out$sdtab)))
  expect_true(all(is.finite(out$sdtab$NPD)))
})

test_that("ferx_calc_npde on a [mixture] fit draws its reference at the fitted class Omega (#672)", {
  mod <- test_path("fixtures", "mixture_override.ferx")
  dat <- test_path("fixtures", "mixture_iv.csv")
  fit <- mixture_override_fit()
  expect_length(fit$mixture_overrides, 1L)
  npd_sd <- function(f) {
    out <- muffle_mixl_lint(ferx_calc_npde(f, nsim = 300L, seed = 4L,
                                           model = mod, data = dat))
    expect_true(all(is.finite(out$sdtab$NPD)))
    stats::sd(out$sdtab$NPD)
  }
  # The reference simulation reads the class-2 override from the fit: widening
  # it twentyfold spreads class 2's reference draws, so the observations fall
  # nearer their middle and NPD narrows. A reference built at the model file's
  # classes reads neither value and gives the same NPD for both.
  wide <- fitrx_edit_fit(fit, function(w) {
    w$mixture_overrides <- lapply(w$mixture_overrides, function(v) 20 * v)
    w
  })
  base_sd <- npd_sd(fit)
  wide_sd <- npd_sd(wide)
  expect_lt(wide_sd, 0.85 * base_sd,
            label = paste(signif(wide_sd, 3), "vs", signif(base_sd, 3)))
})

# --- EPRED / PWRES / EWRES (ferx-core #2871) ---------------------------------

test_that(".ferx_attach_npde carries EPRED/EWRES/PWRES, and only the columns the table has", {
  sd <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1), DV = c(5, 4, 6))
  np <- data.frame(ID = c(1, 1, 2), TIME = c(1, 2, 1),
                   NPDE = c(0.1, 0.2, 0.3), NPD = c(0.4, 0.5, 0.6),
                   EPRED = c(5.5, 4.5, 5.0), EWRES = c(0.7, NA, -0.2),
                   PWRES = c(0.7, NA, -0.2))
  out <- .attach(sd, np)
  expect_identical(out$EPRED, c(5.5, 4.5, 5.0))
  expect_identical(out$PWRES, c(0.7, NA, -0.2))
  expect_identical(out$EWRES, out$PWRES)
  # A table from an engine without them attaches NPDE/NPD and invents nothing.
  old <- .attach(sd, np[, c("ID", "TIME", "NPDE", "NPD")])
  expect_false(any(c("EPRED", "EWRES", "PWRES") %in% names(old)))
  expect_identical(old$NPD, c(0.4, 0.5, 0.6))
  # The Gaussian-row split carries every column, NA on the other rows.
  sd4 <- data.frame(ID = c(1, 1, 2, 2), TIME = c(1, 2, 1, 2))
  np3 <- np[, setdiff(names(np), "DV")]
  np3$TIME <- c(1, 2, 2)
  np3$ID <- c(1, 1, 2)
  out4 <- .attach(sd4, np3, gaussian_rows = c(1L, 2L, 4L))
  expect_identical(out4$EPRED, c(5.5, 4.5, NA, 5.0))
  expect_identical(out4$PWRES, c(0.7, NA, NA, -0.2))
})

test_that("ferx_calc_npde adds EPRED and PWRES/EWRES from the same simulation", {
  fit <- ferx_calc_npde(warfarin_fit(), nsim = 500L, seed = 1L)
  sd <- fit$sdtab
  expect_true(all(c("EPRED", "EWRES", "PWRES") %in% names(sd)))
  # nsim exceeds every subject's observation count, so nothing is masked.
  expect_true(all(is.finite(sd$EPRED)))
  expect_true(all(is.finite(sd$PWRES)))
  # One residual under two names.
  expect_identical(sd$EWRES, sd$PWRES)
  # EPRED is the Monte-Carlo mean, not the eta = 0 PRED: under exp() IIV the two
  # differ, and both are concentrations on the same scale.
  expect_false(isTRUE(all.equal(sd$EPRED, sd$PRED)))
  expect_lt(abs(stats::median(sd$EPRED / sd$PRED) - 1), 0.5)
  expect_lt(abs(mean(sd$PWRES)), 0.5)
})

test_that("a fit with npde_nsim carries EPRED/PWRES that ferx_calc_npde reproduces", {
  ex  <- ferx_example("warfarin")
  # The call-time overrides of the example's [fit_options] warn; they are the point.
  fit <- suppressWarnings(
    ferx_fit(ex$model, ex$data, method = "focei", verbose = FALSE,
             covariance = FALSE,
             settings = list(maxiter = 30L, npde_nsim = 200L, npde_seed = 7L))
  )
  expect_true(all(c("EPRED", "EWRES", "PWRES") %in% names(fit$sdtab)))
  post <- ferx_calc_npde(fit, nsim = 200L, seed = 7L)
  expect_equal(post$sdtab$EPRED, fit$sdtab$EPRED)
  expect_equal(post$sdtab$PWRES, fit$sdtab$PWRES)
})

# --- single-event time-to-event records (ferx-core #1889) --------------------
#
# The oracle is the engine's own fit-time score: a fit with npde_nsim / npde_seed
# writes `TteObsDiagnostic::npde`, and ferx_calc_npde() at the same nsim and seed
# must put the same numbers on the same rows.

pktte_npde_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      ex <- ferx_example("pktte_joint")
      plain <- suppressWarnings(ferx_fit(ex$model, ex$data, covariance = FALSE,
                                         verbose = FALSE))
      at_fit <- suppressWarnings(ferx_fit(
        ex$model, ex$data, covariance = FALSE, verbose = FALSE,
        settings = list(npde_nsim = 300L, npde_seed = 5L)
      ))
      cache <<- list(plain = plain, at_fit = at_fit)
    }
    cache
  }
})

test_that("ferx_calc_npde scores the TTE rows as a fit-time npde run does", {
  fits <- pktte_npde_fits()
  post <- ferx_calc_npde(fits$plain, nsim = 300L, seed = 5L)
  b <- ferx:::.fitrx_sdtab_blocks(post)
  expect_gt(length(b$tte), 0L)
  sd <- post$sdtab
  # Every TTE record scored, NPDE == NPD there, and finite.
  expect_true(all(is.finite(sd$NPDE[b$tte])))
  expect_identical(sd$NPD[b$tte], sd$NPDE[b$tte])
  # A single event has no moment to weight: EPRED / PWRES stay NA there,
  # while the Gaussian rows carry them.
  expect_true(all(is.na(sd$EPRED[b$tte])))
  expect_true(all(is.na(sd$PWRES[b$tte])))
  expect_true(all(is.finite(sd$EPRED[b$gaussian])))
  # The engine's own fit-time scores, row for row.
  ref <- fits$at_fit$sdtab
  bref <- ferx:::.fitrx_sdtab_blocks(fits$at_fit)
  expect_identical(bref$tte, b$tte)
  expect_equal(sd$NPDE[b$tte], ref$NPDE[bref$tte], tolerance = 1e-12)
  expect_equal(sd$NPDE[b$gaussian], ref$NPDE[bref$gaussian], tolerance = 1e-12)
})

test_that("TTE NPDE survives save/load and reaches core's predictions_tte.csv", {
  post <- ferx_calc_npde(pktte_npde_fits()$plain, nsim = 300L, seed = 5L)
  b <- ferx:::.fitrx_sdtab_blocks(post)
  path <- withr::local_tempfile(fileext = ".fitrx")
  ferx_save_fit(post, path)
  back <- ferx_load_fit(path)
  expect_identical(back$sdtab$NPDE, post$sdtab$NPDE)
  # Core's own TTE table (read by core's loader and the CLI) carries the
  # scores too, in TTE-row order.
  staging <- withr::local_tempdir()
  utils::unzip(path, files = "predictions_tte.csv", exdir = staging)
  tte <- utils::read.csv(file.path(staging, "predictions_tte.csv"))
  expect_true("NPDE" %in% names(tte))
  expect_equal(tte$NPDE, post$sdtab$NPDE[b$tte], tolerance = 1e-12)
})

test_that(".ferx_attach_tte_npde puts the scores on the TTE rows, checked by ID/TIME", {
  attach_tte <- getFromNamespace(".ferx_attach_tte_npde", "ferx")
  sd <- data.frame(ID = c(1, 1, 2, 1, 2), TIME = c(1, 2, 1, 7, 9),
                   NPDE = c(0.1, 0.2, 0.3, NA, NA), NPD = c(0.4, 0.5, 0.6, NA, NA))
  tbl <- data.frame(ID = c(1, 2), TIME = c(7, 9), CMT = 3L, NPDE = c(-1, 1.5))
  out <- attach_tte(sd, tbl, 4:5)
  expect_identical(out$NPDE, c(0.1, 0.2, 0.3, -1, 1.5))
  expect_identical(out$NPD, c(0.4, 0.5, 0.6, -1, 1.5))
  # Nothing to attach leaves the table alone.
  expect_identical(attach_tte(sd, NULL, 4:5), sd)
  expect_error(attach_tte(sd, tbl, 5L), "cannot align")
  bad <- tbl
  bad$TIME[2] <- 8
  expect_error(attach_tte(sd, bad, 4:5), "do not line up")
})

test_that("an RTTE endpoint is skipped with W_NPDE_TTE_SKIPPED, Gaussian rows still scored", {
  m <- tempfile(fileext = ".ferx")
  writeLines(c("[parameters]", "  theta TVLAMBDA(0.1, 0.001, 10.0)",
               "  omega ETA_LAMBDA ~ 0.09", "[event_model]", "  cmt    = 2",
               "  type   = rtte", "  family = exponential",
               "  scale  = TVLAMBDA * exp(ETA_LAMBDA)"), m)
  design <- data.frame(ID = 1:30, TIME = 30, DV = 0, EVID = 0L, AMT = 0,
                       CMT = 2L, MDV = 0L)
  sim <- suppressWarnings(ferx_simulate(m, design, n_sim = 1L, seed = 2L,
                                        horizon = 30))
  d <- data.frame(ID = as.integer(sim$ID), TIME = sim$TIME, DV = sim$OBSERVED,
                  EVID = 0L, AMT = 0, CMT = 2L, MDV = 0L)
  fit <- suppressWarnings(ferx_fit(m, d, covariance = FALSE, verbose = FALSE))
  expect_warning(post <- ferx_calc_npde(fit, nsim = 100L, seed = 1L),
                 "W_NPDE_TTE_SKIPPED: the time-to-event records were not scored")
  expect_true(all(is.na(post$sdtab$NPDE)))
})
