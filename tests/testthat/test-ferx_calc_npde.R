
# ---- header from test-npde.R ----
# Tests for ferx_calc_npde(): post-hoc simulation-based NPDE/NPD from a fit.
#
# The argument-guard tests are pure R and run without the compiled engine: each
# error is raised before the FFI call. The end-to-end tests need the real
# backend and the bundled warfarin example (shared `warfarin_fit()` helper).

# --- argument validation (no engine needed) ----------------------------------


# A minimal structurally-valid fit that passes validate_fit_for_params and the
# sdtab gate but carries no model/data path, so the nsim/seed/file guards fire
# before any FFI call.
.stub_fit <- function() {
  list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2), sigma = 0.05,
       sdtab = data.frame(ID = 1, TIME = 1, DV = 5))
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
              sigma = 0.05)  # no sdtab
  expect_error(ferx_calc_npde(bad), "sdtab` is empty")
})
test_that("ferx_calc_npde surfaces a clean error when the engine returns NULL", {
  # The FFI raises on failure (#385), so a NULL table is not expected; should one
  # arrive anyway, ferx_calc_npde must say so rather than crash in the alignment step.
  fit <- list(theta = c(1, 2), omega = matrix(c(0.1, 0, 0, 0.1), 2, 2),
              sigma = 0.05, sdtab = data.frame(ID = 1, TIME = 1, DV = 5),
              model_path = tempfile(fileext = ".ferx"),
              data_path  = tempfile(fileext = ".csv"))
  file.create(fit$model_path, fit$data_path)
  testthat::local_mocked_bindings(ferx_rust_npde_from_fit = function(...) NULL)
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
