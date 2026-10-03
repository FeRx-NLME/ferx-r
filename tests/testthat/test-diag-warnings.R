# Reader findings reach the caller only through ferx-core's one filter
# (ferx-r #426, ferx-core #1645).
#
# Four glue paths used to relay `Population.warnings` raw, past the filter
# (`reader_warning_suppressed`) that `ferx_fit()` and `ferx_simulate()` apply.
# On a compartment-free model that raw list carries `W_CMT_DEFAULTED` and
# `W_NO_DOSES`, both moot there, so `ferx_predict()` warned where fit and
# simulate stayed quiet. Every path now relays the `warnings` of a core
# `*_diag` result instead.

reader_codes <- function(x) {
  w <- attr(x, "simulation_warnings", exact = TRUE)
  expect_type(w, "character")
  sub(":.*$", "", w)
}

# A copy of a bundled model with a theta nothing references: the parser warns
# about it, a finding only core's bundle carries (the data reader never sees
# the model), so it shows the relayed list is the bundle and not
# `Population.warnings`.
with_unused_theta <- function(model) {
  lines <- readLines(model)
  at <- grep("^\\[parameters\\]", lines)
  path <- tempfile(fileext = ".ferx")
  writeLines(append(lines, "  theta UNUSED(1.0, 0.1, 10.0)", after = at), path)
  path
}

has_unused_theta_note <- function(x) {
  any(grepl("theta 'UNUSED' is declared", attr(x, "simulation_warnings", exact = TRUE),
            fixed = TRUE))
}

# A one-compartment model read against the same dose-free data: here no doses
# is a real finding, so the filter must let `W_NO_DOSES` through.
write_one_cpt_twin <- function() {
  path <- tempfile(fileext = ".ferx")
  writeLines(c(
    "[parameters]",
    "  theta TVCL(1.0, 0.01, 100.0)",
    "  theta TVV(10.0, 0.1, 1000.0)",
    "  omega ETA_CL ~ 0.04",
    "  sigma ADD ~ 1.0 (variance)",
    "[individual_parameters]",
    "  CL = TVCL * exp(ETA_CL)",
    "  V  = TVV",
    "[structural_model]",
    "  pk one_cpt_iv(cl=CL, v=V)",
    "[error_model]",
    "  DV ~ additive(ADD)"
  ), path)
  path
}

test_that("ferx_predict on a compartment-free model raises no moot reader finding", {
  ex <- ferx_example("emax_timecourse")
  expect_warning(pred <- ferx_predict(ex$model, ex$data), NA)
  codes <- reader_codes(pred)
  expect_false("W_CMT_DEFAULTED" %in% codes)
  expect_false("W_NO_DOSES" %in% codes)
  expect_gt(nrow(pred), 0L)
})

test_that("ferx_predict(fit = ) on a compartment-free model is quiet too", {
  # A separate glue entry point (`ferx_rust_predict_from_fit`); the common one
  # after a fit, so it must not be the half that still warns.
  ex <- ferx_example("emax_timecourse")
  fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = FALSE)
  expect_warning(pred <- ferx_predict(ex$model, ex$data, fit = fit), NA)
  codes <- reader_codes(pred)
  expect_false("W_CMT_DEFAULTED" %in% codes)
  expect_false("W_NO_DOSES" %in% codes)
})

test_that("a one-compartment twin on the same dose-free data still raises W_NO_DOSES", {
  # The other side of the filter: the finding is withheld because the model
  # makes it moot, not because predict stopped relaying reader findings.
  ex <- ferx_example("emax_timecourse")
  expect_warning(
    pred <- ferx_predict(write_one_cpt_twin(), ex$data),
    "W_NO_DOSES"
  )
  expect_true("W_NO_DOSES" %in% reader_codes(pred))
})

test_that("ferx_simulate_with_uncertainty on a compartment-free model raises no moot reader finding", {
  ex <- ferx_example("emax_timecourse")
  fit <- ferx_fit(ex$model, ex$data, method = "gn", covariance = TRUE)
  expect_warning(
    sims <- ferx_simulate_with_uncertainty(
      ex$model, ex$data, fit,
      n_uncertainty_draws = 2L, n_sim_per_draw = 1L, seed = 1L
    ),
    NA
  )
  codes <- reader_codes(sims)
  expect_false("W_CMT_DEFAULTED" %in% codes)
  expect_false("W_NO_DOSES" %in% codes)
  expect_gt(nrow(sims), 0L)
})

test_that("ferx_predict relays the bundle ferx_simulate relays, not only reader findings", {
  # A parse warning reached `ferx_simulate()` and never `ferx_predict()`,
  # which relayed the data reader's list alone.
  ex <- ferx_example("warfarin")
  model <- with_unused_theta(ex$model)
  sim <- suppressWarnings(ferx_simulate(model, ex$data, n_sim = 1L, seed = 1L))
  expect_true(has_unused_theta_note(sim))
  expect_warning(pred <- ferx_predict(model, ex$data), "UNUSED")
  expect_true(has_unused_theta_note(pred))
})

# The adaptive controller supplies the whole regimen, so a dose-free
# observation grid is its normal input. `DV = .` makes every row a design
# point, so the reader also raises `W_DESIGN_DV` -- a finding that must still
# get through.
write_adaptive_grid <- function() {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(ID = 1L, TIME = seq(0, 96, by = 12), DV = ".", EVID = 0L,
               AMT = ".", CMT = 1L, MDV = 0L),
    path, row.names = FALSE, quote = FALSE
  )
  path
}

test_that("adaptive on a dose-free grid raises no W_NO_DOSES and still relays the reader", {
  ex <- ferx_example("adaptive_tdm")
  res <- suppressWarnings(
    ferx_simulate_adaptive(ex$model, write_adaptive_grid(), n_sim = 1L, seed = 1L)
  )
  codes <- reader_codes(res)
  expect_false("W_NO_DOSES" %in% codes)
  expect_true("W_DESIGN_DV" %in% codes)
  expect_gt(nrow(res$doses), 0L)
})

test_that("adaptive relays the engine's bundle, not only reader findings", {
  ex <- ferx_example("adaptive_tdm")
  res <- suppressWarnings(ferx_simulate_adaptive(
    with_unused_theta(ex$model), write_adaptive_grid(), n_sim = 1L, seed = 1L
  ))
  expect_true(has_unused_theta_note(res))
  expect_false("W_NO_DOSES" %in% reader_codes(res))
})
