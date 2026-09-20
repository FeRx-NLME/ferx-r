# A model or dataset the engine refuses is an R error on every entry point, not
# console text and a NULL (#385). `ferx_fit()` already behaved this way (#367);
# these are the functions that did not. Fixtures and expectations live in
# helper-engine-errors.R.
#
# Each entry point gets its own test per failure, so a regression in one glue
# function's error arm names that function.

# -- Parse stage --------------------------------------------------------------

for (nm in engine_entry_point_names()) {
  local({
    nm <- nm

    test_that(paste(nm, "raises E_UNKNOWN_BLOCK for a block the engine does not know"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(engine_entry_point(nm)(unknown_block_model(), ex$data))
      expect_coded_refusal(probe, "Unknown block `[not_a_block]`", "E_UNKNOWN_BLOCK")
      expect_identical(probe$cond$block, "not_a_block")
    })

    test_that(paste(nm, "raises a malformed [data_selection] clause"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(engine_entry_point(nm)(bad_selection_model(), ex$data))
      expect_coded_refusal(probe, "malformed filter expression 'DV <'", "E_PARSE")
    })

    test_that(paste(nm, "quotes a '%' from the model back as written"), {
      # The message is `Rf_error()`'s format string: unescaped, `5%': r...` was
      # read as a conversion and came back as "5right-hand side '5-2139062144s".
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(engine_entry_point(nm)(percent_selection_model(), ex$data))
      expect_coded_refusal(
        probe,
        "malformed filter expression 'DV < 5%': right-hand side '5%' is not a number",
        "E_PARSE"
      )
    })
  })
}

# -- Data-read stage ----------------------------------------------------------

for (nm in engine_entry_point_names()) {
  local({
    nm <- nm

    test_that(paste(nm, "raises a dataset the reader refuses"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(engine_entry_point(nm)(ex$model, no_time_data()))
      expect_coded_refusal(probe, "Missing TIME column", "E_DATA")
    })
  })
}

# -- A broken precondition of predict / simulate ------------------------------
#
# A panic at the pinned engine, an `Err` after FeRx-NLME/ferx-core#898; the
# caller must not be able to tell. The phrase is the validation diagnostic's
# own, which the panic text carries.

for (nm in setdiff(engine_entry_point_names(), "ferx_calc_npde()")) {
  local({
    nm <- nm

    test_that(paste(nm, "raises E_DOSE_CMT_NOT_INFUSABLE for an infusion into CMT = 0"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(engine_entry_point(nm)(ex$model, infusion_into_cmt0_data()))
      expect_coded_refusal(probe, "infusion into compartment 0",
                           "E_DOSE_CMT_NOT_INFUSABLE")
    })
  })
}

test_that("ferx_calc_npde() raises an infusion into CMT = 0, uncoded while the engine words it differently", {
  # NPDE reaches the dose without the up-front compartment check, so the engine
  # panics deeper, in different words from the validation diagnostic. A code is
  # attached only on a text match (see the mislabelling tests below), so this
  # one is an ordinary error until the engine reports it the way it reports
  # the others.
  ex    <- ferx_example("warfarin")
  probe <- engine_error_probe(
    engine_entry_point("ferx_calc_npde()")(ex$model, infusion_into_cmt0_data())
  )
  expect_uncoded_refusal(probe, "infusion into compartment 0")
})

# -- Simulate / predict stage -------------------------------------------------
#
# The engine's validation pass has no diagnostic for these, so they arrive as
# ordinary errors: the engine's prose, no code. `test-ferx_simulate.R` holds
# two of them from the user's side.

test_that("ferx_simulate() raises the engine's refusal to sample TTE without a horizon", {
  ex    <- ferx_example("pktte_joint")
  probe <- engine_error_probe(ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L))
  expect_uncoded_refusal(probe, "requires a finite, positive administrative horizon")
})

test_that("ferx_simulate(fit = ) raises the engine's refusal to sample TTE without a horizon", {
  ex  <- ferx_example("pktte_joint")
  fit <- ferx_fit(ex$model, ex$data, method = "focei", covariance = FALSE,
                  verbose = FALSE, settings = list(maxiter = 3L))
  probe <- engine_error_probe(
    ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, fit = fit)
  )
  expect_uncoded_refusal(probe, "requires a finite, positive administrative horizon")
})

test_that("the fit = forms raise a fit that cannot drive the model's kappa", {
  ex  <- ferx_example("warfarin_iov")
  fit <- ferx_fit(ex$model, ex$data, covariance = FALSE, verbose = FALSE,
                  settings = list(maxiter = 3L))
  fit$omega_iov <- NULL

  calls <- list(
    function() ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L, fit = fit),
    function() ferx_predict(ex$model, ex$data, fit = fit),
    function() ferx_predict_survival(ex$model, ex$data, times = c(1, 2), fit = fit),
    function() ferx_calc_npde(fit, nsim = 20L, seed = 1L)
  )
  for (call_it in calls) {
    expect_uncoded_refusal(engine_error_probe(call_it()), "the fit carries no omega_iov")
  }
})

# -- A failure is never labelled with some other finding's code ----------------
#
# Each case pairs a failure validation cannot see with a dataset that carries
# exactly one validation error of its own (the CMT = 0 infusion). Labelling a
# failure with "the one error validation found" got every one of these wrong.

test_that("a fit that does not fit the model is not labelled with the data's diagnostic", {
  ex    <- ferx_example("warfarin")
  short <- warfarin_fit_cov()
  short$theta <- short$theta[-1L]
  data  <- infusion_into_cmt0_data()

  probe <- engine_error_probe(ferx_predict(ex$model, data, fit = short))
  expect_uncoded_refusal(probe, "theta length 2 does not match model (3 expected)")

  probe <- engine_error_probe(
    ferx_calc_npde(short, nsim = 20L, seed = 1L, model = ex$model, data = data)
  )
  expect_uncoded_refusal(probe, "theta length 2 does not match model (3 expected)")

  probe <- engine_error_probe(ferx_simulate_with_uncertainty(
    ex$model, data, fit = short, n_uncertainty_draws = 2L, n_sim_per_draw = 1L
  ))
  expect_uncoded_refusal(probe, "theta length 2 does not match model (3 expected)")
})

test_that("a missing horizon is not labelled with the data's diagnostic", {
  # Not a `fit =` form and not a fit-shape message: the label was wrong here too.
  ex    <- ferx_example("pktte_joint")
  probe <- engine_error_probe(
    ferx_simulate(ex$model, infusion_into_cmt0_data("pktte_joint"), n_sim = 1L, seed = 1L)
  )
  expect_uncoded_refusal(probe, "requires a finite, positive administrative horizon")
})

# -- One handler for every entry point ----------------------------------------

test_that("#385's script: a loop over candidate models stops on the refused one", {
  ex     <- ferx_example("warfarin")
  models <- c(ex$model, unknown_block_model())
  codes  <- vapply(models, function(m) {
    tryCatch(
      {
        suppressWarnings(ferx_predict(m, ex$data))
        "ran"
      },
      ferx_engine_error = function(e) e$code
    )
  }, character(1), USE.NAMES = FALSE)
  expect_identical(codes, c("ran", "E_UNKNOWN_BLOCK"))
})
