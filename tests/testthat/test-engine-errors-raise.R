# A model or dataset the engine refuses is an R error on every entry point, not
# console text and a NULL (#385). `ferx_fit()` already behaved this way (#367);
# these are the functions that did not. Fixtures and expectations live in
# helper-engine-errors.R.
#
# Each entry point gets its own test per failure, so a regression in one glue
# function's error arm names that function.

entry_points <- engine_entry_points(warfarin_fit_cov())

# -- Parse stage --------------------------------------------------------------

for (nm in names(entry_points)) {
  local({
    nm <- nm
    call_it <- entry_points[[nm]]

    test_that(paste(nm, "raises E_UNKNOWN_BLOCK for a block the engine does not know"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(call_it(unknown_block_model(), ex$data))
      expect_coded_refusal(probe, "Unknown block `[not_a_block]`", "E_UNKNOWN_BLOCK")
      expect_identical(probe$cond$block, "not_a_block")
    })

    test_that(paste(nm, "raises a malformed [data_selection] clause"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(call_it(bad_selection_model(), ex$data))
      expect_coded_refusal(probe, "malformed filter expression 'DV <'", "E_PARSE")
    })
  })
}

# -- Data-read stage ----------------------------------------------------------

for (nm in names(entry_points)) {
  local({
    nm <- nm
    call_it <- entry_points[[nm]]

    test_that(paste(nm, "raises a dataset the reader refuses"), {
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(call_it(ex$model, no_time_data()))
      expect_coded_refusal(probe, "Missing TIME column", "E_DATA")
    })
  })
}

# -- A broken precondition of predict / simulate ------------------------------

for (nm in names(entry_points)) {
  local({
    nm <- nm
    call_it <- entry_points[[nm]]

    test_that(paste(nm, "raises E_DOSE_CMT_NOT_INFUSABLE for an infusion into CMT = 0"), {
      # A panic at the pinned engine, an `Err` after FeRx-NLME/ferx-core#898;
      # the caller must not be able to tell. The phrase is the validation
      # diagnostic's own, which both the panic text and the `Err` carry.
      ex    <- ferx_example("warfarin")
      probe <- engine_error_probe(call_it(ex$model, infusion_into_cmt0_data()))
      expect_coded_refusal(probe, "infusion into compartment 0",
                           "E_DOSE_CMT_NOT_INFUSABLE")
    })
  })
}

# -- Simulate / predict stage -------------------------------------------------
#
# The engine's validation pass has no diagnostic for these two, so they arrive
# as ordinary errors: the engine's prose, no code. `test-ferx_simulate.R` holds
# the same two cases from the user's side.

test_that("ferx_simulate() raises the engine's refusal to sample TTE without a horizon", {
  ex    <- ferx_example("pktte_joint")
  probe <- engine_error_probe(ferx_simulate(ex$model, ex$data, n_sim = 1L, seed = 1L))
  expect_refusal(probe, "requires a finite, positive administrative horizon")
  expect_false(inherits(probe$cond, "ferx_engine_error"))
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
    probe <- engine_error_probe(call_it())
    expect_refusal(probe, "the fit carries no omega_iov")
    expect_false(inherits(probe$cond, "ferx_engine_error"))
  }
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
