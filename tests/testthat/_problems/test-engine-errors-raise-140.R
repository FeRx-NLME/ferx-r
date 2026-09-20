# Extracted from test-engine-errors-raise.R:140

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "ferx", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
ex    <- ferx_example("warfarin")
short <- warfarin_fit_cov()
short$theta <- short$theta[-1L]
data  <- infusion_into_cmt0_data()
probe <- engine_error_probe(ferx_predict(ex$model, data, fit = short))
expect_uncoded_refusal(probe, "theta length 2 does not match model (3 expected)")
