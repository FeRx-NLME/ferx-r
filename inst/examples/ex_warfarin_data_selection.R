library(ferx)

# Demonstrates the [data_selection] block for excluding records at read
# time without modifying the CSV -- equivalent to NONMEM $DATA IGNORE=.
# The model drops observations below a surrogate LLOQ of 1.0 mg/L.

ex <- ferx_example("warfarin_data_selection")
ferx_model_show(ex$model)

# --- R-side preview before fitting ----------------------------------------
# ferx_apply_selection() applies the same ignore/accept logic in pure R so you
# can inspect which records will be excluded before committing to a fit.
# Passing the model applies its [data_selection] block, so the preview keeps
# the records the fit will.
sel <- ferx_apply_selection(ex$data, model = ex$model)
cat("Retained:", nrow(sel), "  Excluded:", nrow(ferx_apply_selection(sel, excluded = TRUE)), "\n")

# Inspect the excluded rows and the reason they were dropped:
ferx_apply_selection(sel, excluded = TRUE)[, c("ID", "TIME", "DV", ".exclude_reason")]

# --- Fit with model-file data_selection -----------------------------------
fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Exclusion summary is reported in the print output and in fit$exclusions:
fit$exclusions$n_records_total     # total records read
fit$exclusions$n_obs_excluded      # observations excluded by the filter
fit$exclusions$fired_ignore        # which ignore conditions fired

# The same filter can be applied via the R API instead of the model file: the
# plain warfarin model states no [data_selection] block, so the filter comes
# from the ferx_data and is recorded as a call-time setting. (gradient = "fd"
# matches the [fit_options] of the data-selection model.)
plain <- ferx_example("warfarin")$model
fit2 <- ferx_fit(plain, ferx_apply_selection(ex$data, ignore = "DV < 1.0"),
                 gradient = "fd")
fit2$call_settings$ignore          # "DV < 1.0"
isTRUE(abs(fit$ofv - fit2$ofv) < 1e-4)
