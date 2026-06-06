library(ferx)

# Demonstrates the [data_selection] block for excluding records at read
# time without modifying the CSV -- equivalent to NONMEM $DATA IGNORE=.
# The model drops observations below a surrogate LLOQ of 1.0 mg/L.

ex <- ferx_example("warfarin_data_selection")
ferx_model_show(ex$model)

# --- R-side preview before fitting ----------------------------------------
# ferx_selection() applies the same ignore/accept logic in pure R so you
# can inspect which records will be excluded before committing to a fit.
sel <- ferx_selection(ex$data, ignore = "DV < 1.0")
cat("Retained:", nrow(sel), "  Excluded:", nrow(ferx_selection_excluded(sel)), "\n")

# Inspect the excluded rows and the reason they were dropped:
ferx_selection_excluded(sel)[, c("ID", "TIME", "DV", ".exclude_reason")]

# --- Fit with model-file data_selection -----------------------------------
fit <- ferx_fit(ex$model, ex$data)
print(fit)

# Exclusion summary is reported in the print output and in fit$exclusions:
fit$exclusions$n_records_total     # total records read
fit$exclusions$n_obs_excluded      # observations excluded by the filter
fit$exclusions$fired_ignore        # which ignore conditions fired

# The same filter can be applied via the R API instead of the model file:
fit2 <- ferx_fit(ex$model, ferx_selection(ex$data, ignore = "DV < 1.0"))
isTRUE(abs(fit$ofv - fit2$ofv) < 1e-4)
