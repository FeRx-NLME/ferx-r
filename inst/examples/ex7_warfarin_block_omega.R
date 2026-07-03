library(ferx)
ex <- ferx_example("warfarin_block_omega")

fit <- ferx_fit(ex$model, ex$data)
fit

# Show off-diagonal correlations
fit$cor_matrix
