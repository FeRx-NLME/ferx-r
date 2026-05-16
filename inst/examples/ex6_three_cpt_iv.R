library(ferx)
library(ggplot2)
ex <- ferx_example("three_cpt_iv")

fit <- ferx_fit(ex$model, ex$data)
fit

# Goodness-of-fit
sdtab <- fit$sdtab
ggplot(sdtab, aes(IPRED, DV)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = "steelblue") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DV vs IPRED (three-compartment IV)")
